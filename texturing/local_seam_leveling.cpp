/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <math/accum.h>
#include <core/image_io.h>
#include <core/image_tools.h>

#include "progress_counter.h"
#include "texturing.h"
#include "seam_leveling.h"

TEX_NAMESPACE_BEGIN

#define STRIP_SIZE 20

math::Vec3f
mean_color_of_edge_point(std::vector<ProjectedEdgeInfo> projected_edge_infos,
                         std::vector<TexturePatch::ConstPtr> const & texture_patches,
                         float t) {

    assert(0.0f <= t && t <= 1.0f);
    math::Accum<math::Vec3f> color_accum(math::Vec3f(0.0f));

    for (ProjectedEdgeInfo const & projected_edge_info : projected_edge_infos) {
        TexturePatch::ConstPtr texture_patch = texture_patches[projected_edge_info.texture_patch_id];
        if (texture_patch->get_label() == 0) continue;
        math::Vec2f pixel = projected_edge_info.p1 * t + (1.0f - t) * projected_edge_info.p2;
        math::Vec3f color = texture_patch->get_pixel_value(pixel);
        color_accum.add(color, 1.0f);
    }

    math::Vec3f mean_color = color_accum.normalized();
    return mean_color;
}

void
draw_line(math::Vec2f p1, math::Vec2f p2, std::vector<ProjectedEdgeInfo> const & projected_edge_infos,
    std::vector<TexturePatch::ConstPtr> const & texture_patches, TexturePatch::Ptr texture_patch) {
    /* http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm */

    int x0 = std::floor(p1[0] + 0.5f);
    int y0 = std::floor(p1[1] + 0.5f);
    int const x1 = std::floor(p2[0] + 0.5f);
    int const y1 = std::floor(p2[1] + 0.5f);

    float tdx = static_cast<float>(x1 - x0);
    float tdy = static_cast<float>(y1 - y0);
    float length = std::sqrt(tdx * tdx + tdy * tdy);

    int const dx = std::abs(x1 - x0);
    int const dy = std::abs(y1 - y0) ;
    int const sx = x0 < x1 ? 1 : -1;
    int const sy = y0 < y1 ? 1 : -1;
    int err = dx - dy;

    int x = x0;
    int y = y0;
    while (true) {
        math::Vec2i pixel(x, y);

        tdx = static_cast<float>(x1 - x);
        tdy = static_cast<float>(y1 - y);

        /* If the length is zero we sample the midpoint of the projected edge. */
        float t = length != 0.0f ? std::sqrt(tdx * tdx + tdy * tdy) / length : 0.5f;

        // mean color of projections of edges across multiviews
        texture_patch->set_pixel_value(pixel, mean_color_of_edge_point(projected_edge_infos, texture_patches, t));
        if (x == x1 && y == y1)
            break;

        int const e2 = 2 * err;
        if (e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if (e2 < dx) {
            err += dx;
            y += sy;
        }
    }
}

void
local_seam_leveling(UniGraph const & graph,
                    core::TriangleMesh::ConstPtr mesh,
                    VertexProjectionInfos const & vertex_projection_infos,
                    std::vector<TexturePatch::Ptr> * texture_patches) {

    //TODO huge waste of memory...
    std::vector<TexturePatch::ConstPtr> orig_texture_patches(texture_patches->size());
    for (std::size_t i = 0; i < orig_texture_patches.size(); ++i) {
        orig_texture_patches[i] = texture_patches->at(i)->duplicate();
    }

    // number of vertices
    std::size_t const num_vertices = vertex_projection_infos.size();

    {
        /*Find all the seam edges*/
        std::vector<MeshEdge> seam_edges;
        find_seam_edges(graph, mesh, &seam_edges);

        /*draw colors of edges*/
        for (MeshEdge seam_edge : seam_edges){
            std::vector<ProjectedEdgeInfo> projected_edge_infos;

            // find edge projections
            find_mesh_edge_projections(vertex_projection_infos, seam_edge, &projected_edge_infos);

            for (ProjectedEdgeInfo const & projected_edge_info : projected_edge_infos) {
                draw_line(projected_edge_info.p1, projected_edge_info.p2, projected_edge_infos,
                    orig_texture_patches, texture_patches->at(projected_edge_info.texture_patch_id));
            }
        }
    }

    for (std::size_t i = 0; i < num_vertices; ++i) {
        std::vector<VertexProjectionInfo> const & projection_infos = vertex_projection_infos[i];
        if (projection_infos.size() <= 1) continue;

        math::Accum<math::Vec3f> color_accum(math::Vec3f(0.0f));
        for (std::size_t j = 0; j < projection_infos.size(); ++j) {
            VertexProjectionInfo const & projection_info = projection_infos[j];
            TexturePatch::ConstPtr original_texture_patch = orig_texture_patches[projection_info.texture_patch_id];
            if (original_texture_patch->get_label() == 0) continue;
            math::Vec3f color = original_texture_patch->get_pixel_value(projection_info.projection);
            color_accum.add(color, 1.0f);
        }

        // mean color of vertices across multiple views
        math::Vec3f mean_color = color_accum.normalized();

        for (std::size_t j = 0; j < projection_infos.size(); ++j){
            VertexProjectionInfo const & projection_info = projection_infos[j];
            math::Vec2i pixel(projection_info.projection +  math::Vec2f(0.5f, 0.5f));  // while 0.5 ??????????????
            TexturePatch::Ptr texture_patch = texture_patches->at(projection_info.texture_patch_id);
            texture_patch->set_pixel_value(pixel, mean_color);
        }
    }

    // ====================================Blending====================================================//
    ProgressCounter texture_patch_counter("\tBlending texture patches", texture_patches->size());
    #pragma omp parallel for schedule(dynamic)
    for (std::size_t i = 0; i < texture_patches->size(); ++i) {
        TexturePatch::Ptr texture_patch = texture_patches->at(i);
        texture_patch_counter.progress<SIMPLE>();

        /* Only alter a small strip of texture patches originating from input images. */
        if (texture_patch->get_label() != 0) {
            texture_patch->prepare_blending_mask(STRIP_SIZE);
        }

        if(texture_patch->get_faces().size()<10000)continue;
        char image_name_before[255];
        char image_name_after[255];
        char validity_mask_name[255];
        char blending_mask_name[255];
        sprintf(image_name_before,"examples/task7/texture_patches_poisson_blending/texture_patch_before%d.jpg", i);
        sprintf(image_name_after,"examples/task7/texture_patches_poisson_blending/texture_patch_after%d.jpg", i);
        sprintf(blending_mask_name,"examples/task7/texture_patches_poisson_blending/blending_mask%d.jpg", i);
        sprintf(validity_mask_name,"examples/task7/texture_patches_poisson_blending/validity_mask%d.jpg", i);
        core::FloatImage::Ptr image = texture_patch->get_image()->duplicate();
        core::ByteImage::Ptr validity_mask = texture_patch->get_validity_mask()->duplicate();
        core::ByteImage::Ptr blending_mask = texture_patch->get_blending_mask()->duplicate();

        core::ByteImage::Ptr blending_mask_color = core::ByteImage::create(blending_mask->width()
                , blending_mask->height(), 3);
        //blending_mask_color->fill_color({0});

        for (int y = 0; y < blending_mask->height(); ++y) {
            for (int x = 0; x < blending_mask->width(); ++x) {

                blending_mask_color->at(x, y, 0) = 0;
                blending_mask_color->at(x, y, 1) = 255;
                blending_mask_color->at(x, y, 2) = 0;

                if (blending_mask->at(x, y, 0) == 128) {
                    blending_mask_color->at(x, y, 0) = 255;
                    blending_mask_color->at(x, y, 1) = 0;
                    blending_mask_color->at(x, y, 2) = 0;
                }
                if (blending_mask->at(x, y, 0) == 126) {
                    blending_mask_color->at(x, y, 0) = 0;
                    blending_mask_color->at(x, y, 1) = 0;
                    blending_mask_color->at(x, y, 2) = 255;
                }

                if (blending_mask->at(x, y, 0) == 255) {
                    blending_mask_color->at(x, y, 0) = 255;
                    blending_mask_color->at(x, y, 1) = 255;
                    blending_mask_color->at(x, y, 2) = 255;
                }
            }
        }
        core::image::save_file(core::image::float_to_byte_image(image), image_name_before);
        core::image::save_file(validity_mask, validity_mask_name);
        core::image::save_file(blending_mask_color, blending_mask_name);

        // poisson blending
        texture_patch->blend(orig_texture_patches[i]->get_image());


        image = texture_patch->get_image()->duplicate();
        core::image::save_file(core::image::float_to_byte_image(image), image_name_after);


        texture_patch->release_blending_mask();
        texture_patch_counter.inc();
    }
}

TEX_NAMESPACE_END
