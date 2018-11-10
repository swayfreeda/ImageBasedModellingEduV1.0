/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <fstream>
#include <vector>

#include <util/timer.h>
#include <util/system.h>
#include <util/file_system.h>
#include <core/mesh_io_ply.h>
#include <core/image_tools.h>
#include <core/image_drawing.h>
#include "core/image_io.h"

#include "texturing/util.h"
#include "texturing/timer.h"
#include "texturing/debug.h"
#include "texturing/texturing.h"
#include "texturing/progress_counter.h"

#include "arguments.h"

int main(int argc, char **argv) {
#ifdef RESEARCH
    std::cout << "******************************************************************************" << std::endl
              << " Due to use of the -DRESEARCH=ON compile option, this program is licensed "     << std::endl
              << " for research purposes only. Please pay special attention to the gco license."  << std::endl
              << "******************************************************************************" << std::endl;
#endif

    util::system::register_segfault_handler();
    Timer timer;
    timer.measure("Start");
    util::WallTimer wtimer;
    Arguments conf;
    try {
        conf = parse_args(argc, argv);
    } catch (std::invalid_argument & ia) {
        std::cerr << ia.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (!util::fs::dir_exists(util::fs::dirname(conf.out_prefix).c_str())) {
        std::cerr << "Destination directory does not exist!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    //==================================Load Mesh ===================================//
    std::cout << "Load and prepare mesh: " << std::endl;
    core::TriangleMesh::Ptr mesh;
    try {
        mesh = core::geom::load_ply_mesh(conf.in_mesh);
    } catch (std::exception& e) {
        std::cerr << "\tCould not load mesh: "<< e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }


    //=================================Prepare Mesh=================================//
    core::VertexInfoList::Ptr vertex_infos = core::VertexInfoList::create(mesh);
    tex::prepare_mesh(vertex_infos, mesh);


    //=================================Geneatring texture views=====================//
    std::size_t const num_faces = mesh->get_faces().size() / 3;
    std::cout << "Generating texture views: " << std::endl;
    tex::TextureViews texture_views;
    tex::generate_texture_views(conf.in_scene, &texture_views);
    write_string_to_file(conf.out_prefix + ".conf", conf.to_string());
    timer.measure("Loading");


    //===============================Building adjacency graph=======================//
    std::cout << "Building adjacency graph: " << std::endl; // each facet is corresponding a facet
    tex::Graph graph(num_faces);
    tex::build_adjacency_graph(mesh, vertex_infos, &graph);
    wtimer.reset();


    //===============================View Selection ================================//
    // if labeling file does not exist, compute a view label for each facet via MRF
    if (conf.labeling_file.empty()) {
        std::cout << "View selection:" << std::endl;
        std::size_t const num_faces = mesh->get_faces().size() / 3;
        tex::DataCosts data_costs(num_faces, texture_views.size());

        // if data cost file does not exist, compute the data cost
        if (conf.data_cost_file.empty()) {

            /******************* compute data cost  ***********************************/
            tex::calculate_data_costs(mesh, &texture_views, conf.settings, &data_costs);
            /**************************************************************************/

            if (conf.write_intermediate_results) {
                std::cout << "\tWriting data cost file... " << std::flush;
                    ST::save_to_file(data_costs, conf.out_prefix + "_data_costs.spt");
                std::cout << "done." << std::endl;
            }
        } else {  // if the data cost file exists, just load it from the file
            std::cout << "\tLoading data cost file... " << std::flush;
            try {
                ST::load_from_file(conf.data_cost_file, &data_costs);
            } catch (util::FileException e) {
                std::cout << "failed!" << std::endl;
                std::cerr << e.what() << std::endl;
                std::exit(EXIT_FAILURE);
            }
            std::cout << "done." << std::endl;
        }
        timer.measure("Calculating data costs");


        // MRF Optimization for view selection
        tex::view_selection(data_costs, &graph, conf.settings);
        timer.measure("Running MRF optimization");

        /* Write labeling to file. */
        if (conf.write_intermediate_results) {
            std::vector<std::size_t> labeling(graph.num_nodes());
            for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
                labeling[i] = graph.get_label(i);
            }
            vector_to_file(conf.out_prefix + "_labeling.vec", labeling);
        }
    } else {  // if labeling file has existed, just read it from the file
        std::cout << "Loading labeling from file... " << std::flush;

        /* Load labeling from file. */
        std::vector<std::size_t> labeling = vector_from_file<std::size_t>(conf.labeling_file);
        if (labeling.size() != graph.num_nodes()) {
            std::cerr << "Wrong labeling file for this mesh/scene combination... aborting!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        /* Transfer labeling to graph. */
        for (std::size_t i = 0; i < labeling.size(); ++i) {
            const std::size_t label = labeling[i];
            if (label > texture_views.size()){
                std::cerr << "Wrong labeling file for this mesh/scene combination... aborting!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            graph.set_label(i, label);
        }

        std::cout << "done." << std::endl;
    }
    std::cout << "\tTook: " << wtimer.get_elapsed_sec() << "s" << std::endl;
    // todo rendering the mesh with different colors
    {
        std::vector<std::size_t> labeling(graph.num_nodes());
        for (std::size_t i = 0; i < graph.num_nodes(); ++i) {
            labeling[i] = graph.get_label(i);
        }
        std::vector<std::size_t>::iterator max_iter = std::max_element(labeling.begin(), labeling.end());
        int n_labels = *max_iter + 1;
        std::vector<math::Vec3f> colors(n_labels);
        colors[0][0] = 0;
        colors[0][1] = 0;
        colors[0][2] = 0;

        for(int i=1; i< colors.size(); i++){
            colors[i][0] = rand()&255;
            colors[i][1] = rand()&255;
            colors[i][2] = rand()&255;
        }

        std::ofstream out("./examples/task7/view_selection_result.ply");
        assert(out.is_open());
        out<<"ply"<<std::endl;
        out<<"format ascii 1.0"<<std::endl;
        out<<"element vertex "<<mesh->get_vertices().size()<<std::endl;
        out<<"property float x"<<std::endl;
        out<<"property float y"<<std::endl;
        out<<"property float z"<<std::endl;
        out<<"element face "<<mesh->get_faces().size()/3<<std::endl;
        out<<"property list uchar int vertex_indices"<<std::endl;
        out<<"property uchar red"<<std::endl;
        out<<"property uchar green"<<std::endl;
        out<<"property uchar blue"<<std::endl;
        out<<"end_header"<<std::endl;

        // Face List
        core::TriangleMesh::FaceList const & mesh_faces = mesh->get_faces();
        // Vertices
        core::TriangleMesh::VertexList const & vertices = mesh->get_vertices();
        for(int i=0; i< vertices.size(); i++){
            out<<vertices[i][0]<<" "<<vertices[i][1]<<" "<<vertices[i][2]<<std::endl;
        }

        std::cout<<"labeling size: "<<labeling.size()<<" faces size: "<<num_faces<<std::endl;
        for(int i=0; i< labeling.size(); i++){
            int label = labeling[i];
            assert(label>=0 && label <labeling.size());
            int v0 = mesh_faces[3*i + 0];
            int v1 = mesh_faces[3*i + 1];
            int v2 = mesh_faces[3*i + 2];
            int r = colors[label][0];
            int g = colors[label][1];
            int b = colors[label][2];
            out<<"3 "<<v0<<" "<<v1<<" "<<v2<<" "<<r<<" "<<g<<" "<<b<<std::endl;
        }
        out.close();
    }

    //=============================================texture_atlases====================================================//
    tex::TextureAtlases texture_atlases;
    {
        /* Create texture patches and adjust them. */
        tex::TexturePatches texture_patches;
        tex::VertexProjectionInfos vertex_projection_infos;
        std::cout << "Generating texture patches:" << std::endl;
        tex::generate_texture_patches(graph,
                                      mesh,
                                      vertex_infos,
                                      &texture_views,
                                      &vertex_projection_infos,
                                      &texture_patches);

        // todo save texturePatches and the projection of facets
        {
            for(int i=0; i< texture_patches.size(); i++)
            {
                if(texture_patches[i]->get_faces().size()<10000)continue;

                char image_name[255];
                char validity_mask_name[255];
                char blending_mask_name[255];

                sprintf(image_name,"examples/task7/texture_patches_init/texture_patch%d.jpg", i);
                sprintf(blending_mask_name,"examples/task7/texture_patches_init/blending_mask%d.jpg", i);
                sprintf(validity_mask_name,"examples/task7/texture_patches_init/validity_mask%d.jpg", i);

                core::FloatImage::Ptr image = texture_patches[i]->get_image()->duplicate();
                core::ByteImage::Ptr validity_mask = texture_patches[i]->get_validity_mask()->duplicate();
                core::ByteImage::Ptr blending_mask = texture_patches[i]->get_blending_mask()->duplicate();

                float color[3]={255, 0,0};
                std::vector<math::Vec2f> texcoords = texture_patches[i]->get_texcoords();
                for(int i=0; i< texcoords.size(); i+=3){
                    math::Vec2f v0 = texcoords[i+0];
                    math::Vec2f v1 = texcoords[i+1];
                    math::Vec2f v2 = texcoords[i+2];

                    core::image::draw_line<float>(*image, int(v0[0]), int(v0[1]), int(v1[0]), int(v1[1]), color);
                    core::image::draw_line<float>(*image, int(v1[0]), int(v1[1]), int(v2[0]), int(v2[1]), color);
                    core::image::draw_line<float>(*image, int(v0[0]), int(v0[1]), int(v2[0]), int(v2[1]), color);

                }

                core::image::save_file(core::image::float_to_byte_image(image), image_name);
                core::image::save_file(validity_mask, validity_mask_name);
                core::image::save_file(blending_mask, blending_mask_name);

            }
        }

        // Global seam leveling
        if (conf.settings.global_seam_leveling) {
            std::cout << "Running global seam leveling:" << std::endl;
            tex::global_seam_leveling(graph,
                                      mesh,
                                      vertex_infos,
                                      vertex_projection_infos,
                                      &texture_patches);
            timer.measure("Running global seam leveling");

            {
                for(int i=0; i< texture_patches.size(); i++)
                {
                    if(texture_patches[i]->get_faces().size()<10000)continue;

                    char image_name[255];
                    char validity_mask_name[255];
                    char blending_mask_name[255];

                    sprintf(image_name,"examples/task7/texture_pathes_color_adjustment/texture_patch%d.jpg", i);
                    sprintf(blending_mask_name,"examples/task7/texture_pathes_color_adjustment/blending_mask%d.jpg", i);
                    sprintf(validity_mask_name,"examples/task7/texture_pathes_color_adjustment/validity_mask%d.jpg", i);

                    core::FloatImage::Ptr image = texture_patches[i]->get_image()->duplicate();
                    core::ByteImage::Ptr validity_mask = texture_patches[i]->get_validity_mask()->duplicate();
                    core::ByteImage::Ptr blending_mask = texture_patches[i]->get_blending_mask()->duplicate();

                    float color[3]={255, 0,0};
                    std::vector<math::Vec2f> texcoords = texture_patches[i]->get_texcoords();
                    for(int i=0; i< texcoords.size(); i+=3){
                        math::Vec2f v0 = texcoords[i+0];
                        math::Vec2f v1 = texcoords[i+1];
                        math::Vec2f v2 = texcoords[i+2];

                        core::image::draw_line<float>(*image, int(v0[0]), int(v0[1]), int(v1[0]), int(v1[1]), color);
                        core::image::draw_line<float>(*image, int(v1[0]), int(v1[1]), int(v2[0]), int(v2[1]), color);
                        core::image::draw_line<float>(*image, int(v0[0]), int(v0[1]), int(v2[0]), int(v2[1]), color);

                    }

                    core::image::save_file(core::image::float_to_byte_image(image), image_name);
                    core::image::save_file(validity_mask, validity_mask_name);
                    core::image::save_file(blending_mask, blending_mask_name);

                }
            }

        } else {
            ProgressCounter texture_patch_counter("Calculating validity masks for texture patches", texture_patches.size());
            #pragma omp parallel for schedule(dynamic)
            for (std::size_t i = 0; i < texture_patches.size(); ++i) {
                texture_patch_counter.progress<SIMPLE>();
                TexturePatch::Ptr texture_patch = texture_patches[i];
                std::vector<math::Vec3f> patch_adjust_values(texture_patch->get_faces().size() * 3, math::Vec3f(0.0f));
                texture_patch->adjust_colors(patch_adjust_values);
                texture_patch_counter.inc();
            }
            timer.measure("Calculating texture patch validity masks");
        }


        //======================================local seam leveling===========================//
        // Local seam leveling (Poisson Editing???)
        if (conf.settings.local_seam_leveling) {
            std::cout << "Running local seam leveling:" << std::endl;
            tex::local_seam_leveling(graph, mesh, vertex_projection_infos, &texture_patches);
        }
        timer.measure("Running local seam leveling");


        //====================================Generating textgure atlases====================//
        /* Generate texture atlases. */
        std::cout << "Generating texture atlases:" << std::endl;
        tex::generate_texture_atlases(&texture_patches, &texture_atlases);
    }


       //=================================Write Obj model=====================================//
    /* Create and write out obj model. */
    {
        std::cout << "Building objmodel:" << std::endl;
        tex::Model model;
        tex::build_model(mesh, texture_atlases, &model);
        timer.measure("Building OBJ model");

        std::cout << "\tSaving model... " << std::flush;
        tex::Model::save(model, conf.out_prefix);
        std::cout << "done." << std::endl;
        timer.measure("Saving");
    }

    std::cout << "Whole texturing procedure took: " << wtimer.get_elapsed_sec() << "s" << std::endl;
    timer.measure("Total");
    if (conf.write_timings) {
        timer.write_to_file(conf.out_prefix + "_timings.csv");
    }

    if (conf.write_view_selection_model) {
        texture_atlases.clear();
        std::cout << "Generating debug texture patches:" << std::endl;
        {
            tex::TexturePatches texture_patches;
            generate_debug_embeddings(&texture_views);
            tex::VertexProjectionInfos vertex_projection_infos; // Will only be written
            tex::generate_texture_patches(graph, mesh, vertex_infos, &texture_views, &vertex_projection_infos, &texture_patches);
            tex::generate_texture_atlases(&texture_patches, &texture_atlases);
        }

        std::cout << "Building debug objmodel:" << std::endl;
        {
            tex::Model model;
            tex::build_model(mesh, texture_atlases, &model);
            std::cout << "\tSaving model... " << std::flush;
            tex::Model::save(model, conf.out_prefix + "_view_selection");
            std::cout << "done." << std::endl;
        }
    }
}
