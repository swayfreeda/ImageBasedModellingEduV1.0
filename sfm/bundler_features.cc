/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "util/timer.h"
#include "core/image.h"
#include "core/image_exif.h"
#include "core/image_tools.h"
#include "sfm/bundler_common.h"
#include "sfm/extract_focal_length.h"
#include "sfm/bundler_features.h"

SFM_NAMESPACE_BEGIN
SFM_BUNDLER_NAMESPACE_BEGIN

void
Features::compute (core::Scene::Ptr scene, ViewportList* viewports)
{

    if (scene == nullptr)
        throw std::invalid_argument("Null scene given");
    if (viewports == nullptr)
        throw std::invalid_argument("No viewports given");

    core::Scene::ViewList const& views = scene->get_views();

    /* Initialize viewports. */
    viewports->clear();
    viewports->resize(views.size());

    std::size_t num_views = viewports->size();
    std::size_t num_done = 0;
    std::size_t total_features = 0;

    /* Iterate the scene and compute features. */
    for (std::size_t i = 0; i < views.size(); ++i)
    {
        {
            num_done += 1;
            float percent = (num_done * 1000 / num_views) / 10.0f;
            std::cout << "\rDetecting features, view " << num_done << " of "
                << num_views << " (" << percent << "%)..." << std::flush;
        }


        // 获取图像
        if (views[i] == nullptr)
            continue;

        core::View::Ptr view = views[i];
        core::ByteImage::Ptr image = view->get_byte_image
            (this->opts.image_embedding);
        if (image == nullptr)
            continue;

        // 对图像进行降采样，使得其满足尺寸限制
        /* Rescale image until maximum image size is met. */
        util::WallTimer timer;
        while (this->opts.max_image_size > 0
            && image->width() * image->height() > this->opts.max_image_size)
            image = core::image::rescale_half_size<uint8_t>(image);


        // 计算每个视角的特征点
        Viewport* viewport = &viewports->at(i);
        viewport->features.set_options(this->opts.feature_options);
        viewport->features.compute_features(image);
        // 对图像特征点的位置进行初始化
        viewport->features.normalize_feature_positions();

        {
            std::size_t const num_feats = viewport->features.positions.size();
            std::cout << "\rView ID "
                << util::string::get_filled(view->get_id(), 4, '0') << " ("
                << image->width() << "x" << image->height() << "), "
                << util::string::get_filled(num_feats, 5, ' ') << " features"
                << ", took " << timer.get_elapsed() << " ms." << std::endl;
            total_features += viewport->features.positions.size();
        }

        /* Clean up unused embeddings. */
        image.reset();
        view->cache_cleanup();
    }

    std::cout << "\rComputed " << total_features << " features "
        << "for " << num_views << " views (average "
        << (total_features / num_views) << ")." << std::endl;

}

SFM_BUNDLER_NAMESPACE_END
SFM_NAMESPACE_END
