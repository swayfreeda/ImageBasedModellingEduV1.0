/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>
#include <fstream>

#include "util/aligned_memory.h"
#include "core/image.h"
#include "core/image_io.h"
#include "core/image_tools.h"
#include "features/sift.h"
#include "features/matching.h"
#include "sfm/correspondence.h"
#include "sfm/ransac_fundamental.h"
#include "sfm/ransac_homography.h"
#include "sfm/feature_set.h"
#include "visualizer.h"

#define MAX_PIXELS 1000000

void
normalize_correspondences (sfm::Correspondences2D2D * corr, int width1, int height1,
    int width2, int height2)
{
    float const img1_width = static_cast<float>(width1);
    float const img1_height = static_cast<float>(height1);
    float const img1_maxdim = static_cast<float>(std::max(width1, height1));
    float const img2_width = static_cast<float>(width2);
    float const img2_height = static_cast<float>(height2);
    float const img2_maxdim = static_cast<float>(std::max(width2, height2));
    for (std::size_t i = 0; i < corr->size(); ++i)
    {
        sfm::Correspondence2D2D & c = corr->at(i);
        c.p1[0] = (c.p1[0] + 0.5f - img1_width / 2.0f) / img1_maxdim;
        c.p1[1] = (c.p1[1] + 0.5f - img1_height / 2.0f) / img1_maxdim;
        c.p2[0] = (c.p2[0] + 0.5f - img2_width / 2.0f) / img2_maxdim;
        c.p2[1] = (c.p2[1] + 0.5f - img2_height / 2.0f) / img2_maxdim;
    }
}

void
denormalize_correspondences (sfm::Correspondences2D2D * corr, int width1, int height1,
    int width2, int height2)
{
    float const img1_width = static_cast<float>(width1);
    float const img1_height = static_cast<float>(height1);
    float const img1_maxdim = static_cast<float>(std::max(width1, height1));
    float const img2_width = static_cast<float>(width2);
    float const img2_height = static_cast<float>(height2);
    float const img2_maxdim = static_cast<float>(std::max(width2, height2));
    for (std::size_t i = 0; i < corr->size(); ++i)
    {
        sfm::Correspondence2D2D & c = corr->at(i);
        c.p1[0] = c.p1[0] * img1_maxdim + img1_width / 2.0f - 0.5f;
        c.p1[1] = c.p1[1] * img1_maxdim + img1_height / 2.0f - 0.5f;
        c.p2[0] = c.p2[0] * img2_maxdim + img2_width / 2.0f - 0.5f;
        c.p2[1] = c.p2[1] * img2_maxdim + img2_height / 2.0f - 0.5f;
    }
}

void convert_sift_discriptors(features::Sift::Descriptors const&sift_descrs,
                            util::AlignedMemory<math::Vec128f, 16> *aligned_descr)
{
    aligned_descr->resize(sift_descrs.size());
    float * data_ptr=aligned_descr->data()->begin();
    for(int i=0; i<sift_descrs.size(); ++i, data_ptr+=128)
    {
        features::Sift::Descriptor const& descr = sift_descrs[i];
        std::copy(descr.data.begin(), descr.data.end(), data_ptr);
    }

}

features::Matching::Result feature_matching(features::Sift::Descriptors const&sift_discrs1
                    , features::Sift::Descriptors const &sift_discrs2)
{
    // 进行数据转换
    util::AlignedMemory<math::Vec128f, 16> aligned_descrs1, aligned_descrs2;
    convert_sift_discriptors(sift_discrs1, &aligned_descrs1);
    convert_sift_discriptors(sift_discrs2, &aligned_descrs2);


    features::Matching::Options matching_opts;
    matching_opts.descriptor_length = 128;
    matching_opts.distance_threshold = 1.0f;
    matching_opts.lowe_ratio_threshold = 0.8f;

    features::Matching::Result matching_result;
    features::Matching::twoway_match(matching_opts, aligned_descrs1.data()->begin(), sift_discrs1.size()
                                                    ,aligned_descrs2.data()->begin(), sift_discrs2.size(),&matching_result);

    features::Matching::remove_inconsistent_matches(&matching_result);

    return matching_result;
}



int
main (int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "Syntax: " << argv[0] << " <img1> <img2>" << std::endl;
        return 1;
    }

    /* 加载图像. */
    std::string fname1(argv[1]);
    std::string fname2(argv[2]);
    std::cout << "Loading " << fname1 << "..." << std::endl;
    core::ByteImage::Ptr img1 = core::image::load_file(fname1);
    std::cout << "Loading " << fname2 << "..." << std::endl;
    core::ByteImage::Ptr img2 = core::image::load_file(fname2);
    // 控制图像尺寸
    while(img1->get_pixel_amount() > MAX_PIXELS){
        img1=core::image::rescale_half_size<uint8_t>(img1);
    }
    while(img2->get_pixel_amount() > MAX_PIXELS){
        img2=core::image::rescale_half_size<uint8_t>(img2);
    }


    /*计算sift特征点 */
    sfm::FeatureSet::Options feature_set_opts;
    feature_set_opts.feature_types = sfm::FeatureSet::FEATURE_SIFT;
    feature_set_opts.sift_opts.verbose_output = true;

    sfm::FeatureSet feat1(feature_set_opts);
    feat1.compute_features(img1);

    sfm::FeatureSet feat2(feature_set_opts);
    feat2.compute_features(img2);

    std::cout << "Image 1 (" << img1->width() << "x" << img1->height() << ") "
              << feat1.sift_descriptors.size() << " descriptors." << std::endl;
    std::cout << "Image 2 (" << img2->width() << "x" << img2->height() << ") "
              << feat2.sift_descriptors.size() << " descriptors." << std::endl;


    /*特征匹配*/
    features::Matching::Result mathing_result = feature_matching(feat1.sift_descriptors, feat2.sift_descriptors);
    int n_consitent_matches = features::Matching::count_consistent_matches(mathing_result);
    std::cout << "Consistent Sift Matches: "
              << n_consitent_matches
              << std::endl;


    sfm::Correspondences2D2D corr_all;
    std::vector<int> const & m12 = mathing_result.matches_1_2;
    for(int i=0; i<m12.size(); i++)
    {
        if(m12[i]<0)continue;

        sfm::Correspondence2D2D c2d2d;

        c2d2d.p1[0] = feat1.positions[i][0];
        c2d2d.p1[1] = feat1.positions[i][1];
        c2d2d.p2[0] = feat2.positions[m12[i]][0];
        c2d2d.p2[1] = feat2.positions[m12[i]][1];

        corr_all.push_back(c2d2d);
    }



    normalize_correspondences(&corr_all, img1->width(), img1->height(),
                                img2->width(), img1->height());


    std::ofstream out("./examples/task2/correspondences.txt");
    assert(out.is_open());
    out<<corr_all.size()<<std::endl;
    for(int i=0; i< corr_all.size(); i++){
        out<<corr_all[i].p1[0]<<" "<<corr_all[i].p1[1]<<" "<<corr_all[i].p2[0]<<" "<<corr_all[i].p2[1]<<std::endl;
    }

    /* RANSAC 估计本征矩阵, 并对特征匹配对进行筛选*/
    sfm::RansacFundamental::Options ransac_fundamental_opts;
    ransac_fundamental_opts.max_iterations =1000;
    ransac_fundamental_opts.verbose_output = true;
    sfm::RansacFundamental ransac_fundamental(ransac_fundamental_opts);
    sfm::RansacFundamental::Result ransac_fundamental_result;
    ransac_fundamental.estimate(corr_all, &ransac_fundamental_result);
    // 根据估计的Fundamental矩阵对特征匹配对进行筛选
    sfm::Correspondences2D2D corr_f;
    for(int i=0; i<ransac_fundamental_result.inliers.size(); ++i)
    {
        int inlier_id = ransac_fundamental_result.inliers[i];
        corr_f.push_back(corr_all[inlier_id]);
    }

    std::cout<<"total number is "<< corr_all.size()<<std::endl;
    std::cout<<"inlier number is "<< corr_f.size()<<std::endl;
#if 0
    /* RANSAC 估计单应矩阵, 并对特征匹配对进行筛选*/
    sfm::RansacHomography::Options ransac_homography_opts;
    ransac_homography_opts.verbose_output=true;
    ransac_homography_opts.max_iterations=1000;
    sfm::RansacHomography ransac_homography(ransac_homography_opts);
    sfm::RansacHomography::Result ransac_homograph_result;
    ransac_homography.estimate(corr_all, &ransac_homograph_result);
    // 根据估计的homography对匹配对进行筛选
    sfm::Correspondences2D2D corr_h;
    for(int i=0; i<ransac_homograph_result.inliers.size(); ++i)
    {
        int inlier_id = ransac_homograph_result.inliers[i];
        corr_h.push_back(corr_all[inlier_id]);
    }


    /* 匹配特征可视化. */
    denormalize_correspondences(&corr_all, img1->width(), img1->height(),
                                img2->width(), img2->height());
    denormalize_correspondences(&corr_f, img1->width(), img1->height(),
                                img2->width(), img2->height());
    denormalize_correspondences(&corr_h, img1->width(), img1->height(),
                                img2->width(), img2->height());
    core::ByteImage::Ptr vis_img;
    vis_img = features::Visualizer::draw_matches(img1, img2, corr_all);
    core::image::save_png_file(vis_img, "./tmp/matches_unfiltered.png");
    vis_img = features::Visualizer::draw_matches(img1, img2, corr_f);
    core::image::save_png_file(vis_img, "./tmp/matches_fundamental.png");
    vis_img = features::Visualizer::draw_matches(img1, img2, corr_h);
    core::image::save_png_file(vis_img, "./tmp/matches_homography.png");
#endif
    return 0;
}
