/*
 * Copyright (C) 2015, Ronny Klowsky, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <set>
#include <ctime>

#include <core/scene.h>
#include "mvs/defines.h"
#include "math/vector.h"
#include "math/functions.h"
#include "math/octree_tools.h"
#include "core/image.h"
#include "core/image_tools.h"
#include "util/file_system.h"
#include "util/strings.h"
#include "mvs/settings.h"
#include "mvs/dmrecon.h"
#include "mvs/global_view_selection.h"
#include "mvs/progress.h"
#include "mvs/single_view.h"

MVS_NAMESPACE_BEGIN

DMRecon::DMRecon(core::Scene::Ptr _scene, Settings const& _settings)
    : scene(_scene)
    , settings(_settings)
{
    core::Scene::ViewList const& mve_views(scene->get_views());

    /* Check if master image exists */
    if (settings.refViewNr >= mve_views.size())
        throw std::invalid_argument("Master view index out of bounds");

    /* Check for meaningful scale factor */
    if (settings.scale < 0.f)
        throw std::invalid_argument("Invalid scale factor");

    /* Check if image embedding is set. */
    if (settings.imageEmbedding.empty())
        throw std::invalid_argument("Invalid image embedding");
    /* Fetch bundle file. */
    try {
        this->bundle = this->scene->get_bundle();
    }
    catch (std::exception& e) {
        throw std::runtime_error(std::string("Error reading bundle file: ")
              + e.what());
    }

    //为场景中的每一幅图像创建一个Single view
    /* Create list of SingleView pointers from MVE views. */
    views.resize(mve_views.size());
    for (std::size_t i = 0; i < mve_views.size(); ++i) {
        if (mve_views[i] == nullptr || !mve_views[i]->is_camera_valid() ||
            !mve_views[i]->has_image(this->settings.imageEmbedding,
            core::IMAGE_TYPE_UINT8))
            continue;
        views[i] = mvs::SingleView::create(scene, mve_views[i],
            this->settings.imageEmbedding);
    }

    // 提取refrence view 并创建图像金字塔
    SingleView::Ptr refV = views[settings.refViewNr];
    if (refV == nullptr)
        throw std::invalid_argument("Invalid master view");

    /* Prepare reconstruction */
    // 加载图像，并且构建图像金字塔
    refV->loadColorImage(this->settings.scale);
    // 创建相关的depth image, normal image, dz image, conf image深度图像，法向量图像等等
    refV->prepareMasterView(settings.scale);
    core::ByteImage::ConstPtr scaled_img = refV->getScaledImg();
    // 需要重建的尺度的图像的尺寸
    this->width = scaled_img->width();
    this->height = scaled_img->height();

    if (!settings.quiet)
        std::cout << "scaled image size: " << this->width << " x "
                  << this->height << std::endl;
}

void
DMRecon::start()
{
    try
    {
        progress.start_time = std::time(nullptr);

        // 对稀疏特征进行重建
        analyzeFeatures();

        // 全局视角选择
        globalViewSelection();

        // 处理特征，对当前的三维点投影到图像上进行深度值估计
        // 并且将重建的特征点添加到队列中，作为种子点
        processFeatures();

        // 处理队列
        processQueue();


        // 保存图像
        if (progress.cancelled)
        {
            progress.status = RECON_CANCELLED;
            return;
        }

        progress.status = RECON_SAVING;
        SingleView::Ptr refV(views[settings.refViewNr]);
        if (settings.writePlyFile)
        {
            if (!settings.quiet)
                std::cout << "Saving ply file as "
                    << settings.plyPath << "/"
                    << refV->createFileName(settings.scale)
                    << ".ply" << std::endl;
            refV->saveReconAsPly(settings.plyPath, settings.scale);
        }

        // Save images to view
        core::View::Ptr view = refV->getMVEView();

        std::string name("depth-L");
        name += util::string::get(settings.scale);
        view->set_image(refV->depthImg, name);

        if (settings.keepDzMap)
        {
            name = "dz-L";
            name += util::string::get(settings.scale);
            view->set_image(refV->dzImg, name);
        }

        if (settings.keepConfidenceMap)
        {
            name = "conf-L";
            name += util::string::get(settings.scale);
            view->set_image(refV->confImg, name);
        }

        if (settings.scale != 0)
        {
            name = "undist-L";
            name += util::string::get(settings.scale);
            view->set_image(refV->getScaledImg()->duplicate(), name);
        }

        progress.status = RECON_IDLE;

        /* Output percentage of filled pixels */
        {
            int nrPix = this->width * this->height;
            float percent = (float) progress.filled / (float) nrPix;
            if (!settings.quiet)
                std::cout << "Filled " << progress.filled << " pixels, i.e. "
                          << util::string::get_fixed(percent * 100.f, 1)
                          << " %." << std::endl;
        }

        /* Output required time to process the image */
        size_t mvs_time = std::time(nullptr) - progress.start_time;
        if (!settings.quiet)
            std::cout << "MVS took " << mvs_time << " seconds." << std::endl;
    }
    catch (util::Exception e)
    {
        if (!settings.quiet)
            std::cout << "Reconstruction failed: " << e << std::endl;

        progress.status = RECON_CANCELLED;
        return;
    }
}

/*
 * Attach features that are visible in the reference view (according to
 * the bundle) to all other views if inside the frustum.
 */
void
DMRecon::analyzeFeatures()
{
    progress.status = RECON_FEATURES;

    // 获取参考视角
    SingleView::ConstPtr refV = views[settings.refViewNr];
    // 读取所有重建的3D稀疏点
    core::Bundle::Features const& features = bundle->get_features();

    // 对每一个特征点
    for (std::size_t i = 0; i < features.size() && !progress.cancelled; ++i)
    {
        // 三维点在该参考视角中可见
        if (!features[i].contains_view_id(settings.refViewNr))
            continue;

        // 将特征点的三维坐标投影到原始尺度的图像中，判断：1) 相机坐标系中的Z坐标是否为正； 2) 投影到图像上的坐标是否超出了
        // 图像的范围
        math::Vec3f featurePos(features[i].pos);
        if (!refV->pointInFrustum(featurePos))
            continue;

        // 重建的特征点需要在设定的空间范围内( aabbMin 和 aabbMax确定的包围盒内)
        if (!math::geom::point_box_overlap(featurePos,
            this->settings.aabbMin, this->settings.aabbMax))
            continue;

        // 如果该三维点是在参考视角中可见的，则将其添加到所有可见的视角中--每个视角包含哪些可见的特征点
        for (std::size_t j = 0; j < features[i].refs.size(); ++j)
        {
            int view_id = features[i].refs[j].view_id;
            if (view_id < 0 || view_id >= static_cast<int>(views.size())
                || views[view_id] == nullptr)
                continue;

            // 判断视角是否在其它视角中可见
            if (views[view_id]->pointInFrustum(featurePos))
                views[view_id]->addFeature(i);
        }
    }
}

void
DMRecon::globalViewSelection()
{
    progress.status = RECON_GLOBALVS;
    if (progress.cancelled)
        return;

    //执行全局的视角选择
    /* Perform global view selection. */
    GlobalViewSelection globalVS(views, bundle->get_features(), settings);
    globalVS.performVS();
    neighViews = globalVS.getSelectedIDs();

    // 全局的视角选择失败
    if (neighViews.empty())
        throw std::runtime_error("Global View Selection failed");

    /* Print result of global view selection. */
    if (!settings.quiet)
    {
        std::cout << "Global View Selection:";
        for (IndexSet::const_iterator iter = neighViews.begin();
            iter != neighViews.end(); ++iter)
            std::cout << " " << *iter;
        std::cout << std::endl;
    }

    /* Load selected images. */
    if (!settings.quiet)
        std::cout << "Loading color images..." << std::endl;

    // 对全局选择的邻域视角加载图像
    for (IndexSet::const_iterator iter = neighViews.begin();
        iter != neighViews.end() && !progress.cancelled; ++iter)
        views[*iter]->loadColorImage(0);
}

void
DMRecon::processFeatures()
{
    progress.status = RECON_FEATURES;
    if (progress.cancelled)
        return;

    // 参考视角
    SingleView::Ptr refV = views[settings.refViewNr];

    // 所有的三维特征点
    core::Bundle::Features const& features = bundle->get_features();

    if (!settings.quiet)
        std::cout << "Processing " << features.size()
            << " features..." << std::endl;

    std::size_t success = 0;
    std::size_t processed = 0;
    for (std::size_t i = 0; i < features.size() && !progress.cancelled; ++i) {
        /*
         * Use feature if visible in reference view or
         * at least one neighboring view.
         */
        bool useFeature = false;

        // 特征在参考视角中可见
        if (features[i].contains_view_id(settings.refViewNr))
            useFeature = true;

        // 或者在至少一个邻域视角中可见，两者满足其一即可
        for (IndexSet::const_iterator id = neighViews.begin();
            useFeature == false && id != neighViews.end(); ++id){
            if (features[i].contains_view_id(*id))
                useFeature = true;
        }
        if (!useFeature)
            continue;

        // 特征的三维坐标
        math::Vec3f featPos(features[i].pos);

        // 特征要在参考视角中可见
        if (!refV->pointInFrustum(featPos))
            continue;

        /* Check if feature is inside AABB. */
        // 特征要在设定的空间范围内
        if (!math::geom::point_box_overlap(featPos,
            this->settings.aabbMin, this->settings.aabbMax))
            continue;

        /* Start processing the feature. */
        processed += 1;

        // 将三维特征点投影到对应尺度的参考图像上
        math::Vec2f pixPosF = refV->worldToScreenScaled(featPos);
        int const x = math::round(pixPosF[0]);
        int const y = math::round(pixPosF[1]);

        // 初始的深度设置为三维特征点到参考图像的距离
        float initDepth = (featPos - refV->camPos).norm();

        // 对三维点的深度进行优化，同时得到法向量，深度图，以及优化的可信度
        PatchOptimization patch(views, settings, x, y, initDepth,
            0.f, 0.f, neighViews, IndexSet());
        patch.doAutoOptimization();
        float conf = patch.computeConfidence();
        if (conf <= 0.0f)
            continue;

        // 深度优化成功
        /* Feature depth optimization was successful. */
        success += 1;
        // 特征点在图像中的索引
        int const index = y * this->width + x;
        // 获取深度值
        float depth = patch.getDepth();
        // 获取法向量
        math::Vec3f normal = patch.getNormal();
        // 如果当前优化的置信度大于之前的置信度，即更可靠
        if (refV->confImg->at(index) < conf)
        {
            if (refV->confImg->at(index) <= 0)
                ++progress.filled;

            refV->depthImg->at(index) = depth;
            refV->normalImg->at(index, 0) = normal[0];
            refV->normalImg->at(index, 1) = normal[1];
            refV->normalImg->at(index, 2) = normal[2];
            refV->dzImg->at(index, 0) = patch.getDzI();
            refV->dzImg->at(index, 1) = patch.getDzJ();
            refV->confImg->at(index) = conf;

            // 将当前重建的结果添加到队列(prQueue)中
            QueueData tmpData;
            tmpData.confidence = conf;
            tmpData.depth = depth;
            tmpData.dz_i = patch.getDzI();
            tmpData.dz_j = patch.getDzJ();
            tmpData.localViewIDs = patch.getLocalViewIDs();
            tmpData.x = x;
            tmpData.y = y;
            prQueue.push(tmpData);
        }
    }
    if (!settings.quiet)
        std::cout << "Processed " << processed << " features, from which "
                  << success << " succeeded optimization." << std::endl;
}

void
DMRecon::processQueue()
{
    progress.status = RECON_QUEUE;
    if (progress.cancelled)  return;

    // 当前的参考视角
    SingleView::Ptr refV = this->views[settings.refViewNr];

    if (!settings.quiet)
        std::cout << "Process queue ..." << std::endl;

    size_t count = 0, lastStatus = 1;
    progress.queueSize = prQueue.size();

    // 初始状态
    if (!settings.quiet)
        std::cout << "Count: " << std::setw(8) << count
                  << "  filled: " << std::setw(8) << progress.filled  // 填充深度值的像素个数
                  << "  Queue: " << std::setw(8) << progress.queueSize  // 队列中的种子点个数
                  << std::endl;
    lastStatus = progress.filled;

    // 队列不为空
    while (!prQueue.empty() && !progress.cancelled)
    {
        progress.queueSize = prQueue.size();

        // 每填充1000个像素并且填充像素在增加，则打印输出填充结果
        if ((progress.filled % 1000 == 0) && (progress.filled != lastStatus)) {
            if (!settings.quiet)
                std::cout << "Count: " << std::setw(8) << count
                          << "  filled: " << std::setw(8) << progress.filled
                          << "  Queue: " << std::setw(8) << progress.queueSize
                          << std::endl;
            lastStatus = progress.filled;
        }

        // 取出第一个种子点
        QueueData tmpData = prQueue.top();
        // 去除种子点
        prQueue.pop();
        ++count;
        float x = tmpData.x;
        float y = tmpData.y;
        int index = y * this->width + x;

        // 一般情况下这里应该相等 todo ?????
        if (refV->confImg->at(index) > tmpData.confidence) {
            continue ;
        }

        PatchOptimization patch(views, settings, x, y, tmpData.depth,
            tmpData.dz_i, tmpData.dz_j, neighViews, tmpData.localViewIDs);
        patch.doAutoOptimization();
        tmpData.confidence = patch.computeConfidence();
        if (tmpData.confidence == 0) {
            continue;
        }

        float new_depth = patch.getDepth();
        tmpData.depth = new_depth;
        tmpData.dz_i = patch.getDzI();
        tmpData.dz_j = patch.getDzJ();
        math::Vec3f normal = patch.getNormal();
        tmpData.localViewIDs = patch.getLocalViewIDs();

        // 从未重建的点进行了重建
        if (refV->confImg->at(index) <= 0) {
            ++progress.filled;
        }
        // 性能进行了提升
        if (refV->confImg->at(index) < tmpData.confidence) {
            refV->depthImg->at(index) = tmpData.depth;
            refV->normalImg->at(index, 0) = normal[0];
            refV->normalImg->at(index, 1) = normal[1];
            refV->normalImg->at(index, 2) = normal[2];
            refV->dzImg->at(index, 0) = tmpData.dz_i;
            refV->dzImg->at(index, 1) = tmpData.dz_j;
            refV->confImg->at(index) = tmpData.confidence;

            // left
            tmpData.x = x - 1; tmpData.y = y;
            index = tmpData.y * this->width + tmpData.x;
            if (refV->confImg->at(index) < tmpData.confidence - 0.05f ||
                refV->confImg->at(index) == 0.f)
            {
                prQueue.push(tmpData);
            }
            // right
            tmpData.x = x + 1; tmpData.y = y;
            index = tmpData.y * this->width + tmpData.x;
            if (refV->confImg->at(index) < tmpData.confidence - 0.05f ||
                refV->confImg->at(index) == 0.f)
            {
                prQueue.push(tmpData);
            }
            // top
            tmpData.x = x; tmpData.y = y - 1;
            index = tmpData.y * this->width + tmpData.x;
            if (refV->confImg->at(index) < tmpData.confidence - 0.05f ||
                refV->confImg->at(index) == 0.f)
            {
                prQueue.push(tmpData);
            }
            // bottom
            tmpData.x = x; tmpData.y = y + 1;
            index = tmpData.y * this->width + tmpData.x;
            if (refV->confImg->at(index) < tmpData.confidence - 0.05f ||
                refV->confImg->at(index) == 0.f)
            {
                prQueue.push(tmpData);
            }
        }
    }
}

MVS_NAMESPACE_END
