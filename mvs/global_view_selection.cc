/*
 * Copyright (C) 2015, Ronny Klowsky, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "math/vector.h"
#include "mvs/global_view_selection.h"
#include "mvs/mvs_tools.h"
#include "mvs/settings.h"

MVS_NAMESPACE_BEGIN

GlobalViewSelection::GlobalViewSelection(
    std::vector<SingleView::Ptr> const& views,
    core::Bundle::Features const& features,
    Settings const& settings)
    : ViewSelection(settings)
    , views(views)
    , features(features)
{
    available.clear();
    available.resize(views.size(), true);

    // 参考视角设置成unavailable
    available[settings.refViewNr] = false;

    // 如果视角没有加载，则设置成unavailable
    for (std::size_t i = 0; i < views.size(); ++i)
        if (views[i] == nullptr)
            available[i] = false;
}

void
GlobalViewSelection::performVS()
{
    selected.clear();
    bool foundOne = true;
    while (foundOne && (selected.size() < settings.globalVSMax))
    {
        float maxBenefit = 0.f;
        std::size_t maxView = 0;
        foundOne = false;
        for (std::size_t i = 0; i < views.size(); ++i)
        {
            if (!available[i])
                continue;

            float benefit = benefitFromView(i);
            if (benefit > maxBenefit) {
                maxBenefit = benefit;
                maxView = i;
                foundOne = true;
            }
        }
        if (foundOne) {
            selected.insert(maxView);
            available[maxView] = false;
        }
    }
}

float
GlobalViewSelection::benefitFromView(std::size_t i)
{
    // 参考视角
    SingleView::Ptr refV = views[settings.refViewNr];

    // 需要计算的视角
    SingleView::Ptr tmpV = views[i];

    // 共享特征的个数
    std::vector<std::size_t> nFeatIDs = tmpV->getFeatureIndices();

    // 遍历所有的共享特征
    // Go over all features visible in view i and reference view
    float benefit = 0;
    for (std::size_t k = 0; k < nFeatIDs.size(); ++k) {
        float score = 1.f;
        // Parallax with reference view
        // 特征点的三维坐标
        math::Vec3f ftPos(features[nFeatIDs[k]].pos);
        // 计算两个视角的视差(三维点与相机位置确定的射线间的夹角）
        float plx = parallax(ftPos, refV, tmpV);

        // 如果视差小于一定阈值(10度）
        if (plx < settings.minParallax)
            score *= sqr(plx / 10.f);

        // Resolution compared to reference view
        float mfp = refV->footPrintScaled(ftPos);
        float nfp = tmpV->footPrint(ftPos);
        float ratio = mfp / nfp;

        // reference view的分辨率要比tmp view 大2倍以上
        if (ratio > 2.)
            ratio = 2. / ratio;
        else if (ratio > 1.)
            ratio = 1.;
        score *= ratio;

        // 该特征点与已经选定的视角的视差，与任何一个视角的残差小于特定值，都会减少score
        // Parallax with other selected views that see the same feature
        IndexSet::const_iterator citV;
        for (citV = selected.begin(); citV != selected.end(); ++citV) {
            plx = parallax(ftPos, views[*citV], tmpV);
            if (plx < settings.minParallax)
                score *= sqr(plx / 10.f);
        }
        benefit += score;
    }
    return benefit;
}


MVS_NAMESPACE_END
