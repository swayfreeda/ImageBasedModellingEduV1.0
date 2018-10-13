/*
 * Copyright (C) 2015, Ronny Klowsky, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef DMRECON_PATCH_OPTIMIZATION_H
#define DMRECON_PATCH_OPTIMIZATION_H

#include <iostream>

#include "math/vector.h"
#include "mvs/defines.h"
#include "mvs/patch_sampler.h"
#include "mvs/single_view.h"
#include "mvs/local_view_selection.h"

MVS_NAMESPACE_BEGIN

struct Status
{
    std::size_t iterationCount;
    bool converged;
    bool optiSuccess;
};

/**
 * \description 对patch的深度进行优化
 */
class PatchOptimization
{
public:

    // constructor
    PatchOptimization(
        std::vector<SingleView::Ptr> const& _views,
        Settings const& _settings,
        int _x,          // Pixel position
        int _y,
        float _depth,
        float _dzI,
        float _dzJ,
        IndexSet const& _globalViewIDs,
        IndexSet const& _localViewIDs);

    // todo 颜色尺度是什么意思？？？
    void computeColorScale();
    float computeConfidence();
    float derivNorm();
    void doAutoOptimization();
    float getDepth() const;
    float getDzI() const;
    float getDzJ() const;

    // 获取局部视角
    IndexSet const& getLocalViewIDs() const;

    // 获取法向量
    math::Vec3f getNormal() const;
    float objFunValue();

    // 优化深度
    void optimizeDepthOnly();

    // 优化深度和法向量
    void optimizeDepthAndNormal();

private:
    std::vector<SingleView::Ptr> const& views;
    Settings const& settings;
    // initial values and settings
    const int midx;
    const int midy;

    float depth;
    float dzI, dzJ;// represents patch normal
    std::map<std::size_t, math::Vec3f, std::less<std::size_t> > colorScale;
    Status status;

    PatchSampler::Ptr sampler;

    // patch点的x 和 y 坐标
    std::vector<int> ii, jj;

    // 像素权重
    std::vector<float> pixel_weight;
    LocalViewSelection localVS;
};

inline float
PatchOptimization::getDepth() const
{
    return depth;
}

inline float
PatchOptimization::getDzI() const
{
    return dzI;
}

inline float
PatchOptimization::getDzJ() const
{
    return dzJ;
}

inline IndexSet const&
PatchOptimization::getLocalViewIDs() const
{
    return localVS.getSelectedIDs();
}

inline math::Vec3f
PatchOptimization::getNormal() const
{
    return sampler->getPatchNormal();
}

MVS_NAMESPACE_END

#endif
