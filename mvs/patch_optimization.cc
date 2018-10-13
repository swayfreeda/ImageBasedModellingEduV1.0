/*
 * Copyright (C) 2015, Ronny Klowsky, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include "util/strings.h"
#include "math/algo.h"
#include "math/defines.h"
#include "math/matrix.h"
#include "math/matrix_tools.h"
#include "math/vector.h"
#include "mvs/patch_optimization.h"
#include "mvs/settings.h"

MVS_NAMESPACE_BEGIN

PatchOptimization::PatchOptimization(
    std::vector<SingleView::Ptr> const& _views,
    Settings const& _settings,
    int _x,          // Pixel position
    int _y,
    float _depth,
    float _dzI,
    float _dzJ,
    IndexSet const & _globalViewIDs,
    IndexSet const & _localViewIDs) :
    views(_views),
    settings(_settings),
    midx(_x),
    midy(_y),
    depth(_depth),
    dzI(_dzI),
    dzJ(_dzJ),
    sampler(PatchSampler::create(views, settings, midx, midy, depth, dzI, dzJ)),
    ii(sqr(settings.filterWidth)),
    jj(sqr(settings.filterWidth)),
    localVS(views, settings, _globalViewIDs, _localViewIDs, sampler)
{
    status.iterationCount = 0;
    status.optiSuccess = true;
    status.converged = false;

    // 初始化成功，计算出参考视角中的
    if (!sampler->success[settings.refViewNr]) {
        // Sampler could not be initialized properly
        status.optiSuccess = false;
        return;
    }

    std::size_t count = 0;

    // patch 宽度的一半
    int halfFW = (int) settings.filterWidth / 2;
    // 像素权重
    pixel_weight.resize(sampler->getNrSamples());

    // 像素权重初始化
    for (int j = -halfFW; j <= halfFW; ++j)
        for (int i = -halfFW; i <= halfFW; ++i) {
            ii[count] = i;
            jj[count] = j;
            pixel_weight[count] = 1.f;
            ++count;
        }

    // 进行局部视角选择
    localVS.performVS();
    if (!localVS.success) {
        status.optiSuccess = false;
        return;
    }

    // brute force initialize all colorScale entries
    float masterMeanCol = sampler->getMasterMeanColor();
    for (std::size_t idx = 0; idx < views.size(); ++idx) {
        colorScale[idx] = math::Vec3f(1.f / masterMeanCol);
    }

    // 计算每个视角的颜色尺度 ??????????????????
    computeColorScale();
}

void
PatchOptimization::computeColorScale()
{
    if (!settings.useColorScale)
        return;
    // 参考视角中的样本点的颜色
    Samples const & mCol = sampler->getMasterColorSamples();
    // 局部视角
    IndexSet const & neighIDs = localVS.getSelectedIDs();
    IndexSet::const_iterator id;

    // 对于每一个局部视角
    for (id = neighIDs.begin(); id != neighIDs.end(); ++id) {
        // just copied from old mvs:
        // 获取邻域样本点的颜色值
        Samples const & nCol = sampler->getNeighColorSamples(*id);
        if (!sampler->success[*id])
            return;
        // 对于每一个颜色通道
        for (int c = 0; c < 3; ++c) {
            float ab = 0.f;
            float aa = 0.f;
            // 每个样本点
            for (std::size_t i = 0; i < mCol.size(); ++i) {
                ab += (mCol[i][c] - nCol[i][c] * colorScale[*id][c]) * nCol[i][c];
                aa += sqr(nCol[i][c]);
            }
            if (std::abs(aa) > 1e-6) {
                colorScale[*id][c] += ab / aa;
                if (colorScale[*id][c] > 1e3)
                    status.optiSuccess = false;
            }
            else
                status.optiSuccess = false;
        }
    }
}

float
PatchOptimization::computeConfidence()
{
    SingleView::Ptr refV = views[settings.refViewNr];
    if (!status.converged)
        return 0.f;

    /* Compute mean NCC between reference view and local neighbors,
       where each NCC has to be higher than acceptance NCC */
    float meanNCC = 0.f;
    IndexSet const & neighIDs = localVS.getSelectedIDs();
    IndexSet::const_iterator id;
    for (id = neighIDs.begin(); id != neighIDs.end(); ++id) {
        meanNCC += sampler->getFastNCC(*id);
    }
    meanNCC /= neighIDs.size();

    float score = (meanNCC - settings.acceptNCC) /
                  (1.f - settings.acceptNCC);

    /* Compute angle between estimated surface normal and view direction
       and weight current score with dot product */
    math::Vec3f viewDir(refV->viewRayScaled(midx, midy));
    math::Vec3f normal(sampler->getPatchNormal());
    float dotP = - normal.dot(viewDir);
    if (dotP < 0.2f) {
        return 0.f;
    }
    return score;
}

float
PatchOptimization::derivNorm()
{
    IndexSet const & neighIDs = localVS.getSelectedIDs();
    IndexSet::const_iterator id;
    std::size_t nrSamples = sampler->getNrSamples();

    float norm(0);
    for (id = neighIDs.begin(); id != neighIDs.end(); ++id)
    {
        Samples nCol, nDeriv;
        sampler->fastColAndDeriv(*id, nCol, nDeriv);
        if (!sampler->success[*id]) {
            status.optiSuccess = false;
            return -1.f;
        }

        math::Vec3f cs(colorScale[*id]);
        for (std::size_t i = 0; i < nrSamples; ++i) {
            norm += pixel_weight[i] * (cs.cw_mult(nDeriv[i])).square_norm();
        }
    }
    return norm;
}

void
PatchOptimization::doAutoOptimization()
{
    if (!localVS.success || !status.optiSuccess) {
        return;
    }

    // first four iterations only refine depth
    while ((status.iterationCount < 4) && (status.optiSuccess)) {
        optimizeDepthOnly();
        ++status.iterationCount;
    }

    bool viewRemoved = false;
    bool converged = false;
    while (status.iterationCount < settings.maxIterations &&
        localVS.success && status.optiSuccess)
    {
        // 邻域视角
        IndexSet const & neighIDs = localVS.getSelectedIDs();
        IndexSet::const_iterator id;

        // 保存当前每个邻域视角的NCC值
        std::vector<float> oldNCC;
        for (id = neighIDs.begin(); id != neighIDs.end(); ++id) {
            // 快速的计算参考视角和当前视角的NCC值
            oldNCC.push_back(sampler->getFastNCC(*id));
        }

        // 当有视角移除或者优化一定次数之后，重新进行深度和法向量的优化，以及颜色尺度的计算
        // 否则仅仅优化深度值
        status.optiSuccess = false;
        if (status.iterationCount % 5 == 4 || viewRemoved) {
            optimizeDepthAndNormal();
            computeColorScale();
            viewRemoved = false;
        }
        else {
            optimizeDepthOnly();
        }

        if (!status.optiSuccess)
            return;

        // 检查每个视角
        converged = true;
        std::size_t count = 0;
        IndexSet toBeReplaced;
        for (id = neighIDs.begin(); id != neighIDs.end(); ++id, ++count) {
            // 快速计算NCC
            float ncc = sampler->getFastNCC(*id);

            // 如果当前ncc的值比前一次迭代ncc的值大于一定的阈值，则未收敛
            if (std::abs(ncc - oldNCC[count]) > settings.minRefineDiff)
                converged = false;

            // ncc的值小于可接受的范围或者优化次数达到一定值但是尚未收敛
            if ((ncc < settings.acceptNCC) || (status.iterationCount == 14 &&
                 std::abs(ncc - oldNCC[count]) > settings.minRefineDiff)) {
                toBeReplaced.insert(*id);
                viewRemoved = true;
            }
        }

        // 如果有视角移除，重新进行局部视角选择和重新计算颜色尺度
        if (viewRemoved) {
            localVS.replaceViews(toBeReplaced);
            if (!localVS.success) {
                return;
            }
            computeColorScale();
        }
        else if (!status.optiSuccess) {
            // all views valid but optimization failed -> give up
            return;
        }
        else if (converged) {
            status.converged = true;
            return;
        }
        ++status.iterationCount;
    }
}

float
PatchOptimization::objFunValue()
{
    Samples const & mCol = sampler->getMasterColorSamples();
    std::size_t nrSamples = sampler->getNrSamples();
    IndexSet const & neighIDs = localVS.getSelectedIDs();
    IndexSet::const_iterator id;
    float obj = 0.f;
    for (id = neighIDs.begin(); id != neighIDs.end(); ++id) {
        Samples const & nCol = sampler->getNeighColorSamples(*id);
        if (!sampler->success[*id])
            return -1.f;
        math::Vec3f cs(colorScale[*id]);
        for (std::size_t i = 0; i < nrSamples; ++i) {
            obj += pixel_weight[i] * (mCol[i] - cs.cw_mult(nCol[i])).square_norm();
        }
    }
    return obj;
}

void
PatchOptimization::optimizeDepthOnly()
{
    float numerator = 0.f;
    float denom = 0.f;
    // 参考视角的样本颜色值
    Samples const & mCol = sampler->getMasterColorSamples();
    // 所有的局部视角
    IndexSet const & neighIDs = localVS.getSelectedIDs();
    // patch三维点的个数
    IndexSet::const_iterator id;
    std::size_t nrSamples = sampler->getNrSamples();

    // 对于每一个邻域视角
    for (id = neighIDs.begin(); id != neighIDs.end(); ++id)
    {
        // 计算邻域视角采样点的颜色和梯度
        Samples nCol, nDeriv;
        sampler->fastColAndDeriv(*id, nCol, nDeriv);
        if (!sampler->success[*id]) {
            status.optiSuccess = false;
            return;
        }

        // 颜色尺度
        math::Vec3f cs(colorScale[*id]);

        // 对于每一个样本点 比较核心！！！！！！！！
        for (std::size_t i = 0; i < nrSamples; ++i) {
            numerator += pixel_weight[i] * (cs.cw_mult(nDeriv[i])).dot
                (mCol[i] - cs.cw_mult(nCol[i]));
            denom += pixel_weight[i] * (cs.cw_mult(nDeriv[i])).square_norm();
        }
    }

    if (denom > 0) {
        depth += numerator / denom;
        sampler->update(depth, dzI, dzJ);
        if (sampler->success[settings.refViewNr])
            status.optiSuccess = true;
        else
            status.optiSuccess = false;
    }
}

void
PatchOptimization::optimizeDepthAndNormal()
{
    if (!localVS.success) {
        return;
    }
    IndexSet const & neighIDs = localVS.getSelectedIDs();
    std::size_t nrSamples = sampler->getNrSamples();

    // Solve linear system A*x = b using Moore-Penrose pseudoinverse
    // Fill matrix ATA and vector ATb:
    math::Matrix3d ATA(0.f);
    math::Vec3d ATb(0.f);
    Samples const & mCol = sampler->getMasterColorSamples();
    IndexSet::const_iterator id;
    std::size_t row = 0;
    for (id = neighIDs.begin(); id != neighIDs.end(); ++id)
    {
        Samples nCol, nDeriv;
        sampler->fastColAndDeriv(*id, nCol, nDeriv);
        if (!sampler->success[*id]) {
            status.optiSuccess = false;
            return;
        }
        math::Vec3f cs(colorScale[*id]);
        for (std::size_t i = 0; i < nrSamples; ++i) {
            for (int c = 0; c < 3; ++c) {
                math::Vec3f a_i;
                a_i[0] = pixel_weight[i] *         cs[c] * nDeriv[i][c];
                a_i[1] = pixel_weight[i] * ii[i] * cs[c] * nDeriv[i][c];
                a_i[2] = pixel_weight[i] * jj[i] * cs[c] * nDeriv[i][c];
                float b_i = pixel_weight[i] * (mCol[i][c] - cs[c] * nCol[i][c]);
                assert(!MATH_ISINF(a_i[0]));
                assert(!MATH_ISINF(a_i[1]));
                assert(!MATH_ISINF(a_i[2]));
                assert(!MATH_ISINF(b_i));
                ATA(0,0) += a_i[0] * a_i[0];
                ATA(0,1) += a_i[0] * a_i[1];
                ATA(0,2) += a_i[0] * a_i[2];
                ATA(1,1) += a_i[1] * a_i[1];
                ATA(1,2) += a_i[1] * a_i[2];
                ATA(2,2) += a_i[2] * a_i[2];
                ATb += a_i * b_i;
                ++row;
            }
        }
    }
    ATA(1,0) = ATA(0,1); ATA(2,0) = ATA(0,2); ATA(2,1) = ATA(1,2);
    double detATA = math::matrix_determinant(ATA);
    if (detATA == 0.f) {
        status.optiSuccess = false;
        return;
    }
    math::Matrix3d ATAinv = math::matrix_inverse(ATA, detATA);
    math::Vec3f X = ATAinv * ATb;

    // update depth and normal
    dzI += X[1];
    dzJ += X[2];
    depth += X[0];
    sampler->update(depth, dzI, dzJ);
    if (sampler->success[settings.refViewNr])
        status.optiSuccess = true;
    else
        status.optiSuccess = false;
}


MVS_NAMESPACE_END
