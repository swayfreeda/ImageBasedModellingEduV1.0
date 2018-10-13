/*
 * Copyright (C) 2015, Ronny Klowsky, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef DMRECON_DMRECON_H
#define DMRECON_DMRECON_H

#include <fstream>
#include <string>
#include <vector>
#include <queue>

#include "core/bundle.h"
#include "core/image.h"
#include "core/scene.h"
#include "mvs/defines.h"
#include "mvs/patch_optimization.h"
#include "mvs/single_view.h"
#include "mvs/progress.h"

MVS_NAMESPACE_BEGIN

// 优先级队列
struct QueueData
{
    // 像素坐标
    int x;
    int y;
    // 置信度
    float confidence;
    // 深度值
    float depth;
    //x,y方向上三维空间中的深度变化率
    float dz_i, dz_j;

    // 局部视角，用于进行depeth, dz_i, dz_j的优化
    IndexSet localViewIDs;

    bool operator< (const QueueData& rhs) const;
};

class DMRecon
{
public:
    DMRecon(core::Scene::Ptr scene, Settings const& settings);

    std::size_t getRefViewNr() const;
    Progress const& getProgress() const;
    Progress& getProgress();
    void start();

private:
    core::Scene::Ptr scene;
    core::Bundle::ConstPtr bundle;
    std::vector<SingleView::Ptr> views;

    Settings settings;
    std::priority_queue<QueueData> prQueue;
    IndexSet neighViews;
    std::vector<SingleView::Ptr> imgNeighbors;
    int width;
    int height;
    Progress progress;

    void analyzeFeatures();
    void globalViewSelection();
    void processFeatures();
    void processQueue();
    void refillQueueFromLowRes();
};

/* ------------------------- Implementation ----------------------- */

inline bool
QueueData::operator< (const QueueData& rhs) const
{
    return (confidence < rhs.confidence);
}

inline const Progress&
DMRecon::getProgress() const
{
    return progress;
}

inline Progress&
DMRecon::getProgress()
{
    return progress;
}

inline std::size_t
DMRecon::getRefViewNr() const
{
    return settings.refViewNr;
}

MVS_NAMESPACE_END

#endif
