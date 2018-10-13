/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef DMRECON_LOCAL_VIEW_SELECTION_H
#define DMRECON_LOCAL_VIEW_SELECTION_H


#include "mvs/defines.h"
#include "mvs/settings.h"
#include "mvs/view_selection.h"
#include "mvs/patch_sampler.h"
#include "mvs/single_view.h"


MVS_NAMESPACE_BEGIN

/**
 * 局部视角选择是从全局视角中选取一些视角(最多4个)，局部视角是针对每个patch的
 * 因此更加灵活，也就是每个patch的局部视角可以是不同的，但是全局视角都是相同的。
 */
class LocalViewSelection : public ViewSelection
{
public:
    LocalViewSelection(
        std::vector<SingleView::Ptr> const& views,
        Settings const& settings,
        IndexSet const& globalViews,
        IndexSet const& propagated,
        PatchSampler::Ptr sampler);
    void performVS();
    void replaceViews(IndexSet const& toBeReplaced);

    bool success;

private:
    std::vector<SingleView::Ptr> const& views;
    PatchSampler::Ptr sampler;
};

MVS_NAMESPACE_END

#endif
