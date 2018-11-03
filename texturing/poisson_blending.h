/*
 * Copyright (C) 2015, Nils Moehrle
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef TEX_POISSONBLENDING_HEADER
#define TEX_POISSONBLENDING_HEADER

#include "core/image.h"

void
poisson_blend(core::FloatImage::ConstPtr src, core::ByteImage::ConstPtr mask,
    core::FloatImage::Ptr dest, float alpha);


#endif /* TEX_POISSONBLENDING_HEADER */
