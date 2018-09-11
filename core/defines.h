/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef CORE_DEFINES_HEADER
#define CORE_DEFINES_HEADER

#define CORE_NAMESPACE_BEGIN namespace core {
#define CORE_NAMESPACE_END }

#define CORE_IMAGE_NAMESPACE_BEGIN namespace image {
#define CORE_IMAGE_NAMESPACE_END }

#define CORE_GEOM_NAMESPACE_BEGIN namespace geom {
#define CORE_GEOM_NAMESPACE_END }

#ifndef STD_NAMESPACE_BEGIN
#   define STD_NAMESPACE_BEGIN namespace std {
#   define STD_NAMESPACE_END }
#endif

/** Multi-View Environment library. */
CORE_NAMESPACE_BEGIN
/** Image tools, loading and processing functions. */
CORE_IMAGE_NAMESPACE_BEGIN CORE_IMAGE_NAMESPACE_END
/** Geometric tools, loading and processing functions. */
CORE_GEOM_NAMESPACE_BEGIN CORE_GEOM_NAMESPACE_END
CORE_NAMESPACE_END

#endif /* MVE_DEFINES_HEADER */
