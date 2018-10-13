/*
 * Copyright (C) 2015, Simon Fuhrmann
 * TU Darmstadt - Graphics, Capture and Massively Parallel Computing
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef MVE_PBRTFILE_HEADER
#define MVE_PBRTFILE_HEADER

#include "core/mesh.h"

CORE_NAMESPACE_BEGIN
CORE_GEOM_NAMESPACE_BEGIN

/**
 * Saves a PBRT compatible mesh from a triangle mesh.
 */
void
save_pbrt_mesh (TriangleMesh::ConstPtr mesh, std::string const& filename);

CORE_GEOM_NAMESPACE_END
CORE_NAMESPACE_END

#endif /* MVE_PBRTFILE_HEADER */
