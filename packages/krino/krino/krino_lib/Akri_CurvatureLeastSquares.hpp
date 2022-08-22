// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_CURVATURELEASTSQUARES_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_CURVATURELEASTSQUARES_HPP_
#include <array>
#include <vector>
#include "Akri_Vec.hpp"

namespace krino {

class Quaternion;

void set_rotation_matrix_for_rotating_normal_to_zDir(std::array<std::array<double,3>,3> & m, const Vector3d & normalDir);
Vector3d rotate_3d_vector(const std::array<std::array<double,3>,3> & m, const Vector3d & v);
Vector3d reverse_rotate_3d_vector(const std::array<std::array<double,3>,3> & m, const Vector3d & v);

Vector3d compute_patch_normal(const std::vector<Vector3d> & haloNodeLocs, const std::vector<std::array<int,2>> & haloSegments);
Vector3d compute_least_squares_curvature_times_normal(const std::vector<Vector3d> & haloNodeLocs, const std::vector<std::array<int,2>> & haloSegments);
Vector3d compute_least_squares_curvature_times_normal(const Vector3d & approximateNormal, const std::vector<Vector3d> & neighborNodeLocs);
Vector3d compute_least_squares_normal(const Vector3d & approximateNormal, const std::vector<Vector3d> & neighborNodeLocs);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_CURVATURELEASTSQUARES_HPP_ */
