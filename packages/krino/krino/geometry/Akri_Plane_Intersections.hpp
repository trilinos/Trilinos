// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_PLANE_INTERSECTIONS_H_
#define KRINO_INCLUDE_AKRI_PLANE_INTERSECTIONS_H_
#include <stk_math/StkPlane.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

bool find_intersection_of_three_planes(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2, stk::math::Vector3d & point);

bool find_intersection_of_three_planes_within_tet(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2, stk::math::Vector3d & point);

bool find_intersection_of_two_planes_and_side_of_tet(const int side, const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, stk::math::Vector3d & intersectionPoint);

bool find_intersection_of_two_2D_planes(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, stk::math::Vector3d & intersectionPoint);

bool find_intersection_of_two_2D_planes_within_tri(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, stk::math::Vector3d & intersectionPoint);

stk::math::Vector3d triangle_parametric_coordinates_of_projected_point(const std::array<stk::math::Vector3d,3> & triCoords, const stk::math::Vector3d & pt);

bool within_tri_bounds(const stk::math::Vector3d & triangleParamCoords);
}

#endif /* KRINO_INCLUDE_AKRI_PLANE_INTERSECTIONS_H_ */
