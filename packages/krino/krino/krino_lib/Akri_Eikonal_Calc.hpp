// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_EIKONAL_CALC_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_EIKONAL_CALC_HPP_
#include <stk_mesh/base/Entity.hpp>
#include <Akri_FieldRef.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

double calculate_gradient_magnitude_triangle(const std::array<stk::math::Vector3d,3> & x, const std::array<double,3> & d);

double calculate_gradient_magnitude_tetrahedron(const std::array<stk::math::Vector3d,4> & x, const std::array<double,4> & d);

double calculate_gradient_magnitude(const int npe,
    const stk::mesh::Entity * elem_nodes,
    const FieldRef dRef,
    const std::function<const stk::math::Vector3d &(stk::mesh::Entity)> & get_coordinates);

double calculate_gradient_magnitude(const int npe,
    const stk::mesh::Entity * elem_nodes,
    const FieldRef dRef,
    const FieldRef xRef);

double eikonal_solve_triangle(const std::array<stk::math::Vector3d,3> & x,
    const std::array<double,2> & d,
    const int sign,
    const double far,
    const double speed);

double eikonal_solve_triangle(const std::array<stk::math::Vector3d,3> & x,
    const std::array<double,2> & d,
    const double far,
    const double speed);

double eikonal_solve_tetrahedron(const std::array<stk::math::Vector3d,4> & x,
    const std::array<double,3> & dSigned,
    const int sign,
    const double far,
    const double speed);

double eikonal_solve_tetrahedron(const std::array<stk::math::Vector3d,4> & x,
    const std::array<double,3> & d,
    const double far,
    const double speed);

std::pair<double,double> eikonal_solve_triangle_for_distance_and_extension_speed(const std::array<stk::math::Vector3d,3> & x,
    const std::array<double,2> & d,
    const std::array<double,2> & extSpeed,
    const double far);

std::pair<double,double> eikonal_solve_triangle_for_distance_and_extension_speed(const std::array<stk::math::Vector3d,3> & x,
    const std::array<double,2> & dSigned,
    const std::array<double,2> & extSpeed,
    const int sign,
    const double far);

std::pair<double,double> eikonal_solve_tetrahedron_for_distance_and_extension_speed(const std::array<stk::math::Vector3d,4> & x,
    const std::array<double,3> & d,
    const std::array<double,3> & extSpeed,
    const double far);

std::pair<double,double> eikonal_solve_tetrahedron_for_distance_and_extension_speed(const std::array<stk::math::Vector3d,4> & x,
    const std::array<double,3> & dSigned,
    const std::array<double,3> & extSpeed,
    const int sign,
    const double far);


std::array<int,3> get_oriented_nodes_triangle(const int nodeToUpdate);
std::array<int,4> get_oriented_nodes_tetrahedron(const int nodeToUpdate);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_EIKONAL_CALC_HPP_ */
