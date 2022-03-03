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
#include <Akri_Vec.hpp>

namespace krino {

double calculate_gradient_magnitude_triangle(const std::array<Vector3d,3> & x, const std::array<double,3> & d);

double calculate_gradient_magnitude_tetrahedron(const std::array<Vector3d,4> & x, const std::array<double,4> & d);

double calculate_gradient_magnitude(const int npe,
    const stk::mesh::Entity * elem_nodes,
    const FieldRef dRef,
    const std::function<const Vector3d &(stk::mesh::Entity)> & get_coordinates);

double eikonal_solve_triangle(const std::array<Vector3d,3> & x,
    const std::array<double,3> & d,
    const int sign,
    const int dim,
    const double speed);

double eikonal_solve_tetrahedron(const std::array<Vector3d,4> & x,
    const std::array<double,4> & d,
    const int sign,
    const double speed);


std::array<int,3> get_oriented_nodes_triangle(const int nodeToUpdate);
std::array<int,4> get_oriented_nodes_tetrahedron(const int nodeToUpdate);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_EIKONAL_CALC_HPP_ */
