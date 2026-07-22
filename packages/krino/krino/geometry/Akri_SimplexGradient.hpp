/*
 * Akri_SimplexCalcs.hpp
 *
 *  Created on: Sep 11, 2025
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_GEOMETRY_AKRI_SIMPLEXGRADIENT_HPP_
#define KRINO_KRINO_GEOMETRY_AKRI_SIMPLEXGRADIENT_HPP_

#include <array>
#include <stk_math/StkVector.hpp>

namespace krino {

void calculate_triangle2d_gradient_and_area(const std::array<stk::math::Vector3d,3> & x, const std::array<double,3> & d, stk::math::Vector3d & grad, double & area);
stk::math::Vector3d calculate_triangle2d_gradient(const std::array<stk::math::Vector3d,3> & x, const std::array<double,3> & d);
double calculate_triangle2d_gradient_magnitude(const std::array<stk::math::Vector3d,3> & x, const std::array<double,3> & d);

void calculate_tetrahedron_gradient_and_volume(const std::array<stk::math::Vector3d,4> & x, const std::array<double,4> & d, stk::math::Vector3d & grad, double & vol);
stk::math::Vector3d calculate_tetrahedron_gradient(const std::array<stk::math::Vector3d,4> & x, const std::array<double,4> & d);
double calculate_tetrahedron_gradient_magnitude(const std::array<stk::math::Vector3d,4> & x, const std::array<double,4> & d);

}



#endif /* KRINO_KRINO_GEOMETRY_AKRI_SIMPLEXGRADIENT_HPP_ */
