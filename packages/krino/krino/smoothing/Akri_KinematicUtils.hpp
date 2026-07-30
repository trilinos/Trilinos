// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_KINEMATICUTILS_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_KINEMATICUTILS_HPP_
#include <Akri_SethHillConfig.hpp>
#include <MiniTensor.h>

namespace krino {
namespace kinematicUtils {

// Precompute ideal reference configurations for all elements from initial coords
template<typename ELEMCOORDS>
double compute_ideal_tet4_h(const ELEMCOORDS & coords,
    SmoothReference mode);
std::array<stk::math::Vector3d, 4> compute_ideal_tet4_reference(
    const double h);
std::array<stk::math::Vector3d, 4> compute_ideal_tet4_reference(
    const std::array<stk::math::Vector3d, 4> & coords,
    SmoothReference smoothRef);
std::array<stk::math::Vector2d, 3> compute_ideal_tri3_2d_reference(
    const double h);

template<typename ELEMCOORDS>
minitensor::Tensor<double, 3> compute_deformation_gradient_tet4(
    const double refSize,
    const ELEMCOORDS & cur);

minitensor::Tensor<double, 3> compute_deformation_gradient_tet4(
    const std::array<stk::math::Vector3d, 4> & ref,
    const std::array<stk::math::Vector3d, 4> & cur);

double tet4_volume(const double size);
double tet4_volume(const std::array<stk::math::Vector3d, 4> & ref);
double tri3_area(const double size);

template<typename ELEMCOORDS>
minitensor::Tensor<double, 2> compute_deformation_gradient_tri3_2d(
    const double refSize,
    const ELEMCOORDS & cur);

minitensor::Tensor<double, 2> compute_deformation_gradient_tri3_2d(
    const std::array<stk::math::Vector2d, 3> & ref,
    const std::array<stk::math::Vector2d, 3> & cur);

std::array<stk::math::Vector3d, 4>
tet4_forces_from_stress(const double refSize,
    const minitensor::Tensor<double, 3> & P);

std::array<stk::math::Vector3d, 4>
tet4_forces_from_stress(const std::array<stk::math::Vector3d, 4> & ref,
    const minitensor::Tensor<double, 3> & P);

std::array<double, 6>
tri3_2d_forces_from_stress(const double refSize, const minitensor::Tensor<double, 2> & P);

std::array<double, 6>
tri3_2d_forces_from_stress(const std::array<stk::math::Vector2d, 3> & ref,
    const minitensor::Tensor<double, 2> & P);

inline double pow_p23(const double x)
{
  return std::cbrt(x*x);
}

inline double pow_m23(const double x)
{
  return 1. / pow_p23(x);
}

template<minitensor::Index N>
bool is_invertible(const minitensor::Tensor<double, N> & A)
{
  const double detA = det(A);
  return detA > 0 && std::isfinite(detA);
}

template<minitensor::Index N>
std::pair<bool, minitensor::Tensor<double, N>> safe_inverse(const minitensor::Tensor<double, N> & A)
{
  if (is_invertible(A))
    return std::make_pair(true, minitensor::inverse(A));
  return std::make_pair(false, A);
}

}
}



#endif /* KRINO_KRINO_SMOOTHING_AKRI_KINEMATICUTILS_HPP_ */
