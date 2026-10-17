// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_SMOOTHING_AKRI_KINEMATICELEMENT_HPP_
#define KRINO_KRINO_SMOOTHING_AKRI_KINEMATICELEMENT_HPP_

#include <array>
#include <cmath>
#include <limits>

#include <Akri_KinematicUtils.hpp>
#include <stk_math/StkVector.hpp>

namespace krino {

// Element energy and nodal forces for any kinematic model.  A MODEL provides energy(F, J) and
// stress(F, J), the energy density and its first Piola-Kirchhoff stress P = dW/dF, for the 3D
// Tensor and the 2D Tensor2d.  The reference element is either an ideal element of size refSize or
// explicit reference coordinates.  A non-invertible element has infinite energy and zero forces.

template<class MODEL, typename ELEMCOORDS>
double tet4_element_energy(const MODEL & model, const double refSize, const ELEMCOORDS & cur)
{
  const auto & [F, J] = kinematicUtils::compute_deformation_gradient_and_jacobian_tet4(refSize, cur);
  if (kinematicUtils::is_invertible(J))
    return kinematicUtils::element_energy(model.energy(F, J), kinematicUtils::tet4_volume(refSize));
  return std::numeric_limits<double>::max();
}

template<class MODEL>
double tet4_element_energy(const MODEL & model, const std::array<stk::math::Vector3d, 4> & ref, const std::array<stk::math::Vector3d, 4> & cur)
{
  const auto & [F, J] = kinematicUtils::compute_deformation_gradient_and_jacobian_tet4(ref, cur);
  if (kinematicUtils::is_invertible(J))
    return kinematicUtils::element_energy(model.energy(F, J), kinematicUtils::tet4_volume(ref));
  return std::numeric_limits<double>::max();
}

template<class MODEL, typename ELEMCOORDS>
std::array<double, 12> tet4_element_forces(const MODEL & model, const double refSize, const ELEMCOORDS & cur)
{
  const auto & [F, J] = kinematicUtils::compute_deformation_gradient_and_jacobian_tet4(refSize, cur);
  if (kinematicUtils::is_invertible(J))
    return kinematicUtils::tet4_forces_from_stress(refSize, model.stress(F, J));
  return {};
}

template<class MODEL>
std::array<double, 12> tet4_element_forces(const MODEL & model, const std::array<stk::math::Vector3d, 4> & ref, const std::array<stk::math::Vector3d, 4> & cur)
{
  const auto & [F, J] = kinematicUtils::compute_deformation_gradient_and_jacobian_tet4(ref, cur);
  if (kinematicUtils::is_invertible(J))
    return kinematicUtils::tet4_forces_from_stress(ref, model.stress(F, J));
  return {};
}

template<class MODEL, typename ELEMCOORDS>
double tri3_2d_element_energy(const MODEL & model, const double refSize, const ELEMCOORDS & cur)
{
  const auto & [F, J] = kinematicUtils::compute_deformation_gradient_and_jacobian_tri3_2d(refSize, cur);
  if (kinematicUtils::is_invertible(J))
    return kinematicUtils::element_energy(model.energy(F, J), kinematicUtils::tri3_area(refSize));
  return std::numeric_limits<double>::max();
}

template<class MODEL, typename ELEMCOORDS>
std::array<double, 6> tri3_2d_element_forces(const MODEL & model, const double refSize, const ELEMCOORDS & cur)
{
  const auto & [F, J] = kinematicUtils::compute_deformation_gradient_and_jacobian_tri3_2d(refSize, cur);
  if (kinematicUtils::is_invertible(J))
    return kinematicUtils::tri3_2d_forces_from_stress(refSize, model.stress(F, J));
  return {};
}

template<class MODEL>
std::array<double, 6> tri3_2d_element_forces(const MODEL & model, const std::array<stk::math::Vector2d, 3> & ref, const std::array<stk::math::Vector2d, 3> & cur)
{
  const auto & [F, J] = kinematicUtils::compute_deformation_gradient_and_jacobian_tri3_2d(ref, cur);
  if (kinematicUtils::is_invertible(J))
    return kinematicUtils::tri3_2d_forces_from_stress(ref, model.stress(F, J));
  return {};
}

// Anisotropic tri3 and tet4 forms.  The metric M = F_M^T F_M, with metric mapping F_M = M_loc^1/2 Q^T
// from the spectral decomposition M = Q M_loc Q^T (metric::compute_metric_mapping and
// compute_metric_mapping_2d), maps the element's ideal reference configuration X to the unit
// configuration Y = F_M X, where Y is the unit equilateral triangle or tet.  The model is evaluated on
// the deformation gradient in unit space, F_U = F_M F F_M^-1 = F_M dx/dY, with F = dx/dX.  The energy
// is integrated over the reference element X, of area or volume vol(Y)/|det(F_M)|, and the forces are
// integrated over X with the stress pushed back to the reference configuration, P = F_M^T P_U F_M^-T.
// Since dN/dX = F_M^T dN/dY, these are the unit-element forces of |det(F_M)|^-1 F_M^T P_U.
//
// F_M must have det(F_M) > 0, as compute_metric_mapping and compute_metric_mapping_2d ensure: one with
// det(F_M) < 0 reflects the ideal element, so that elements of positive orientation appear inverted.
// Otherwise F_M is not unique (the eigenvectors' order and signs are arbitrary), but any two choices
// differ by a rotation, and the energy and forces depend on it only through M and
// det(F_M) = sqrt(det(M)), so they are the same for any choice.

template<class MODEL, typename ELEMCOORDS>
double tri3_2d_element_energy(const MODEL & model, const double * metricMapping, const ELEMCOORDS & cur)
{
  const minitensor::Tensor<double, 2> Fm(metricMapping);
  const minitensor::Tensor<double, 2> Fu = kinematicUtils::compute_anisotropic_deformation_gradient_tri3_2d(Fm, cur);
  const double Ju = minitensor::det(Fu);
  if (kinematicUtils::is_invertible(Ju))
    return kinematicUtils::element_energy(model.energy(Fu, Ju), kinematicUtils::tri3_area(1.) / std::abs(minitensor::det(Fm)));
  return std::numeric_limits<double>::max();
}

template<class MODEL, typename ELEMCOORDS>
std::array<double, 6> tri3_2d_element_forces(const MODEL & model, const double * metricMapping, const ELEMCOORDS & cur)
{
  const minitensor::Tensor<double, 2> Fm(metricMapping);
  const minitensor::Tensor<double, 2> Fu = kinematicUtils::compute_anisotropic_deformation_gradient_tri3_2d(Fm, cur);
  const double Ju = minitensor::det(Fu);
  if (kinematicUtils::is_invertible(Ju))
    return kinematicUtils::tri3_2d_forces_from_stress(1., (1./std::abs(minitensor::det(Fm))) * (minitensor::transpose(Fm) * model.stress(Fu, Ju)));
  return {};
}

template<class MODEL, typename ELEMCOORDS>
double tet4_element_energy(const MODEL & model, const double * metricMapping, const ELEMCOORDS & cur)
{
  const minitensor::Tensor<double, 3> Fm(metricMapping);
  const minitensor::Tensor<double, 3> Fu = kinematicUtils::compute_anisotropic_deformation_gradient_tet4(Fm, cur);
  const double Ju = minitensor::det(Fu);
  if (kinematicUtils::is_invertible(Ju))
    return kinematicUtils::element_energy(model.energy(Fu, Ju), kinematicUtils::tet4_volume(1.) / std::abs(minitensor::det(Fm)));
  return std::numeric_limits<double>::max();
}

template<class MODEL, typename ELEMCOORDS>
std::array<double, 12> tet4_element_forces(const MODEL & model, const double * metricMapping, const ELEMCOORDS & cur)
{
  const minitensor::Tensor<double, 3> Fm(metricMapping);
  const minitensor::Tensor<double, 3> Fu = kinematicUtils::compute_anisotropic_deformation_gradient_tet4(Fm, cur);
  const double Ju = minitensor::det(Fu);
  if (kinematicUtils::is_invertible(Ju))
    return kinematicUtils::tet4_forces_from_stress(1., (1./std::abs(minitensor::det(Fm))) * (minitensor::transpose(Fm) * model.stress(Fu, Ju)));
  return {};
}

} // namespace krino

#endif /* KRINO_KRINO_SMOOTHING_AKRI_KINEMATICELEMENT_HPP_ */
