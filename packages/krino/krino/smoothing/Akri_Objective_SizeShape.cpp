// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <Akri_Objective_SizeShape.hpp>

#include <Akri_DistributedVector.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_KinematicUtils.hpp>
#include <MiniTensor.h>

namespace krino {

using Tensor = minitensor::Tensor<double, 3>;
using Tensor2d = minitensor::Tensor<double, 2>;

static double sizeShape_element_pseudo_energy(const Tensor & F, const double J)
{
  const Tensor C = minitensor::transpose(F) * F;
  const double Jm23 = kinematicUtils::pow_m23(J);
  const Tensor Cbar = Jm23 * C;
  const double trCbar   = minitensor::trace(Cbar);

  const double etaShape = 1.0/3.0 * trCbar;
  const double etaSize = kinematicUtils::pow_p23(0.5*(J + 1./J));
  return etaShape*etaSize - 1.;
}

static double sizeShape_element_pseudo_energy_2d(const Tensor2d & F, const double J)
{
  const Tensor2d C = minitensor::transpose(F) * F;
  const double Jinv = 1./J;
  const Tensor2d Cbar = Jinv * C;
  const double trCbar   = minitensor::trace(Cbar);

  const double etaShape = 0.5 * trCbar;
  const double etaSize = 0.5*(J + Jinv);
  return etaShape*etaSize - 1.;
}

static Tensor sizeShape_element_pseudo_stress(const Tensor & F, const double J)
{
  const Tensor C = minitensor::transpose(F) * F;
  const Tensor FinvT = minitensor::transpose(minitensor::inverse(F));

  const double Jm23 = kinematicUtils::pow_m23(J);
  const Tensor Cbar = Jm23 * C;
  const double trCbar   = minitensor::trace(Cbar);

  const double etaShape = 1.0/3.0 * trCbar;
  const double etaSize = kinematicUtils::pow_p23(0.5*(J + 1./J));

  const Tensor dEtaShape = 2./3. * Jm23 * F - 2./9.*trCbar * FinvT;
  const Tensor dEtaSize = 1./3. * (J - 1./J) / std::cbrt(0.5*(J + 1./J)) * FinvT;

  return etaSize*dEtaShape + etaShape*dEtaSize;
}

static Tensor2d sizeShape_element_pseudo_stress_2d(const Tensor2d & F, const double J)
{
  const Tensor2d C = minitensor::transpose(F) * F;
  const Tensor2d FinvT = minitensor::transpose(minitensor::inverse(F));
  const double Jinv = 1./J;
  const Tensor2d Cbar = Jinv * C;
  const double trCbar = minitensor::trace(Cbar);

  const double etaShape = 0.5 * trCbar;
  const double etaSize = 0.5*(J + Jinv);

  const Tensor2d dEtaShape = Jinv * F - 0.5*trCbar * FinvT;
  const Tensor2d dEtaSize = 0.5 * (J - Jinv) * FinvT;

  return etaSize*dEtaShape + etaShape*dEtaSize;
}

std::array<double, 12>
SizeShapeObjective::compute_tet4_element_pseudo_forces(
    const std::array<stk::math::Vector3d, 4> & ref,
    const std::array<stk::math::Vector3d, 4> & current_coords)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor P = sizeShape_element_pseudo_stress(F, J);
  return kinematicUtils::tet4_forces_from_stress(ref, P);
}

template<typename ELEMCOORDS>
std::array<double, 12>
SizeShapeObjective::compute_tet4_element_pseudo_forces(const double refSize,
    const ELEMCOORDS & current_coords)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor P = sizeShape_element_pseudo_stress(F, J);
  return kinematicUtils::tet4_forces_from_stress(refSize, P);
}

double SizeShapeObjective::compute_tet4_element_pseudo_energy(
    const std::array<stk::math::Vector3d, 4> & ref,
    const std::array<stk::math::Vector3d, 4> & current_coords)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = sizeShape_element_pseudo_energy(F, J);
  return W * kinematicUtils::tet4_volume(ref);
}

template<typename ELEMCOORDS>
double SizeShapeObjective::compute_tet4_element_pseudo_energy(const double refSize,
    const ELEMCOORDS & current_coords)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = sizeShape_element_pseudo_energy(F, J);
  return W * kinematicUtils::tet4_volume(refSize);
}

template<typename ELEMCOORDS>
double SizeShapeObjective::compute_tri3_2d_element_pseudo_energy(const double refSize,
    const ELEMCOORDS & current_coords)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = sizeShape_element_pseudo_energy_2d(F, J);
  return W * kinematicUtils::tri3_area(refSize);
}

template<typename ELEMCOORDS>
std::array<double, 6>
SizeShapeObjective::compute_tri3_2d_element_pseudo_forces(const double refSize,
    const ELEMCOORDS & current_coords)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor2d P = sizeShape_element_pseudo_stress_2d(F, J);
  const std::array<double, 6> forces = kinematicUtils::tri3_2d_forces_from_stress(refSize, P);

  return forces;
}

std::array<double, 6>
SizeShapeObjective::compute_tri3_2d_element_pseudo_forces(const std::array<stk::math::Vector2d, 3> & ref,
    const std::array<stk::math::Vector2d, 3> & current_coords)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor2d P = sizeShape_element_pseudo_stress_2d(F, J);
  const std::array<double, 6> forces = kinematicUtils::tri3_2d_forces_from_stress(ref, P);

  return forces;
}

double SizeShapeObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const double refSize,
    const double * const * elemNodeCoords) const
{
  STK_ThrowAssert((spatialDim == 2 && (npe == 3 || npe == 6)) || (spatialDim == 3 && (npe == 4 || npe == 10)));
  if (spatialDim == 2)
    return compute_tri3_2d_element_pseudo_energy(refSize, elemNodeCoords);
  return compute_tet4_element_pseudo_energy(refSize, elemNodeCoords);
}

void SizeShapeObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const double refSize,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  STK_ThrowAssert((spatialDim == 2 && (npe == 3 || npe == 6)) || (spatialDim == 3 && (npe == 4 || npe == 10)));
  if (spatialDim == 2)
    kinematicUtils::fill_element_sensitivity(compute_tri3_2d_element_pseudo_forces(refSize, elemNodeCoords), elemGradContrib);
  else
    kinematicUtils::fill_element_sensitivity(compute_tet4_element_pseudo_forces(refSize, elemNodeCoords), elemGradContrib);
}

// explicit template instantiation
using CONSTPTRCONSTPTR = double const * const *;
template double SizeShapeObjective::compute_tet4_element_pseudo_energy(const double refSize, const std::array<stk::math::Vector3d, 4> & current_coords);
template double SizeShapeObjective::compute_tet4_element_pseudo_energy(const double refSize, const CONSTPTRCONSTPTR & current_coords);
template std::array<double, 12> SizeShapeObjective::compute_tet4_element_pseudo_forces(const double refSize, const std::array<stk::math::Vector3d, 4> & current_coords);
template std::array<double, 12> SizeShapeObjective::compute_tet4_element_pseudo_forces(const double refSize, const CONSTPTRCONSTPTR & current_coords);
template double SizeShapeObjective::compute_tri3_2d_element_pseudo_energy(const double refSize, const std::array<stk::math::Vector2d, 3> & current_coords);
template double SizeShapeObjective::compute_tri3_2d_element_pseudo_energy(const double refSize, const std::array<stk::math::Vector3d, 3> & current_coords);
template double SizeShapeObjective::compute_tri3_2d_element_pseudo_energy(const double refSize, const CONSTPTRCONSTPTR & current_coords);
template std::array<double, 6> SizeShapeObjective::compute_tri3_2d_element_pseudo_forces(const double refSize, const std::array<stk::math::Vector2d, 3> & current_coords);
template std::array<double, 6> SizeShapeObjective::compute_tri3_2d_element_pseudo_forces(const double refSize, const std::array<stk::math::Vector3d, 3> & current_coords);
template std::array<double, 6> SizeShapeObjective::compute_tri3_2d_element_pseudo_forces(const double refSize, const CONSTPTRCONSTPTR & current_coords);

}


