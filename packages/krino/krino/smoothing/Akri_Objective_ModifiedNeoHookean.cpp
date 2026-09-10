// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <Akri_Objective_ModifiedNeoHookean.hpp>

#include <Akri_DistributedVector.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_KinematicUtils.hpp>
#include <MiniTensor.h>

using Tensor = minitensor::Tensor<double, 3>;
using Tensor2d = minitensor::Tensor<double, 2>;

namespace krino {

static double modifiedNeoHookean_element_energy(const Tensor & F, const double J)
{
  const Tensor C = minitensor::transpose(F) * F;
  const double Jm23 = kinematicUtils::pow_m23(J);
  const Tensor Cbar = Jm23 * C;

  const double WShape = minitensor::trace(Cbar) - 3.;

  const double Jn = 1.0 / J;
  const double WSize  = 1./4. * ((J - 1.0) * (J - 1.0) + (Jn - 1.0) * (Jn - 1.0));

  return WShape + WSize;
}

static double modifiedNeoHookean_element_energy_2d(const Tensor2d & F, const double J)
{
  const Tensor2d C = minitensor::transpose(F) * F;
  const double Jinv = 1./J;
  const Tensor2d Cbar = Jinv * C;

  const double WShape = minitensor::trace(Cbar) - 2.;

  const double Jn = 1.0 / J;
  const double WSize  = 1./4. * ((J - 1.0) * (J - 1.0) + (Jn - 1.0) * (Jn - 1.0));
  return WShape + WSize;
}

static Tensor modifiedNeoHookean_element_stress(const Tensor & F, const double J)
{
  const Tensor C = minitensor::transpose(F) * F;
  const Tensor FinvT = minitensor::transpose(minitensor::inverse(F));

  const double Jm23 = kinematicUtils::pow_m23(J);
  const Tensor Cbar = Jm23 * C;
  const double trCbar   = minitensor::trace(Cbar);

  const Tensor dWShape = 2. * Jm23 * F - 2./3. * trCbar * FinvT;

  const double Jn = 1.0 / J;
  const double J2  = J * J;
  const double Jn2 = 1.0 / J2;

  const Tensor dWSize = 1./2. * (J2 - J - Jn2 + Jn) * FinvT;

  return dWShape + dWSize;
}

static Tensor2d modifiedNeoHookean_element_stress_2d(const Tensor2d & F, const double J)
{
  const Tensor2d C = minitensor::transpose(F) * F;
  const Tensor2d FinvT = minitensor::transpose(minitensor::inverse(F));

  const double Jinv = 1./J;
  const Tensor2d Cbar = Jinv * C;
  const double trCbar = minitensor::trace(Cbar);

  const Tensor2d dWShape = 2.0 * Jinv * F - trCbar * FinvT;

  const double Jn = 1.0 / J;
  const double J2  = J * J;
  const double Jn2 = 1.0 / J2;

  const Tensor2d dWSize = 1./2. * (J2 - J - Jn2 + Jn) * FinvT;

  return dWShape + dWSize;
}

std::array<double, 12>
ModifiedNeoHookeanObjective::compute_tet4_element_forces(
    const std::array<stk::math::Vector3d, 4> & ref,
    const std::array<stk::math::Vector3d, 4> & current_coords)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor P = modifiedNeoHookean_element_stress(F, J);
  return kinematicUtils::tet4_forces_from_stress(ref, P);
}

template<typename ELEMCOORDS>
std::array<double, 12>
ModifiedNeoHookeanObjective::compute_tet4_element_forces(const double refSize,
    const ELEMCOORDS & current_coords)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor P = modifiedNeoHookean_element_stress(F, J);
  return kinematicUtils::tet4_forces_from_stress(refSize, P);
}

double ModifiedNeoHookeanObjective::compute_tet4_element_energy(
    const std::array<stk::math::Vector3d, 4> & ref,
    const std::array<stk::math::Vector3d, 4> & current_coords)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = modifiedNeoHookean_element_energy(F, J);
  return W * kinematicUtils::tet4_volume(ref);
}

template<typename ELEMCOORDS>
double ModifiedNeoHookeanObjective::compute_tet4_element_energy(const double refSize,
    const ELEMCOORDS & current_coords)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = modifiedNeoHookean_element_energy(F, J);
  return W * kinematicUtils::tet4_volume(refSize);
}

template<typename ELEMCOORDS>
double ModifiedNeoHookeanObjective::compute_tri3_2d_element_energy(const double refSize,
    const ELEMCOORDS & current_coords)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = modifiedNeoHookean_element_energy_2d(F, J);
  return W * kinematicUtils::tri3_area(refSize);
}

template<typename ELEMCOORDS>
std::array<double, 6>
ModifiedNeoHookeanObjective::compute_tri3_2d_element_forces(const double refSize,
    const ELEMCOORDS & current_coords)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor2d P = modifiedNeoHookean_element_stress_2d(F, J);
  const std::array<double, 6> forces = kinematicUtils::tri3_2d_forces_from_stress(refSize, P);

  return forces;
}

std::array<double, 6>
ModifiedNeoHookeanObjective::compute_tri3_2d_element_forces(const std::array<stk::math::Vector2d, 3> & ref,
    const std::array<stk::math::Vector2d, 3> & current_coords)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor2d P = modifiedNeoHookean_element_stress_2d(F, J);
  const std::array<double, 6> forces = kinematicUtils::tri3_2d_forces_from_stress(ref, P);

  return forces;
}

double ModifiedNeoHookeanObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const double refSize,
    const double * const * elemNodeCoords) const
{
  STK_ThrowAssert((spatialDim == 2 && (npe == 3 || npe == 6)) || (spatialDim == 3 && (npe == 4 || npe == 10)));
  if (spatialDim == 2)
    return compute_tri3_2d_element_energy(refSize, elemNodeCoords);
  return compute_tet4_element_energy(refSize, elemNodeCoords);
}

void ModifiedNeoHookeanObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const double refSize,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  STK_ThrowAssert((spatialDim == 2 && (npe == 3 || npe == 6)) || (spatialDim == 3 && (npe == 4 || npe == 10)));
  if (spatialDim == 2)
    kinematicUtils::fill_element_sensitivity(compute_tri3_2d_element_forces(refSize, elemNodeCoords), elemGradContrib);
  else
    kinematicUtils::fill_element_sensitivity(compute_tet4_element_forces(refSize, elemNodeCoords), elemGradContrib);
}

// explicit template instantiation
using CONSTPTRCONSTPTR = double const * const *;
template double ModifiedNeoHookeanObjective::compute_tet4_element_energy(const double refSize, const std::array<stk::math::Vector3d, 4> & current_coords);
template double ModifiedNeoHookeanObjective::compute_tet4_element_energy(const double refSize, const CONSTPTRCONSTPTR & current_coords);
template std::array<double, 12> ModifiedNeoHookeanObjective::compute_tet4_element_forces(const double refSize, const std::array<stk::math::Vector3d, 4> & current_coords);
template std::array<double, 12> ModifiedNeoHookeanObjective::compute_tet4_element_forces(const double refSize, const CONSTPTRCONSTPTR & current_coords);
template double ModifiedNeoHookeanObjective::compute_tri3_2d_element_energy(const double refSize, const std::array<stk::math::Vector2d, 3> & current_coords);
template double ModifiedNeoHookeanObjective::compute_tri3_2d_element_energy(const double refSize, const std::array<stk::math::Vector3d, 3> & current_coords);
template double ModifiedNeoHookeanObjective::compute_tri3_2d_element_energy(const double refSize, const CONSTPTRCONSTPTR & current_coords);
template std::array<double, 6> ModifiedNeoHookeanObjective::compute_tri3_2d_element_forces(const double refSize, const std::array<stk::math::Vector2d, 3> & current_coords);
template std::array<double, 6> ModifiedNeoHookeanObjective::compute_tri3_2d_element_forces(const double refSize, const std::array<stk::math::Vector3d, 3> & current_coords);
template std::array<double, 6> ModifiedNeoHookeanObjective::compute_tri3_2d_element_forces(const double refSize, const CONSTPTRCONSTPTR & current_coords);

}
