// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// Ported from Norma.jl (https://github.com/sandialabs/Norma.jl)
// by Alejandro Mota

#include <Akri_Objective_SethHill.hpp>
#include <Akri_AllReduce.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <stk_math/StkVector.hpp>
#include <Akri_KinematicUtils.hpp>
#include <MiniTensor.h>

namespace krino {

namespace {

using Tensor = minitensor::Tensor<double, 3>;
using Tensor2d = minitensor::Tensor<double, 2>;

// ==========================================================================
// MiniTensor helpers
// ==========================================================================

// Matrix power A^p for integer p.
// Uses binary_powering for |p| and inverse for negative p.
template<minitensor::Index N>
minitensor::Tensor<double, N> power(const minitensor::Tensor<double, N> & A, int p)
{
  if (p == 0)
    return minitensor::identity<double, N>();
  if (p == 1)
    return A;
  if (p == -1)
    return minitensor::inverse(A);
  if (p > 0)
    return minitensor::binary_powering(A, static_cast<minitensor::Index>(p));
  // p < -1: inverse of A^|p|
  return minitensor::inverse(
      minitensor::binary_powering(A, static_cast<minitensor::Index>(-p)));
}

// Matrix power A^p for integer p.
// Uses binary_powering for |p| and inverse for negative p.
template<minitensor::Index N>
std::pair<bool,minitensor::Tensor<double, N>> safe_power(const minitensor::Tensor<double, N> & A, int p)
{
  if (p == 0)
    return std::make_pair(true, minitensor::identity<double, N>());
  if (p == 1)
    return std::make_pair(true,A);
  if (p == -1)
  {
    return kinematicUtils::safe_inverse(A);
  }
  if (p > 0)
    return std::make_pair(true, minitensor::binary_powering(A, static_cast<minitensor::Index>(p)));
  // p < -1: inverse of A^|p|
  return kinematicUtils::safe_inverse(minitensor::binary_powering(A, static_cast<minitensor::Index>(-p)));
}

// ==========================================================================
// Seth-Hill constitutive model
// Ported from Norma.jl src/constitutive.jl lines 308-453
// ==========================================================================

// Strain energy W(F) — Norma.jl constitutive.jl:308-326
double seth_hill_energy(const Tensor & F, const double J, const SethHillParams & p)
{
  const Tensor C = minitensor::transpose(F) * F;
  const double Jm23 = kinematicUtils::pow_m23(J);
  const Tensor Cbar = Jm23 * C;

  const double Jm  = std::pow(J, p.m);
  const double Jnm = 1.0 / Jm;

  const auto [successCbar_n, Cbar_n] = safe_power(Cbar, p.n);
  const auto [successCbar_nn, Cbar_nn] = safe_power(Cbar, -p.n);
  if (!successCbar_n || !successCbar_nn)
    return std::numeric_limits<double>::max();

  const Tensor Cbar_2n  = Cbar_n * Cbar_n;
  const Tensor Cbar_n2n = Cbar_nn * Cbar_nn;

  const double trCbar_n   = minitensor::trace(Cbar_n);
  const double trCbar_nn  = minitensor::trace(Cbar_nn);
  const double trCbar_2n  = minitensor::trace(Cbar_2n);
  const double trCbar_n2n = minitensor::trace(Cbar_n2n);

  const double m2 = static_cast<double>(p.m * p.m);
  const double n2 = static_cast<double>(p.n * p.n);

  const double Wbulk  = p.kappa / (4.0 * m2) *
      ((Jm - 1.0) * (Jm - 1.0) + (Jnm - 1.0) * (Jnm - 1.0));
  const double Wshear = p.mu / (4.0 * n2) *
      (trCbar_2n + trCbar_n2n - 2.0 * trCbar_n - 2.0 * trCbar_nn + 6.0);

  return Wbulk + Wshear;
}

// 1st Piola-Kirchhoff stress P = dW/dF — Norma.jl constitutive.jl:424-453
Tensor seth_hill_stress(const Tensor & F, const double J, const SethHillParams & p)
{
  const Tensor C = minitensor::transpose(F) * F;
  const Tensor Finv = minitensor::inverse(F);
  const Tensor FinvT = minitensor::transpose(Finv);

  const double Jm   = std::pow(J, p.m);
  const double Jnm  = 1.0 / Jm;
  const double J2m  = Jm * Jm;
  const double Jn2m = 1.0 / J2m;

  const double Jm23 = kinematicUtils::pow_m23(J);
  const double Jp23 = kinematicUtils::pow_p23(J);
  const Tensor Cbar     = Jm23 * C;
  const Tensor Cbar_inv = Jp23 * (Finv * FinvT);

  const Tensor Cbar_n   = power(Cbar, p.n);
  const Tensor Cbar_nn  = power(Cbar_inv, p.n);
  const Tensor Cbar_2n  = Cbar_n * Cbar_n;
  const Tensor Cbar_n2n = Cbar_nn * Cbar_nn;

  const double trCbar_n   = minitensor::trace(Cbar_n);
  const double trCbar_nn  = minitensor::trace(Cbar_nn);
  const double trCbar_2n  = minitensor::trace(Cbar_2n);
  const double trCbar_n2n = minitensor::trace(Cbar_n2n);

  // Pbulk = kappa/(2m) * (J^2m - J^m - J^-2m + J^-m) * F^-T
  const double bulk_coeff = p.kappa / (2.0 * p.m) *
      (J2m - Jm - Jn2m + Jnm);
  const Tensor Pbulk = bulk_coeff * FinvT;

  // Pshear = mu/n * [1/3*(-trCbar2n + trCbarn + trCbar-2n - trCbar-n)*F^-T
  //                   + F^-T * (Cbar2n - Cbarn - Cbar-2n + Cbar-n)]
  const double iso_coeff = (1.0 / 3.0) *
      (-trCbar_2n + trCbar_n + trCbar_n2n - trCbar_nn);
  const Tensor dev_part = (Cbar_2n - Cbar_n) - (Cbar_n2n - Cbar_nn);
  const Tensor Pshear = (p.mu / p.n) *
      (iso_coeff * FinvT + FinvT * dev_part);

  return Pbulk + Pshear;
}

double seth_hill_energy_2d(const Tensor2d & F, const double J, const SethHillParams & p)
{
  const Tensor2d C = minitensor::transpose(F) * F;
  const double Jinv = 1.0 / J;
  const Tensor2d Cbar = Jinv * C;

  const double Jm  = std::pow(J, p.m);
  const double Jnm = 1.0 / Jm;

  const auto [successCbar_n,  Cbar_n]  = safe_power(Cbar,  p.n);
  const auto [successCbar_nn, Cbar_nn] = safe_power(Cbar, -p.n);
  if (!successCbar_n || !successCbar_nn)
    return std::numeric_limits<double>::max();

  const Tensor2d Cbar_2n  = Cbar_n * Cbar_n;
  const Tensor2d Cbar_n2n = Cbar_nn * Cbar_nn;

  const double trCbar_n   = minitensor::trace(Cbar_n);
  const double trCbar_nn  = minitensor::trace(Cbar_nn);
  const double trCbar_2n  = minitensor::trace(Cbar_2n);
  const double trCbar_n2n = minitensor::trace(Cbar_n2n);

  const double m2 = static_cast<double>(p.m * p.m);
  const double n2 = static_cast<double>(p.n * p.n);

  const double Wbulk  = p.kappa / (4.0 * m2) *
      ((Jm - 1.0) * (Jm - 1.0) + (Jnm - 1.0) * (Jnm - 1.0));

  const double Wshear = p.mu / (4.0 * n2) *
      (trCbar_2n + trCbar_n2n - 2.0 * trCbar_n - 2.0 * trCbar_nn + 4.0);

  return Wbulk + Wshear;
}

Tensor2d seth_hill_stress_2d(const Tensor2d & F, const double J, const SethHillParams & p)
{
  const Tensor2d C = minitensor::transpose(F) * F;
  const Tensor2d Finv = minitensor::inverse(F);
  const Tensor2d FinvT = minitensor::transpose(Finv);

  const double Jm   = std::pow(J, p.m);
  const double Jnm  = 1.0 / Jm;
  const double J2m  = Jm * Jm;
  const double Jn2m = 1.0 / J2m;

  const double Jinv = 1.0 / J;
  const Tensor2d Cbar     = Jinv * C;
  const Tensor2d Cbar_inv = J * (Finv * FinvT);

  const Tensor2d Cbar_n   = power(Cbar, p.n);
  const Tensor2d Cbar_nn  = power(Cbar_inv, p.n);
  const Tensor2d Cbar_2n  = Cbar_n * Cbar_n;
  const Tensor2d Cbar_n2n = Cbar_nn * Cbar_nn;

  const double trCbar_n   = minitensor::trace(Cbar_n);
  const double trCbar_nn  = minitensor::trace(Cbar_nn);
  const double trCbar_2n  = minitensor::trace(Cbar_2n);
  const double trCbar_n2n = minitensor::trace(Cbar_n2n);

  // Pbulk = kappa/(2m) * (J^2m - J^m - J^-2m + J^-m) * F^-T
  const double bulk_coeff = p.kappa / (2.0 * p.m) *
      (J2m - Jm - Jn2m + Jnm);
  const Tensor2d Pbulk = bulk_coeff * FinvT;

  // 2D:
  // Pshear = mu/n * [1/2*(-trCbar2n + trCbarn + trCbar-2n - trCbar-n)*F^-T
  //                  + F^-T * (Cbar2n - Cbarn - Cbar-2n + Cbar-n)]
  const double iso_coeff = 0.5 *
      (-trCbar_2n + trCbar_n + trCbar_n2n - trCbar_nn);
  const Tensor2d dev_part = (Cbar_2n - Cbar_n) - (Cbar_n2n - Cbar_nn);
  const Tensor2d Pshear = (p.mu / p.n) *
      (iso_coeff * FinvT + FinvT * dev_part);

  return Pbulk + Pshear;
}

}

template<typename ELEMCOORDS>
std::array<double, 12>
SethHillSmoothingObjective::compute_tet4_element_forces(
    const double refSize,
    const ELEMCOORDS & current_coords,
    const SethHillParams & params)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor P = seth_hill_stress(F, J, params);
  return kinematicUtils::tet4_forces_from_stress(refSize, P);
}

std::array<double, 12>
SethHillSmoothingObjective::compute_tet4_element_forces(
    const std::array<stk::math::Vector3d, 4> & ref,
    const std::array<stk::math::Vector3d, 4> & current_coords,
    const SethHillParams & params)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor P = seth_hill_stress(F, J, params);
  return kinematicUtils::tet4_forces_from_stress(ref, P);
}

template<typename ELEMCOORDS>
double SethHillSmoothingObjective::compute_tet4_element_energy(
    const double refSize,
    const ELEMCOORDS & current_coords,
    const SethHillParams & params)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = seth_hill_energy(F, J, params);
  if (W == std::numeric_limits<double>::max())
    return std::numeric_limits<double>::max();

  return W * kinematicUtils::tet4_volume(refSize);
}

double SethHillSmoothingObjective::compute_tet4_element_energy(
    const std::array<stk::math::Vector3d, 4> & ref,
    const std::array<stk::math::Vector3d, 4> & current_coords,
    const SethHillParams & params)
{
  const Tensor F = kinematicUtils::compute_deformation_gradient_tet4(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = seth_hill_energy(F, J, params);

  return W * kinematicUtils::tet4_volume(ref);
}

template<typename ELEMCOORDS>
double SethHillSmoothingObjective::compute_tri3_2d_element_energy(const double refSize,
    const ELEMCOORDS & current_coords,
    const SethHillParams & params)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return std::numeric_limits<double>::max();

  const double W = seth_hill_energy_2d(F, J, params);
  return W * kinematicUtils::tri3_area(refSize);
}

template<typename ELEMCOORDS>
std::array<double, 6>
SethHillSmoothingObjective::compute_tri3_2d_element_forces(const double refSize,
    const ELEMCOORDS & current_coords,
    const SethHillParams & params)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(refSize, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor2d P = seth_hill_stress_2d(F, J, params);
  const std::array<double, 6> forces = kinematicUtils::tri3_2d_forces_from_stress(refSize, P);

  return forces;
}

std::array<double, 6>
SethHillSmoothingObjective::compute_tri3_2d_element_forces(const std::array<stk::math::Vector2d, 3> & ref,
    const std::array<stk::math::Vector2d, 3> & current_coords,
    const SethHillParams & params)
{
  const Tensor2d F = kinematicUtils::compute_deformation_gradient_tri3_2d(ref, current_coords);
  const double J = minitensor::det(F);

  if (J <= 0.0 || !std::isfinite(J))
    return {};

  const Tensor2d P = seth_hill_stress_2d(F, J, params);
  const std::array<double, 6> forces = kinematicUtils::tri3_2d_forces_from_stress(ref, P);

  return forces;
}

double SethHillSmoothingObjective::compute_element_objective(
    const unsigned spatialDim,
    const unsigned npe,
    const double refSize,
    const double * const * elemNodeCoords) const
{
  STK_ThrowAssert((spatialDim == 2 && (npe == 3 || npe == 6)) || (spatialDim == 3 && (npe == 4 || npe == 10)));
  if (spatialDim == 2)
    return compute_tri3_2d_element_energy(refSize, elemNodeCoords, mySethHillParams);
  return compute_tet4_element_energy(refSize, elemNodeCoords, mySethHillParams);
}

void SethHillSmoothingObjective::fill_element_sensitivity(
    const unsigned spatialDim,
    const unsigned npe,
    const double refSize,
    const double * const * elemNodeCoords,
    std::vector<double> & elemGradContrib) const
{
  STK_ThrowAssert((spatialDim == 2 && (npe == 3 || npe == 6)) || (spatialDim == 3 && (npe == 4 || npe == 10)));
  if (spatialDim == 2)
    kinematicUtils::fill_element_sensitivity(compute_tri3_2d_element_forces(refSize, elemNodeCoords, mySethHillParams), elemGradContrib);
  else
    kinematicUtils::fill_element_sensitivity(compute_tet4_element_forces(refSize, elemNodeCoords, mySethHillParams), elemGradContrib);
}

// explicit template instantiation
using CONSTPTRCONSTPTR = double const * const *;
template double SethHillSmoothingObjective::compute_tet4_element_energy(const double refSize, const std::array<stk::math::Vector3d, 4> & current_coords, const SethHillParams & params);
template double SethHillSmoothingObjective::compute_tet4_element_energy(const double refSize, const CONSTPTRCONSTPTR & current_coords, const SethHillParams & params);
template std::array<double, 12> SethHillSmoothingObjective::compute_tet4_element_forces(const double refSize, const std::array<stk::math::Vector3d, 4> & current_coords, const SethHillParams & params);
template std::array<double, 12> SethHillSmoothingObjective::compute_tet4_element_forces(const double refSize, const CONSTPTRCONSTPTR & current_coords, const SethHillParams & params);
template double SethHillSmoothingObjective::compute_tri3_2d_element_energy(const double refSize, const std::array<stk::math::Vector2d, 3> & current_coords, const SethHillParams & params);
template double SethHillSmoothingObjective::compute_tri3_2d_element_energy(const double refSize, const std::array<stk::math::Vector3d, 3> & current_coords, const SethHillParams & params);
template double SethHillSmoothingObjective::compute_tri3_2d_element_energy(const double refSize, const CONSTPTRCONSTPTR & current_coords, const SethHillParams & params);
template std::array<double, 6> SethHillSmoothingObjective::compute_tri3_2d_element_forces(const double refSize, const std::array<stk::math::Vector2d, 3> & current_coords, const SethHillParams & params);
template std::array<double, 6> SethHillSmoothingObjective::compute_tri3_2d_element_forces(const double refSize, const std::array<stk::math::Vector3d, 3> & current_coords, const SethHillParams & params);
template std::array<double, 6> SethHillSmoothingObjective::compute_tri3_2d_element_forces(const double refSize, const CONSTPTRCONSTPTR & current_coords, const SethHillParams & params);

}


