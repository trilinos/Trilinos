// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_KinematicUtils.hpp>

#include <Akri_SethHillConfig.hpp>
#include <Akri_AllReduce.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_Smoothing_Utils.hpp>
#include <stk_math/StkVector.hpp>
#include <MiniTensor.h>

namespace krino {
namespace kinematicUtils {

using Tensor = minitensor::Tensor<double, 3>;
using Tensor2d = minitensor::Tensor<double, 2>;

// ==========================================================================
// Ideal tet4 reference configuration
// Ported from Norma.jl src/model.jl lines 170-206
// ==========================================================================

double equal_volume_tet_h(const stk::math::Vector3d & u,
    const stk::math::Vector3d & v, const stk::math::Vector3d & w)
{
  // No abs — matches Norma.jl exactly (allows inverted element handling)
  return std::cbrt(std::sqrt(2.0) * Dot(u, Cross(v, w)));
}

double avg_edge_length_tet_h(const stk::math::Vector3d & u,
    const stk::math::Vector3d & v, const stk::math::Vector3d & w)
{
  const double lu = std::sqrt(Dot(u, u));
  const double lv = std::sqrt(Dot(v, v));
  const double lw = std::sqrt(Dot(w, w));
  const stk::math::Vector3d uv = u - v;
  const stk::math::Vector3d uw = u - w;
  const stk::math::Vector3d vw = v - w;
  const double luv = std::sqrt(Dot(uv, uv));
  const double luw = std::sqrt(Dot(uw, uw));
  const double lvw = std::sqrt(Dot(vw, vw));
  return (lu + lv + lw + luv + luw + lvw) / 6.0;
}

template<typename ELEMCOORDS>
double compute_ideal_tet4_h(const ELEMCOORDS & coords,
                            SmoothReference mode)
{
  const stk::math::Vector3d u(coords[1][0]-coords[0][0], coords[1][1]-coords[0][1], coords[1][2]-coords[0][2]);
  const stk::math::Vector3d v(coords[2][0]-coords[0][0], coords[2][1]-coords[0][1], coords[2][2]-coords[0][2]);
  const stk::math::Vector3d w(coords[3][0]-coords[0][0], coords[3][1]-coords[0][1], coords[3][2]-coords[0][2]);

  double h = 0.0;
  switch (mode)
  {
    case SmoothReference::EqualVolume:
      h = equal_volume_tet_h(u, v, w);
      break;
    case SmoothReference::AvgEdgeLength:
      h = avg_edge_length_tet_h(u, v, w);
      break;
    case SmoothReference::Max:
      h = std::max(equal_volume_tet_h(u, v, w), avg_edge_length_tet_h(u, v, w));
      break;
  }
  return h;
}

std::array<stk::math::Vector3d, 4> compute_ideal_tet4_reference(
    const double h)
{
  const double c = h * 0.5 / std::sqrt(2.0);
  return {{
    stk::math::Vector3d( c,  c,  c),
    stk::math::Vector3d(-c, -c,  c),
    stk::math::Vector3d(-c,  c, -c),
    stk::math::Vector3d( c, -c, -c)
  }};
}

std::array<stk::math::Vector2d, 3> compute_ideal_tri3_2d_reference(
    const double h)
{
  const double cx = h*0.5;
  const double cy = h*0.5*std::sqrt(3);
  return {{
    stk::math::Vector2d( 0, cy),
    stk::math::Vector2d(-cx, 0),
    stk::math::Vector2d( cx, 0)
  }};
}

std::array<stk::math::Vector3d, 4> compute_ideal_tet4_reference(
    const std::array<stk::math::Vector3d, 4> & coords,
    SmoothReference mode)
{
  const double h = compute_ideal_tet4_h(coords, mode);
  return compute_ideal_tet4_reference(h);
}

// ==========================================================================
// Element-level energy and force computation
// Ported from Norma.jl src/model.jl lines 649-676
// ==========================================================================

Tensor compute_deformation_gradient_tet4(
    const std::array<stk::math::Vector3d, 4> & ref,
    const std::array<stk::math::Vector3d, 4> & cur)
{
  // dXdxi(i, j) = X_j(node i+1) - X_j(node 0) = ∂X_j/∂ξ_i  (transposed Jacobian)
  // dxdxi(i, j) = ∂x_j/∂ξ_i
  // (inverse(dXdxi))(i, j) = ∂ξ_j/∂X_i
  // So (inverse(dXdxi) * dxdxi)(i, j) = ∂x_j/∂X_i — the transpose of F.
  // Take one more transpose to return the standard convention F(i, j) = ∂x_i/∂X_j.
  Tensor dXdxi, dxdxi;
  for (int j = 0; j < 3; ++j)
  {
    dXdxi(0, j) = ref[1][j] - ref[0][j];
    dXdxi(1, j) = ref[2][j] - ref[0][j];
    dXdxi(2, j) = ref[3][j] - ref[0][j];

    dxdxi(0, j) = cur[1][j] - cur[0][j];
    dxdxi(1, j) = cur[2][j] - cur[0][j];
    dxdxi(2, j) = cur[3][j] - cur[0][j];
  }

  return minitensor::transpose(minitensor::inverse(dXdxi) * dxdxi);
}

Tensor2d compute_deformation_gradient_tri3_2d(
    const std::array<stk::math::Vector2d, 3> & ref,
    const std::array<stk::math::Vector2d, 3> & cur)
{
  // dXdxi(i, j) = X_j(node i+1) - X_j(node 0) = ∂X_j/∂ξ_i  (transposed Jacobian)
  // dxdxi(i, j) = ∂x_j/∂ξ_i
  // (inverse(dXdxi))(i, j) = ∂ξ_j/∂X_i
  // So (inverse(dXdxi) * dxdxi)(i, j) = ∂x_j/∂X_i — the transpose of F.
  // Take one more transpose to return the standard convention F(i, j) = ∂x_i/∂X_j.
  Tensor2d dXdxi, dxdxi;
  for (int j = 0; j < 2; ++j)
  {
    dXdxi(0, j) = ref[1][j] - ref[0][j];
    dXdxi(1, j) = ref[2][j] - ref[0][j];

    dxdxi(0, j) = cur[1][j] - cur[0][j];
    dxdxi(1, j) = cur[2][j] - cur[0][j];
  }

  return minitensor::transpose(minitensor::inverse(dXdxi) * dxdxi);
}

template<typename ELEMCOORDS>
Tensor2d compute_deformation_gradient_tri3_2d(
    const double refSize,
    const ELEMCOORDS & cur)
{
  const double invRefSize = 1./refSize;
  const double c = invRefSize/std::sqrt(3.0);
  Tensor2d F;
  for (int i = 0; i < 2; ++i)
  {
    F(i,0) = invRefSize *(cur[2][i] - cur[1][i]);
    F(i,1) =  c * (2.0 * cur[0][i] - cur[1][i] - cur[2][i]);
  }

  return F;
}

template<typename ELEMCOORDS>
Tensor compute_deformation_gradient_tet4(
    const double refSize,
    const ELEMCOORDS & cur)
{
  const double c = 0.5 *std::sqrt(2.0) / refSize;
  Tensor F;
  for (int i = 0; i < 3; ++i)
  {
    F(i,0) = c * (cur[0][i] - cur[1][i] - cur[2][i] + cur[3][i]);
    F(i,1) = c * (cur[0][i] - cur[1][i] + cur[2][i] - cur[3][i]);
    F(i,2) = c * (cur[0][i] + cur[1][i] - cur[2][i] - cur[3][i]);
  }

  return F;
}

double tri3_area(const double size)
{
  return std::sqrt(3.0)/4. * size*size;
}

double tet4_volume(const double size)
{
  return std::sqrt(2.0)/12. * size*size*size;
}

double tet4_volume(const std::array<stk::math::Vector3d, 4> & ref)
{
  // Build dXdxi and its inverse for force assembly
  Tensor dXdxi;
  for (int j = 0; j < 3; ++j)
  {
    dXdxi(0, j) = ref[1][j] - ref[0][j];
    dXdxi(1, j) = ref[2][j] - ref[0][j];
    dXdxi(2, j) = ref[3][j] - ref[0][j];
  }
  const double tet_reference_domain_volume = 1.0 / 6.0;
  const double Vref = minitensor::det(dXdxi) * tet_reference_domain_volume;

  return Vref;
}

std::array<stk::math::Vector3d, 4>
tet4_forces_from_stress(const std::array<stk::math::Vector3d, 4> & ref,
    const Tensor & P)
{
  // Build dXdxi and its inverse for force assembly
  Tensor dXdxi;
  for (int j = 0; j < 3; ++j)
  {
    dXdxi(0, j) = ref[1][j] - ref[0][j];
    dXdxi(1, j) = ref[2][j] - ref[0][j];
    dXdxi(2, j) = ref[3][j] - ref[0][j];
  }
  const double tet_reference_domain_volume = 1.0 / 6.0;
  const double Vref = minitensor::det(dXdxi) * tet_reference_domain_volume;
  const Tensor dXdxi_inv = minitensor::inverse(dXdxi);

  // dN_a/dX_k  (a = node index, k = reference spatial direction)
  std::array<std::array<double, 3>, 4> dNdX;
  for (int k = 0; k < 3; ++k)
  {
    dNdX[0][k] = -(dXdxi_inv(k, 0) + dXdxi_inv(k, 1) + dXdxi_inv(k, 2));
    dNdX[1][k] = dXdxi_inv(k, 0);
    dNdX[2][k] = dXdxi_inv(k, 1);
    dNdX[3][k] = dXdxi_inv(k, 2);
  }

  // Weak form: f_a[j] = sum_k P(j, k) * dNdX[a][k] * Vref
  std::array<stk::math::Vector3d, 4> forces;
  for (int a = 0; a < 4; ++a)
    for (int j = 0; j < 3; ++j)
    {
      double val = 0.0;
      for (int k = 0; k < 3; ++k)
        val += P(j, k) * dNdX[a][k];
      forces[a][j] = val * Vref;
    }

  return forces;
}

std::array<stk::math::Vector3d, 4>
tet4_forces_from_stress(const double refSize, const Tensor & P)
{
  const double c = 0.5 *std::sqrt(2.0) / refSize;
  const double cVref = c*tet4_volume(refSize);

  // dN_a/dX_k  (a = node index, k = reference spatial direction)
  // Weak form: f_a[j] = sum_k P(j, k) * dNdX[a][k] * Vref
  std::array<stk::math::Vector3d, 4> forces;
  for (int j = 0; j < 3; ++j)
  {
    forces[0][j] = cVref*( P(j,0)+P(j,1)+P(j,2));
    forces[1][j] = cVref*(-P(j,0)-P(j,1)+P(j,2));
    forces[2][j] = cVref*(-P(j,0)+P(j,1)-P(j,2));
    forces[3][j] = cVref*( P(j,0)-P(j,1)-P(j,2));
  }

  return forces;
}

std::array<double, 6>
tri3_2d_forces_from_stress(const std::array<stk::math::Vector2d, 3> & ref,
    const Tensor2d & P)
{
  // Build dXdxi and its inverse for force assembly
  Tensor2d dXdxi;
  for (int j = 0; j < 2; ++j)
  {
    dXdxi(0, j) = ref[1][j] - ref[0][j];
    dXdxi(1, j) = ref[2][j] - ref[0][j];
  }
  const double tri_reference_domain_volume = 0.5;
  const double Aref = minitensor::det(dXdxi) * tri_reference_domain_volume;
  const Tensor2d dXdxi_inv = minitensor::inverse(dXdxi);

  // dN_a/dX_k  (a = node index, k = reference spatial direction)
  std::array<std::array<double, 2>, 3> dNdX;
  for (int k = 0; k < 2; ++k)
  {
    dNdX[0][k] = -(dXdxi_inv(k, 0) + dXdxi_inv(k, 1));
    dNdX[1][k] = dXdxi_inv(k, 0);
    dNdX[2][k] = dXdxi_inv(k, 1);
  }

  // Weak form: f_a[j] = sum_k P(j, k) * dNdX[a][k] * Vref
  std::array<double, 6> forces;
  for (int a = 0; a < 3; ++a)
    for (int j = 0; j < 2; ++j)
    {
      double val = 0.0;
      for (int k = 0; k < 2; ++k)
        val += P(j, k) * dNdX[a][k];
      forces[a*2+j] = val * Aref;
    }

  return forces;
}

std::array<double, 6>
tri3_2d_forces_from_stress(const double refSize, const Tensor2d& P)
{
  const double c = 0.25 * refSize;
  const double s = std::sqrt(3.0) * c;

  std::array<double, 6> forces;
  for (int j = 0; j < 2; ++j)
  {
    forces[0*2+j] = 2.0*c * P(j,1);
    forces[1*2+j] = -s*P(j,0)-c*P(j,1);
    forces[2*2+j] =  s*P(j,0)-c*P(j,1);
  }

  return forces;
}

// explicit template instantiation
using CONSTPTRCONSTPTR = double const * const *;
template double compute_ideal_tet4_h(const std::array<stk::math::Vector3d, 4> & coords, SmoothReference mode);
template double compute_ideal_tet4_h(const CONSTPTRCONSTPTR & coords, SmoothReference mode);
template Tensor compute_deformation_gradient_tet4(const double refSize, const std::array<stk::math::Vector3d, 4> & cur);
template Tensor compute_deformation_gradient_tet4(const double refSize, const CONSTPTRCONSTPTR & cur);
template Tensor2d compute_deformation_gradient_tri3_2d(const double refSize, const std::array<stk::math::Vector3d, 3> & cur);
template Tensor2d compute_deformation_gradient_tri3_2d(const double refSize, const std::array<stk::math::Vector2d, 3> & cur);
template Tensor2d compute_deformation_gradient_tri3_2d(const double refSize, const CONSTPTRCONSTPTR & cur);

}
}

