// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_QualityMetricWithSensitivities.hpp>
#include <Akri_QualityMetric.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include "Sacado.hpp"

namespace krino
{

using FAD_Type3 = Sacado::Fad::SFad<double,3>;
using Vector3FAD = stk::math::Vec<FAD_Type3,3>;

using FAD_Type12 = Sacado::Fad::SFad<double,12>;
using Vector3FAD12 = stk::math::Vec<FAD_Type12,3>;

using FAD_Type6 = Sacado::Fad::SFad<double,6>;
using Vector2FAD6 = stk::math::Vec<FAD_Type6,2>;

static FAD_Type3 Dot(const Vector3FAD & a, const stk::math::Vector3d& b)
{
  FAD_Type3 dot = 0.0;
  for (unsigned i=0; i<3; ++i) dot += a[i]*b[i];
  return dot;
}

static FAD_Type3 tet_scaled_jacobian_and_sensitivity_to_last_vertex_using_FAD(const stk::math::Vector3d & v0, const stk::math::Vector3d & v1, const stk::math::Vector3d & v2, const stk::math::Vector3d & v3)
{
  const stk::math::Vector3d edge0 = v1 - v0;
  const stk::math::Vector3d edge1 = v2 - v1;
  const stk::math::Vector3d edge2 = v0 - v2;
  const Vector3FAD edge3(FAD_Type3(3, 0, v3[0]-v0[0]), FAD_Type3(3, 1, v3[1]-v0[1]), FAD_Type3(3, 2, v3[2]-v0[2]));
  const Vector3FAD edge4(FAD_Type3(3, 0, v3[0]-v1[0]), FAD_Type3(3, 1, v3[1]-v1[1]), FAD_Type3(3, 2, v3[2]-v1[2]));
  const Vector3FAD edge5(FAD_Type3(3, 0, v3[0]-v2[0]), FAD_Type3(3, 1, v3[1]-v2[1]), FAD_Type3(3, 2, v3[2]-v2[2]));

  const FAD_Type3 jacobi = Dot(edge3, Cross(edge2, edge0));

  // products of lengths squared of each edge attached to a node.
  const double edge0_length_squared = edge0.length_squared();
  const double edge1_length_squared = edge1.length_squared();
  const double edge2_length_squared = edge2.length_squared();
  const FAD_Type3 edge3_length_squared = edge3.length_squared();
  const FAD_Type3 edge4_length_squared = edge4.length_squared();
  const FAD_Type3 edge5_length_squared = edge5.length_squared();

  const std::array<FAD_Type3,4> length_squared = {{
    edge0_length_squared * edge2_length_squared * edge3_length_squared,
    edge0_length_squared * edge1_length_squared * edge4_length_squared,
    edge1_length_squared * edge2_length_squared * edge5_length_squared,
    edge3_length_squared * edge4_length_squared * edge5_length_squared
  }};
  int which_node = 0;
  if(length_squared[1] > length_squared[which_node])
    which_node = 1;
  if(length_squared[2] > length_squared[which_node])
    which_node = 2;
  if(length_squared[3] > length_squared[which_node])
    which_node = 3;

  FAD_Type3 length_product = std::sqrt( length_squared[which_node] );
  if(length_product < std::abs(jacobi))
    length_product = std::abs(jacobi);

  FAD_Type3 result = 0.0;
  const double lengthMin = 1.0E-30;
  if( length_product > lengthMin )
  {
    static const double root_of_2 = std::sqrt(2.0);
    result = root_of_2 * jacobi / length_product;
  }

  return result;
}

static FAD_Type3 tet_mean_ratio_and_sensitivity_to_last_vertex_using_FAD(const stk::math::Vector3d & v0, const stk::math::Vector3d & v1, const stk::math::Vector3d & v2, const stk::math::Vector3d & v3)
{
  const stk::math::Vector3d edge0 = v1 - v0;
  const stk::math::Vector3d edge1 = v2 - v1;
  const stk::math::Vector3d edge2 = v0 - v2;
  const Vector3FAD edge3(FAD_Type3(3, 0, v3[0]-v0[0]), FAD_Type3(3, 1, v3[1]-v0[1]), FAD_Type3(3, 2, v3[2]-v0[2]));
  const Vector3FAD edge4(FAD_Type3(3, 0, v3[0]-v1[0]), FAD_Type3(3, 1, v3[1]-v1[1]), FAD_Type3(3, 2, v3[2]-v1[2]));
  const Vector3FAD edge5(FAD_Type3(3, 0, v3[0]-v2[0]), FAD_Type3(3, 1, v3[1]-v2[1]), FAD_Type3(3, 2, v3[2]-v2[2]));

  static const double sqrt2 = std::sqrt(2.);
  const FAD_Type3 tetVolumeTimes6 = Dot(edge3, Cross(edge2, edge0));
  const FAD_Type3 sum = (edge0.length_squared() + edge1.length_squared() + edge2.length_squared() +
      edge3.length_squared() + edge4.length_squared() + edge5.length_squared()) / 6;
  return sqrt2 * tetVolumeTimes6 / (sum*std::sqrt(sum));
}

static FAD_Type12 tet_mean_ratio_and_sensitivities_using_FAD(const stk::math::Vector3d & x0, const stk::math::Vector3d & x1, const stk::math::Vector3d & x2, const stk::math::Vector3d & x3)
{
  const Vector3FAD12 v0(FAD_Type12(12, 0, x0[0]), FAD_Type12(12, 1, x0[1]), FAD_Type12(12, 2, x0[2]));
  const Vector3FAD12 v1(FAD_Type12(12, 3, x1[0]), FAD_Type12(12, 4, x1[1]), FAD_Type12(12, 5, x1[2]));
  const Vector3FAD12 v2(FAD_Type12(12, 6, x2[0]), FAD_Type12(12, 7, x2[1]), FAD_Type12(12, 8, x2[2]));
  const Vector3FAD12 v3(FAD_Type12(12, 9, x3[0]), FAD_Type12(12,10, x3[1]), FAD_Type12(12,11, x3[2]));
  const Vector3FAD12 edge0 = v1 - v0;
  const Vector3FAD12 edge1 = v2 - v1;
  const Vector3FAD12 edge2 = v0 - v2;
  const Vector3FAD12 edge3 = v3 - v0;
  const Vector3FAD12 edge4 = v3 - v1;
  const Vector3FAD12 edge5 = v3 - v2;

  static const double sqrt2 = std::sqrt(2.);
  const FAD_Type12 tetVolumeTimes6 = Dot(edge3, Cross(edge2, edge0));
  const FAD_Type12 sum = (edge0.length_squared() + edge1.length_squared() + edge2.length_squared() +
      edge3.length_squared() + edge4.length_squared() + edge5.length_squared()) / 6;

  return sqrt2 * tetVolumeTimes6 / (sum*std::sqrt(sum));
}

static FAD_Type6 tri2d_mean_ratio_and_sensitivities_using_FAD(const double * x0, const double * x1, const double * x2)
{
  const Vector2FAD6 v0(FAD_Type6(6, 0, x0[0]), FAD_Type6(6, 1, x0[1]));
  const Vector2FAD6 v1(FAD_Type6(6, 2, x1[0]), FAD_Type6(6, 3, x1[1]));
  const Vector2FAD6 v2(FAD_Type6(6, 4, x2[0]), FAD_Type6(6, 5, x2[1]));
  const Vector2FAD6 edge0 = v1 - v0;
  const Vector2FAD6 edge1 = v2 - v1;
  const Vector2FAD6 edge2 = v0 - v2;

  static const double coeff = 2./std::sqrt(3.);
  const FAD_Type6 triAreaTimes2 = (edge2[0]*edge0[1]-edge0[0]*edge2[1]);
  const FAD_Type6 sumSq = (edge0.length_squared() + edge1.length_squared() + edge2.length_squared()) / 3;
  return coeff * triAreaTimes2 / sumSq;
}

static double tet_mean_ratio_and_sensitivities(
    const stk::math::Vector3d & v0,
    const stk::math::Vector3d & v1,
    const stk::math::Vector3d & v2,
    const stk::math::Vector3d & v3,
    stk::math::Vector3d & grad0,
    stk::math::Vector3d & grad1,
    stk::math::Vector3d & grad2,
    stk::math::Vector3d & grad3)
{
  static const double sqrt2 = std::sqrt(2.0);
  static constexpr double oneThird = 1.0 / 3.0;

  // Relative-to-v0 edge vectors
  const stk::math::Vector3d a = v1 - v0;
  const stk::math::Vector3d b = v2 - v0;
  const stk::math::Vector3d c = v3 - v0;

  const stk::math::Vector3d axb = Cross(a,b);

  // T = 6 * signed volume
  const double T = Dot(c, axb);

  const double aa = a.length_squared();
  const double bb = b.length_squared();
  const double cc = c.length_squared();
  const double ab = Dot(a, b);
  const double ac = Dot(a, c);
  const double bc = Dot(b, c);

  // s = (1/6) * sum of squared edge lengths
  // Use ||u-v||^2 = ||u||^2 + ||v||^2 - 2 u.v
  const double s = 0.5 * (aa + bb + cc) - oneThird * (ab + ac + bc);

  if (s <= 0.0)
  {
    grad0 = stk::math::Vector3d::ZERO;
    grad1 = stk::math::Vector3d::ZERO;
    grad2 = stk::math::Vector3d::ZERO;
    grad3 = stk::math::Vector3d::ZERO;
    return 0.0;
  }

  const double sqrtS = std::sqrt(s);
  const double invS = 1.0 / s;
  const double invS32 = 1.0 / (s * sqrtS);
  const double value = sqrt2 * T * invS32;

  // grad(T)
  const stk::math::Vector3d dTdv1 = Cross(b, c);
  const stk::math::Vector3d dTdv2 = Cross(c, a);
  const stk::math::Vector3d dTdv3 = axb;
  const stk::math::Vector3d dTdv0 = -(dTdv1 + dTdv2 + dTdv3);

  // grad(s)
  const stk::math::Vector3d dsdv1 = a - oneThird * (b + c);
  const stk::math::Vector3d dsdv2 = b - oneThird * (a + c);
  const stk::math::Vector3d dsdv3 = c - oneThird * (a + b);
  const stk::math::Vector3d dsdv0 = -(dsdv1 + dsdv2 + dsdv3);

  const double factor = sqrt2 * invS32;
  const double alpha = 1.5 * T * invS;

  grad0 = factor * (dTdv0 - alpha * dsdv0);
  grad1 = factor * (dTdv1 - alpha * dsdv1);
  grad2 = factor * (dTdv2 - alpha * dsdv2);
  grad3 = factor * (dTdv3 - alpha * dsdv3);

  return value;
}

static double tri3_2d_mean_ratio_and_sensitivities(
    const double * v0,
    const double * v1,
    const double * v2,
    double * grad)
{
  static const double coeff = 2.0 / std::sqrt(3.0);
  static constexpr double twoThirds = 2.0 / 3.0;

  const stk::math::Vector2d a(v1[0] - v0[0], v1[1] - v0[1]);
  const stk::math::Vector2d b(v2[0] - v0[0], v2[1] - v0[1]);

  const double area2 = a[0]*b[1] - a[1]*b[0];

  const double aa = a[0]*a[0] + a[1]*a[1];
  const double bb = b[0]*b[0] + b[1]*b[1];
  const double ab = a[0]*b[0] + a[1]*b[1];
  const double s = twoThirds * (aa + bb - ab);

  if (s <= 0.0) {
    for (int i=0; i<6; ++i)
      grad[i] = 0.;
    return 0.0;
  }

  const double value = coeff * area2 / s;

  const stk::math::Vector2d dAdv1( b[1], -b[0]);
  const stk::math::Vector2d dAdv2(-a[1],  a[0]);
  const stk::math::Vector2d dAdv0 = -(dAdv1 + dAdv2);

  const stk::math::Vector2d dsdv1 = twoThirds * (2.0*a - b);
  const stk::math::Vector2d dsdv2 = twoThirds * (2.0*b - a);
  const stk::math::Vector2d dsdv0 = -(dsdv1 + dsdv2);

  const double factor = coeff / s;
  const double alpha = area2 / s;

  const stk::math::Vector2d g0 = factor * (dAdv0 - alpha * dsdv0);
  const stk::math::Vector2d g1 = factor * (dAdv1 - alpha * dsdv1);
  const stk::math::Vector2d g2 = factor * (dAdv2 - alpha * dsdv2);

  for (int i=0; i<2; ++i)
  {
    grad[2*0+i] = g0[i];
    grad[2*1+i] = g1[i];
    grad[2*2+i] = g2[i];
  }

  return value;
}

static std::tuple<double, std::array<double,3>> tet_scaled_jacobian_and_sensitivity_to_last_vertex(const stk::math::Vector3d & v0, const stk::math::Vector3d & v1, const stk::math::Vector3d & v2, const stk::math::Vector3d & v3)
{
  FAD_Type3 sj = tet_scaled_jacobian_and_sensitivity_to_last_vertex_using_FAD(v0, v1, v2, v3);
  return std::make_tuple(sj.val(), std::array<double,3>{sj.dx(0), sj.dx(1), sj.dx(2)});
}

static std::tuple<double, std::array<double,3>> tet_mean_ratio_and_sensitivity_to_last_vertex(const stk::math::Vector3d & v0, const stk::math::Vector3d & v1, const stk::math::Vector3d & v2, const stk::math::Vector3d & v3)
{
  FAD_Type3 mr = tet_mean_ratio_and_sensitivity_to_last_vertex_using_FAD(v0, v1, v2, v3);
  return std::make_tuple(mr.val(), std::array<double,3>{mr.dx(0), mr.dx(1), mr.dx(2)});
}

static std::tuple<double, std::array<double,3>> tet_volume_and_sensitivity_to_last_vertex(const stk::math::Vector3d & v0, const stk::math::Vector3d & v1, const stk::math::Vector3d & v2, const stk::math::Vector3d & v3)
{
  const stk::math::Vector3d edge0 = v1 - v0;
  const stk::math::Vector3d edge2 = v0 - v2;
  const stk::math::Vector3d edge3 = v3 - v0;

  static constexpr double oneSixth = 1./6.;
  const stk::math::Vector3d dvol_dv3 = oneSixth*Cross(edge2, edge0);
  const double vol = Dot(edge3, dvol_dv3);
  const std::array<double,3> sens = {dvol_dv3[0], dvol_dv3[1], dvol_dv3[2]};

  return {vol, sens};
}

static const std::array<int,4> & get_tet_permutation_with_nth_vertex_last(const int n)
{
  static std::array< std::array<int,4>, 4> perms {{
    {{1,3,2,0}},
    {{2,3,0,1}},
    {{0,3,1,2}},
    {{0,1,2,3}}
  }};
  return perms[n];
}

std::tuple<double, std::array<double,3>> ScaledJacobianQualityMetricWithSensitivities::tet_quality_and_sensitivity_to_nth_vertex(const std::array<stk::math::Vector3d, 4> & verts, const int n)
{
  const std::array<int,4> & perm = get_tet_permutation_with_nth_vertex_last(n);
  return tet_scaled_jacobian_and_sensitivity_to_last_vertex(verts[perm[0]], verts[perm[1]], verts[perm[2]], verts[perm[3]]);
}

std::tuple<double, std::array<double,3>> MeanRatioQualityMetricWithSensitivities::tet_quality_and_sensitivity_to_nth_vertex(const std::array<stk::math::Vector3d, 4> & verts, const int n)
{
  const std::array<int,4> & perm = get_tet_permutation_with_nth_vertex_last(n);
  return tet_mean_ratio_and_sensitivity_to_last_vertex(verts[perm[0]], verts[perm[1]], verts[perm[2]], verts[perm[3]]);
}

std::tuple<double, std::array<double,3>> tet_volume_and_sensitivity_to_nth_vertex(const std::array<stk::math::Vector3d, 4> & verts, const int n)
{
  const std::array<int,4> & perm = get_tet_permutation_with_nth_vertex_last(n);
  return tet_volume_and_sensitivity_to_last_vertex(verts[perm[0]], verts[perm[1]], verts[perm[2]], verts[perm[3]]);
}

std::tuple<double, std::array<stk::math::Vector3d, 4>> MeanRatioQualityMetricWithSensitivities::tet_quality_and_sensitivities_using_FAD(const std::array<stk::math::Vector3d, 4> & verts)
{
  const FAD_Type12 mr = tet_mean_ratio_and_sensitivities_using_FAD(verts[0], verts[1], verts[2], verts[3]);
  const std::array<stk::math::Vector3d, 4> sens = { stk::math::Vector3d(mr.dx()), stk::math::Vector3d(mr.dx()+3), stk::math::Vector3d(mr.dx()+6), stk::math::Vector3d(mr.dx()+9) };
  return std::make_tuple(mr.val(), sens);
}

std::tuple<double, std::array<stk::math::Vector3d, 4>> MeanRatioQualityMetricWithSensitivities::tet_quality_and_sensitivities(const std::array<stk::math::Vector3d, 4> & verts)
{
  std::array<stk::math::Vector3d, 4> sens;
  const double mr = tet_mean_ratio_and_sensitivities(verts[0], verts[1], verts[2], verts[3], sens[0], sens[1], sens[2], sens[3]);
  return std::make_tuple(mr, sens);
}

std::tuple<double, std::array<double, 6>> MeanRatioQualityMetricWithSensitivities::tri2d_quality_and_sensitivities_using_FAD(const double * v0, const double * v1, const double * v2)
{
  const FAD_Type6 mr = tri2d_mean_ratio_and_sensitivities_using_FAD(v0, v1, v2);
  std::array<double, 6> sens;
  for (unsigned i=0; i<6; ++i)
    sens[i] = mr.fastAccessDx(i);
  return std::make_tuple(mr.val(), sens);
}

std::tuple<double, std::array<double, 6>> MeanRatioQualityMetricWithSensitivities::tri2d_quality_and_sensitivities(const double * v0, const double * v1, const double * v2)
{
  std::array<double, 6> sens;
  const double mr = tri3_2d_mean_ratio_and_sensitivities(v0, v1, v2, sens.data());
  return std::make_tuple(mr, sens);
}

}
