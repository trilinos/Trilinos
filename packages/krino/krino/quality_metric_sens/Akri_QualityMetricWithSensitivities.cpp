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
  const stk::math::Vector3d edge2 = v0 - v2;
  const Vector3FAD edge3(FAD_Type3(3, 0, v3[0]-v0[0]), FAD_Type3(3, 1, v3[1]-v0[1]), FAD_Type3(3, 2, v3[2]-v0[2]));

  const double volumeMin = 1.0E-30;
  static constexpr double oneSixth = 1./6.;
  const FAD_Type3 tetVolume = oneSixth*Dot(edge3, Cross(edge2, edge0));
  if( std::abs( tetVolume.val() ) < volumeMin )
      return 0.0;

  const stk::math::Vector3d edge1 = v2 - v1;
  const Vector3FAD edge4(FAD_Type3(3, 0, v3[0]-v1[0]), FAD_Type3(3, 1, v3[1]-v1[1]), FAD_Type3(3, 2, v3[2]-v1[2]));
  const Vector3FAD edge5(FAD_Type3(3, 0, v3[0]-v2[0]), FAD_Type3(3, 1, v3[1]-v2[1]), FAD_Type3(3, 2, v3[2]-v2[2]));

  const double edge0_length_squared = edge0.length_squared();
  const double edge1_length_squared = edge1.length_squared();
  const double edge2_length_squared = edge2.length_squared();
  const FAD_Type3 edge3_length_squared = edge3.length_squared();
  const FAD_Type3 edge4_length_squared = edge4.length_squared();
  const FAD_Type3 edge5_length_squared = edge5.length_squared();

  const int sign = tetVolume < 0. ? -1 : 1;
  return sign * 12. * std::pow(3.*std::abs(tetVolume), 2./3.) / (edge0_length_squared + edge1_length_squared + edge2_length_squared + edge3_length_squared + edge4_length_squared + edge5_length_squared);
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

}
