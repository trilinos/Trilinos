// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_QualityMetric.hpp>
#include <stk_util/util/ReportHandler.hpp>

namespace krino
{

double calculate_tet_volume_using_sides(const stk::math::Vector3d &side0, const stk::math::Vector3d &side2, const stk::math::Vector3d &side3)
{
    return Dot(side3, Cross(side2, side0))/6.0;
}

double MeanRatioQualityMetric::tet_mean_ratio(const std::vector<stk::math::Vector3d> &nodeLocations)
{
  ThrowAssert(nodeLocations.size() == 4 || nodeLocations.size() == 10);

  const stk::math::Vector3d side0 = nodeLocations[1] - nodeLocations[0];
  const stk::math::Vector3d side2 = nodeLocations[0] - nodeLocations[2];
  const stk::math::Vector3d side3 = nodeLocations[3] - nodeLocations[0];

  const double volumeMin = 1.0E-30;
  const double tetVolume = calculate_tet_volume_using_sides(side0, side2, side3);
  if( std::abs( tetVolume ) < volumeMin )
      return 0.0;

  const stk::math::Vector3d side1 = nodeLocations[2] - nodeLocations[1];
  const stk::math::Vector3d side4 = nodeLocations[3] - nodeLocations[1];
  const stk::math::Vector3d side5 = nodeLocations[3] - nodeLocations[2];

  const double side0_length_squared = side0.length_squared();
  const double side1_length_squared = side1.length_squared();
  const double side2_length_squared = side2.length_squared();
  const double side3_length_squared = side3.length_squared();
  const double side4_length_squared = side4.length_squared();
  const double side5_length_squared = side5.length_squared();

  const int sign = tetVolume < 0. ? -1 : 1;
  return sign * 12. * std::pow(3.*std::abs(tetVolume), 2./3.) / (side0_length_squared + side1_length_squared + side2_length_squared + side3_length_squared + side4_length_squared + side5_length_squared);
}

double ScaledJacobianQualityMetric::tet_scaled_jacobian(const std::vector<stk::math::Vector3d> &nodeLocations)
{
  ThrowAssert(nodeLocations.size() == 4 || nodeLocations.size() == 10);

  const stk::math::Vector3d side0 = nodeLocations[1] - nodeLocations[0];
  const stk::math::Vector3d side1 = nodeLocations[2] - nodeLocations[1];
  const stk::math::Vector3d side2 = nodeLocations[0] - nodeLocations[2];
  const stk::math::Vector3d side3 = nodeLocations[3] - nodeLocations[0];
  const stk::math::Vector3d side4 = nodeLocations[3] - nodeLocations[1];
  const stk::math::Vector3d side5 = nodeLocations[3] - nodeLocations[2];

  const double jacobi = Dot(side3, Cross(side2, side0));

  // products of lengths squared of each edge attached to a node.
  const double side0_length_squared = side0.length_squared();
  const double side1_length_squared = side1.length_squared();
  const double side2_length_squared = side2.length_squared();
  const double side3_length_squared = side3.length_squared();
  const double side4_length_squared = side4.length_squared();
  const double side5_length_squared = side5.length_squared();

  const double length_squared[4] = {
    side0_length_squared * side2_length_squared * side3_length_squared,
    side0_length_squared * side1_length_squared * side4_length_squared,
    side1_length_squared * side2_length_squared * side5_length_squared,
    side3_length_squared * side4_length_squared * side5_length_squared
  };
  int which_node = 0;
  if(length_squared[1] > length_squared[which_node])
    which_node = 1;
  if(length_squared[2] > length_squared[which_node])
    which_node = 2;
  if(length_squared[3] > length_squared[which_node])
    which_node = 3;

  double length_product = std::sqrt( length_squared[which_node] );
  if(length_product < std::abs(jacobi))
    length_product = std::abs(jacobi);

  const double lengthMin = 1.0E-30;
  if( length_product < lengthMin )
    return 0.0;

  static const double root_of_2 = std::sqrt(2.0);
  return root_of_2 * jacobi / length_product;
}

double ScaledJacobianQualityMetric::tri2d_scaled_jacobian(const std::vector<stk::math::Vector3d> &nodeLocations)
{
  const double absScaledJacobian = tri3d_scaled_jacobian(nodeLocations);
  const double normalZ =
        (nodeLocations[1][0]-nodeLocations[0][0])*(nodeLocations[2][1]-nodeLocations[0][1]) -
        (nodeLocations[1][1]-nodeLocations[0][1])*(nodeLocations[2][0]-nodeLocations[0][0]);
  return (normalZ > 0.) ? absScaledJacobian : -absScaledJacobian;
}

double ScaledJacobianQualityMetric::tri3d_scaled_jacobian(const std::vector<stk::math::Vector3d> &nodeLocations)
{
  ThrowAssert(nodeLocations.size() == 3 || nodeLocations.size() == 6);

  const stk::math::Vector3d edge0 = nodeLocations[1] - nodeLocations[0];
  const stk::math::Vector3d edge1 = nodeLocations[2] - nodeLocations[0];
  const stk::math::Vector3d edge2 = nodeLocations[2] - nodeLocations[1];

  const double lenSqr0 = edge0.length_squared();
  const double lenSqr1 = edge1.length_squared();
  const double lenSqr2 = edge2.length_squared();

  const double maxEdgeLengthProduct = std::sqrt( std::max(lenSqr0*lenSqr1, std::max(lenSqr1*lenSqr2, lenSqr0*lenSqr2)) );

  const double lengthMin = 1.0E-30;
  if( maxEdgeLengthProduct < lengthMin )
    return 0.0;

  static const double two_over_root_of_3 = 2./sqrt(3.0);
  const double jacobian = Cross(edge0, edge1).length();
  return jacobian*two_over_root_of_3/maxEdgeLengthProduct;
}

}
