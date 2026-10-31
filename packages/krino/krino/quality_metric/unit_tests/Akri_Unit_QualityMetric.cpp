// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_QualityMetric.hpp>

#include <stk_math/StkVector.hpp>

#include <cmath>
#include <vector>

namespace krino {

// Pure-geometry coverage for the element quality metrics.  Both the
// MeanRatioQualityMetric and the ScaledJacobianQualityMetric are normalized so
// that an ideal (regular) element scores exactly 1.0; degenerate (zero-volume)
// elements score 0.  All gold values below are computed analytically.  Also
// exercised: the get_element_quality_metric dim/size dispatch, the "higher is
// better" comparison, and the quad-diagonal angle helpers.

namespace {
const double tol = 1.e-9;

// Regular tetrahedron with all edge lengths equal to 1.
std::vector<stk::math::Vector3d> regular_tet()
{
  const double s3 = std::sqrt(3.);
  const double s6 = std::sqrt(6.);
  return {
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0.5, s3/2., 0.),
    stk::math::Vector3d(0.5, s3/6., s6/3.)};
}

// Right tetrahedron on the coordinate axes (unit legs).
std::vector<stk::math::Vector3d> right_tet()
{
  return {
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0., 1., 0.),
    stk::math::Vector3d(0., 0., 1.)};
}
}

//--------------------------------------------------------------------------
// MeanRatioQualityMetric
//--------------------------------------------------------------------------

TEST(MeanRatioQualityMetric, regularTet_scoresOne)
{
  EXPECT_NEAR(1., MeanRatioQualityMetric::tet_mean_ratio(regular_tet()), tol);
}

TEST(MeanRatioQualityMetric, rightTet_scoresKnownValueLessThanOne)
{
  // metric = sqrt(2) * (6V) / (meanEdgeLenSq)^1.5, with 6V=1 and
  // meanEdgeLenSq = (1+2+1+1+2+2)/6 = 1.5.
  const double gold = std::sqrt(2.) / std::pow(1.5, 1.5);
  EXPECT_NEAR(gold, MeanRatioQualityMetric::tet_mean_ratio(right_tet()), tol);
  EXPECT_LT(MeanRatioQualityMetric::tet_mean_ratio(right_tet()), 1.);
}

TEST(MeanRatioQualityMetric, flatTet_scoresZero)
{
  const std::vector<stk::math::Vector3d> flatTet{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0., 1., 0.),
    stk::math::Vector3d(1., 1., 0.)};
  EXPECT_NEAR(0., MeanRatioQualityMetric::tet_mean_ratio(flatTet), tol);
}

TEST(MeanRatioQualityMetric, equilateralTri2d_scoresOne)
{
  const std::vector<stk::math::Vector3d> equilateral{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0.5, std::sqrt(3.)/2., 0.)};
  EXPECT_NEAR(1., MeanRatioQualityMetric::tri2d_mean_ratio(equilateral), tol);
}

TEST(MeanRatioQualityMetric, rightTri2d_scoresKnownValue)
{
  // coeff=2/sqrt(3), 2*area=1, meanEdgeLenSq=(1+2+1)/3=4/3 -> sqrt(3)/2.
  const std::vector<stk::math::Vector3d> rightTri{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0., 1., 0.)};
  EXPECT_NEAR(std::sqrt(3.)/2., MeanRatioQualityMetric::tri2d_mean_ratio(rightTri), tol);
}

TEST(MeanRatioQualityMetric, getElementQualityMetric_dispatchesOnDim)
{
  MeanRatioQualityMetric metric;

  const std::vector<stk::math::Vector3d> tri{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0., 1., 0.)};
  EXPECT_NEAR(MeanRatioQualityMetric::tri2d_mean_ratio(tri), metric.get_element_quality_metric(2, tri), tol);

  const std::vector<stk::math::Vector3d> tet = right_tet();
  EXPECT_NEAR(MeanRatioQualityMetric::tet_mean_ratio(tet), metric.get_element_quality_metric(3, tet), tol);
}

TEST(MeanRatioQualityMetric, bestAndComparison_higherIsBetter)
{
  MeanRatioQualityMetric metric;
  EXPECT_NEAR(1., metric.get_best_value_for_metric(), tol);
  EXPECT_NEAR(0.2, metric.get_acceptable_value_for_metric(), tol);
  EXPECT_TRUE(metric.is_first_quality_metric_better_than_second(0.9, 0.5));
  EXPECT_FALSE(metric.is_first_quality_metric_better_than_second(0.5, 0.9));
}

//--------------------------------------------------------------------------
// ScaledJacobianQualityMetric
//--------------------------------------------------------------------------

TEST(ScaledJacobianQualityMetric, regularTet_scoresOne)
{
  EXPECT_NEAR(1., ScaledJacobianQualityMetric::tet_scaled_jacobian(regular_tet()), tol);
}

TEST(ScaledJacobianQualityMetric, rightTet_scoresKnownValue)
{
  // jacobi=1, longest edge-length-product at a node = sqrt(2*2)=2 -> sqrt(2)*1/2.
  EXPECT_NEAR(std::sqrt(2.)/2., ScaledJacobianQualityMetric::tet_scaled_jacobian(right_tet()), tol);
}

TEST(ScaledJacobianQualityMetric, flatTet_scoresZero)
{
  const std::vector<stk::math::Vector3d> flatTet{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0., 1., 0.),
    stk::math::Vector3d(1., 1., 0.)};
  EXPECT_NEAR(0., ScaledJacobianQualityMetric::tet_scaled_jacobian(flatTet), tol);
}

TEST(ScaledJacobianQualityMetric, equilateralTri3d_scoresOne)
{
  const std::vector<stk::math::Vector3d> equilateral{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0.5, std::sqrt(3.)/2., 0.)};
  EXPECT_NEAR(1., ScaledJacobianQualityMetric::tri3d_scaled_jacobian(equilateral), tol);
}

TEST(ScaledJacobianQualityMetric, tri2d_signFollowsOrientation)
{
  // 2d scaled jacobian carries the sign of the (CCW-positive) signed area.
  const double gold = 2./std::sqrt(6.); // 2/sqrt(3) * 1 / sqrt(2) for the right triangle
  const std::vector<stk::math::Vector3d> ccw{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0., 1., 0.)};
  EXPECT_NEAR(gold, ScaledJacobianQualityMetric::tri2d_scaled_jacobian(ccw), tol);

  const std::vector<stk::math::Vector3d> cw{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(0., 1., 0.),
    stk::math::Vector3d(1., 0., 0.)};
  EXPECT_NEAR(-gold, ScaledJacobianQualityMetric::tri2d_scaled_jacobian(cw), tol);
}

TEST(ScaledJacobianQualityMetric, getElementQualityMetric_dispatchesOnDimAndSize)
{
  ScaledJacobianQualityMetric metric;

  const std::vector<stk::math::Vector3d> tet = right_tet();
  EXPECT_NEAR(ScaledJacobianQualityMetric::tet_scaled_jacobian(tet), metric.get_element_quality_metric(3, tet), tol);

  // 3 nodes in 3d -> triangle path (not tet).
  const std::vector<stk::math::Vector3d> tri3d{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0.5, std::sqrt(3.)/2., 0.)};
  EXPECT_NEAR(ScaledJacobianQualityMetric::tri3d_scaled_jacobian(tri3d), metric.get_element_quality_metric(3, tri3d), tol);

  const std::vector<stk::math::Vector3d> tri2d{
    stk::math::Vector3d(0., 0., 0.),
    stk::math::Vector3d(1., 0., 0.),
    stk::math::Vector3d(0., 1., 0.)};
  EXPECT_NEAR(ScaledJacobianQualityMetric::tri2d_scaled_jacobian(tri2d), metric.get_element_quality_metric(2, tri2d), tol);
}

TEST(ScaledJacobianQualityMetric, bestAndComparison_higherIsBetter)
{
  ScaledJacobianQualityMetric metric;
  EXPECT_NEAR(1., metric.get_best_value_for_metric(), tol);
  EXPECT_NEAR(0.1, metric.get_acceptable_value_for_metric(), tol);
  EXPECT_TRUE(metric.is_first_quality_metric_better_than_second(0.8, 0.3));
  EXPECT_FALSE(metric.is_first_quality_metric_better_than_second(0.3, 0.8));
}

//--------------------------------------------------------------------------
// lexicographic comparison and quad-diagonal helpers
//--------------------------------------------------------------------------

TEST(QualityMetricHelpers, isLessThanInXThenYThenZ)
{
  EXPECT_TRUE(is_less_than_in_x_then_y_then_z(stk::math::Vector3d(0,5,5), stk::math::Vector3d(1,0,0)));
  EXPECT_TRUE(is_less_than_in_x_then_y_then_z(stk::math::Vector3d(1,0,5), stk::math::Vector3d(1,1,0)));
  EXPECT_TRUE(is_less_than_in_x_then_y_then_z(stk::math::Vector3d(1,1,0), stk::math::Vector3d(1,1,1)));
  EXPECT_FALSE(is_less_than_in_x_then_y_then_z(stk::math::Vector3d(1,1,1), stk::math::Vector3d(1,1,1)));
  EXPECT_FALSE(is_less_than_in_x_then_y_then_z(stk::math::Vector3d(2,0,0), stk::math::Vector3d(1,9,9)));
}

TEST(QualityMetricHelpers, quadDiagonal_symmetricSquareIsAmbiguous)
{
  // A unit square has four equal right angles: neither diagonal cuts a larger
  // angle, so the result is 0 and cutting 0-2 is not preferred.
  const stk::math::Vector3d x0(0,0,0), x1(1,0,0), x2(1,1,0), x3(0,1,0);
  EXPECT_EQ(0, determine_diagonal_of_quad_that_cuts_largest_angle(x0,x1,x2,x3));
  EXPECT_FALSE(will_cutting_quad_from_0to2_cut_largest_angle(x0,x1,x2,x3));
}

TEST(QualityMetricHelpers, quadDiagonal_largestAnglesAtNodes0And2)
{
  // A rhombus pinched along the 0-2 axis: the obtuse angles sit at nodes 0 and
  // 2, so cutting the 0-2 diagonal splits the largest angles (-> -1, true).
  const stk::math::Vector3d x0(0,0.2,0), x1(1,0,0), x2(0,-0.2,0), x3(-1,0,0);
  EXPECT_EQ(-1, determine_diagonal_of_quad_that_cuts_largest_angle(x0,x1,x2,x3));
  EXPECT_TRUE(will_cutting_quad_from_0to2_cut_largest_angle(x0,x1,x2,x3));
}

TEST(QualityMetricHelpers, quadDiagonal_largestAnglesAtNodes1And3)
{
  // The complementary rhombus: obtuse angles at nodes 1 and 3, so the 1-3
  // diagonal is preferred (-> +1) and cutting 0-2 is not.
  const stk::math::Vector3d x0(-1,0,0), x1(0,0.2,0), x2(1,0,0), x3(0,-0.2,0);
  EXPECT_EQ(1, determine_diagonal_of_quad_that_cuts_largest_angle(x0,x1,x2,x3));
  EXPECT_FALSE(will_cutting_quad_from_0to2_cut_largest_angle(x0,x1,x2,x3));
}

}
