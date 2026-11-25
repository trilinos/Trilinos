// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_SharpFeatureCalcs.hpp>

#include <cmath>
#include <numbers>

namespace krino
{

// Triangles sharing the edge along x, with their off-edge vertices flat, at a right angle, at 45
// degrees and folded onto each other.
TEST(SharpFeatureCalcs, cosineOfAngleBetweenTrianglesAtEdge_flatRightAngleAndFolded)
{
  const stk::math::Vector3d edgePt0(0., 0., 0.);
  const stk::math::Vector3d edgePt1(1., 0., 0.);
  EXPECT_EQ(-1.,
      cosine_of_angle_between_triangles_at_edge(edgePt0, edgePt1, {0.5, 1., 0.}, {0.5, -1., 0.}));
  EXPECT_EQ(0.,
      cosine_of_angle_between_triangles_at_edge(edgePt0, edgePt1, {0.5, 1., 0.}, {0.5, 0., 1.}));
  EXPECT_NEAR(std::sqrt(0.5),
      cosine_of_angle_between_triangles_at_edge(edgePt0, edgePt1, {0.5, 1., 0.}, {0.5, 1., 1.}),
      1.e-15);
  EXPECT_EQ(1.,
      cosine_of_angle_between_triangles_at_edge(edgePt0, edgePt1, {0.5, 1., 0.}, {0.2, 2., 0.}));
}

TEST(SharpFeatureCalcs, cosineOfAngleBetweenSegmentsAtVertex_straightRightAngleAndFolded)
{
  const stk::math::Vector3d vertex(1., 1., 0.);
  EXPECT_EQ(-1., cosine_of_angle_between_segments_at_vertex(vertex, {3., 1., 0.}, {0., 1., 0.}));
  EXPECT_EQ(0., cosine_of_angle_between_segments_at_vertex(vertex, {3., 1., 0.}, {1., 2., 0.}));
  EXPECT_EQ(1., cosine_of_angle_between_segments_at_vertex(vertex, {3., 1., 0.}, {2., 1., 0.}));
}

// A feature is sharp when its angle is less than the feature angle.
TEST(SharpFeatureCalcs, isAngleSharperThanFeatureAngle_onlyAnglesLessThanTheFeatureAngle)
{
  const double cosFeatureAngle = std::cos(135. * std::numbers::pi / 180.);
  EXPECT_TRUE(is_angle_sharper_than_feature_angle(0., cosFeatureAngle));   // 90 degrees
  EXPECT_FALSE(is_angle_sharper_than_feature_angle(-1., cosFeatureAngle)); // flat
  EXPECT_FALSE(is_angle_sharper_than_feature_angle(
      cosFeatureAngle, cosFeatureAngle)); // exactly the feature angle
}

} // namespace krino
