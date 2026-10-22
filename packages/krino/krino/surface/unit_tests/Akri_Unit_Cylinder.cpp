// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_AnalyticSurf.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkVector.hpp>

namespace krino {

// Signed-distance and closest-point-normal coverage for the finite Cylinder
// (a capped cylinder / "can").  point_signed_distance has three regimes --
// alongside the barrel (inside vs. outside the radius), and beyond an end cap
// (over the flat face vs. off the circular rim); closest_point_normal
// classifies bottom cap / top cap / lateral surface / rim.  Reference cylinder
// runs along +z from (0,0,0) to (0,0,4) with radius 1.

namespace {
const stk::math::Vector3d endPt1(0.,0.,0.);
const stk::math::Vector3d endPt2(0.,0.,4.);
const double radius = 1.;
const double tol = 1.e-9;
}

//--------------------------------------------------------------------------
// point_signed_distance
//--------------------------------------------------------------------------

TEST(Cylinder, signedDistance_alongsideBarrel)
{
  Cylinder cylinder(endPt1, endPt2, radius);

  // On the axis at the midpoint: closest surface is the barrel, distance -radius.
  EXPECT_NEAR(-1., cylinder.point_signed_distance(stk::math::Vector3d(0.,0.,2.)), tol);
  // Inside, partway to the barrel.
  EXPECT_NEAR(-0.5, cylinder.point_signed_distance(stk::math::Vector3d(0.5,0.,2.)), tol);
  // On the lateral surface.
  EXPECT_NEAR(0., cylinder.point_signed_distance(stk::math::Vector3d(1.,0.,2.)), tol);
  // Outside radially.
  EXPECT_NEAR(1., cylinder.point_signed_distance(stk::math::Vector3d(2.,0.,2.)), tol);
}

TEST(Cylinder, signedDistance_insideNearCapIsDistanceToCap)
{
  Cylinder cylinder(endPt1, endPt2, radius);

  // On the axis 0.25 above the bottom face: closest surface is the cap, not the barrel.
  EXPECT_NEAR(-0.25, cylinder.point_signed_distance(stk::math::Vector3d(0.,0.,0.25)), tol);
  // Symmetric near the top face.
  EXPECT_NEAR(-0.25, cylinder.point_signed_distance(stk::math::Vector3d(0.,0.,3.75)), tol);
}

TEST(Cylinder, signedDistance_beyondCapOverFlatFace)
{
  Cylinder cylinder(endPt1, endPt2, radius);

  // Directly beyond the top face (still within the radius): distance is axial gap.
  EXPECT_NEAR(2., cylinder.point_signed_distance(stk::math::Vector3d(0.,0.,6.)), tol);
  // Directly beyond the bottom face.
  EXPECT_NEAR(3., cylinder.point_signed_distance(stk::math::Vector3d(0.,0.,-3.)), tol);
}

TEST(Cylinder, signedDistance_beyondCapOffRimIsEuclidean)
{
  Cylinder cylinder(endPt1, endPt2, radius);

  // (4,0,8): 4 beyond the top face axially and radius+3 from the axis, so the
  // closest point is the top rim -> sqrt(4^2 + 3^2) = 5.
  EXPECT_NEAR(5., cylinder.point_signed_distance(stk::math::Vector3d(4.,0.,8.)), tol);
}

TEST(Cylinder, negativeSign_flipsDistance)
{
  Cylinder cylinder(endPt1, endPt2, radius, -1);

  EXPECT_NEAR(1., cylinder.point_signed_distance(stk::math::Vector3d(0.,0.,2.)), tol);  // inside now positive
  EXPECT_NEAR(-1., cylinder.point_signed_distance(stk::math::Vector3d(2.,0.,2.)), tol); // outside now negative
}

TEST(Cylinder, hasEdges)
{
  Cylinder cylinder(endPt1, endPt2, radius);
  EXPECT_TRUE(cylinder.has_edges());
}

//--------------------------------------------------------------------------
// closest_point_normal
//--------------------------------------------------------------------------

TEST(Cylinder, closestPointNormal_lateralSurfacePointsRadiallyOutward)
{
  Cylinder cylinder(endPt1, endPt2, radius);
  EXPECT_TRUE(cylinder.can_approximate_closest_point_normal());

  // Outside the barrel at the midpoint.
  expect_eq(stk::math::Vector3d(1.,0.,0.), cylinder.closest_point_normal(stk::math::Vector3d(2.,0.,2.)));
  // Inside near the barrel at the midpoint still resolves to the lateral normal.
  expect_eq(stk::math::Vector3d(1.,0.,0.), cylinder.closest_point_normal(stk::math::Vector3d(0.5,0.,2.)));
}

TEST(Cylinder, closestPointNormal_capsPointAlongAxis)
{
  Cylinder cylinder(endPt1, endPt2, radius);

  // Near the bottom face, within the radius: normal is -axis.
  expect_eq(stk::math::Vector3d(0.,0.,-1.), cylinder.closest_point_normal(stk::math::Vector3d(0.2,0.,-1.)));
  // Near the top face, within the radius: normal is +axis.
  expect_eq(stk::math::Vector3d(0.,0.,1.), cylinder.closest_point_normal(stk::math::Vector3d(0.2,0.,5.)));
}

//--------------------------------------------------------------------------
// does_intersect / insert_into (padded axis-aligned bounding box)
//--------------------------------------------------------------------------

TEST(Cylinder, doesIntersect_trueWhenBoxOverlapsPaddedBoxElseFalse)
{
  Cylinder cylinder(endPt1, endPt2, radius);

  const BoundingBox overlapping(stk::math::Vector3d(0.5,0.5,0.5), stk::math::Vector3d(2.,2.,2.));
  EXPECT_TRUE(cylinder.does_intersect(overlapping));

  const BoundingBox faraway(stk::math::Vector3d(10.,10.,10.), stk::math::Vector3d(11.,11.,11.));
  EXPECT_FALSE(cylinder.does_intersect(faraway));
}

TEST(Cylinder, insertInto_accommodatesPaddedExtent)
{
  Cylinder cylinder(endPt1, endPt2, radius);

  BoundingBox bbox;
  cylinder.insert_into(bbox);
  // Endpoints padded by the radius: z spans [-1,5], x/y span [-1,1].
  EXPECT_TRUE(bbox.contains(stk::math::Vector3d(1.,1.,5.)));
  EXPECT_TRUE(bbox.contains(stk::math::Vector3d(-1.,-1.,-1.)));
  EXPECT_FALSE(bbox.contains(stk::math::Vector3d(0.,0.,6.)));
}

}
