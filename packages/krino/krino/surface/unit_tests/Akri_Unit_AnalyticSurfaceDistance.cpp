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
#include <Akri_Transformation.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkVector.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace krino {

// These tests lock down the signed-distance / closest-point contracts of the
// analytic surfaces (Sphere, Plane, Point, Ellipsoid).  The companion file
// Akri_Unit_AnalyticSurface.cpp already covers Cuboid distances thoroughly and
// one Sphere edge intersection; here we exercise the point_signed_distance,
// closest_point_normal and does_intersect behavior those tests leave untested,
// including the sign flag, the (surprising) squared-distance return of Point,
// and Ellipsoid rotation.

//--------------------------------------------------------------------------
// Sphere: signed distance = sign*(|x-center| - radius).
//--------------------------------------------------------------------------

TEST(Sphere, signedDistance_insideZeroOutside)
{
  const stk::math::Vector3d center(1.,2.,3.);
  const double radius = 2.;
  Sphere sphere(center, radius);

  const double tol = 1.e-9;
  EXPECT_NEAR(-radius, sphere.point_signed_distance(center), tol);              // at center
  EXPECT_NEAR(0., sphere.point_signed_distance(center + stk::math::Vector3d(radius,0,0)), tol); // on surface
  EXPECT_NEAR(3., sphere.point_signed_distance(center + stk::math::Vector3d(5,0,0)), tol);      // outside
  EXPECT_NEAR(-1., sphere.point_signed_distance(center + stk::math::Vector3d(0,1,0)), tol);     // inside
}

TEST(Sphere, negativeSign_flipsDistanceOnly)
{
  const stk::math::Vector3d center = stk::math::Vector3d::ZERO;
  const double radius = 2.;
  Sphere sphere(center, radius, -1);

  const double tol = 1.e-9;
  EXPECT_NEAR(radius, sphere.point_signed_distance(center), tol);        // inside now positive
  EXPECT_NEAR(0., sphere.point_signed_distance(stk::math::Vector3d(radius,0,0)), tol);
  EXPECT_NEAR(-3., sphere.point_signed_distance(stk::math::Vector3d(5,0,0)), tol); // outside now negative
}

TEST(Sphere, closestPointNormal_pointsRadiallyOutward)
{
  const stk::math::Vector3d center(1.,2.,3.);
  Sphere sphere(center, 2.);

  EXPECT_TRUE(sphere.can_approximate_closest_point_normal());
  expect_eq(stk::math::Vector3d(1,0,0), sphere.closest_point_normal(center + stk::math::Vector3d(5,0,0)));
  expect_eq(stk::math::Vector3d(0,1,0), sphere.closest_point_normal(center + stk::math::Vector3d(0,0.5,0)));
  const stk::math::Vector3d diagonal(1,1,1);
  expect_eq(diagonal.unit_vector(), sphere.closest_point_normal(center + diagonal));
}

TEST(Sphere, negativeSign_closestPointNormalStillOutward)
{
  // The sign flag flips the distance, not the normal direction.
  const stk::math::Vector3d center = stk::math::Vector3d::ZERO;
  Sphere sphere(center, 2., -1);
  expect_eq(stk::math::Vector3d(1,0,0), sphere.closest_point_normal(stk::math::Vector3d(5,0,0)));
}

TEST(Sphere, doesIntersect_trueWhenBoxOverlapsSphereBoxElseFalse)
{
  Sphere sphere(stk::math::Vector3d::ZERO, 1.);

  const BoundingBox overlapping(stk::math::Vector3d(0.5,0.5,0.5), stk::math::Vector3d(2,2,2));
  EXPECT_TRUE(sphere.does_intersect(overlapping));

  const BoundingBox faraway(stk::math::Vector3d(10,10,10), stk::math::Vector3d(11,11,11));
  EXPECT_FALSE(sphere.does_intersect(faraway));
}

//--------------------------------------------------------------------------
// Plane: signed distance = multiplier*(unit_normal . x + offset).
//--------------------------------------------------------------------------

TEST(Plane, signedDistance_normalizesNormalAndUsesOffset)
{
  // Non-unit normal along +x; plane x = 2 expressed as normal.x + offset = 0
  // with offset = -2 (after the normal is unitized in the constructor).
  Plane plane(stk::math::Vector3d(5,0,0), -2.);

  const double tol = 1.e-9;
  EXPECT_NEAR(0., plane.point_signed_distance(stk::math::Vector3d(2,0,0)), tol);   // on plane
  EXPECT_NEAR(1., plane.point_signed_distance(stk::math::Vector3d(3,7,-4)), tol);  // +x side, y/z irrelevant
  EXPECT_NEAR(-1.5, plane.point_signed_distance(stk::math::Vector3d(0.5,0,0)), tol); // -x side
}

TEST(Plane, multiplier_scalesAndFlipsSign)
{
  Plane plane(stk::math::Vector3d(0,1,0), 0., -2.);

  const double tol = 1.e-9;
  EXPECT_NEAR(0., plane.point_signed_distance(stk::math::Vector3d(0,0,0)), tol);
  EXPECT_NEAR(-2., plane.point_signed_distance(stk::math::Vector3d(0,1,0)), tol);  // multiplier scales & flips
  EXPECT_NEAR(6., plane.point_signed_distance(stk::math::Vector3d(0,-3,0)), tol);
}

TEST(Plane, pointConstructor_planePassesThroughPoint)
{
  // Plane through (1,2,3) with non-unit normal along +y: the plane y = 2.
  Plane plane(stk::math::Vector3d(0,4,0), stk::math::Vector3d(1,2,3), 1.);

  const double tol = 1.e-9;
  EXPECT_NEAR(0., plane.point_signed_distance(stk::math::Vector3d(1,2,3)), tol);
  EXPECT_NEAR(0., plane.point_signed_distance(stk::math::Vector3d(-5,2,8)), tol);   // anywhere on y = 2
  EXPECT_NEAR(1.5, plane.point_signed_distance(stk::math::Vector3d(0,3.5,0)), tol);
  EXPECT_NEAR(-2., plane.point_signed_distance(stk::math::Vector3d(0,0,0)), tol);
}

TEST(Plane, pointConstructor_matchesOffsetConstructor)
{
  const stk::math::Vector3d normal(1,-2,2);
  const stk::math::Vector3d point(0.5,0.25,-1.);
  const double multiplier = -3.;
  // Offset of the plane through point, for the unit normal.
  const double offset = -Dot(normal, point)/normal.length();

  const Plane fromPoint(normal, point, multiplier);
  const Plane fromOffset(normal, offset, multiplier);

  const double tol = 1.e-9;
  for (const stk::math::Vector3d & x : {stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,2,3), stk::math::Vector3d(-4,0.5,7)})
    EXPECT_NEAR(fromOffset.point_signed_distance(x), fromPoint.point_signed_distance(x), tol);
}

TEST(Plane, fractionOfExtent_offsetSetFromExtentAlongUnitNormal)
{
  // Non-unit normal along (1,1,0); extent of normal.x over a mesh given along the unit normal.
  const stk::math::Vector3d normal(2,2,0);
  const stk::math::Vector3d unitNormal = normal.unit_vector();
  const double minDot = -1.;
  const double maxDot = 3.;

  const double tol = 1.e-9;
  for (const double fraction : {0., 0.25, 1., -0.5, 1.5})
  {
    auto plane = Plane::at_fraction_of_extent(normal, fraction, 1.);
    EXPECT_TRUE(plane->needs_offset_from_extent());
    EXPECT_DOUBLE_EQ(fraction, plane->fraction_of_extent());
    EXPECT_NEAR(0., (plane->normal() - unitNormal).length(), tol);

    plane->set_offset_from_extent(minDot, maxDot);
    EXPECT_FALSE(plane->needs_offset_from_extent());

    const double planeDot = minDot + fraction*(maxDot - minDot);
    EXPECT_NEAR(0., plane->point_signed_distance(planeDot*unitNormal), tol) << fraction;
    EXPECT_NEAR(0., plane->point_signed_distance(planeDot*unitNormal + stk::math::Vector3d(-4,4,7)), tol) << fraction; // along the plane
    EXPECT_NEAR(2., plane->point_signed_distance((planeDot+2.)*unitNormal), tol) << fraction;
  }
}

TEST(Plane, fractionOfExtent_multiplierScalesAndFlipsSign)
{
  auto plane = Plane::at_fraction_of_extent(stk::math::Vector3d(0,1,0), 0.5, -2.);
  plane->set_offset_from_extent(0., 4.);

  const double tol = 1.e-9;
  EXPECT_NEAR(0., plane->point_signed_distance(stk::math::Vector3d(0,2,0)), tol);
  EXPECT_NEAR(-2., plane->point_signed_distance(stk::math::Vector3d(0,3,0)), tol);
}

TEST(Plane, doesIntersect_alwaysTrue_andInsertsEntireDomain)
{
  Plane plane(stk::math::Vector3d(0,0,1), 0.);

  const BoundingBox anywhere(stk::math::Vector3d(10,10,10), stk::math::Vector3d(11,11,11));
  EXPECT_TRUE(plane.does_intersect(anywhere));

  BoundingBox bbox;
  plane.insert_into(bbox);
  EXPECT_TRUE(bbox.contains(stk::math::Vector3d(1e30,-1e30,0.)));
}

//--------------------------------------------------------------------------
// Point: point_signed_distance returns the SQUARED distance (never negative).
//--------------------------------------------------------------------------

TEST(Point, signedDistance_returnsSquaredDistance_nonNegative)
{
  const stk::math::Vector3d loc(1.,2.,3.);
  Point point(loc);

  const double tol = 1.e-9;
  EXPECT_NEAR(0., point.point_signed_distance(loc), tol);
  EXPECT_NEAR(9., point.point_signed_distance(loc + stk::math::Vector3d(3,0,0)), tol);   // 3^2
  EXPECT_NEAR(25., point.point_signed_distance(loc + stk::math::Vector3d(0,-3,4)), tol); // 3^2+4^2
}

TEST(Point, doesIntersect_trueOnlyWhenBoxContainsPoint)
{
  Point point(stk::math::Vector3d(1.,1.,1.));

  const BoundingBox containing(stk::math::Vector3d(0,0,0), stk::math::Vector3d(2,2,2));
  EXPECT_TRUE(point.does_intersect(containing));

  const BoundingBox notContaining(stk::math::Vector3d(0,0,0), stk::math::Vector3d(0.5,0.5,0.5));
  EXPECT_FALSE(point.does_intersect(notContaining));
}

//--------------------------------------------------------------------------
// Ellipsoid: approximate signed distance, exact only when it reduces to a
// sphere.  distance = |semiAxes|*(sqrt(sum((d_i/a_i)^2)) - 1).
//--------------------------------------------------------------------------

TEST(Ellipsoid, sphereEquivalent_matchesSphereDistanceExactly)
{
  // When all semi-axes are equal the ellipsoid is a sphere, and the
  // characteristic-radius scaling makes the signed distance exact.
  const std::vector<double> center{1.,2.,3.};
  const std::vector<double> semiAxes{2.,2.,2.};
  Ellipsoid ellipsoid(center, semiAxes, {}, 1);
  Sphere sphere(stk::math::Vector3d(1,2,3), 2.);

  const double tol = 1.e-9;
  const stk::math::Vector3d c(1,2,3);
  for (const stk::math::Vector3d & x : {c, c+stk::math::Vector3d(2,0,0), c+stk::math::Vector3d(5,0,0), c+stk::math::Vector3d(0,1,0)})
    EXPECT_NEAR(sphere.point_signed_distance(x), ellipsoid.point_signed_distance(x), tol);
}

TEST(Ellipsoid, unequalAxes_zeroOnSurfaceNegativeInsidePositiveOutside)
{
  const std::vector<double> center{0.,0.,0.};
  const std::vector<double> semiAxes{4.,2.,1.};
  Ellipsoid ellipsoid(center, semiAxes, {}, 1);

  const double tol = 1.e-9;
  // Zero at each semi-axis tip.
  EXPECT_NEAR(0., ellipsoid.point_signed_distance(stk::math::Vector3d(4,0,0)), tol);
  EXPECT_NEAR(0., ellipsoid.point_signed_distance(stk::math::Vector3d(0,2,0)), tol);
  EXPECT_NEAR(0., ellipsoid.point_signed_distance(stk::math::Vector3d(0,0,1)), tol);
  // Negative at center, positive well outside.
  EXPECT_LT(ellipsoid.point_signed_distance(stk::math::Vector3d::ZERO), 0.);
  EXPECT_GT(ellipsoid.point_signed_distance(stk::math::Vector3d(8,0,0)), 0.);
}

TEST(Ellipsoid, negativeSign_flipsInsideOutside)
{
  const std::vector<double> center{0.,0.,0.};
  const std::vector<double> semiAxes{4.,2.,1.};
  Ellipsoid ellipsoid(center, semiAxes, {}, -1);

  EXPECT_GT(ellipsoid.point_signed_distance(stk::math::Vector3d::ZERO), 0.);
  EXPECT_LT(ellipsoid.point_signed_distance(stk::math::Vector3d(8,0,0)), 0.);
}

TEST(Ellipsoid, rotated_matchesUnrotatedAtRotatedQueryPoint)
{
  // A query point rotated by the ellipsoid's rotation should yield the same
  // distance from the rotated ellipsoid as the unrotated point does from the
  // unrotated ellipsoid.
  const std::vector<double> center{0.,0.,0.};
  const std::vector<double> semiAxes{4.,2.,1.};
  const std::vector<double> rotationVec{0.3,-0.5,0.8};

  Ellipsoid unrotated(center, semiAxes, {}, 1);
  Ellipsoid rotated(center, semiAxes, rotationVec, 1);

  Quaternion q;
  q.set_from_rotation_vector(stk::math::Vector3d(rotationVec[0], rotationVec[1], rotationVec[2]));

  const double tol = 1.e-9;
  for (const stk::math::Vector3d & x : {stk::math::Vector3d(3,0,0), stk::math::Vector3d(0,1.5,0.5), stk::math::Vector3d(1,-1,0.7)})
  {
    const stk::math::Vector3d xRotated = q.rotate_3d_vector(x);
    EXPECT_NEAR(unrotated.point_signed_distance(x), rotated.point_signed_distance(xRotated), tol);
  }
}

// A quarter turn about z of an ellipsoid off the origin swaps the x and y extents of its bounding box
// about its center.
TEST(Ellipsoid, rotatedAndTranslated_boundingBoxOfTheTurnedEllipsoid)
{
  const stk::math::Vector3d center(1.,-2.,3.);
  Ellipsoid ellipsoid({center[0], center[1], center[2]}, {1.,2.,3.}, {0.,0.,0.5*M_PI}, 1);
  BoundingBox bbox;
  ellipsoid.insert_into(bbox);
  const stk::math::Vector3d halfExtents(2.,1.,3.);
  const double bboxTol = 1.e-5; // BoundingBox is single precision
  for (int i=0; i<3; ++i)
  {
    EXPECT_NEAR((center-halfExtents)[i], bbox.get_min()[i], bboxTol);
    EXPECT_NEAR((center+halfExtents)[i], bbox.get_max()[i], bboxTol);
  }
}

// For a general rotation of an ellipsoid off the origin, the bounding box contains the points on its
// surface in many directions from its center.
TEST(Ellipsoid, rotatedAndTranslated_boundingBoxContainsPointsOnTheSurface)
{
  const stk::math::Vector3d center(1.,-2.,3.);
  const stk::math::Vector3d semiAxes(1.,2.,3.);
  Ellipsoid ellipsoid({center[0], center[1], center[2]}, {semiAxes[0], semiAxes[1], semiAxes[2]}, {0.5,0.25,0.1}, 1);
  BoundingBox bbox;
  ellipsoid.insert_into(bbox);

  // The points where the distance is zero, found by bisection, so the test doesn't depend on the
  // rotation convention.
  for (const stk::math::Vector3d & dir : {stk::math::Vector3d(1.,0.,0.), stk::math::Vector3d(0.,1.,0.), stk::math::Vector3d(0.,0.,1.),
         stk::math::Vector3d(1.,1.,0.), stk::math::Vector3d(0.,1.,1.), stk::math::Vector3d(1.,0.,1.), stk::math::Vector3d(1.,1.,1.)})
  {
    for (const double sign : {-1., 1.})
    {
      const stk::math::Vector3d unitDir = sign*dir.unit_vector();
      double lo = 0., hi = 4.;
      for (int iter=0; iter<60; ++iter)
      {
        const double mid = 0.5*(lo+hi);
        if (ellipsoid.point_signed_distance(center + mid*unitDir) < 0.) lo = mid; else hi = mid;
      }
      const stk::math::Vector3d onSurface = center + lo*unitDir;
      const double bboxTol = 1.e-5; // BoundingBox is single precision
      for (int i=0; i<3; ++i)
      {
        EXPECT_LE(bbox.get_min()[i]-bboxTol, onSurface[i]);
        EXPECT_GE(bbox.get_max()[i]+bboxTol, onSurface[i]);
      }
    }
  }
}

// A rotation vector turns an ellipsoid and a cuboid the same way, about their centers, by the right
// hand rule: a quarter turn about z turns their x axes to y.  So the ends of the ellipsoid's axes
// are at the middles of the faces of the cuboid around it.
TEST(Ellipsoid, rotation_sameAsCuboid)
{
  const stk::math::Vector3d center(1.,-2.,3.);
  const stk::math::Vector3d semiAxes(1.,2.,3.);
  for (const stk::math::Vector3d & rotation : {stk::math::Vector3d(0.,0.,0.5*M_PI), stk::math::Vector3d(0.5,0.25,0.1)})
  {
    Ellipsoid ellipsoid({center[0], center[1], center[2]}, {semiAxes[0], semiAxes[1], semiAxes[2]}, {rotation[0], rotation[1], rotation[2]}, 1);
    const Cuboid cuboid(center, 2.*semiAxes, rotation);
    Quaternion q;
    q.set_from_rotation_vector(rotation);
    for (int i=0; i<3; ++i)
    {
      stk::math::Vector3d axisEnd = stk::math::Vector3d::ZERO;
      axisEnd[i] = semiAxes[i];
      const stk::math::Vector3d end = center + q.rotate_3d_vector(axisEnd);
      EXPECT_NEAR(0., ellipsoid.point_signed_distance(end), 1.e-12) << "axis " << i << ", rotation " << rotation;
      EXPECT_NEAR(0., cuboid.point_signed_distance(end), 1.e-12) << "axis " << i << ", rotation " << rotation;
    }
  }
  Quaternion quarterTurnAboutZ;
  quarterTurnAboutZ.set_from_rotation_vector(stk::math::Vector3d(0.,0.,0.5*M_PI));
  const stk::math::Vector3d xTurned = quarterTurnAboutZ.rotate_3d_vector(stk::math::Vector3d(1.,0.,0.));
  EXPECT_NEAR(0., (xTurned-stk::math::Vector3d(0.,1.,0.)).length(), 1.e-12);
}

//--------------------------------------------------------------------------
// Random: each query gives the next value of a sequence set by the seed,
// uniform in [-amplitude, amplitude], wherever the query point is.
//--------------------------------------------------------------------------

namespace {
std::vector<double> next_random_values(const Random & random, const unsigned num)
{
  std::vector<double> values;
  for (unsigned i=0; i<num; ++i)
    values.push_back(random.point_signed_distance(stk::math::Vector3d(i, 0., 0.)));
  return values;
}
}

TEST(Random, sameSeed_sameSequence_otherSeed_otherSequence)
{
  const Random random(123);
  const Random sameSeed(123);
  const Random otherSeed(124);
  const std::vector<double> values = next_random_values(random, 10);
  EXPECT_EQ(values, next_random_values(sameSeed, 10));
  EXPECT_NE(values, next_random_values(otherSeed, 10));
}

TEST(Random, values_withinAmplitudeSetByTruncationLength)
{
  const BoundingBox anyBox(stk::math::Vector3d(-1.,-1.,-1.), stk::math::Vector3d(1.,1.,1.));
  Random random(7);
  for (const double amplitude : {1., 0.25})
  {
    if (amplitude != 1.)
      random.prepare_to_compute(0., anyBox, amplitude);
    const std::vector<double> values = next_random_values(random, 100);
    const auto [minValue, maxValue] = std::minmax_element(values.begin(), values.end());
    EXPECT_LE(-amplitude, *minValue);
    EXPECT_GE(amplitude, *maxValue);
    EXPECT_LT(*minValue, *maxValue);
    EXPECT_LT(0.5*amplitude, *maxValue - *minValue); // spread across the interval
  }
  // A truncation length of zero leaves the amplitude.
  random.prepare_to_compute(0., anyBox, 0.);
  for (const double value : next_random_values(random, 100))
  {
    EXPECT_GE(0.25, std::abs(value));
  }
}

TEST(Random, everywhere)
{
  const Random random(7);
  BoundingBox bbox;
  random.insert_into(bbox);
  EXPECT_TRUE(bbox.contains(BoundingBox::ENTIRE_DOMAIN));
  EXPECT_TRUE(random.does_intersect(BoundingBox(stk::math::Vector3d(100.,100.,100.), stk::math::Vector3d(101.,101.,101.))));
  EXPECT_FALSE(random.truncated_distance_may_have_wrong_sign());
}

//--------------------------------------------------------------------------
// LevelSet_String_Function: the expression of x, y, z and t, at the time
// given to prepare_to_compute (0 before that).
//--------------------------------------------------------------------------

TEST(LevelSet_String_Function, evaluatesTheExpression)
{
  const LevelSet_String_Function surface("x*x + y*y + z*z - 1.0");
  EXPECT_DOUBLE_EQ(-1., surface.point_signed_distance(stk::math::Vector3d(0.,0.,0.)));
  EXPECT_DOUBLE_EQ(2., surface.point_signed_distance(stk::math::Vector3d(1.,1.,1.)));
  EXPECT_DOUBLE_EQ(3., surface.point_signed_distance(stk::math::Vector3d(0.,-2.,0.)));
}

TEST(LevelSet_String_Function, expressionOfTime_atTheTimeGivenToPrepareToCompute)
{
  const BoundingBox anyBox(stk::math::Vector3d(-1.,-1.,-1.), stk::math::Vector3d(1.,1.,1.));
  LevelSet_String_Function surface("x - t");
  const stk::math::Vector3d x(3.,0.,0.);
  EXPECT_DOUBLE_EQ(3., surface.point_signed_distance(x));
  for (const double time : {2., 0.5})
  {
    surface.prepare_to_compute(time, anyBox, 0.);
    EXPECT_DOUBLE_EQ(3.-time, surface.point_signed_distance(x));
  }
}

TEST(LevelSet_String_Function, changesWithTime_onlyIfTheExpressionUsesTime)
{
  EXPECT_TRUE(LevelSet_String_Function("x - t").changes_with_time());
  EXPECT_FALSE(LevelSet_String_Function("x - 0.5").changes_with_time());
  EXPECT_FALSE(Sphere(stk::math::Vector3d::ZERO, 1.).changes_with_time());
}

TEST(LevelSet_String_Function, boundingBox_entireDomainUnlessSet)
{
  LevelSet_String_Function surface("x");
  BoundingBox bbox;
  surface.insert_into(bbox);
  EXPECT_TRUE(bbox.contains(BoundingBox::ENTIRE_DOMAIN));
  const BoundingBox farBox(stk::math::Vector3d(100.,100.,100.), stk::math::Vector3d(101.,101.,101.));
  EXPECT_TRUE(surface.does_intersect(farBox));

  const BoundingBox givenBox(stk::math::Vector3d(-1.,-1.,-1.), stk::math::Vector3d(1.,1.,1.));
  surface.set_bounding_box(givenBox);
  BoundingBox setBbox;
  surface.insert_into(setBbox);
  EXPECT_TRUE(setBbox.contains(givenBox));
  EXPECT_TRUE(givenBox.contains(setBbox));
  EXPECT_FALSE(surface.does_intersect(farBox));
  EXPECT_TRUE(surface.does_intersect(BoundingBox(stk::math::Vector3d(0.5,0.5,0.5), stk::math::Vector3d(2.,2.,2.))));
}

} // namespace krino
