// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_Segment.hpp>
#include <Akri_Triangle.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkVector.hpp>

#include <array>

namespace krino {

// CalcSegment3 / CalcTriangle3 are pure closest-point geometry kernels. The
// existing Akri_Unit_Geometry tests cover CalcTriangle3's normal/area helpers;
// these tests target the untested closest-point projection: for the segment,
// the three clamp regions of the parametric location; for the triangle, the
// seven Voronoi regions (three vertices, three edges, the interior face) that
// closest_point_projection classifies.

using Seg = CalcSegment3<double>;
using Tri = CalcTriangle3<double>;

//--------------------------------------------------------------------------
// CalcSegment3: closest point / parametric location on a segment.
// Reference segment p0=(0,0,0), p1=(2,0,0) lies along +x.
//--------------------------------------------------------------------------

TEST(CalcSegment3, closestParametricLocation_beforeStart_clampsToZero)
{
  const stk::math::Vector3d p0(0,0,0), p1(2,0,0);
  EXPECT_DOUBLE_EQ(0.0, Seg::closest_parametric_location(p0, p1, stk::math::Vector3d(-1,5,0)));
}

TEST(CalcSegment3, closestParametricLocation_pastEnd_clampsToOne)
{
  const stk::math::Vector3d p0(0,0,0), p1(2,0,0);
  EXPECT_DOUBLE_EQ(1.0, Seg::closest_parametric_location(p0, p1, stk::math::Vector3d(3,5,0)));
}

TEST(CalcSegment3, closestParametricLocation_interior_projectsOntoSegment)
{
  const stk::math::Vector3d p0(0,0,0), p1(2,0,0);
  // (0.5,3,0) projects to x=0.5 -> param 0.5/2 = 0.25.
  EXPECT_DOUBLE_EQ(0.25, Seg::closest_parametric_location(p0, p1, stk::math::Vector3d(0.5,3,0)));
}

TEST(CalcSegment3, closestParametricLocation_degenerateZeroLength_returnsZero)
{
  const stk::math::Vector3d p0(1,1,1);
  EXPECT_DOUBLE_EQ(0.0, Seg::closest_parametric_location(p0, p0, stk::math::Vector3d(5,5,5)));
}

TEST(CalcSegment3, closestPoint_interiorAndClamped)
{
  const stk::math::Vector3d p0(0,0,0), p1(2,0,0);
  expect_eq(stk::math::Vector3d(0.5,0,0), Seg::closest_point(p0, p1, stk::math::Vector3d(0.5,3,0)));
  expect_eq(p0, Seg::closest_point(p0, p1, stk::math::Vector3d(-1,5,0)));
  expect_eq(p1, Seg::closest_point(p0, p1, stk::math::Vector3d(3,5,0)));
}

TEST(CalcSegment3, closestPointAndParametricCoord_reportsBoth)
{
  const stk::math::Vector3d p0(0,0,0), p1(2,0,0);
  stk::math::Vector3d closestPt;
  double paramLoc;
  Seg::closest_point_and_parametric_coord(p0, p1, stk::math::Vector3d(0.5,3,0), closestPt, paramLoc);
  expect_eq(stk::math::Vector3d(0.5,0,0), closestPt);
  EXPECT_DOUBLE_EQ(0.25, paramLoc);
}

TEST(CalcSegment3, distanceSquared_perpendicularDrop)
{
  const stk::math::Vector3d p0(0,0,0), p1(2,0,0);
  // Interior drop of length 3 -> distance^2 = 9.
  EXPECT_DOUBLE_EQ(9.0, Seg::distance_squared(p0, p1, stk::math::Vector3d(0.5,3,0)));
  // Off the end: closest is p1, distance^2 = 1^2 + 5^2 = 26.
  EXPECT_DOUBLE_EQ(26.0, Seg::distance_squared(p0, p1, stk::math::Vector3d(3,5,0)));
}

TEST(CalcSegment3, isProjectionInsideSegment_interiorTrueEndsFalse)
{
  const stk::math::Vector3d p0(0,0,0), p1(2,0,0);
  EXPECT_TRUE(Seg::is_projection_of_point_inside_segment(p0, p1, stk::math::Vector3d(0.5,3,0)));
  EXPECT_FALSE(Seg::is_projection_of_point_inside_segment(p0, p1, stk::math::Vector3d(-1,0,0)));
  EXPECT_FALSE(Seg::is_projection_of_point_inside_segment(p0, p1, stk::math::Vector3d(3,0,0)));
}

TEST(CalcSegment3, length_matchesEuclidean)
{
  const stk::math::Vector3d p0(0,0,0), p1(3,4,0);
  EXPECT_DOUBLE_EQ(5.0, Seg::length(p0, p1));
  EXPECT_DOUBLE_EQ(25.0, Seg::length_squared(p0, p1));
}

//--------------------------------------------------------------------------
// CalcTriangle3::closest_point_projection: the seven Voronoi regions.
// Reference triangle a=(0,0,0), b=(1,0,0), c=(0,1,0) in the z=0 plane;
// parametric coords (r,s) with a->(0,0), b->(1,0), c->(0,1).
//--------------------------------------------------------------------------

namespace {
const stk::math::Vector3d a(0,0,0), b(1,0,0), c(0,1,0);

void expect_closest(const stk::math::Vector3d & queryPt,
    const stk::math::Vector3d & goldClosest,
    const stk::math::Vector3d & goldParam)
{
  stk::math::Vector3d closestPt;
  stk::math::Vec<double,2> paramPt;
  Tri::closest_point_and_parametric_coords(a, b, c, queryPt, closestPt, paramPt);
  expect_eq(goldClosest, closestPt);
  EXPECT_NEAR(goldParam[0], paramPt[0], 1.e-9);
  EXPECT_NEAR(goldParam[1], paramPt[1], 1.e-9);
}
}

TEST(CalcTriangle3, closestPoint_vertexRegionA)
{
  expect_closest(stk::math::Vector3d(-1,-1,0), a, stk::math::Vector3d(0,0,0));
}

TEST(CalcTriangle3, closestPoint_vertexRegionB)
{
  expect_closest(stk::math::Vector3d(2,-1,0), b, stk::math::Vector3d(1,0,0));
}

TEST(CalcTriangle3, closestPoint_vertexRegionC)
{
  expect_closest(stk::math::Vector3d(-1,2,0), c, stk::math::Vector3d(0,1,0));
}

TEST(CalcTriangle3, closestPoint_edgeAB)
{
  // Drops perpendicular onto edge AB (the y=0 edge) at x=0.5.
  expect_closest(stk::math::Vector3d(0.5,-1,0), stk::math::Vector3d(0.5,0,0), stk::math::Vector3d(0.5,0,0));
}

TEST(CalcTriangle3, closestPoint_edgeAC)
{
  // Drops perpendicular onto edge AC (the x=0 edge) at y=0.5.
  expect_closest(stk::math::Vector3d(-1,0.5,0), stk::math::Vector3d(0,0.5,0), stk::math::Vector3d(0,0.5,0));
}

TEST(CalcTriangle3, closestPoint_edgeBC)
{
  // (1,1,0) drops onto the hypotenuse midpoint (0.5,0.5,0).
  expect_closest(stk::math::Vector3d(1,1,0), stk::math::Vector3d(0.5,0.5,0), stk::math::Vector3d(0.5,0.5,0));
}

TEST(CalcTriangle3, closestPoint_faceInterior_projectsThroughPlane)
{
  // A point lifted off the triangle plane projects straight down into the face.
  expect_closest(stk::math::Vector3d(0.25,0.25,5), stk::math::Vector3d(0.25,0.25,0), stk::math::Vector3d(0.25,0.25,0));
}

TEST(CalcTriangle3, distanceSquared_faceInterior)
{
  // (0.25,0.25,5) is 5 above the face -> distance^2 = 25.
  EXPECT_DOUBLE_EQ(25.0, Tri::distance_squared(a, b, c, stk::math::Vector3d(0.25,0.25,5)));
}

TEST(CalcTriangle3, isProjectionInsideTriangle_faceTrueEdgeVertexFalse)
{
  EXPECT_TRUE(Tri::is_projection_of_point_inside_triangle(a, b, c, stk::math::Vector3d(0.25,0.25,5)));
  EXPECT_FALSE(Tri::is_projection_of_point_inside_triangle(a, b, c, stk::math::Vector3d(0.5,-1,0)));  // edge AB
  EXPECT_FALSE(Tri::is_projection_of_point_inside_triangle(a, b, c, stk::math::Vector3d(-1,-1,0)));   // vertex A
}

TEST(CalcTriangle3, parametricToRealCoords_roundTripsWithProjection)
{
  const stk::math::Vec<double,2> param(0.25, 0.25);
  expect_eq(stk::math::Vector3d(0.25,0.25,0), Tri::parametric_to_real_coords(a, b, c, param));
}

}
