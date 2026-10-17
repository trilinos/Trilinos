// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_Plane_Intersections.hpp>
#include <Akri_Triangle.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkPlane.hpp>
#include <stk_math/StkVector.hpp>

#include <array>

namespace krino {

// Akri_Plane_Intersections is pure geometry: it solves for the point common to
// three planes (3D) or two planes plus a coordinate plane (2D / tet-side), and
// clips against the reference tet/tri. Planes are stk::math::Plane3d, built here
// from a normal and a point on the plane (the normal is unitized internally).

namespace {
// The plane {x = c}: normal +x, passing through (c,0,0).
stk::math::Plane3d plane_x(const double c) { return stk::math::Plane3d(stk::math::Vector3d(1,0,0), stk::math::Vector3d(c,0,0)); }
stk::math::Plane3d plane_y(const double c) { return stk::math::Plane3d(stk::math::Vector3d(0,1,0), stk::math::Vector3d(0,c,0)); }
stk::math::Plane3d plane_z(const double c) { return stk::math::Plane3d(stk::math::Vector3d(0,0,1), stk::math::Vector3d(0,0,c)); }
}

//--------------------------------------------------------------------------
// find_intersection_of_three_planes
//--------------------------------------------------------------------------

// Three axis-aligned planes x=a, y=b, z=c meet at the single point (a,b,c).
TEST(PlaneIntersections, threePlanes_axisAligned_meetAtExpectedPoint)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_three_planes(plane_x(0.2), plane_y(0.3), plane_z(0.4), pt);
  ASSERT_TRUE(found);
  expect_eq(stk::math::Vector3d(0.2, 0.3, 0.4), pt);
}

// Two parallel planes make the system singular (det == 0): no unique point.
TEST(PlaneIntersections, threePlanes_twoParallel_returnsFalse)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_three_planes(plane_x(0.2), plane_x(0.5), plane_y(0.3), pt);
  EXPECT_FALSE(found);
}

//--------------------------------------------------------------------------
// find_intersection_of_three_planes_within_tet: the three-plane solve, then
// clipped to the reference tet (x>=0, y>=0, z>=0, x+y+z<=1).
//--------------------------------------------------------------------------

TEST(PlaneIntersections, threePlanesWithinTet_insidePoint_returnsTrue)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_three_planes_within_tet(plane_x(0.2), plane_y(0.2), plane_z(0.2), pt);
  ASSERT_TRUE(found);
  expect_eq(stk::math::Vector3d(0.2, 0.2, 0.2), pt);
}

TEST(PlaneIntersections, threePlanesWithinTet_outsidePoint_returnsFalse)
{
  stk::math::Vector3d pt;
  // (0.6,0.6,0.6): x+y+z = 1.8 > 1, so outside the tet even though the planes meet.
  const bool found = find_intersection_of_three_planes_within_tet(plane_x(0.6), plane_y(0.6), plane_z(0.6), pt);
  EXPECT_FALSE(found);
}

//--------------------------------------------------------------------------
// find_intersection_of_two_planes_and_side_of_tet: intersect two planes with the
// plane of a given tet side, and clip to the tet. Sides 0/2/3 are the coordinate
// planes y=0 / x=0 / z=0; side 1 is the slanted face x+y+z=1.
//--------------------------------------------------------------------------

// Side 0 is the y=0 face. Intersect x=0.3 and z=0.4 with it -> (0.3, 0, 0.4).
TEST(PlaneIntersections, twoPlanesAndSideOfTet_side0_yZeroFace)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_two_planes_and_side_of_tet(0, plane_x(0.3), plane_z(0.4), pt);
  ASSERT_TRUE(found);
  expect_eq(stk::math::Vector3d(0.3, 0.0, 0.4), pt);
}

// Side 2 is the x=0 face. Intersect y=0.3 and z=0.4 with it -> (0, 0.3, 0.4).
TEST(PlaneIntersections, twoPlanesAndSideOfTet_side2_xZeroFace)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_two_planes_and_side_of_tet(2, plane_y(0.3), plane_z(0.4), pt);
  ASSERT_TRUE(found);
  expect_eq(stk::math::Vector3d(0.0, 0.3, 0.4), pt);
}

// Side 3 is the z=0 face. Intersect x=0.3 and y=0.4 with it -> (0.3, 0.4, 0).
TEST(PlaneIntersections, twoPlanesAndSideOfTet_side3_zZeroFace)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_two_planes_and_side_of_tet(3, plane_x(0.3), plane_y(0.4), pt);
  ASSERT_TRUE(found);
  expect_eq(stk::math::Vector3d(0.3, 0.4, 0.0), pt);
}

// Side 1 is the slanted face x+y+z=1. Intersect x=0.25 and y=0.25 with it ->
// (0.25, 0.25, 0.5), which lies on the face and inside the tet.
TEST(PlaneIntersections, twoPlanesAndSideOfTet_side1_slantedFace)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_two_planes_and_side_of_tet(1, plane_x(0.25), plane_y(0.25), pt);
  ASSERT_TRUE(found);
  expect_eq(stk::math::Vector3d(0.25, 0.25, 0.5), pt);
}

// The intersection point lands outside the tet -> clipped away, returns false.
TEST(PlaneIntersections, twoPlanesAndSideOfTet_outsideTet_returnsFalse)
{
  stk::math::Vector3d pt;
  // On the y=0 face but x=1.3 pushes it past x+y+z<=1.
  const bool found = find_intersection_of_two_planes_and_side_of_tet(0, plane_x(1.3), plane_z(0.1), pt);
  EXPECT_FALSE(found);
}

//--------------------------------------------------------------------------
// find_intersection_of_two_2D_planes: two planes whose normals lie in z=0,
// intersected in the z=0 plane -> the point where the two lines cross.
//--------------------------------------------------------------------------

TEST(PlaneIntersections, two2DPlanes_perpendicularLines_crossAtExpectedPoint)
{
  stk::math::Vector3d pt;
  // Lines x=0.3 and y=0.7 in the z=0 plane cross at (0.3, 0.7, 0).
  const bool found = find_intersection_of_two_2D_planes(plane_x(0.3), plane_y(0.7), pt);
  ASSERT_TRUE(found);
  expect_eq(stk::math::Vector3d(0.3, 0.7, 0.0), pt);
}

TEST(PlaneIntersections, two2DPlanes_parallelLines_returnsFalse)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_two_2D_planes(plane_x(0.3), plane_x(0.6), pt);
  EXPECT_FALSE(found);
}

//--------------------------------------------------------------------------
// find_intersection_of_two_2D_planes_within_tri: as above, clipped to the
// reference triangle (x>=0, y>=0, x+y<=1).
//--------------------------------------------------------------------------

TEST(PlaneIntersections, two2DPlanesWithinTri_insidePoint_returnsTrue)
{
  stk::math::Vector3d pt;
  const bool found = find_intersection_of_two_2D_planes_within_tri(plane_x(0.25), plane_y(0.25), pt);
  ASSERT_TRUE(found);
  expect_eq(stk::math::Vector3d(0.25, 0.25, 0.0), pt);
}

TEST(PlaneIntersections, two2DPlanesWithinTri_outsidePoint_returnsFalse)
{
  stk::math::Vector3d pt;
  // (0.7,0.7): x+y = 1.4 > 1, outside the triangle.
  const bool found = find_intersection_of_two_2D_planes_within_tri(plane_x(0.7), plane_y(0.7), pt);
  EXPECT_FALSE(found);
}

//--------------------------------------------------------------------------
// triangle_parametric_coordinates_of_projected_point: parametric (r,s) of the
// orthogonal projection of a point onto a triangle's plane. Vertex 0 -> (0,0),
// vertex 1 -> (1,0), vertex 2 -> (0,1).
//--------------------------------------------------------------------------

TEST(PlaneIntersections, triangleParametricCoords_verticesMapToUnitCorners)
{
  const std::array<stk::math::Vector3d,3> tri = {{ {0,0,0}, {2,0,0}, {0,4,0} }};

  expect_eq(stk::math::Vector3d(0,0,0), triangle_parametric_coordinates_of_projected_point(tri, tri[0]));
  expect_eq(stk::math::Vector3d(1,0,0), triangle_parametric_coordinates_of_projected_point(tri, tri[1]));
  expect_eq(stk::math::Vector3d(0,1,0), triangle_parametric_coordinates_of_projected_point(tri, tri[2]));
}

TEST(PlaneIntersections, triangleParametricCoords_centroidAndProjection)
{
  const std::array<stk::math::Vector3d,3> tri = {{ {0,0,0}, {2,0,0}, {0,4,0} }};

  const stk::math::Vector3d centroid = 1./3.*(tri[0]+tri[1]+tri[2]);
  expect_eq(stk::math::Vector3d(1./3., 1./3., 0.0), triangle_parametric_coordinates_of_projected_point(tri, centroid));

  // A point lifted off the triangle plane projects to the same parametric coords
  // as its in-plane footprint: (1,2,0) is the midpoint-ish interior point.
  const stk::math::Vector3d lifted = stk::math::Vector3d(1.0, 2.0, 5.0);
  expect_eq(stk::math::Vector3d(0.5, 0.5, 0.0), triangle_parametric_coordinates_of_projected_point(tri, lifted));
}

//--------------------------------------------------------------------------
// within_tri_bounds: barycentric point is inside the reference triangle.
//--------------------------------------------------------------------------

TEST(PlaneIntersections, withinTriBounds_insideAndOnBoundary_true)
{
  EXPECT_TRUE(within_tri_bounds(stk::math::Vector3d(0.25, 0.25, 0.0)));
  EXPECT_TRUE(within_tri_bounds(stk::math::Vector3d(0.0, 0.0, 0.0)));  // vertex
  EXPECT_TRUE(within_tri_bounds(stk::math::Vector3d(0.5, 0.5, 0.0)));  // on hypotenuse
}

TEST(PlaneIntersections, withinTriBounds_outside_false)
{
  EXPECT_FALSE(within_tri_bounds(stk::math::Vector3d(0.7, 0.7, 0.0)));   // r+s > 1
  EXPECT_FALSE(within_tri_bounds(stk::math::Vector3d(-0.1, 0.5, 0.0)));  // r < 0
}

//==========================================================================
// The suites below drive the same intersection functions but from the mesh
// perspective: planes built from the faces of the reference tet / tri (three
// points) rather than the axis-aligned analytic planes above. They were moved
// here from the krino_lib Akri_Unit_Geometry monolith, which tested geometry-
// library functions. The two styles are complementary coverage, not
// duplicates: these exercise the four tet faces / three tri edges as a set.
//==========================================================================

namespace {
void expect_intersection(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2, const stk::math::Vector3d & goldIntersection)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_three_planes(plane0, plane1, plane2, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_intersection_with_side_of_tet(const int side, const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Vector3d & goldIntersection)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_two_planes_and_side_of_tet(side, plane0, plane1, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_no_intersection(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_three_planes(plane0, plane1, plane2, intersectionPoint));
}

void expect_no_intersection_within_tet(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Plane3d & plane2)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_three_planes_within_tet(plane0, plane1, plane2, intersectionPoint));
}

void expect_no_intersection_with_side_of_tet(const int side, const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_planes_and_side_of_tet(side, plane0, plane1, intersectionPoint));
}
}

struct IntersectPlanes : public ::testing::Test
{
  const stk::math::Vector3d pt0{0., 0., 0.};
  const stk::math::Vector3d pt1{1., 0., 0.};
  const stk::math::Vector3d pt2{0., 1., 0.};
  const stk::math::Vector3d pt3{0., 0., 1.};
  const stk::math::Plane3d plane0{pt0, pt1, pt3};
  const stk::math::Plane3d plane1{pt1, pt2, pt3};
  const stk::math::Plane3d plane2{pt2, pt0, pt3};
  const stk::math::Plane3d plane3{pt0, pt2, pt1};

  const stk::math::Vector3d offsetPt0{-0.1, -0.1, -0.1};
  const stk::math::Vector3d offsetPt1{1.1, -0.1, -0.1};
  const stk::math::Vector3d offsetPt2{-0.1, 1.1, -0.1};
  const stk::math::Vector3d offsetPt3{-0.1, -0.1, 1.1};
  const stk::math::Plane3d offsetPlane0{offsetPt0, offsetPt1, offsetPt3};
  const stk::math::Plane3d offsetPlane1{offsetPt1, offsetPt2, offsetPt3};
  const stk::math::Plane3d offsetPlane2{offsetPt2, offsetPt0, offsetPt3};
  const stk::math::Plane3d offsetPlane3{offsetPt0, offsetPt2, offsetPt1};

};

TEST_F(IntersectPlanes, given3Planes_FindCorrectIntersection)
{
  expect_intersection(plane0, plane1, plane2, stk::math::Vector3d{0., 0., 1.});
  expect_intersection(plane0, plane1, plane3, stk::math::Vector3d{1., 0., 0.});
  expect_intersection(plane1, plane2, plane3, stk::math::Vector3d{0., 1., 0.});
  expect_intersection(plane0, plane2, plane3, stk::math::Vector3d{0., 0., 0.});

  expect_no_intersection(plane0, plane0, plane1);
}

TEST_F(IntersectPlanes, given3PlanesThatDoNotIntersectWithinTetBounds_FindNoIntersection)
{
  expect_no_intersection_within_tet(offsetPlane0, offsetPlane1, offsetPlane2);
  expect_no_intersection_within_tet(offsetPlane0, offsetPlane1, offsetPlane3);
  expect_no_intersection_within_tet(offsetPlane1, offsetPlane2, offsetPlane3);
  expect_no_intersection_within_tet(offsetPlane0, offsetPlane2, offsetPlane3);
}

TEST_F(IntersectPlanes, given2PlanesThatIntersectSideOfTet_FindCorrectIntersection)
{
  expect_intersection_with_side_of_tet(0, plane1, plane2, stk::math::Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(0, plane1, plane3, stk::math::Vector3d{1., 0., 0.});
  expect_intersection_with_side_of_tet(0, plane2, plane3, stk::math::Vector3d{0., 0., 0.});

  expect_intersection_with_side_of_tet(1, plane2, plane3, stk::math::Vector3d{0., 1., 0.});
  expect_intersection_with_side_of_tet(1, plane0, plane2, stk::math::Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(1, plane0, plane3, stk::math::Vector3d{1., 0., 0.});

  expect_intersection_with_side_of_tet(2, plane0, plane1, stk::math::Vector3d{0., 0., 1.});
  expect_intersection_with_side_of_tet(2, plane0, plane3, stk::math::Vector3d{0., 0., 0.});
  expect_intersection_with_side_of_tet(2, plane1, plane3, stk::math::Vector3d{0., 1., 0.});

  expect_intersection_with_side_of_tet(3, plane0, plane1, stk::math::Vector3d{1., 0., 0.});
  expect_intersection_with_side_of_tet(3, plane0, plane2, stk::math::Vector3d{0., 0., 0.});
  expect_intersection_with_side_of_tet(3, plane1, plane2, stk::math::Vector3d{0., 1., 0.});
}

TEST_F(IntersectPlanes, given2PlanesThatDoNotIntersectSideOfTetAtAll_FindNoIntersection)
{
  expect_no_intersection_with_side_of_tet(0, plane0, plane2);
  expect_no_intersection_with_side_of_tet(1, plane1, plane2);
  expect_no_intersection_with_side_of_tet(2, plane2, plane1);
  expect_no_intersection_with_side_of_tet(3, plane3, plane2);
}

TEST_F(IntersectPlanes, given2PlanesThatDoNotIntersectSideOfTetWithinTetBounds_FindNoIntersection)
{
  expect_no_intersection_with_side_of_tet(0, offsetPlane1, offsetPlane2);
  expect_no_intersection_with_side_of_tet(0, offsetPlane1, offsetPlane3);
  expect_no_intersection_with_side_of_tet(0, offsetPlane2, offsetPlane3);

  expect_no_intersection_with_side_of_tet(1, offsetPlane2, offsetPlane3);
  expect_no_intersection_with_side_of_tet(1, offsetPlane0, offsetPlane2);
  expect_no_intersection_with_side_of_tet(1, offsetPlane0, offsetPlane3);

  expect_no_intersection_with_side_of_tet(2, offsetPlane0, offsetPlane1);
  expect_no_intersection_with_side_of_tet(2, offsetPlane0, offsetPlane3);
  expect_no_intersection_with_side_of_tet(2, offsetPlane1, offsetPlane3);

  expect_no_intersection_with_side_of_tet(3, offsetPlane0, offsetPlane1);
  expect_no_intersection_with_side_of_tet(3, offsetPlane0, offsetPlane2);
  expect_no_intersection_with_side_of_tet(3, offsetPlane1, offsetPlane2);
}

namespace {
stk::math::Vector3d compute_2d_plane_direction(const stk::math::Vector3d & pt0, const stk::math::Vector3d & pt1)
{
  return stk::math::Vector3d(pt1[1]-pt0[1],pt0[0]-pt1[0],0.);
}

stk::math::Plane3d build_2d_plane(const stk::math::Vector3d & pt0, const stk::math::Vector3d & pt1)
{
  return stk::math::Plane3d(compute_2d_plane_direction(pt0,pt1),pt0);
}
}

struct Intersect2DPlanes : public ::testing::Test
{
  Intersect2DPlanes()
  : plane0{build_2d_plane(pt0, pt1)},
    plane1{build_2d_plane(pt1, pt2)},
    plane2{build_2d_plane(pt0, pt2)},
    offsetPlane0{build_2d_plane(offsetPt0, offsetPt1)},
    offsetPlane1{build_2d_plane(offsetPt1, offsetPt2)},
    offsetPlane2{build_2d_plane(offsetPt0, offsetPt2)}
  {
  }

  const stk::math::Vector3d pt0{0., 0., 0.};
  const stk::math::Vector3d pt1{1., 0., 0.};
  const stk::math::Vector3d pt2{0., 1., 0.};

  const stk::math::Vector3d offsetPt0{-0.1, -0.1, 0.};
  const stk::math::Vector3d offsetPt1{1.1, -0.1, 0.};
  const stk::math::Vector3d offsetPt2{-0.1, 1.1, 0.};

  stk::math::Plane3d plane0;
  stk::math::Plane3d plane1;
  stk::math::Plane3d plane2;

  stk::math::Plane3d offsetPlane0;
  stk::math::Plane3d offsetPlane1;
  stk::math::Plane3d offsetPlane2;
};

namespace {
void expect_2d_intersection(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1, const stk::math::Vector3d & goldIntersection)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_TRUE(find_intersection_of_two_2D_planes(plane0, plane1, intersectionPoint));
  expect_eq(goldIntersection, intersectionPoint);
}

void expect_no_2d_intersection_within_tri(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_2D_planes_within_tri(plane0, plane1, intersectionPoint));
}

void expect_no_2d_intersection(const stk::math::Plane3d & plane0, const stk::math::Plane3d & plane1)
{
  stk::math::Vector3d intersectionPoint;
  EXPECT_FALSE(find_intersection_of_two_2D_planes(plane0, plane1, intersectionPoint));
}
}

TEST_F(Intersect2DPlanes, given2Plane2Ds_FindCorrectIntersection)
{
  expect_2d_intersection(plane0, plane1, stk::math::Vector3d{1., 0., 0.});
  expect_2d_intersection(plane1, plane2, stk::math::Vector3d{0., 1., 0.});
  expect_2d_intersection(plane0, plane2, stk::math::Vector3d{0., 0., 0.});

  expect_no_2d_intersection(plane0, plane0);
}

TEST_F(Intersect2DPlanes, given2Plane2DsThatDoNotIntersectWithinTriBounds_FindNoIntersection)
{
  expect_no_2d_intersection_within_tri(offsetPlane0, offsetPlane1);
  expect_no_2d_intersection_within_tri(offsetPlane1, offsetPlane2);
  expect_no_2d_intersection_within_tri(offsetPlane0, offsetPlane2);
}

TEST(ProjectionOf3DPointsIntoTriangle, givenTriangleCheckProjectionOfPointOnAndOffPlane)
{
  const stk::math::Vector3d pt0{0., 0., 0.};
  const stk::math::Vector3d pt1{1., 0., 0.};
  const stk::math::Vector3d pt2{0., 1., 1.};
  const std::array<stk::math::Vector3d,3> triCoords{{ pt0, pt1, pt2 }};

  expect_eq(stk::math::Vector3d{0.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt0));
  expect_eq(stk::math::Vector3d{1.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt1));
  expect_eq(stk::math::Vector3d{0.,1.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt2));

  const stk::math::Vector3d normal = CalcTriangle3<double>::normal(pt0, pt1, pt2);
  expect_eq(stk::math::Vector3d{0.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt0+normal));
  expect_eq(stk::math::Vector3d{1.,0.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt1+normal));
  expect_eq(stk::math::Vector3d{0.,1.,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, pt2+normal));

  const double oneThird = 1./3.;
  const stk::math::Vector3d midPt = oneThird*(pt0+pt1+pt2);
  expect_eq(stk::math::Vector3d{oneThird,oneThird,0.}, triangle_parametric_coordinates_of_projected_point(triCoords, midPt+normal));
}

}
