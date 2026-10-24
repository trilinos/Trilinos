// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <Akri_IntersectionUtils.hpp>
#include <Akri_BoundingBox.hpp>
#include <Akri_Triangle.hpp>
#include <Akri_UnitTestUtils.hpp>

#include <stk_math/StkVector.hpp>

#include <array>
#include <vector>

namespace krino {

// Coverage for the free functions in Akri_IntersectionUtils:
//  - does_bounding_box_intersect_{tetrahedron,triangle_3d,edge_2d}: cheap
//    separating-axis tests (per axis-aligned component; ndim=3 for tet/tri,
//    ndim=2 for edge) that report a POSSIBLE overlap, not a true geometric one.
//  - append_tetrahedron_intersection_parametric_coordinates: appends a
//    candidate's parametric coordinates iff it lies within the tet.
//  - append_triangle_edge_intersection_parametric_coordinates: appends the
//    (r,s,0) parametric coordinates of an edge/triangle-plane crossing iff the
//    crossing lands inside the triangle.

//--------------------------------------------------------------------------
// does_bounding_box_intersect_tetrahedron
//--------------------------------------------------------------------------

TEST(IntersectionUtils, boundingBoxTetrahedron_overlapVsSeparated)
{
  const std::array<stk::math::Vector3d,4> tet{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0),
    stk::math::Vector3d(0,1,0), stk::math::Vector3d(0,0,1)};

  const BoundingBox overlapping(stk::math::Vector3d(0.1,0.1,0.1), stk::math::Vector3d(0.2,0.2,0.2));
  EXPECT_TRUE(does_bounding_box_intersect_tetrahedron(overlapping, tet));

  // Every tet node is strictly below the box minimum on each component.
  const BoundingBox separated(stk::math::Vector3d(2,2,2), stk::math::Vector3d(3,3,3));
  EXPECT_FALSE(does_bounding_box_intersect_tetrahedron(separated, tet));
}

TEST(IntersectionUtils, boundingBoxTetrahedron_separatedOnSingleAxisAbove)
{
  const std::array<stk::math::Vector3d,4> tet{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0),
    stk::math::Vector3d(0,1,0), stk::math::Vector3d(0,0,1)};

  // Overlaps in x and y, but every tet node is strictly above the box max in z.
  const BoundingBox separatedInZ(stk::math::Vector3d(-1,-1,-3), stk::math::Vector3d(1,1,-2));
  EXPECT_FALSE(does_bounding_box_intersect_tetrahedron(separatedInZ, tet));
}

//--------------------------------------------------------------------------
// does_bounding_box_intersect_triangle_3d
//--------------------------------------------------------------------------

TEST(IntersectionUtils, boundingBoxTriangle3d_overlapVsSeparated)
{
  const std::array<stk::math::Vector3d,3> tri{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(2,0,0), stk::math::Vector3d(0,2,0)};

  const BoundingBox overlapping(stk::math::Vector3d(-1,-1,-1), stk::math::Vector3d(1,1,1));
  EXPECT_TRUE(does_bounding_box_intersect_triangle_3d(overlapping, tri));

  const BoundingBox separated(stk::math::Vector3d(5,5,5), stk::math::Vector3d(6,6,6));
  EXPECT_FALSE(does_bounding_box_intersect_triangle_3d(separated, tri));
}

//--------------------------------------------------------------------------
// does_bounding_box_intersect_edge_2d (only x,y components are considered)
//--------------------------------------------------------------------------

TEST(IntersectionUtils, boundingBoxEdge2d_overlapVsSeparated)
{
  const std::array<stk::math::Vector3d,2> edge{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0)};

  const BoundingBox overlapping(stk::math::Vector3d(-1,-1,-1), stk::math::Vector3d(2,2,2));
  EXPECT_TRUE(does_bounding_box_intersect_edge_2d(overlapping, edge));

  const BoundingBox separated(stk::math::Vector3d(5,5,5), stk::math::Vector3d(6,6,6));
  EXPECT_FALSE(does_bounding_box_intersect_edge_2d(separated, edge));
}

TEST(IntersectionUtils, boundingBoxEdge2d_ignoresZComponent)
{
  // The edge is far from the box in z, but the 2d test only checks x and y.
  const std::array<stk::math::Vector3d,2> edge{
    stk::math::Vector3d(0,0,100), stk::math::Vector3d(1,0,100)};

  const BoundingBox overlappingInXY(stk::math::Vector3d(-1,-1,-1), stk::math::Vector3d(2,2,2));
  EXPECT_TRUE(does_bounding_box_intersect_edge_2d(overlappingInXY, edge));
}

//--------------------------------------------------------------------------
// append_tetrahedron_intersection_parametric_coordinates
//--------------------------------------------------------------------------

TEST(IntersectionUtils, appendTetIntersection_insideAppendsParametricCoords)
{
  // Reference unit tet: parametric coordinates of a point equal its (x,y,z).
  const std::array<stk::math::Vector3d,4> tet{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0),
    stk::math::Vector3d(0,1,0), stk::math::Vector3d(0,0,1)};

  std::vector<stk::math::Vector3d> intParamCoords;
  append_tetrahedron_intersection_parametric_coordinates(tet, stk::math::Vector3d(0.1,0.2,0.3), intParamCoords);

  ASSERT_EQ(1u, intParamCoords.size());
  expect_eq(stk::math::Vector3d(0.1,0.2,0.3), intParamCoords[0]);
}

TEST(IntersectionUtils, appendTetIntersection_outsideAppendsNothing)
{
  const std::array<stk::math::Vector3d,4> tet{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0),
    stk::math::Vector3d(0,1,0), stk::math::Vector3d(0,0,1)};

  std::vector<stk::math::Vector3d> intParamCoords;
  // Sum of parametric coordinates exceeds one (beyond the slanted face).
  append_tetrahedron_intersection_parametric_coordinates(tet, stk::math::Vector3d(0.5,0.5,0.5), intParamCoords);
  // Negative parametric coordinate (behind a face).
  append_tetrahedron_intersection_parametric_coordinates(tet, stk::math::Vector3d(-0.1,0.2,0.2), intParamCoords);

  EXPECT_TRUE(intParamCoords.empty());
}

//--------------------------------------------------------------------------
// append_triangle_edge_intersection_parametric_coordinates
//--------------------------------------------------------------------------

TEST(IntersectionUtils, appendTriangleEdgeIntersection_crossingInsideAppendsCoords)
{
  // Reference triangle in z=0; parametric coords of an in-plane point (px,py,0)
  // reduce exactly to (px,py,0).
  const std::array<stk::math::Vector3d,3> tri{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0)};
  const auto [triArea, triNormal] = CalcTriangle3<double>::area_and_normal(tri);

  std::vector<stk::math::Vector3d> intParamCoords;
  // Edge pierces the plane at (0.25,0.25,0), which is inside the triangle.
  append_triangle_edge_intersection_parametric_coordinates(tri, triArea, triNormal,
      stk::math::Vector3d(0.25,0.25,-1.), stk::math::Vector3d(0.25,0.25,1.), intParamCoords);

  ASSERT_EQ(1u, intParamCoords.size());
  expect_eq(stk::math::Vector3d(0.25,0.25,0.), intParamCoords[0]);
}

TEST(IntersectionUtils, appendTriangleEdgeIntersection_crossingOutsideAppendsNothing)
{
  const std::array<stk::math::Vector3d,3> tri{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0)};
  const auto [triArea, triNormal] = CalcTriangle3<double>::area_and_normal(tri);

  std::vector<stk::math::Vector3d> intParamCoords;
  // Crosses the plane at (0.8,0.8,0): in-plane but outside the triangle (sum>1).
  append_triangle_edge_intersection_parametric_coordinates(tri, triArea, triNormal,
      stk::math::Vector3d(0.8,0.8,-1.), stk::math::Vector3d(0.8,0.8,1.), intParamCoords);

  EXPECT_TRUE(intParamCoords.empty());
}

TEST(IntersectionUtils, appendTriangleEdgeIntersection_noPlaneCrossingAppendsNothing)
{
  const std::array<stk::math::Vector3d,3> tri{
    stk::math::Vector3d(0,0,0), stk::math::Vector3d(1,0,0), stk::math::Vector3d(0,1,0)};
  const auto [triArea, triNormal] = CalcTriangle3<double>::area_and_normal(tri);

  std::vector<stk::math::Vector3d> intParamCoords;
  // Both endpoints on the +z side of the plane: no sign change, no crossing.
  append_triangle_edge_intersection_parametric_coordinates(tri, triArea, triNormal,
      stk::math::Vector3d(0.25,0.25,1.), stk::math::Vector3d(0.25,0.25,2.), intParamCoords);

  EXPECT_TRUE(intParamCoords.empty());
}

}
