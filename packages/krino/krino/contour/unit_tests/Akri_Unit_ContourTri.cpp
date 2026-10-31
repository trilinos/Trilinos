// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_ContourTri.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_UnitTestUtils.hpp>

namespace krino {

// 2D analog of ContourTetFixture: contour a single right triangle whose nodes
// sit at (0,0), (1,0) and (0,1). ContourTri::append_facets_for_tri3 is pure
// geometry -- it takes the three vertex coordinates plus a per-node distance
// and appends the zero-level line segment(s), so no BulkData is required.
class ContourTriFixture : public ::testing::Test
{
public:
  std::array<stk::math::Vector3d,3> get_coords() const
  {
    return {{ stk::math::Vector3d(meshSpec.nodeLocs[0][0], meshSpec.nodeLocs[0][1], 0.0),
              stk::math::Vector3d(meshSpec.nodeLocs[1][0], meshSpec.nodeLocs[1][1], 0.0),
              stk::math::Vector3d(meshSpec.nodeLocs[2][0], meshSpec.nodeLocs[2][1], 0.0) }};
  }
  const Facet2d & get_facet(unsigned i) const { return facets.get_facets()[i]; }

protected:
  RightTri meshSpec;
  Faceted_Surface<Facet2d> facets;
};

TEST_F(ContourTriFixture, AppendFacetsTri3_AllNegative_AppendsNothing)
{
  const std::array<double,3> dist = {{-1.0, -2.0, -3.0}};
  ContourTri::append_facets_for_tri3(get_coords(), dist, facets);
  EXPECT_TRUE(facets.get_facets().empty());
}

TEST_F(ContourTriFixture, AppendFacetsTri3_AllPositive_AppendsNothing)
{
  const std::array<double,3> dist = {{1.0, 2.0, 3.0}};
  ContourTri::append_facets_for_tri3(get_coords(), dist, facets);
  EXPECT_TRUE(facets.get_facets().empty());
}

// One node positive, two negative (caseId 2): the interface crosses edges 0-1
// and 2-0. The segment is emitted from the edge-2-0 crossing to the edge-0-1
// crossing.
TEST_F(ContourTriFixture, AppendFacetsTri3_Case2_AppendsOneSegmentWithExpectedVertices)
{
  const std::array<double,3> dist = {{1.0, -1.0, -1.0}};
  ContourTri::append_facets_for_tri3(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  const stk::math::Vector3d x3{0.5, 0.0, 0.0}; // crossing on edge 0-1
  const stk::math::Vector3d x5{0.0, 0.5, 0.0}; // crossing on edge 2-0

  expect_eq(get_facet(0).facet_vertex(0), x5);
  expect_eq(get_facet(0).facet_vertex(1), x3);
}

// The complementary sign pattern (caseId 24) emits the same two crossings with
// the opposite winding.
TEST_F(ContourTriFixture, AppendFacetsTri3_Case24_AppendsOneSegmentWithOppositeWinding)
{
  const std::array<double,3> dist = {{-1.0, 1.0, 1.0}};
  ContourTri::append_facets_for_tri3(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  const stk::math::Vector3d x3{0.5, 0.0, 0.0}; // crossing on edge 0-1
  const stk::math::Vector3d x5{0.0, 0.5, 0.0}; // crossing on edge 2-0

  expect_eq(get_facet(0).facet_vertex(0), x3);
  expect_eq(get_facet(0).facet_vertex(1), x5);
}

// Two nodes exactly zero and the third negative (caseId 4): the interface lies
// along the 0-1 edge itself, so the segment connects the two on-contour nodes.
TEST_F(ContourTriFixture, AppendFacetsTri3_Case4_AppendsSegmentAlongZeroEdge)
{
  const std::array<double,3> dist = {{0.0, 0.0, -1.0}};
  ContourTri::append_facets_for_tri3(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  expect_eq(get_facet(0).facet_vertex(0), get_coords()[0]);
  expect_eq(get_facet(0).facet_vertex(1), get_coords()[1]);
}

// One node positive, one exactly zero, one negative (caseId 5): the segment
// runs from the edge-2-0 crossing to the on-contour node 1.
TEST_F(ContourTriFixture, AppendFacetsTri3_Case5_AppendsSegmentUsingZeroNode)
{
  const std::array<double,3> dist = {{1.0, 0.0, -1.0}};
  ContourTri::append_facets_for_tri3(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  const stk::math::Vector3d x5{0.0, 0.5, 0.0}; // crossing on edge 2-0

  expect_eq(get_facet(0).facet_vertex(0), x5);
  expect_eq(get_facet(0).facet_vertex(1), get_coords()[1]);
}

}
