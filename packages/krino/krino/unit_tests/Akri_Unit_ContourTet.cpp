// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>
#include <Akri_ContourTet.hpp>
#include <Akri_MeshSpecs.hpp>
#include <Akri_Faceted_Surface.hpp>
#include <Akri_UnitTestUtils.hpp>

namespace krino {


class ContourTetFixture : public ::testing::Test
{
public:
  std::array<stk::math::Vector3d,4> get_coords() const { return {{ meshSpec.nodeLocs[0], meshSpec.nodeLocs[1], meshSpec.nodeLocs[2], meshSpec.nodeLocs[3] }}; }
  const Facet3d & get_facet(unsigned i) const { return facets.get_facets()[i]; }

protected:
  RightTet meshSpec;
  Faceted_Surface<Facet3d> facets;
};


TEST_F(ContourTetFixture, AppendFacetsTet4_AllNegative_AppendsNothing)
{
  const std::array<double,4> dist = {{-1.0, -2.0, -3.0, -4.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);
  EXPECT_TRUE(facets.get_facets().empty());
}

TEST_F(ContourTetFixture, AppendFacetsTet4_AllPositive_AppendsNothing)
{
  const std::array<double,4> dist = {{1.0, 2.0, 3.0, 4.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);
  EXPECT_TRUE(facets.get_facets().empty());
}

TEST_F(ContourTetFixture, AppendFacetsTet4_Case2_AppendsOneTriangleWithExpectedVertices)
{
  const std::array<double,4> dist = {{1.0, -1.0, -1.0, -1.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  const stk::math::Vector3d x4{0.5, 0.0, 0.0};
  const stk::math::Vector3d x6{0.0, 0.5, 0.0};
  const stk::math::Vector3d x7{0.0, 0.0, 0.5};

  expect_eq(get_facet(0).facet_vertex(0), x6);
  expect_eq(get_facet(0).facet_vertex(1), x4);
  expect_eq(get_facet(0).facet_vertex(2), x7);
}

TEST_F(ContourTetFixture, AppendFacetsTet4_Case78_AppendsOneTriangleWithOppositeWinding)
{
  const std::array<double,4> dist = {{-1.0, 1.0, 1.0, 1.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  const stk::math::Vector3d x4{0.5, 0.0, 0.0};
  const stk::math::Vector3d x6{0.0, 0.5, 0.0};
  const stk::math::Vector3d x7{0.0, 0.0, 0.5};

  expect_eq(get_facet(0).facet_vertex(0), x4);
  expect_eq(get_facet(0).facet_vertex(1), x6);
  expect_eq(get_facet(0).facet_vertex(2), x7);
}

TEST_F(ContourTetFixture, AppendFacetsTet4_Case5_AppendsTriangleUsingZeroNode)
{
  const std::array<double,4> dist = {{1.0, 0.0, -1.0, -1.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  const stk::math::Vector3d x6{0.0, 0.5, 0.0};
  const stk::math::Vector3d x7{0.0, 0.0, 0.5};

  expect_eq(get_facet(0).facet_vertex(0), x6);
  expect_eq(get_facet(0).facet_vertex(1), get_coords()[1]);
  expect_eq(get_facet(0).facet_vertex(2), x7);
}

TEST_F(ContourTetFixture, AppendFacetsTet4_Case75_AppendsTriangleUsingZeroNodeOppositeWinding)
{
  const std::array<double,4> dist = {{-1.0, 0.0, 1.0, 1.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  const stk::math::Vector3d x6{0.0, 0.5, 0.0};
  const stk::math::Vector3d x7{0.0, 0.0, 0.5};

  expect_eq(get_facet(0).facet_vertex(0), x7);
  expect_eq(get_facet(0).facet_vertex(1), get_coords()[1]);
  expect_eq(get_facet(0).facet_vertex(2), x6);
}

TEST_F(ContourTetFixture, AppendFacetsTet4_Case14_AppendsTriangleUsingTwoZeroNodes)
{
  const std::array<double,4> dist = {{1.0, 0.0, 0.0, -1.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);

  ASSERT_EQ(facets.size(), 1u);

  const stk::math::Vector3d x7{0.0, 0.0, 0.5};

  expect_eq(get_facet(0).facet_vertex(0), get_coords()[1]);
  expect_eq(get_facet(0).facet_vertex(1), x7);
  expect_eq(get_facet(0).facet_vertex(2), get_coords()[2]);
}

TEST_F(ContourTetFixture, AppendFacetsTet4_Case8_Face4False_AppendsExpectedTwoTriangles)
{
  const std::array<double,4> dist = {{1.0, 1.0, -1.0, -1.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);

  // will_cutting_quad_from_0to2_cut_largest_angle(...) returns false.

  ASSERT_EQ(facets.size(), 2u);

  const stk::math::Vector3d x5{0.5, 0.5, 0.0}; // edge 1-2
  const stk::math::Vector3d x6{0.0, 0.5, 0.0}; // edge 0-2
  const stk::math::Vector3d x7{0.0, 0.0, 0.5}; // edge 0-3
  const stk::math::Vector3d x8{0.5, 0.0, 0.5}; // edge 1-3

  expect_eq(get_facet(0).facet_vertex(0), x5);
  expect_eq(get_facet(0).facet_vertex(1), x8);
  expect_eq(get_facet(0).facet_vertex(2), x7);

  expect_eq(get_facet(1).facet_vertex(0), x5);
  expect_eq(get_facet(1).facet_vertex(1), x7);
  expect_eq(get_facet(1).facet_vertex(2), x6);
}

TEST_F(ContourTetFixture, AppendFacetsTet4_Case8_Face4True_AppendsExpectedTwoTriangles)
{
  const std::array<double,4> dist = {{1.0, 1.0, -3.0, -1.0}};
  ContourTet::append_facets_for_tet4(get_coords(), dist, facets);

  // will_cutting_quad_from_0to2_cut_largest_angle(...) returns true.

  ASSERT_EQ(facets.size(), 2u);

  const stk::math::Vector3d x5{0.75, 0.25, 0.0}; // edge 1-2
  const stk::math::Vector3d x6{0.0,  0.25, 0.0}; // edge 0-2
  const stk::math::Vector3d x7{0.0,  0.0,  0.5}; // edge 0-3
  const stk::math::Vector3d x8{0.5,  0.0,  0.5}; // edge 1-3

  expect_eq(get_facet(0).facet_vertex(0), x8);
  expect_eq(get_facet(0).facet_vertex(1), x7);
  expect_eq(get_facet(0).facet_vertex(2), x6);

  expect_eq(get_facet(1).facet_vertex(0), x8);
  expect_eq(get_facet(1).facet_vertex(1), x6);
  expect_eq(get_facet(1).facet_vertex(2), x5);
}

}


