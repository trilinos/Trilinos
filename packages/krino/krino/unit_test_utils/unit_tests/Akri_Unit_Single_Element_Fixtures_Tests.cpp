// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Unit_Single_Element_Fixtures.hpp>

#include <gtest/gtest.h>

#include <stk_topology/topology.hpp>

namespace krino {

TEST(SingleElementFixture, tri3)
{
  stk::topology tri3 = stk::topology::TRIANGLE_3_2D;
  SingleElementFixture test_fixture(tri3);

  test_fixture.generate_mesh();

  EXPECT_EQ(2u, test_fixture.stk_fixture.meta_data().spatial_dimension());

  const stk::mesh::BucketVector & elem_buckets = test_fixture.stk_fixture.bulk_data().buckets(stk::topology::ELEMENT_RANK);
  ASSERT_EQ(1u, elem_buckets.size());
  EXPECT_EQ(1u, (*elem_buckets.begin())->size());
  EXPECT_EQ(3u, test_fixture.stk_fixture.bulk_data().num_nodes((*elem_buckets[0])[0]));
  const stk::mesh::BucketVector & node_buckets = test_fixture.stk_fixture.bulk_data().buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(1u, node_buckets.size());
  EXPECT_EQ(3u, (*node_buckets.begin())->size());
}

TEST(SingleElementFixture, tet4)
{
  stk::topology tet4 = stk::topology::TETRAHEDRON_4;
  SingleElementFixture test_fixture(tet4);

  test_fixture.generate_mesh();

  EXPECT_EQ(3u, test_fixture.stk_fixture.meta_data().spatial_dimension());

  const stk::mesh::BucketVector & elem_buckets = test_fixture.stk_fixture.bulk_data().buckets(stk::topology::ELEMENT_RANK);
  ASSERT_EQ(1u, elem_buckets.size());
  EXPECT_EQ(1u, (*elem_buckets.begin())->size());
  EXPECT_EQ(4u, test_fixture.stk_fixture.bulk_data().num_nodes((*elem_buckets[0])[0]));
  const stk::mesh::BucketVector & node_buckets = test_fixture.stk_fixture.bulk_data().buckets(stk::topology::NODE_RANK);
  ASSERT_EQ(1u, node_buckets.size());
  EXPECT_EQ(4u, (*node_buckets.begin())->size());
}

}
