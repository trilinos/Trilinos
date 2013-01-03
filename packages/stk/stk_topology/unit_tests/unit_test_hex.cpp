#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_topology/topology.hpp>

STKUNIT_UNIT_TEST( stk_topology, hex8)
{
  stk::topology hex8 = stk::topology::HEX_8;

  STKUNIT_EXPECT_TRUE(hex8.is_valid());
  STKUNIT_EXPECT_TRUE(hex8.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(hex8.is_shell());

  STKUNIT_EXPECT_EQ(hex8.rank(),stk::topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(hex8.side_rank(),stk::topology::FACE_RANK);

  const int num_nodes = 8;
  const int num_edges = 12;
  const int num_faces = 6;

  STKUNIT_EXPECT_EQ(hex8.num_nodes(),num_nodes);
  STKUNIT_EXPECT_EQ(hex8.num_vertices(),num_nodes);
  STKUNIT_EXPECT_EQ(hex8.num_edges(),num_edges);
  STKUNIT_EXPECT_EQ(hex8.num_faces(),num_faces);

  STKUNIT_EXPECT_FALSE(hex8.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(hex8.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(hex8.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(hex8.base(),stk::topology::HEX_8);

}

