#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_topology/topology.hpp>

STKUNIT_UNIT_TEST( stk_topology, node)
{
  using stk::topology;

  topology t = topology::NODE;

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::NODE_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::INVALID_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),0);
  STKUNIT_EXPECT_EQ(t.num_nodes(),0);
  STKUNIT_EXPECT_EQ(t.num_vertices(),0);
  STKUNIT_EXPECT_EQ(t.num_edges(),0);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_permutations(),0);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),0);

  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::NODE);

}

