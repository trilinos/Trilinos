#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>

TEST( stk_topology, invalid_topology)
{
  using stk::topology;

  topology t = topology::INVALID_TOPOLOGY;


  EXPECT_FALSE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::INVALID_RANK);
  EXPECT_EQ(t.side_rank(),topology::INVALID_RANK);


  EXPECT_EQ(t.dimension(),0u);
  EXPECT_EQ(t.num_nodes(),0u);
  EXPECT_EQ(t.num_vertices(),0u);
  EXPECT_EQ(t.num_edges(),0u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),0u);
  EXPECT_EQ(t.num_positive_permutations(),0u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::INVALID_TOPOLOGY);

}

