#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>

TEST( stk_topology, particle)
{
  using stk::topology;

  topology t = topology::PARTICLE;

  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_FALSE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),topology::NODE_RANK);


  EXPECT_EQ(t.dimension(),1u);
  EXPECT_EQ(t.num_nodes(),1u);
  EXPECT_EQ(t.num_vertices(),1u);
  EXPECT_EQ(t.num_edges(),0u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),1u);
  EXPECT_EQ(t.num_positive_permutations(),1u);

  EXPECT_TRUE(t.defined_on_spatial_dimension(1));
  EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::PARTICLE);

}


