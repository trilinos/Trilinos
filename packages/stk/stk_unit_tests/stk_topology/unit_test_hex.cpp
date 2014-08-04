#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>

TEST( stk_topology, hex8)
{
  stk::topology hex8 = stk::topology::HEX_8;

  EXPECT_TRUE(hex8.is_valid());
  EXPECT_TRUE(hex8.has_homogeneous_faces());
  EXPECT_FALSE(hex8.is_shell());

  EXPECT_EQ(hex8.rank(),stk::topology::ELEMENT_RANK);
  EXPECT_EQ(hex8.side_rank(),stk::topology::FACE_RANK);

  const unsigned num_nodes = 8;
  const unsigned num_edges = 12;
  const unsigned num_faces = 6;

  EXPECT_EQ(hex8.num_nodes(),num_nodes);
  EXPECT_EQ(hex8.num_vertices(),num_nodes);
  EXPECT_EQ(hex8.num_edges(),num_edges);
  EXPECT_EQ(hex8.num_faces(),num_faces);

  EXPECT_FALSE(hex8.defined_on_spatial_dimension(1));
  EXPECT_FALSE(hex8.defined_on_spatial_dimension(2));
  EXPECT_TRUE(hex8.defined_on_spatial_dimension(3));

  EXPECT_EQ(hex8.base(),stk::topology::HEX_8);

}

