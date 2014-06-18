#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>

TEST( stk_topology, shell_line_2)
{
  using stk::topology;

  topology t = topology::SHELL_LINE_2;


  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),2u);
  EXPECT_EQ(t.num_vertices(),2u);
  EXPECT_EQ(t.num_edges(),2u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),2u);
  EXPECT_EQ(t.num_positive_permutations(),1u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::SHELL_LINE_2);

  const char a[] = "ab";
  {
    const char b[] = "ab";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "ba";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    char edge[2];
    t.edge_nodes(a,0,edge);
    t.edge_topology().equivalent(edge,"ab");

    t.edge_nodes(a,1,edge);
    t.edge_topology().equivalent(edge,"ba");
  }
}

TEST( stk_topology, shell_line_3)
{
  using stk::topology;

  topology t = topology::SHELL_LINE_3;


  EXPECT_TRUE(t.is_valid());
  EXPECT_FALSE(t.has_homogeneous_faces());
  EXPECT_TRUE(t.is_shell());

  EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  EXPECT_EQ(t.dimension(),2u);
  EXPECT_EQ(t.num_nodes(),3u);
  EXPECT_EQ(t.num_vertices(),2u);
  EXPECT_EQ(t.num_edges(),2u);
  EXPECT_EQ(t.num_faces(),0u);
  EXPECT_EQ(t.num_permutations(),2u);
  EXPECT_EQ(t.num_positive_permutations(),1u);

  EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  EXPECT_EQ(t.base(),topology::SHELL_LINE_2);

  const char a[] = "abc";
  {
    const char b[] = "abc";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "bac";
    EXPECT_TRUE(t.equivalent(a,b).first);
    EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    char edge[3];
    t.edge_nodes(a,0,edge);
    t.edge_topology().equivalent(edge,"abc");

    t.edge_nodes(a,1,edge);
    t.edge_topology().equivalent(edge,"bac");
  }
}

