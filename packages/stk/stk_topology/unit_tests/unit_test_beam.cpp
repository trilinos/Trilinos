#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_topology/topology.hpp>

STKUNIT_UNIT_TEST( stk_topology, beam_2)
{
  using stk::topology;

  topology t = topology::BEAM_2;


  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),2u);
  STKUNIT_EXPECT_EQ(t.num_nodes(),2u);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2u);
  STKUNIT_EXPECT_EQ(t.num_edges(),1u);
  STKUNIT_EXPECT_EQ(t.num_faces(),0u);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2u);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1u);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1u));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2u));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3u));

  STKUNIT_EXPECT_EQ(t.base(),topology::BEAM_2);

  const char a[] = "ab";
  {
    const char b[] = "ab";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "ba";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    char edge[2];
    t.edge_nodes(a,0,edge);
    t.edge_topology().equivalent(edge,"ab");
  }
}

STKUNIT_UNIT_TEST( stk_topology, beam_3)
{
  using stk::topology;

  topology t = topology::BEAM_3;


  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),2u);
  STKUNIT_EXPECT_EQ(t.num_nodes(),3u);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2u);
  STKUNIT_EXPECT_EQ(t.num_edges(),1u);
  STKUNIT_EXPECT_EQ(t.num_faces(),0u);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2u);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1u);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::BEAM_2);

  const char a[] = "abc";
  {
    const char b[] = "abc";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "bac";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    char edge[3];
    t.edge_nodes(a,0,edge);
    t.edge_topology().equivalent(edge,"abc");
  }
}


