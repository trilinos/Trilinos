#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_topology/topology.hpp>


STKUNIT_UNIT_TEST( stk_topology, line_2)
{
  using stk::topology;

  topology t = topology::LINE_2;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"LINE_2") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::EDGE_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::NODE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),1);
  STKUNIT_EXPECT_EQ(t.num_nodes(),2);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2);
  STKUNIT_EXPECT_EQ(t.num_edges(),0);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),2);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::LINE_2);

  const char a[] = "ab";
  {
    const char b[] = "ab";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0);
  }

  {
    const char b[] = "ba";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1);
  }
}

STKUNIT_UNIT_TEST( stk_topology, line_3)
{
  using stk::topology;

  topology t = topology::LINE_3;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"LINE_3") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::EDGE_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::NODE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),1);
  STKUNIT_EXPECT_EQ(t.num_nodes(),3);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2);
  STKUNIT_EXPECT_EQ(t.num_edges(),0);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),2);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::LINE_2);

  const char a[] = "abc";
  {
    const char b[] = "abc";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0);
  }

  {
    const char b[] = "bac";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1);
  }
}


STKUNIT_UNIT_TEST( stk_topology, line_2_1d)
{
  using stk::topology;

  topology t = topology::LINE_2_1D;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"LINE_2_1D") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::NODE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),1);
  STKUNIT_EXPECT_EQ(t.num_nodes(),2);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2);
  STKUNIT_EXPECT_EQ(t.num_edges(),0);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),2);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::LINE_2_1D);

  const char a[] = "ab";
  {
    const char b[] = "ab";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0);
  }

  {
    const char b[] = "ba";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1);
  }
}

STKUNIT_UNIT_TEST( stk_topology, line_3_1d)
{
  using stk::topology;

  topology t = topology::LINE_3_1D;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"LINE_3_1D") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::NODE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),1);
  STKUNIT_EXPECT_EQ(t.num_nodes(),3);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2);
  STKUNIT_EXPECT_EQ(t.num_edges(),0);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),2);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::LINE_2_1D);

  const char a[] = "abc";
  {
    const char b[] = "abc";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0);
  }

  {
    const char b[] = "bac";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1);
  }
}
