#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_topology/topology.hpp>

STKUNIT_UNIT_TEST( stk_topology, tri_3)
{
  using stk::topology;

  topology t = topology::TRI_3;

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::FACE_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);

  STKUNIT_EXPECT_EQ(t.dimension(),2);
  STKUNIT_EXPECT_EQ(t.num_nodes(),3);
  STKUNIT_EXPECT_EQ(t.num_vertices(),3);
  STKUNIT_EXPECT_EQ(t.num_edges(),3);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_permutations(),6);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),3);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::TRI_3);

  const char a[] = "abc";

  {
    const char b[] = "abc";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0);
  }

  {
    const char b[] = "cab";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1);
  }

  {
    const char b[] = "bca";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,2);
  }

  {
    const char b[] = "acb";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,3);
  }

  {
    const char b[] = "cba";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,4);
  }

  {
    const char b[] = "bac";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,5);
  }
}

STKUNIT_UNIT_TEST( stk_topology, tri_4)
{
  using stk::topology;

  topology t = topology::TRI_4;

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::FACE_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);

  STKUNIT_EXPECT_EQ(t.dimension(),2);
  STKUNIT_EXPECT_EQ(t.num_nodes(),4);
  STKUNIT_EXPECT_EQ(t.num_vertices(),3);
  STKUNIT_EXPECT_EQ(t.num_edges(),3);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_permutations(),6);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),3);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::TRI_3);
  const char a[] = "abcd";

  {
    const char b[] = "abcd";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0);
  }

  {
    const char b[] = "cabd";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1);
  }

  {
    const char b[] = "bcad";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,2);
  }

  {
    const char b[] = "acbd";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,3);
  }

  {
    const char b[] = "cbad";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,4);
  }

  {
    const char b[] = "bacd";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,5);
  }
}

STKUNIT_UNIT_TEST( stk_topology, tri_6)
{
  using stk::topology;

  topology t = topology::TRI_6;

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::FACE_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);

  STKUNIT_EXPECT_EQ(t.dimension(),2);
  STKUNIT_EXPECT_EQ(t.num_nodes(),6);
  STKUNIT_EXPECT_EQ(t.num_vertices(),3);
  STKUNIT_EXPECT_EQ(t.num_edges(),3);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_permutations(),6);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),3);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  const char a[] = "abc012";

  {
    const char b[] = "abc012";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0);
  }

  {
    const char b[] = "cab201";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1);
  }

  {
    const char b[] = "bca120";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,2);
  }

  {
    const char b[] = "acb210";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,3);
  }

  {
    const char b[] = "cba102";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,4);
  }

  {
    const char b[] = "bac021";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,5);
  }
}

