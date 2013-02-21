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

  STKUNIT_EXPECT_EQ(t.dimension(),2u);
  STKUNIT_EXPECT_EQ(t.num_nodes(),3u);
  STKUNIT_EXPECT_EQ(t.num_vertices(),3u);
  STKUNIT_EXPECT_EQ(t.num_edges(),3u);
  STKUNIT_EXPECT_EQ(t.num_faces(),0u);
  STKUNIT_EXPECT_EQ(t.num_permutations(),6u);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),3u);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::TRI_3);

  const char a[] = "abc";

  {
    const char b[] = "abc";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "cab";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    const char b[] = "bca";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,2u);
  }

  {
    const char b[] = "acb";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,3u);
  }

  {
    const char b[] = "cba";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,4u);
  }

  {
    const char b[] = "bac";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,5u);
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

  STKUNIT_EXPECT_EQ(t.dimension(),2u);
  STKUNIT_EXPECT_EQ(t.num_nodes(),4u);
  STKUNIT_EXPECT_EQ(t.num_vertices(),3u);
  STKUNIT_EXPECT_EQ(t.num_edges(),3u);
  STKUNIT_EXPECT_EQ(t.num_faces(),0u);
  STKUNIT_EXPECT_EQ(t.num_permutations(),6u);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),3u);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::TRI_3);
  const char a[] = "abcd";

  {
    const char b[] = "abcd";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "cabd";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    const char b[] = "bcad";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,2u);
  }

  {
    const char b[] = "acbd";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,3u);
  }

  {
    const char b[] = "cbad";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,4u);
  }

  {
    const char b[] = "bacd";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,5u);
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

  STKUNIT_EXPECT_EQ(t.dimension(),2u);
  STKUNIT_EXPECT_EQ(t.num_nodes(),6u);
  STKUNIT_EXPECT_EQ(t.num_vertices(),3u);
  STKUNIT_EXPECT_EQ(t.num_edges(),3u);
  STKUNIT_EXPECT_EQ(t.num_faces(),0u);
  STKUNIT_EXPECT_EQ(t.num_permutations(),6u);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),3u);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  const char a[] = "abc012";

  {
    const char b[] = "abc012";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,0u);
  }

  {
    const char b[] = "cab201";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,1u);
  }

  {
    const char b[] = "bca120";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,2u);
  }

  {
    const char b[] = "acb210";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,3u);
  }

  {
    const char b[] = "cba102";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,4u);
  }

  {
    const char b[] = "bac021";
    STKUNIT_EXPECT_TRUE(t.equivalent(a,b).first);
    STKUNIT_EXPECT_EQ(t.equivalent(a,b).second,5u);
  }
}

