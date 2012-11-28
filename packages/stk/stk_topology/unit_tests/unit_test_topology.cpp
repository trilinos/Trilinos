#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_topology/topology.hpp>
#include <stk_topology/pretty_print.hpp>

#include <cstring>

#include <iostream>

STKUNIT_MAIN(argc,argv)

STKUNIT_UNIT_TEST( stk_topology, invalid_topology)
{
  using stk::topology;

  topology t = topology::INVALID_TOPOLOGY;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"INVALID_TOPOLOGY") == 0 );

  STKUNIT_EXPECT_FALSE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::INVALID_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::INVALID_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),0);
  STKUNIT_EXPECT_EQ(t.num_nodes(),0);
  STKUNIT_EXPECT_EQ(t.num_vertices(),0);
  STKUNIT_EXPECT_EQ(t.num_edges(),0);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),0);
  STKUNIT_EXPECT_EQ(t.num_permutations(),0);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),0);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::INVALID_TOPOLOGY);

}

STKUNIT_UNIT_TEST( stk_topology, node)
{
  using stk::topology;

  topology t = topology::NODE;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"NODE") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::NODE_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::INVALID_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),0);
  STKUNIT_EXPECT_EQ(t.num_nodes(),0);
  STKUNIT_EXPECT_EQ(t.num_vertices(),0);
  STKUNIT_EXPECT_EQ(t.num_edges(),0);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),0);
  STKUNIT_EXPECT_EQ(t.num_permutations(),0);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),0);

  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::NODE);

}

STKUNIT_UNIT_TEST( stk_topology, particle)
{
  using stk::topology;

  topology t = topology::PARTICLE;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"PARTICLE") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::NODE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),1);
  STKUNIT_EXPECT_EQ(t.num_nodes(),1);
  STKUNIT_EXPECT_EQ(t.num_vertices(),1);
  STKUNIT_EXPECT_EQ(t.num_edges(),0);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),1);
  STKUNIT_EXPECT_EQ(t.num_permutations(),1);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::PARTICLE);

}

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


STKUNIT_UNIT_TEST( stk_topology, beam_2)
{
  using stk::topology;

  topology t = topology::BEAM_2;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"BEAM_2") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),2);
  STKUNIT_EXPECT_EQ(t.num_nodes(),2);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2);
  STKUNIT_EXPECT_EQ(t.num_edges(),1);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),1);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::BEAM_2);

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

  {
    char edge[2];
    t.edge_nodes(a,0,edge);
    t.edge_topology(0).equivalent(edge,"ab");
  }
}

STKUNIT_UNIT_TEST( stk_topology, beam_3)
{
  using stk::topology;

  topology t = topology::BEAM_3;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"BEAM_3") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_FALSE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),2);
  STKUNIT_EXPECT_EQ(t.num_nodes(),3);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2);
  STKUNIT_EXPECT_EQ(t.num_edges(),1);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),1);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::BEAM_2);

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

  {
    char edge[3];
    t.edge_nodes(a,0,edge);
    t.edge_topology(0).equivalent(edge,"abc");
  }
}


STKUNIT_UNIT_TEST( stk_topology, shell_line_2)
{
  using stk::topology;

  topology t = topology::SHELL_LINE_2;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"SHELL_LINE_2") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_TRUE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),2);
  STKUNIT_EXPECT_EQ(t.num_nodes(),2);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2);
  STKUNIT_EXPECT_EQ(t.num_edges(),2);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),2);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::SHELL_LINE_2);

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

  {
    char edge[2];
    t.edge_nodes(a,0,edge);
    t.edge_topology(0).equivalent(edge,"ab");

    t.edge_nodes(a,1,edge);
    t.edge_topology(1).equivalent(edge,"ba");
  }
}

STKUNIT_UNIT_TEST( stk_topology, shell_line_3)
{
  using stk::topology;

  topology t = topology::SHELL_LINE_3;

  STKUNIT_EXPECT_TRUE( strcmp(t.name(),"SHELL_LINE_3") == 0 );

  STKUNIT_EXPECT_TRUE(t.is_valid());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_edges());
  STKUNIT_EXPECT_FALSE(t.has_homogeneous_faces());
  STKUNIT_EXPECT_TRUE(t.has_homogeneous_sides());
  STKUNIT_EXPECT_TRUE(t.is_shell());

  STKUNIT_EXPECT_EQ(t.rank(),topology::ELEMENT_RANK);
  STKUNIT_EXPECT_EQ(t.side_rank(),topology::EDGE_RANK);


  STKUNIT_EXPECT_EQ(t.dimension(),2);
  STKUNIT_EXPECT_EQ(t.num_nodes(),3);
  STKUNIT_EXPECT_EQ(t.num_vertices(),2);
  STKUNIT_EXPECT_EQ(t.num_edges(),2);
  STKUNIT_EXPECT_EQ(t.num_faces(),0);
  STKUNIT_EXPECT_EQ(t.num_sides(),2);
  STKUNIT_EXPECT_EQ(t.num_permutations(),2);
  STKUNIT_EXPECT_EQ(t.num_positive_permutations(),1);

  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(1));
  STKUNIT_EXPECT_TRUE(t.defined_on_spatial_dimension(2));
  STKUNIT_EXPECT_FALSE(t.defined_on_spatial_dimension(3));

  STKUNIT_EXPECT_EQ(t.base(),topology::SHELL_LINE_2);

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

  {
    char edge[3];
    t.edge_nodes(a,0,edge);
    t.edge_topology(0).equivalent(edge,"abc");

    t.edge_nodes(a,1,edge);
    t.edge_topology(1).equivalent(edge,"bac");
  }
}






















STKUNIT_UNIT_TEST( stk_topology, equivalent_tri_3)
{
  using stk::topology;

  topology t = topology::TRI_3;
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

STKUNIT_UNIT_TEST( stk_topology, equivalent_tri_4)
{
  using stk::topology;

  topology t = topology::TRI_4;
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

STKUNIT_UNIT_TEST( stk_topology, equivalent_tri_6)
{
  using stk::topology;

  topology t = topology::TRI_6;
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
