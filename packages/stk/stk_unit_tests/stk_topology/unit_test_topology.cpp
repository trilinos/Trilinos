#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>


TEST( stk_topology, lexicographical_smallest_permutation)
{
  using stk::topology;

  topology t = topology::TRI_3;

  const char nodes[]="bac";

  char permutation_nodes[4] = "";

  unsigned permutation_index = t.lexicographical_smallest_permutation(nodes);
  t.permutation_nodes(nodes,permutation_index,permutation_nodes);

  EXPECT_EQ( std::string("abc"), std::string(permutation_nodes));

  permutation_index = t.lexicographical_smallest_permutation(nodes,true); // only consider positive permutations (true means this)
  t.permutation_nodes(nodes,permutation_index,permutation_nodes);

  EXPECT_EQ( std::string("acb"), std::string(permutation_nodes));
}

TEST( stk_topology, side_node_ordinals)
{
  using stk::topology;

  const char nodes[] = "12345678";

  {
    topology t = topology::QUAD_4_2D;
    std::cout << "QUAD_4_2D side_nodes\n";
    for (unsigned s=0; s<t.num_sides(); ++s) {
      char side_nodes[9] = {};
      t.side_nodes( nodes, s, side_nodes );
      std::cout << "  " << side_nodes << std::endl;
    }
  }

  {
    topology t = topology::HEX_8;
    std::cout << "HEX_8 side_nodes\n";
    for (unsigned s=0; s<t.num_sides(); ++s) {
      char side_nodes[9] = {};
      t.side_nodes( nodes, s, side_nodes );
      std::cout << "  " << side_nodes << std::endl;
    }
  }

}

TEST( stk_topology, superelement_topology )
{
  using stk::topology;

  topology t = stk::create_superelement_topology(6);

  EXPECT_EQ( t.num_nodes(), 6u);
  EXPECT_EQ( t.rank(), topology::ELEMENT_RANK);

  EXPECT_EQ( true, t.is_superelement());
  {
    std::ostringstream name;
    name << t ;
    std::string goldName("SUPERELEMENT_TOPOLOGY_6");
    EXPECT_EQ( goldName, name.str() );
  }

  topology notSuper = topology::HEX_8;
  EXPECT_FALSE( notSuper.is_superelement());

  topology newT = stk::create_superelement_topology(8);

  EXPECT_NE( newT.num_nodes(), 6u);
  EXPECT_EQ( newT.rank(), topology::ELEMENT_RANK);

  EXPECT_EQ( true, newT.is_superelement());
  {
    std::ostringstream name;
    name << newT ;
    std::string goldName("SUPERELEMENT_TOPOLOGY_6");
    EXPECT_NE( goldName, name.str() );
  }

  topology anotherT = stk::create_superelement_topology(6);
  EXPECT_EQ(t, anotherT);

  topology badT = stk::create_superelement_topology(-2);

  EXPECT_EQ( badT.rank(), topology::INVALID_RANK);
}

TEST( stk_topology, arrayMesh )
{
  using stk::topology;

  const int nodes[] = {0,1,2,3,4,5,6,7};
  topology t = topology::HEX_8;
  int side_nodes[4] = {};
  t.side_nodes( nodes, 0, side_nodes );
  EXPECT_EQ( 0, side_nodes[0] );
  EXPECT_EQ( 1, side_nodes[1] );
  EXPECT_EQ( 5, side_nodes[2] );
  EXPECT_EQ( 4, side_nodes[3] );
}
