// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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
