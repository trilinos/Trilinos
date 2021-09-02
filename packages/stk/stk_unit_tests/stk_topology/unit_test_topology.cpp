// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include "gtest/gtest.h"              // for Message, TestPartResult, Test, EXPECT_EQ, Assertion...
#include "stk_topology/topology.hpp"  // for topology, operator<<, create_superedge_topology
#include <cstddef>                    // for size_t
#include <array>                      // for array
#include <iostream>                   // for operator<<, ostringstream, ostream, basic_ostream
#include <string>                     // for string, char_traits



TEST( stk_topology, lexicographical_smallest_permutation)
{
  using stk::topology;

  topology t = topology::TRI_3;

  const char nodes[]="bac";

  char permutation_nodes[4] = "";

  unsigned permutation_index = t.lexicographical_smallest_permutation((char*)nodes);
  t.permutation_nodes((char*)nodes,permutation_index,(char*)permutation_nodes);

  EXPECT_EQ( std::string("abc"), std::string(permutation_nodes));

  permutation_index = t.lexicographical_smallest_permutation((char*)nodes,true); // only consider positive permutations (true means this)
  t.permutation_nodes((char*)nodes,permutation_index,(char*)permutation_nodes);

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
      t.side_nodes( (char*)nodes, s, (char*)side_nodes );
      std::cout << "  " << side_nodes << std::endl;
    }
  }

  {
    topology t = topology::HEX_8;
    std::cout << "HEX_8 side_nodes\n";
    for (unsigned s=0; s<t.num_sides(); ++s) {
      char side_nodes[9] = {};
      t.side_nodes( (char*)nodes, s, (char*)side_nodes );
      std::cout << "  " << side_nodes << std::endl;
    }
  }

}

TEST( stk_topology, superedge_topology)
{
    stk::topology t = stk::create_superedge_topology(4);
    EXPECT_EQ( t.num_nodes(), 4u);
    EXPECT_EQ( t.rank(), stk::topology::EDGE_RANK);
    EXPECT_EQ( 0u, t.num_edges());
    EXPECT_EQ( 0u, t.num_permutations());

    EXPECT_EQ( true, t.is_superedge());
    {
      std::ostringstream name;
      name << t ;
      std::string goldName("SUPEREDGE_TOPOLOGY_4");
      EXPECT_EQ( goldName, name.str() );
    }

    stk::topology notSuper = stk::topology::LINE_2;
    EXPECT_FALSE( notSuper.is_superedge());

    stk::topology newT = stk::create_superedge_topology(8);

    EXPECT_NE( newT.num_nodes(), 6u);
    EXPECT_EQ( newT.rank(), stk::topology::EDGE_RANK);

    EXPECT_EQ( true, newT.is_superedge());
    {
      std::ostringstream name;
      name << newT ;
      std::string goldName("SUPEREDGE_TOPOLOGY_6");
      EXPECT_NE( goldName, name.str() );
    }

    stk::topology anotherT = stk::create_superedge_topology(4);
    EXPECT_EQ(t, anotherT);

    stk::topology badT = stk::create_superedge_topology(-2);
    EXPECT_TRUE(!badT.is_valid());
    EXPECT_EQ(badT, stk::topology::INVALID_TOPOLOGY);
    EXPECT_EQ( badT.rank(), stk::topology::INVALID_RANK);
}

TEST( stk_topology, superface_topology)
{
    stk::topology t = stk::create_superface_topology(4);
    EXPECT_EQ( t.num_nodes(), 4u);
    EXPECT_EQ( t.rank(), stk::topology::FACE_RANK);
    EXPECT_EQ( 0u, t.num_edges());
    EXPECT_EQ( 0u, t.num_permutations());

    EXPECT_EQ( true, t.is_superface());
    {
      std::ostringstream name;
      name << t ;
      std::string goldName("SUPERFACE_TOPOLOGY_4");
      EXPECT_EQ( goldName, name.str() );
    }

    stk::topology notSuper = stk::topology::QUAD_4;
    EXPECT_FALSE( notSuper.is_superface());

    stk::topology newT = stk::create_superface_topology(8);

    EXPECT_NE( newT.num_nodes(), 6u);
    EXPECT_EQ( newT.rank(), stk::topology::FACE_RANK);

    EXPECT_EQ( true, newT.is_superface());
    {
      std::ostringstream name;
      name << newT ;
      std::string goldName("SUPERFACE_TOPOLOGY_6");
      EXPECT_NE( goldName, name.str() );
    }

    stk::topology anotherT = stk::create_superface_topology(4);
    EXPECT_EQ(t, anotherT);

    stk::topology badT = stk::create_superface_topology(-2);
    EXPECT_TRUE(!badT.is_valid());
    EXPECT_EQ(badT, stk::topology::INVALID_TOPOLOGY);
    EXPECT_EQ( badT.rank(), stk::topology::INVALID_RANK);
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
  t.side_nodes( (int*)nodes, 0, (int*)side_nodes );
  EXPECT_EQ( 0, side_nodes[0] );
  EXPECT_EQ( 1, side_nodes[1] );
  EXPECT_EQ( 5, side_nodes[2] );
  EXPECT_EQ( 4, side_nodes[3] );
}

TEST( stk_topology, positive_polarity)
{
    std::array<stk::topology, 2> topologiesToTest = {{stk::topology::QUAD_4, stk::topology::HEX_8}};

    for (size_t i = 0; i < topologiesToTest.size(); ++i)
    {
        unsigned numNodePermutations = topologiesToTest[i].num_permutations();
        unsigned numPositiveNodePermutations = topologiesToTest[i].num_positive_permutations();

        for (unsigned permIndex = 0; permIndex < numPositiveNodePermutations; ++permIndex)
        {
            EXPECT_TRUE(topologiesToTest[i].is_positive_polarity(permIndex));
        }
        for (unsigned permIndex = numPositiveNodePermutations; permIndex < numNodePermutations; ++permIndex)
        {
            EXPECT_TRUE(!topologiesToTest[i].is_positive_polarity(permIndex));
        }
    }
}
