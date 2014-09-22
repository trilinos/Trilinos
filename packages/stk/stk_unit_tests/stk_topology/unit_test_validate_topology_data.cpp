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

#include <map>
#include <vector>

#include <boost/mpl/assert.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/at.hpp>

using namespace stk;
using namespace stk::topology_detail;
using namespace boost;

namespace {

template <typename TopologyData, typename PermutationVector, unsigned NumNodes, unsigned Permutation = 0, unsigned NumPermutations = mpl::size<PermutationVector>::value>
struct check_permutations
{
  static const bool value =    (NumNodes == mpl::size< typename mpl::at_c<PermutationVector,Permutation>::type>::value)
                            && check_permutations<TopologyData,PermutationVector,NumNodes,Permutation+1>::value;

  BOOST_MPL_ASSERT_MSG(   value
                        , PERMUTATION_NODE_ORDINALS_NUM_NODES_ERROR
                        , (TopologyData)
                      );

};

template <typename TopologyData, typename PermutationVector, unsigned NumNodes, unsigned Permutation>
struct check_permutations<TopologyData,PermutationVector,NumNodes,Permutation,Permutation>
{
  static const bool value = true;
};



template <typename TopologyData, unsigned Face = 0, unsigned NumFaces = TopologyData::num_faces >
struct check_faces
{

  typedef typename mpl::at_c<typename TopologyData::face_topology_vector, Face>::type face_topology_;
  typedef topology_data< face_topology_::value > face_topology;

  typedef typename mpl::at_c<typename TopologyData::face_node_ordinals_vector, Face>::type face_nodes;

  static const bool value =    (face_topology::num_nodes == mpl::size< face_nodes >::value)
                            && check_faces<TopologyData,Face+1>::value;

  BOOST_MPL_ASSERT_MSG(   value
                        , FACE_NODE_ORDINALS_NUM_NODES_ERROR
                        , (TopologyData)
                      );

};

template <typename TopologyData, unsigned Face>
struct check_faces<TopologyData,Face,Face>
{
  static const bool value = true;
};

template <typename TopologyData>
struct validate_topology_data
{
  static const bool edge_node_ordinals = TopologyData::num_edges == mpl::size<typename TopologyData::edge_node_ordinals_vector>::value;
  BOOST_MPL_ASSERT_MSG(   edge_node_ordinals
                        , EDGE_NODE_ORDINALS_SIZE_ERROR
                        , (TopologyData)
                      );

  static const bool face_topology = TopologyData::num_faces == mpl::size<typename TopologyData::face_topology_vector>::value;
  BOOST_MPL_ASSERT_MSG(   face_topology
                        , FACE_TOPOLOGY_SIZE_ERROR
                        , (TopologyData)
                      );

  static const bool face_node_ordinals = TopologyData::num_faces == mpl::size<typename TopologyData::face_node_ordinals_vector>::value;
  BOOST_MPL_ASSERT_MSG(   face_node_ordinals
                        , FACE_NODE_ORDINALS_SIZE_ERROR
                        , (TopologyData)
                      );

  static const bool permutation_node_ordinals = TopologyData::num_permutations == mpl::size<typename TopologyData::permutation_node_ordinals_vector>::value;
  BOOST_MPL_ASSERT_MSG(   permutation_node_ordinals
                        , PERMUTATION_NODE_ORDINALS_SIZE_ERROR
                        , (TopologyData)
                      );

  static const bool check_permutations_nodes = check_permutations<TopologyData, typename TopologyData::permutation_node_ordinals_vector, TopologyData::num_nodes>::value;

  static const bool check_face_nodes = check_faces<TopologyData>::value;

  static const bool value =    edge_node_ordinals
                            && face_topology
                            && face_node_ordinals
                            && permutation_node_ordinals
                            && check_permutations_nodes
                            && check_face_nodes
                           ;

};

} //unnamed namespace

TEST( stk_topology, validate_topology)
{
  EXPECT_TRUE( validate_topology_data< topology_data< topology::NODE          > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::LINE_2        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::LINE_3        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_3         > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_4         > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_6         > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_4        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_8        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_9        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::PARTICLE      > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::LINE_2_1D     > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::LINE_3_1D     > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::BEAM_2        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::BEAM_3        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_LINE_2  > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_LINE_3  > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_3_2D      > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_4_2D      > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TRI_6_2D      > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_4_2D     > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_8_2D     > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::QUAD_9_2D     > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_TRI_3   > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_TRI_4   > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_TRI_6   > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_QUAD_4  > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_QUAD_8  > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::SHELL_QUAD_9  > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TET_4         > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TET_8         > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TET_10        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::TET_11        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::PYRAMID_5     > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::PYRAMID_13    > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::PYRAMID_14    > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::WEDGE_6       > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::WEDGE_15      > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::WEDGE_18      > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::HEX_8         > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::HEX_20        > >::value );
  EXPECT_TRUE( validate_topology_data< topology_data< topology::HEX_27        > >::value );

  // check that the permutations define the same sides
  for (stk::topology topo = stk::topology::BEGIN_TOPOLOGY; topo < stk::topology::END_TOPOLOGY; ++topo) {

    if (topo.num_permutations() > 1u && topo.side_rank() > stk::topology::NODE_RANK ) {
      const unsigned num_permutations = topo.num_permutations();
      const unsigned num_sides = topo.num_sides();

      std::map< std::vector<unsigned>, unsigned> side_map;
      for (unsigned side=0; side < num_sides; ++side) {
        stk::topology side_topo = topo.side_topology(side);
        std::vector<unsigned> tmp_side_nodes(side_topo.num_nodes());
        std::vector<unsigned> side_nodes(tmp_side_nodes.size());

        topo.side_node_ordinals(side, tmp_side_nodes.begin());
        unsigned side_perm = side_topo.lexicographical_smallest_permutation(tmp_side_nodes);
        side_topo.permutation_nodes(tmp_side_nodes, side_perm, side_nodes.begin());

        side_map[side_nodes] = 0;
      }

      std::vector<unsigned> nodes(topo.num_nodes());
      for (unsigned perm = 0; perm < num_permutations; ++perm) {
        topo.permutation_node_ordinals(perm, nodes.begin());

        for (unsigned side=0; side < num_sides; ++side) {
          stk::topology side_topo = topo.side_topology(side);
          std::vector<unsigned> tmp_side_nodes(side_topo.num_nodes());
          std::vector<unsigned> side_nodes(tmp_side_nodes.size());

          topo.side_nodes(nodes, side, tmp_side_nodes.begin());
          unsigned side_perm = side_topo.lexicographical_smallest_permutation(tmp_side_nodes);
          side_topo.permutation_nodes(tmp_side_nodes, side_perm, side_nodes.begin());

          side_map[side_nodes] += 1;
        }
      }

      if (!topo.is_shell()) {
        EXPECT_EQ(side_map.size(), num_sides);
      }
      else {
        EXPECT_EQ(side_map.size(), num_sides/2u);
      }

      for (std::map< std::vector<unsigned>, unsigned>::const_iterator itr = side_map.begin(); itr != side_map.end(); ++itr) {
        // Expect that the side has been touched for each permutation
        if (!topo.is_shell()) {
          EXPECT_EQ( itr->second, num_permutations);
        }
        else {
          EXPECT_EQ( itr->second/2u, num_permutations);
        }
      }
    }
  }
}


