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

#include "gtest/gtest.h"              // for AssertionResult, Message, TestPartResult, EXPECT_TRUE
#include "stk_topology/topology.hpp"  // for topology, topology::topology_t, topology_data, oper...
#include <map>                        // for map, _Rb_tree_const_iterator, map<>::const_iterator
#include <type_traits>                // for enable_if
#include <utility>                    // for pair
#include <vector>                     // for vector, operator<, operator==


using namespace stk;
using namespace stk::topology_detail;

namespace {

template <stk::topology::topology_t Topology>
typename std::enable_if<(topology_data<Topology>::num_edges > 0), bool>::type check_edge_node_offsets()
{
  using TopologyData = topology_data<Topology>;
  bool edge_node_ordinals = true;

  unsigned totalEdgeNodeOffset = 0;
  EXPECT_EQ(totalEdgeNodeOffset, TopologyData::edge_node_ordinals_offsets[0]);
  edge_node_ordinals &= (0 == TopologyData::edge_node_ordinals_offsets[0]);

  for (unsigned edge = 0; edge < TopologyData::num_edges; ++edge) {
    stk::topology edgeTopo = stk::topology(Topology).edge_topology(edge);
    unsigned numEdgeNodes = edgeTopo.num_nodes();
    totalEdgeNodeOffset += numEdgeNodes;

    EXPECT_EQ(totalEdgeNodeOffset, TopologyData::edge_node_ordinals_offsets[edge+1]);
    edge_node_ordinals &= (totalEdgeNodeOffset == TopologyData::edge_node_ordinals_offsets[edge+1]);
  }

  unsigned ordinalVectorLength = sizeof(TopologyData::edge_node_ordinals_vector) / sizeof(*TopologyData::edge_node_ordinals_vector);
  EXPECT_EQ(totalEdgeNodeOffset, ordinalVectorLength);
  edge_node_ordinals &= (totalEdgeNodeOffset == ordinalVectorLength);

  return edge_node_ordinals;
}

template <stk::topology::topology_t Topology>
typename std::enable_if<(topology_data<Topology>::num_edges == 0), bool>::type check_edge_node_offsets()
{
  return true;
}


template <stk::topology::topology_t Topology>
typename std::enable_if<(topology_data<Topology>::num_faces > 0), bool>::type check_face_topology()
{
  using TopologyData = topology_data<Topology>;
  bool face_topology = true;

  unsigned faceTopologyVectorLength = sizeof(TopologyData::face_topology_vector) / sizeof(*TopologyData::face_topology_vector);
  face_topology = (stk::topology(Topology).num_faces() == faceTopologyVectorLength);
  EXPECT_EQ(stk::topology(Topology).num_faces(), faceTopologyVectorLength);

  return face_topology;
}

template <stk::topology::topology_t Topology>
typename std::enable_if<(topology_data<Topology>::num_faces == 0), bool>::type check_face_topology()
{
  return true;
}


template <stk::topology::topology_t Topology>
typename std::enable_if<(topology_data<Topology>::num_faces > 0), bool>::type check_face_node_offsets()
{
  using TopologyData = topology_data<Topology>;
  bool face_node_ordinals = true;
  unsigned totalFaceNodeOffset = 0;

  EXPECT_EQ(totalFaceNodeOffset, TopologyData::face_node_ordinals_offsets[0]);
  face_node_ordinals &= (0 == TopologyData::face_node_ordinals_offsets[0]);

  for (unsigned face = 0; face < TopologyData::num_faces; ++face) {
    stk::topology faceTopo = stk::topology(Topology).face_topology(face);
    unsigned numFaceNodes = faceTopo.num_nodes();
    totalFaceNodeOffset += numFaceNodes;

    EXPECT_EQ(totalFaceNodeOffset, TopologyData::face_node_ordinals_offsets[face+1]);
    face_node_ordinals &= (totalFaceNodeOffset == TopologyData::face_node_ordinals_offsets[face+1]);
  }

  unsigned ordinalVectorLength = sizeof(TopologyData::face_node_ordinals_vector) / sizeof(*TopologyData::face_node_ordinals_vector);
  EXPECT_EQ(totalFaceNodeOffset, ordinalVectorLength);
  face_node_ordinals &= (totalFaceNodeOffset == ordinalVectorLength);

  return face_node_ordinals;
}

template <stk::topology::topology_t Topology>
typename std::enable_if<(topology_data<Topology>::num_faces == 0), bool>::type check_face_node_offsets()
{
  return true;
}


template <stk::topology::topology_t Topology>
typename std::enable_if<((topology_data<Topology>::num_faces > 0) && (topology_data<Topology>::num_nodes > 0)), bool>::type check_permutation_node_offsets()
{
  using TopologyData = topology_data<Topology>;

  unsigned permutationVectorLength = sizeof(TopologyData::permutation_node_ordinals_vector) / sizeof(*TopologyData::permutation_node_ordinals_vector);
  unsigned numPermutations = permutationVectorLength / TopologyData::num_nodes;
  bool permutation_node_ordinals = (TopologyData::num_permutations == numPermutations);
  EXPECT_EQ(stk::topology(Topology).num_permutations(), numPermutations);

  return permutation_node_ordinals;
}

template <stk::topology::topology_t Topology>
typename std::enable_if<((topology_data<Topology>::num_faces == 0) || (topology_data<Topology>::num_nodes == 0)), bool>::type check_permutation_node_offsets()
{
  return true;
}


template <stk::topology::topology_t Topology>
bool validate_topology_data()
{
  return  (   check_edge_node_offsets<Topology>()
           && check_face_topology<Topology>()
           && check_face_node_offsets<Topology>()
           && check_permutation_node_offsets<Topology>()
           );
}

} //unnamed namespace

TEST( stk_topology, validate_topology)
{
  EXPECT_TRUE( validate_topology_data< topology::NODE                        >() );
  EXPECT_TRUE( validate_topology_data< topology::LINE_2                      >() );
  EXPECT_TRUE( validate_topology_data< topology::LINE_3                      >() );
  EXPECT_TRUE( validate_topology_data< topology::TRI_3                       >() );
  EXPECT_TRUE( validate_topology_data< topology::TRI_4                       >() );
  EXPECT_TRUE( validate_topology_data< topology::TRI_6                       >() );
  EXPECT_TRUE( validate_topology_data< topology::QUAD_4                      >() );
  EXPECT_TRUE( validate_topology_data< topology::QUAD_6                      >() );
  EXPECT_TRUE( validate_topology_data< topology::QUAD_8                      >() );
  EXPECT_TRUE( validate_topology_data< topology::QUAD_9                      >() );
  EXPECT_TRUE( validate_topology_data< topology::PARTICLE                    >() );
  EXPECT_TRUE( validate_topology_data< topology::LINE_2_1D                   >() );
  EXPECT_TRUE( validate_topology_data< topology::LINE_3_1D                   >() );
  EXPECT_TRUE( validate_topology_data< topology::BEAM_2                      >() );
  EXPECT_TRUE( validate_topology_data< topology::BEAM_3                      >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_LINE_2                >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_LINE_3                >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_SIDE_BEAM_2           >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_SIDE_BEAM_3           >() );
  EXPECT_TRUE( validate_topology_data< topology::SPRING_2                    >() );
  EXPECT_TRUE( validate_topology_data< topology::SPRING_3                    >() );
  EXPECT_TRUE( validate_topology_data< topology::TRI_3_2D                    >() );
  EXPECT_TRUE( validate_topology_data< topology::TRI_4_2D                    >() );
  EXPECT_TRUE( validate_topology_data< topology::TRI_6_2D                    >() );
  EXPECT_TRUE( validate_topology_data< topology::QUAD_4_2D                   >() );
  EXPECT_TRUE( validate_topology_data< topology::QUAD_8_2D                   >() );
  EXPECT_TRUE( validate_topology_data< topology::QUAD_9_2D                   >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_TRI_3                 >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_TRI_4                 >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_TRI_6                 >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_TRI_3_ALL_FACE_SIDES  >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_TRI_4_ALL_FACE_SIDES  >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_TRI_6_ALL_FACE_SIDES  >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_QUAD_4                >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_QUAD_8                >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_QUAD_9                >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_QUAD_4_ALL_FACE_SIDES >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_QUAD_8_ALL_FACE_SIDES >() );
  EXPECT_TRUE( validate_topology_data< topology::SHELL_QUAD_9_ALL_FACE_SIDES >() );
  EXPECT_TRUE( validate_topology_data< topology::TET_4                       >() );
  EXPECT_TRUE( validate_topology_data< topology::TET_8                       >() );
  EXPECT_TRUE( validate_topology_data< topology::TET_10                      >() );
  EXPECT_TRUE( validate_topology_data< topology::TET_11                      >() );
  EXPECT_TRUE( validate_topology_data< topology::PYRAMID_5                   >() );
  EXPECT_TRUE( validate_topology_data< topology::PYRAMID_13                  >() );
  EXPECT_TRUE( validate_topology_data< topology::PYRAMID_14                  >() );
  EXPECT_TRUE( validate_topology_data< topology::WEDGE_6                     >() );
  EXPECT_TRUE( validate_topology_data< topology::WEDGE_12                    >() );
  EXPECT_TRUE( validate_topology_data< topology::WEDGE_15                    >() );
  EXPECT_TRUE( validate_topology_data< topology::WEDGE_18                    >() );
  EXPECT_TRUE( validate_topology_data< topology::HEX_8                       >() );
  EXPECT_TRUE( validate_topology_data< topology::HEX_20                      >() );
  EXPECT_TRUE( validate_topology_data< topology::HEX_27                      >() );

  // check that the permutations define the same sides
  for (stk::topology topo = stk::topology::BEGIN_TOPOLOGY; topo < stk::topology::END_TOPOLOGY; ++topo) {

    if(topo == stk::topology::QUAD_6 || topo == stk::topology::WEDGE_12) { continue; }
    if(topo.is_shell() && topo.side_rank() == stk::topology::FACE_RANK) { continue; }

    if (topo.num_permutations() > 1u && topo.side_rank() > stk::topology::NODE_RANK ) {
      const unsigned num_permutations = topo.num_permutations();
      const unsigned num_sides = topo.num_sides();

      std::map< std::vector<unsigned>, unsigned> side_map;
      for (unsigned side=0; side < num_sides; ++side) {
        stk::topology side_topo = topo.side_topology(side);
        std::vector<unsigned> tmp_side_nodes(side_topo.num_nodes());
        std::vector<unsigned> side_nodes(tmp_side_nodes.size());

        topo.side_node_ordinals(side, tmp_side_nodes.data());
        unsigned side_perm = side_topo.lexicographical_smallest_permutation(tmp_side_nodes.data());
        side_topo.permutation_nodes(tmp_side_nodes.data(), side_perm, side_nodes.data());

        side_map[side_nodes] = 0;
      }

      std::vector<unsigned> nodes(topo.num_nodes());
      for (unsigned perm = 0; perm < num_permutations; ++perm) {
        topo.permutation_node_ordinals(perm, nodes.data());

        for (unsigned side=0; side < num_sides; ++side) {
          stk::topology side_topo = topo.side_topology(side);
          std::vector<unsigned> tmp_side_nodes(side_topo.num_nodes());
          std::vector<unsigned> side_nodes(tmp_side_nodes.size());

          topo.side_nodes(nodes.data(), side, tmp_side_nodes.data());
          unsigned side_perm = side_topo.lexicographical_smallest_permutation(tmp_side_nodes.data());
          side_topo.permutation_nodes(tmp_side_nodes.data(), side_perm, side_nodes.data());

          side_map[side_nodes] += 1;
        }
      }

      if (!topo.is_shell()) {
        EXPECT_EQ(side_map.size(), num_sides);
      }
      else if(topo.rank() != stk::topology::FACE_RANK) {
        EXPECT_EQ(side_map.size(), num_sides/2u);
      }

      for (std::map< std::vector<unsigned>, unsigned>::const_iterator itr = side_map.begin(); itr != side_map.end(); ++itr) {
        // Expect that the side has been touched for each permutation
        if (!topo.is_shell()) {
          EXPECT_EQ( itr->second, num_permutations);
        }
        else if(topo.rank() != stk::topology::FACE_RANK) {
          EXPECT_EQ( itr->second/2u, num_permutations);
        }
      }
    }
  }
}

TEST( stk_topology, verify_3D_topology)
{
    std::vector<stk::topology::topology_t> goldValues = {stk::topology::TET_4, stk::topology::TET_8, stk::topology::TET_10,
                                                         stk::topology::TET_11, stk::topology::PYRAMID_5,
                                                         stk::topology::PYRAMID_13, stk::topology::PYRAMID_14,
                                                         stk::topology::WEDGE_6, stk::topology::WEDGE_12, stk::topology::WEDGE_15, 
                                                         stk::topology::WEDGE_18, stk::topology::HEX_8, stk::topology::HEX_20, stk::topology::HEX_27};

    std::vector<stk::topology::topology_t> results;
    for(stk::topology::topology_t i = stk::topology::BEGIN_TOPOLOGY; i < stk::topology::END_TOPOLOGY; ++i)
    {
        stk::topology t(i);
        if(stk::is_solid_element(t))
            results.push_back(t());
    }

    EXPECT_TRUE(goldValues==results);
}
