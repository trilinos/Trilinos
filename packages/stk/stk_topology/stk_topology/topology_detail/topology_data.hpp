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

#ifndef STKTOPOLOGY_DETAIL_TOPOLOGY_DATA_HPP
#define STKTOPOLOGY_DETAIL_TOPOLOGY_DATA_HPP

#include <stk_topology/topology.hpp>
#include <type_traits>

//TODO implement permutations for tets, pyramids, wedges and hexes
//TODO implement permutations polarity

namespace stk {
namespace topology_detail {

template <topology::topology_t Topology, typename Enable = void>
struct topology_data;


//***************************************************************************
// topology::INVALID -- topology::INVALID_TOPOLOGY
//***************************************************************************

template <>
struct topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::INVALID_TOPOLOGY;
  static constexpr topology::topology_t base = topology::INVALID_TOPOLOGY;

  static constexpr bool is_valid = false;
  static constexpr topology::rank_t rank = topology::INVALID_RANK;
  static constexpr topology::rank_t side_rank = topology::INVALID_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::INVALID_TOPOLOGY};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 0;
  static constexpr unsigned num_nodes = 0;
  static constexpr unsigned num_vertices = 0;
  static constexpr unsigned num_edges = 0;
  static constexpr unsigned num_faces = 0;
  static constexpr unsigned num_permutations = 0;
  static constexpr unsigned num_positive_permutations = 0;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       false,  // 2d
                                                       false}; // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::INVALID_TOPOLOGY};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0};
  static constexpr unsigned edge_node_ordinals_vector[] = {0};

  static constexpr unsigned face_node_ordinals_offsets[] = {0};
  static constexpr unsigned face_node_ordinals_vector[] = {0};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0};
};

//***************************************************************************
// topology::NODE -- topology::NODE_RANK
//***************************************************************************

template <>
struct topology_data<topology::NODE>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::NODE;
  static constexpr topology::topology_t base = topology::NODE;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::NODE_RANK;
  static constexpr topology::rank_t side_rank = topology::INVALID_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::INVALID_TOPOLOGY};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 0;
  static constexpr unsigned num_nodes = 0;
  static constexpr unsigned num_vertices = 0;
  static constexpr unsigned num_edges = 0;
  static constexpr unsigned num_faces = 0;
  static constexpr unsigned num_permutations = 0;
  static constexpr unsigned num_positive_permutations = 0;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       true,   // 1d
                                                       true,   // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::INVALID_TOPOLOGY};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0};
  static constexpr unsigned edge_node_ordinals_vector[] = {0};

  static constexpr unsigned face_node_ordinals_offsets[] = {0};
  static constexpr unsigned face_node_ordinals_vector[] = {0};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0};
};


//***************************************************************************
// PARTICLE -- topology::ELEMENT_RANK
// one node
//***************************************************************************

template <>
struct topology_data<topology::PARTICLE>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::PARTICLE;
  static constexpr topology::topology_t base = topology::PARTICLE;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::INVALID_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::INVALID_TOPOLOGY};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 1;
  static constexpr unsigned num_nodes = 1;
  static constexpr unsigned num_vertices = 1;
  static constexpr unsigned num_edges = 0;
  static constexpr unsigned num_faces = 0;
  static constexpr unsigned num_permutations = 1;
  static constexpr unsigned num_positive_permutations = 1;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       true,   // 1d
                                                       true,   // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::INVALID_TOPOLOGY};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0};
  static constexpr unsigned edge_node_ordinals_vector[] = {0};

  static constexpr unsigned face_node_ordinals_offsets[] = {0};
  static constexpr unsigned face_node_ordinals_vector[] = {0};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0};
};


//***************************************************************************
// topology::LINE -- topology::EDGE_RANK
// 2 or 3 nodes
//
//  o------o------o
//  0      2      1
//
//***************************************************************************

template <>
struct topology_data<topology::LINE_2>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::LINE_2;
  static constexpr topology::topology_t base = topology::LINE_2;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::EDGE_RANK;
  static constexpr topology::rank_t side_rank = topology::NODE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::INVALID_TOPOLOGY};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 1;
  static constexpr unsigned num_nodes = 2;
  static constexpr unsigned num_vertices = 2;
  static constexpr unsigned num_edges = 0;
  static constexpr unsigned num_faces = 0;
  static constexpr unsigned num_permutations = 2;
  static constexpr unsigned num_positive_permutations = 1;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       true,   // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::INVALID_TOPOLOGY};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0};
  static constexpr unsigned edge_node_ordinals_vector[] = {0};

  static constexpr unsigned face_node_ordinals_offsets[] = {0};
  static constexpr unsigned face_node_ordinals_vector[] = {0};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1,
                                                                  1, 0};
};


template <>
struct topology_data<topology::LINE_3>
  : public topology_data<topology::LINE_2>
{
  static constexpr topology::topology_t value = topology::LINE_3;
  static constexpr unsigned num_nodes = 3;

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1,  2,
                                                                  1, 0,  2};
};

//***************************************************************************
// topology::LINE 1D -- topology::ELEMENT_RANK
// only defined on 1d problems
//
//  o------o------o
//  0      2      1
//
//***************************************************************************

template <>
struct topology_data<topology::LINE_2_1D>
  : public topology_data<topology::LINE_2>
{
  static constexpr topology::topology_t value = topology::LINE_2_1D;
  static constexpr topology::topology_t base = topology::LINE_2_1D;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       true,    // 1d
                                                       false,   // 2d
                                                       false};  // 3d
};

template <>
struct topology_data<topology::LINE_3_1D>
  : public topology_data<topology::LINE_3>
{
  static constexpr topology::topology_t value = topology::LINE_3_1D;
  static constexpr topology::topology_t base = topology::LINE_2_1D;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       true,    // 1d
                                                       false,   // 2d
                                                       false};  // 3d
};

//***************************************************************************
// topology::BEAM -- topology::ELEMENT_RANK
// 2 or 3 nodes with a single edge
//
//  o------o------o
//  0      2      1
//       Edge 0
//***************************************************************************

template <>
struct topology_data<topology::BEAM_2>
  : public topology_data<topology::LINE_2>
{
  static constexpr topology::topology_t value = topology::BEAM_2;
  static constexpr topology::topology_t base = topology::BEAM_2;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::EDGE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_2};

  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 2;
  static constexpr unsigned num_edges = 1;

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 2};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1};
};

template <>
struct topology_data<topology::BEAM_3>
  : public topology_data<topology::LINE_3>
{
  static constexpr topology::topology_t value = topology::BEAM_3;
  static constexpr topology::topology_t base = topology::BEAM_2;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::EDGE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3};

  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 2;
  static constexpr unsigned num_edges = 1;

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  2};
};

//***************************************************************************
// topology::SHELL_LINE -- topology::ELEMENT_RANK
// only defined on 2d problems
// 2 or 3 nodes with a two edge
//
//       Edge 1: (1,0,2)
//
//  o------o------o
//  0      2      1
//
//       Edge 0: (0,1,2)
//***************************************************************************

template <>
struct topology_data<topology::SHELL_LINE_2>
  : public topology_data<topology::LINE_2>
{
  static constexpr topology::topology_t value = topology::SHELL_LINE_2;
  static constexpr topology::topology_t base = topology::SHELL_LINE_2;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::EDGE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_2,
                                                                  topology::LINE_2};

  static constexpr bool is_shell = true;
  static constexpr unsigned dimension = 2;
  static constexpr unsigned num_edges = 2;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       false,   // 1d
                                                       true,    // 2d
                                                       false};  // 3d

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 2, 4};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,
                                                           1, 0};
};

template <>
struct topology_data<topology::SHELL_LINE_3>
  : public topology_data<topology::LINE_3>
{
  static constexpr topology::topology_t value = topology::SHELL_LINE_3;
  static constexpr topology::topology_t base = topology::SHELL_LINE_2;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::EDGE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_3};

  static constexpr bool is_shell = true;
  static constexpr unsigned dimension = 2;
  static constexpr unsigned num_edges = 2;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       false,   // 1d
                                                       true,    // 2d
                                                       false};  // 3d

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 6};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  2,
                                                           1, 0,  2};
};

//***************************************************************************
// topology::SPRING -- topology::ELEM_RANK
// 2 or 3 nodes
//
//  o------o------o
//  0      2      1
//
//***************************************************************************

template <>
struct topology_data<topology::SPRING_2>
  : public topology_data<topology::LINE_2>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::SPRING_2;
  static constexpr topology::topology_t base = topology::SPRING_2;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::NODE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::INVALID_TOPOLOGY};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 2;
  static constexpr unsigned num_nodes = 2;
  static constexpr unsigned num_vertices = 2;
  static constexpr unsigned num_edges = 0;
  static constexpr unsigned num_faces = 0;
  static constexpr unsigned num_permutations = 2;
  static constexpr unsigned num_positive_permutations = 2;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       true,   // 1d
                                                       true,   // 2d
                                                       true};  // 3d
};


template <>
struct topology_data<topology::SPRING_3>
  : public topology_data<topology::LINE_3>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::SPRING_3;
  static constexpr topology::topology_t base = topology::SPRING_2;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::NODE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::INVALID_TOPOLOGY};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 2;
  static constexpr unsigned num_nodes = 3;
  static constexpr unsigned num_vertices = 2;
  static constexpr unsigned num_edges = 0;
  static constexpr unsigned num_faces = 0;
  static constexpr unsigned num_permutations = 2;
  static constexpr unsigned num_positive_permutations = 2;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       true,   // 1d
                                                       true,   // 2d
                                                       true};  // 3d
};

//***************************************************************************
// topology::TRIANGLE -- topology::FACE_RANK
// defined on spatial dimension 3d
// 3, 4, or 6 nodes with 3 edges
/*
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2    5 o       o 4   Edge #1
//               /         \
//              /           \
//             /             \
//            o-------o-------o
//           0        3        1
//
//                  Edge #0
//
//
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2      /       \     Edge #1
//               /    4    \
//              /     o     \
//             /             \
//            o---------------o
//           0                 1
//
//                  Edge #0
*/
//***************************************************************************

template <>
struct topology_data<topology::TRI_3>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::TRI_3;
  static constexpr topology::topology_t base = topology::TRI_3;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::FACE_RANK;
  static constexpr topology::rank_t side_rank = topology::EDGE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 2;
  static constexpr unsigned num_nodes = 3;
  static constexpr unsigned num_vertices = 3;
  static constexpr unsigned num_edges = 3;
  static constexpr unsigned num_faces = 0;
  static constexpr unsigned num_permutations = 6;
  static constexpr unsigned num_positive_permutations = 3;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       false,  // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::INVALID_TOPOLOGY};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 2, 4, 6};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,
                                                           1, 2,
                                                           2, 0};

  static constexpr unsigned face_node_ordinals_offsets[] = {0};
  static constexpr unsigned face_node_ordinals_vector[] = {0};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2,
                                                                  2, 0, 1,
                                                                  1, 2, 0,
                                                                  0, 2, 1,
                                                                  2, 1, 0,
                                                                  1, 0, 2};
};

template <>
struct topology_data<topology::TRI_4>
  : public topology_data<topology::TRI_3>
{
  static constexpr topology::topology_t value = topology::TRI_4;
  static constexpr unsigned num_nodes = 4;

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2,  3,
                                                                  2, 0, 1,  3,
                                                                  1, 2, 0,  3,
                                                                  0, 2, 1,  3,
                                                                  2, 1, 0,  3,
                                                                  1, 0, 2,  3};
};

template <>
struct topology_data<topology::TRI_6>
  : public topology_data<topology::TRI_3>
{
  static constexpr topology::topology_t value = topology::TRI_6;
  static constexpr unsigned num_nodes = 6;

  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 6, 9};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  3,
                                                           1, 2,  4,
                                                           2, 0,  5};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2,  3, 4, 5,
                                                                  2, 0, 1,  5, 3, 4,
                                                                  1, 2, 0,  4, 5, 3,
                                                                  0, 2, 1,  5, 4, 3,
                                                                  2, 1, 0,  4, 3, 5,
                                                                  1, 0, 2,  3, 5, 4};
};

//***************************************************************************
// topology::TRIANGLE 2D -- topology::ELEMENT_RANK
// defined on spatial dimension 2d
// 3, 4, or 6 nodes with 3 edges
/*
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2    5 o       o 4   Edge #1
//               /         \
//              /           \
//             /             \
//            o-------o-------o
//           0        3        1
//
//                  Edge #0
//
//
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2      /       \     Edge #1
//               /    4    \
//              /     o     \
//             /             \
//            o---------------o
//           0                 1
//
//                  Edge #0
*/
//***************************************************************************

template <>
struct topology_data<topology::TRI_3_2D>
  : public topology_data<topology::TRI_3>
{
  static constexpr topology::topology_t value = topology::TRI_3_2D;
  static constexpr topology::topology_t base = topology::TRI_3_2D;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       false,   // 1d
                                                       true,    // 2d
                                                       false};  // 3d
};

template <>
struct topology_data<topology::TRI_4_2D>
  : public topology_data<topology::TRI_4>
{
  static constexpr topology::topology_t value = topology::TRI_4_2D;
  static constexpr topology::topology_t base = topology::TRI_3_2D;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       false,   // 1d
                                                       true,    // 2d
                                                       false};  // 3d
};

template <>
struct topology_data<topology::TRI_6_2D>
  : public topology_data<topology::TRI_6>
{
  static constexpr topology::topology_t value = topology::TRI_6_2D;
  static constexpr topology::topology_t base = topology::TRI_3_2D;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       false,   // 1d
                                                       true,    // 2d
                                                       false};  // 3d
};

//***************************************************************************
// topology::SHELL topology::TRIANGLE -- topology::ELEMENT_RANK
// defined on spatial dimension 3d
// 3, 4, or 6 nodes with 3 edges and 2 faces
/*
//
//
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2    5 o       o 4   Edge #1
//               /         \
//              /           \
//             /             \
//            o-------o-------o
//           0        3        1
//
//                  Edge #0
//
//   Face #0 (0, 1, 2,   3, 4, 5)
//   Face #1 (0, 2, 1,   5, 4, 3)
//
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2      /       \     Edge #1
//               /    4    \
//              /     o     \
//             /             \
//            o---------------o
//           0                 1
//
//                  Edge #0
//
//   Face #0 (0, 1, 2,   3)
//   Face #1 (0, 2, 1,   3)
*/
//***************************************************************************

template <>
struct topology_data<topology::SHELL_TRI_3>
  : public topology_data<topology::TRI_3>
{
  static constexpr topology::topology_t value = topology::SHELL_TRI_3;
  static constexpr topology::topology_t base = topology::SHELL_TRI_3;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr bool is_shell = true;
  static constexpr bool has_homogeneous_faces = true;

  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_faces = 2;

  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_3,
                                                                  topology::TRI_3};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 3, 6};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 2,
                                                           0, 2, 1};
};

template <>
struct topology_data<topology::SHELL_TRI_4>
  : public topology_data<topology::TRI_4>
{
  static constexpr topology::topology_t value = topology::SHELL_TRI_4;
  static constexpr topology::topology_t base = topology::SHELL_TRI_3;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr bool is_shell = true;
  static constexpr bool has_homogeneous_faces = true;

  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_faces = 2;

  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_4,
                                                                  topology::TRI_4};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 4, 8};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 2,  3,
                                                           0, 2, 1,  3};
};

template <>
struct topology_data<topology::SHELL_TRI_6>
  : public topology_data<topology::TRI_6>
{
  static constexpr topology::topology_t value = topology::SHELL_TRI_6;
  static constexpr topology::topology_t base = topology::SHELL_TRI_3;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr bool is_shell = true;
  static constexpr bool has_homogeneous_faces = true;

  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_faces = 2;

  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_6,
                                                                  topology::TRI_6};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 6, 12};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 2,  3, 4, 5,
                                                           0, 2, 1,  5, 4, 3};
};

//***************************************************************************
// topology::QUADRILATERAL -- topology::FACE_RANK
// defined on spatial dimension 3d
// 4, 8, or 9 nodes with 4 edges
// for quad_6, extra node is on edge0 and edge2, so that node5 is labeled node6
// in diagram below
//
//                 Edge #2
//
//            3        6        2
//             o-------o-------o
//             |               |
//             |               |
//             |       8       |
//  Edge #3  7 o       o       o 5  Edge #1
//             |               |
//             |               |
//             |               |
//             o-------o-------o
//            0        4        1
//
//                  Edge #0
//
//***************************************************************************

template <>
struct topology_data<topology::QUAD_4>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::QUAD_4;
  static constexpr topology::topology_t base = topology::QUAD_4;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::FACE_RANK;
  static constexpr topology::rank_t side_rank = topology::EDGE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_2, 
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 2;
  static constexpr unsigned num_nodes = 4;
  static constexpr unsigned num_vertices = 4;
  static constexpr unsigned num_edges = 4;
  static constexpr unsigned num_faces = 0;
  static constexpr unsigned num_permutations = 8;
  static constexpr unsigned num_positive_permutations = 4;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       false,  // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::INVALID_TOPOLOGY};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 2, 4, 6, 8};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,
                                                           1, 2,
                                                           2, 3,
                                                           3, 0};

  static constexpr unsigned face_node_ordinals_offsets[] = {0};
  static constexpr unsigned face_node_ordinals_vector[] = {0};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2, 3,
                                                                  3, 0, 1, 2,
                                                                  2, 3, 0, 1,
                                                                  1, 2, 3, 0,
                                                                  0, 3, 2, 1,
                                                                  3, 2, 1, 0,
                                                                  2, 1, 0, 3,
                                                                  1, 0, 3, 2};
};

template <>
struct topology_data<topology::QUAD_6>
  : public topology_data<topology::QUAD_4>
{
  static constexpr topology::topology_t value = topology::QUAD_6;
  static constexpr topology::rank_t rank = topology::FACE_RANK;
  static constexpr unsigned num_nodes = 6;

  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_2,
                                                                  topology::LINE_3,
                                                                  topology::LINE_2};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 5, 8, 10};
  static constexpr unsigned edge_node_ordinals_vector[] = { 0, 1, 4,
                                                            1, 2,
                                                            2, 3, 5,
                                                            3, 0 };

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3,  4, 5,
    3, 0, 1, 2,  4, 5,
    2, 3, 0, 1,  5, 4,
    1, 2, 3, 0,  5, 4,
    0, 3, 2, 1,  5, 4,
    3, 2, 1, 0,  5, 4,
    2, 1, 0, 3,  4, 5,
    1, 0, 3, 2,  4, 5 
  };
};

template <>
struct topology_data<topology::QUAD_8>
  : public topology_data<topology::QUAD_4>
{
  static constexpr topology::topology_t value = topology::QUAD_8;
  static constexpr unsigned num_nodes = 8;

  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 6, 9, 12};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  4,
                                                           1, 2,  5,
                                                           2, 3,  6,
                                                           3, 0,  7};

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3,  4, 5, 6, 7,
    3, 0, 1, 2,  7, 4, 5, 6,
    2, 3, 0, 1,  6, 7, 4, 5,
    1, 2, 3, 0,  5, 6, 7, 4,
    0, 3, 2, 1,  7, 6, 5, 4,
    3, 2, 1, 0,  6, 5, 4, 7,
    2, 1, 0, 3,  5, 4, 7, 6,
    1, 0, 3, 2,  4, 7, 6, 5
  };
};

template <>
struct topology_data<topology::QUAD_9>
  : public topology_data<topology::QUAD_8>
{
  static constexpr topology::topology_t value = topology::QUAD_9;
  static constexpr unsigned num_nodes = 9;

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3,  4, 5, 6, 7,  8,
    3, 0, 1, 2,  7, 4, 5, 6,  8,
    2, 3, 0, 1,  6, 7, 4, 5,  8,
    1, 2, 3, 0,  5, 6, 7, 4,  8,
    0, 3, 2, 1,  7, 6, 5, 4,  8,
    3, 2, 1, 0,  6, 5, 4, 7,  8,
    2, 1, 0, 3,  5, 4, 7, 6,  8,
    1, 0, 3, 2,  4, 7, 6, 5,  8
  };
};

//***************************************************************************
// topology::QUADRILATERAL 2D -- topology::ELEMENT_RANK
// defined on spatial dimension 2d
// 4, 8, or 9 nodes with 4 edges
//
//                 Edge #2
//
//            3        6        2
//             o-------o-------o
//             |               |
//             |               |
//             |       8       |
//  Edge #3  7 o       o       o 5  Edge #1
//             |               |
//             |               |
//             |               |
//             o-------o-------o
//            0        4        1
//
//                  Edge #0
//
//***************************************************************************

template <>
struct topology_data<topology::QUAD_4_2D>
  : public topology_data<topology::QUAD_4>
{
  static constexpr topology::topology_t value = topology::QUAD_4_2D;
  static constexpr topology::topology_t base = topology::QUAD_4_2D;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       false,   // 1d
                                                       true,    // 2d
                                                       false};  // 3d
};

template <>
struct topology_data<topology::QUAD_8_2D>
  : public topology_data<topology::QUAD_8>
{
  static constexpr topology::topology_t value = topology::QUAD_8_2D;
  static constexpr topology::topology_t base = topology::QUAD_4_2D;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       false,   // 1d
                                                       true,    // 2d
                                                       false};  // 3d
};

template <>
struct topology_data<topology::QUAD_9_2D>
  : public topology_data<topology::QUAD_9>
{
  static constexpr topology::topology_t value = topology::QUAD_9_2D;
  static constexpr topology::topology_t base = topology::QUAD_4_2D;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;

  static constexpr bool spatial_dimension_vector[4] = {false,   // 0d
                                                       false,   // 1d
                                                       true,    // 2d
                                                       false};  // 3d
};

//***************************************************************************
//  topology::SHELL topology::QUADRILATERAL -- topology::ELEMENT_RANK
// defined on spatial dimension 3d
// 4, 8, or 9 nodes with 4 edges and 2 faces
//
//                 Edge #2
//
//            3        6        2
//             o-------o-------o
//             |               |
//             |               |
//             |       8       |
//  Edge #3  7 o       o       o 5  Edge #1
//             |               |
//             |               |
//             |               |
//             o-------o-------o
//            0        4        1
//
//                  Edge #0
//
//  Face #0 (0, 1, 2, 3,   4, 5, 6, 7,   8)
//  Face #1 (0, 3, 2, 1,   7, 6, 5, 4,   8)
//
//***************************************************************************

template <>
struct topology_data<topology::SHELL_QUAD_4>
  : public topology_data<topology::QUAD_4>
{
  static constexpr topology::topology_t value = topology::SHELL_QUAD_4;
  static constexpr topology::topology_t base = topology::SHELL_QUAD_4;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr bool is_shell = true;
  static constexpr bool has_homogeneous_faces = true;

  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_faces = 2;

  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_4,
                                                                  topology::QUAD_4};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 4, 8};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 2, 3,
                                                           0, 3, 2, 1};
};

template <>
struct topology_data<topology::SHELL_QUAD_8>
  : public topology_data<topology::QUAD_8>
{
  static constexpr topology::topology_t value = topology::SHELL_QUAD_8;
  static constexpr topology::topology_t base = topology::SHELL_QUAD_4;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr bool is_shell = true;
  static constexpr bool has_homogeneous_faces = true;

  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_faces = 2;

  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_8,
                                                                  topology::QUAD_8};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 8, 16};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 2, 3,  4, 5, 6, 7,
                                                           0, 3, 2, 1,  7, 6, 5, 4};
};

template <>
struct topology_data<topology::SHELL_QUAD_9>
  : public topology_data<topology::QUAD_9>
{
  static constexpr topology::topology_t value = topology::SHELL_QUAD_9;
  static constexpr topology::topology_t base = topology::SHELL_QUAD_4;

  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr bool is_shell = true;
  static constexpr bool has_homogeneous_faces = true;

  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_faces = 2;

  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_9,
                                                                  topology::QUAD_9};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 9, 18};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 2, 3,  4, 5, 6, 7,  8,
                                                           0, 3, 2, 1,  7, 6, 5, 4,  8};
};

//***************************************************************************
// topology::TETRAHEDRON
//***************************************************************************

template <>
struct topology_data<topology::TET_4>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::TET_4;
  static constexpr topology::topology_t base = topology::TET_4;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2};
  static constexpr bool has_homogeneous_faces = true;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_nodes = 4;
  static constexpr unsigned num_vertices = 4;
  static constexpr unsigned num_edges = 6;
  static constexpr unsigned num_faces = 4;
  static constexpr unsigned num_permutations = 12;
  static constexpr unsigned num_positive_permutations = 12;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       false,  // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_3,
                                                                  topology::TRI_3,
                                                                  topology::TRI_3,
                                                                  topology::TRI_3};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 2, 4, 6, 8, 10, 12};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,
                                                           1, 2,
                                                           2, 0,
                                                           0, 3,
                                                           1, 3,
                                                           2, 3};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 3, 6, 9, 12};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 3,
                                                           1, 2, 3,
                                                           0, 3, 2,
                                                           0, 2, 1};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2, 3,
                                                                  1, 2, 0, 3,
                                                                  2, 0, 1, 3,
                                                                  0, 3, 1, 2,
                                                                  3, 1, 0, 2,
                                                                  1, 0, 3, 2,
                                                                  0, 2, 3, 1,
                                                                  2, 3, 0, 1,
                                                                  3, 0, 2, 1,
                                                                  1, 3, 2, 0,
                                                                  3, 2, 1, 0,
                                                                  2, 1, 3, 0};
};

//TODO: Delete TET_8
template <>
struct topology_data<topology::TET_8>
  : public topology_data<topology::TET_4>
{
  static constexpr topology::topology_t value = topology::TET_8;
  static constexpr unsigned num_nodes = 8;

  static constexpr unsigned num_permutations = 1;
  static constexpr unsigned num_positive_permutations = 1;

  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_4,
                                                                  topology::TRI_4,
                                                                  topology::TRI_4,
                                                                  topology::TRI_4};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 4, 8, 12, 16};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 3,  4,
                                                           1, 2, 3,  5,
                                                           0, 3, 2,  7,
                                                           0, 2, 1,  6};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2, 3,  4, 5, 6, 7};
};

template <>
struct topology_data<topology::TET_10>
  : public topology_data<topology::TET_4>
{
  static constexpr topology::topology_t value = topology::TET_10;
  static constexpr unsigned num_nodes = 10;

  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3};

  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_6,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 6, 9, 12, 15, 18};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  4,
                                                           1, 2,  5,
                                                           2, 0,  6,
                                                           0, 3,  7,
                                                           1, 3,  8,
                                                           2, 3,  9};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 6, 12, 18, 24};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 3,  4, 8, 7,
                                                           1, 2, 3,  5, 9, 8,
                                                           0, 3, 2,  7, 9, 6,
                                                           0, 2, 1,  6, 5, 4};

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3,  4, 5, 6, 7, 8, 9,
    1, 2, 0, 3,  5, 6, 4, 8, 9, 7,
    2, 0, 1, 3,  6, 4, 5, 9, 7, 8,
    0, 3, 1, 2,  7, 8, 4, 6, 9, 5,
    3, 1, 0, 2,  8, 4, 7, 9, 5, 6,
    1, 0, 3, 2,  4, 7, 8, 5, 6, 9,
    0, 2, 3, 1,  6, 9, 7, 4, 5, 8,
    2, 3, 0, 1,  9, 7, 6, 5, 8, 4,
    3, 0, 2, 1,  7, 6, 9, 8, 4, 5,
    1, 3, 2, 0,  8, 9, 5, 4, 7, 6,
    3, 2, 1, 0,  9, 5, 8, 7, 6, 4,
    2, 1, 3, 0,  5, 8, 9, 6, 4, 7
  };
};

template <>
struct topology_data<topology::TET_11>
  : public topology_data<topology::TET_10>
{
  static constexpr topology::topology_t value = topology::TET_11;
  static constexpr unsigned num_nodes = 11;

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3,  4, 5, 6, 7, 8, 9,  10,
    1, 2, 0, 3,  5, 6, 4, 8, 9, 7,  10,
    2, 0, 1, 3,  6, 4, 5, 9, 7, 8,  10,
    0, 3, 1, 2,  7, 8, 4, 6, 9, 5,  10,
    3, 1, 0, 2,  8, 4, 7, 9, 5, 6,  10,
    1, 0, 3, 2,  4, 7, 8, 5, 6, 9,  10,
    0, 2, 3, 1,  6, 9, 7, 4, 5, 8,  10,
    2, 3, 0, 1,  9, 7, 6, 5, 8, 4,  10,
    3, 0, 2, 1,  7, 6, 9, 8, 4, 5,  10,
    1, 3, 2, 0,  8, 9, 5, 4, 7, 6,  10,
    3, 2, 1, 0,  9, 5, 8, 7, 6, 4,  10,
    2, 1, 3, 0,  5, 8, 9, 6, 4, 7,  10
  };
};

//***************************************************************************
// topology::PYRAMID
//***************************************************************************
template <>
struct topology_data<topology::PYRAMID_5>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::PYRAMID_5;
  static constexpr topology::topology_t base = topology::PYRAMID_5;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_nodes = 5;
  static constexpr unsigned num_vertices = 5;
  static constexpr unsigned num_edges = 8;
  static constexpr unsigned num_faces = 5;
  static constexpr unsigned num_permutations = 4;
  static constexpr unsigned num_positive_permutations = 4;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       false,  // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_3,
                                                                  topology::TRI_3,
                                                                  topology::TRI_3,
                                                                  topology::TRI_3,
                                                                  topology::QUAD_4};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 2, 4, 6, 8, 10, 12, 14, 16};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,
                                                           1, 2,
                                                           2, 3,
                                                           3, 0,
                                                           0, 4,
                                                           1, 4,
                                                           2, 4,
                                                           3, 4};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 3, 6, 9, 12, 16};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 4,
                                                           1, 2, 4,
                                                           2, 3, 4,
                                                           0, 4, 3,
                                                           0, 3, 2, 1};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2, 3, 4,
                                                                  1, 2, 3, 0, 4,
                                                                  2, 3, 0, 1, 4,
                                                                  3, 0, 1, 2, 4};
};

template <>
struct topology_data<topology::PYRAMID_13>
  : public topology_data<topology::PYRAMID_5>
{
  static constexpr topology::topology_t value = topology::PYRAMID_13;
  static constexpr unsigned num_nodes = 13;

  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3};


  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_6,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6,
                                                                  topology::QUAD_8};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 6, 9, 12, 15, 18, 21, 24};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  5,
                                                           1, 2,  6,
                                                           2, 3,  7,
                                                           3, 0,  8,
                                                           0, 4,  9,
                                                           1, 4,  10,
                                                           2, 4,  11,
                                                           3, 4,  12};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 6, 12, 18, 24, 32};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 4,  5, 10,  9,
                                                           1, 2, 4,  6, 11, 10,
                                                           2, 3, 4,  7, 12, 11,
                                                           3, 0, 4,  8,  9, 12,
                                                           0, 3, 2, 1,  8, 7, 6, 5};

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3, 4,  5, 6, 7, 8,   9, 10, 11, 12,
    1, 2, 3, 0, 4,  6, 7, 8, 5,  10, 11, 12,  9,
    2, 3, 0, 1, 4,  7, 8, 5, 6,  11, 12,  9, 10,
    3, 0, 1, 2, 4,  8, 5, 6, 7,  12,  9, 10, 11
  };
};

template <>
struct topology_data<topology::PYRAMID_14>
  : public topology_data<topology::PYRAMID_13>
{
  static constexpr topology::topology_t value = topology::PYRAMID_14;
  static constexpr unsigned num_nodes = 14;

  static constexpr topology::topology_t face_topology_vector[] = {topology::TRI_6,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6,
                                                                  topology::QUAD_9};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 6, 12, 18, 24, 33};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 4,  5, 10,  9,
                                                           1, 2, 4,  6, 11, 10,
                                                           2, 3, 4,  7, 12, 11,
                                                           3, 0, 4,  8,  9, 12,
                                                           0, 3, 2, 1,  8, 7, 6, 5,  13};


  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3, 4,  5, 6, 7, 8,   9, 10, 11, 12,  13,
    1, 2, 3, 0, 4,  6, 7, 8, 5,  10, 11, 12,  9,  13,
    2, 3, 0, 1, 4,  7, 8, 5, 6,  11, 12,  9, 10,  13,
    3, 0, 1, 2, 4,  8, 5, 6, 7,  12,  9, 10, 11,  13
  };
};

//***************************************************************************
// topology::WEDGE
//***************************************************************************
template <>
struct topology_data<topology::WEDGE_6>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::WEDGE_6;
  static constexpr topology::topology_t base = topology::WEDGE_6;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2};
  static constexpr bool has_homogeneous_faces = false;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_nodes = 6;
  static constexpr unsigned num_vertices = 6;
  static constexpr unsigned num_edges = 9;
  static constexpr unsigned num_faces = 5;
  static constexpr unsigned num_permutations = 6;
  static constexpr unsigned num_positive_permutations = 6;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       false,  // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_4,
                                                                  topology::QUAD_4,
                                                                  topology::QUAD_4,
                                                                  topology::TRI_3,
                                                                  topology::TRI_3};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,
                                                           1, 2,
                                                           2, 0,
                                                           3, 4,
                                                           4, 5,
                                                           5, 3,
                                                           0, 3,
                                                           1, 4,
                                                           2, 5};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 4, 8, 12, 15, 18};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 4, 3,
                                                           1, 2, 5, 4,
                                                           0, 3, 5, 2,
                                                           0, 2, 1,
                                                           3, 4, 5};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2, 3, 4, 5,
                                                                  1, 2, 0, 4, 5, 3,
                                                                  2, 0, 1, 5, 3, 4,
                                                                  3, 5, 4, 0, 2, 1,
                                                                  5, 4, 3, 2, 1, 0,
                                                                  4, 3, 5, 1, 0, 2};
};

template <>
struct topology_data<topology::WEDGE_12>
  : public topology_data<topology::WEDGE_6>
{
  static constexpr topology::topology_t value = topology::WEDGE_12;
  static constexpr unsigned num_nodes = 12;

  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2};


  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_6,
                                                                  topology::QUAD_6,
                                                                  topology::QUAD_6,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 6, 9, 12, 15, 18, 20, 22, 24};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  6,
                                                           1, 2,  7,
                                                           2, 0,  8,
                                                           3, 4,  9,
                                                           4, 5,  10,
                                                           5, 3,  11,
                                                           0, 3,
                                                           1, 4,
                                                           2, 5};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 6, 12, 18, 24, 30};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 4, 3,  6, 9,
                                                           1, 2, 5, 4,  7, 10,
                                                           0, 3, 5, 2,  8, 11,
                                                           0, 2, 1,   8, 7, 6,
                                                           3, 4, 5,  9, 10, 11};

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11,
    1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9,
    2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10,
    3, 5, 4, 0, 2, 1,  9, 11, 10,  8,  7,  6,
    5, 4, 3, 2, 1, 0, 11, 10,  9,  7,  6,  8,
    4, 3, 5, 1, 0, 2, 10,  9, 11,  6,  8,  7
  };
};

template <>
struct topology_data<topology::WEDGE_15>
  : public topology_data<topology::WEDGE_6>
{
  static constexpr topology::topology_t value = topology::WEDGE_15;
  static constexpr unsigned num_nodes = 15;

  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3};


  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_8,
                                                                  topology::QUAD_8,
                                                                  topology::QUAD_8,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  6,
                                                           1, 2,  7,
                                                           2, 0,  8,
                                                           3, 4,  12,
                                                           4, 5,  13,
                                                           5, 3,  14,
                                                           0, 3,  9,
                                                           1, 4,  10,
                                                           2, 5,  11};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 8, 16, 24, 30, 36};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 4, 3,  6, 10, 12,  9,
                                                           1, 2, 5, 4,  7, 11, 13, 10,
                                                           0, 3, 5, 2,  9, 14, 11,  8,
                                                           0, 2, 1,   8,  7,  6,
                                                           3, 4, 5,  12, 13, 14};

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9, 13, 14, 12,
    2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10, 14, 12, 13,
    3, 5, 4, 0, 2, 1, 14, 13, 12,  9, 11, 10,  8,  7,  6,
    5, 4, 3, 2, 1, 0, 13, 12, 14, 11, 10,  9,  7,  6,  8,
    4, 3, 5, 1, 0, 2, 12, 14, 13, 10,  9, 11,  6,  8,  7
  };
};

template <>
struct topology_data<topology::WEDGE_18>
  : public topology_data<topology::WEDGE_15>
{
  static constexpr topology::topology_t value = topology::WEDGE_18;
  static constexpr unsigned num_nodes = 18;

  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_9,
                                                                  topology::QUAD_9,
                                                                  topology::QUAD_9,
                                                                  topology::TRI_6,
                                                                  topology::TRI_6};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 9, 18, 27, 33, 39};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 4, 3,  6, 10, 12,  9,  15,
                                                           1, 2, 5, 4,  7, 11, 13, 10,  16,
                                                           0, 3, 5, 2,  9, 14, 11,  8,  17,
                                                           0, 2, 1,   8,  7,  6,
                                                           3, 4, 5,  12, 13, 14};

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
    1, 2, 0, 4, 5, 3,  7,  8,  6, 10, 11,  9, 13, 14, 12, 16, 17, 15,
    2, 0, 1, 5, 3, 4,  8,  6,  7, 11,  9, 10, 14, 12, 13, 17, 15, 16,
    3, 5, 4, 0, 2, 1, 14, 13, 12,  9, 11, 10,  8,  7,  6, 17, 16, 15,
    5, 4, 3, 2, 1, 0, 13, 12, 14, 11, 10,  9,  7,  6,  8, 16, 15, 17,
    4, 3, 5, 1, 0, 2, 12, 14, 13, 10,  9, 11,  6,  8,  7, 15, 17, 16
  };
};

//***************************************************************************
// topology::HEXAHEDRON -- topology::ELEMENT_RANK
// defined on 3d
/*
//   Linear 8-Node Hexahedron node locations.
//
//          7                    6
//           o------------------o
//          /|                 /|
//         / |                / |
//        /  |               /  |
//       /   |              /   |
//      /    |             /    |
//     /     |            /     |
//  4 /      |         5 /      |
//   o------------------o       |
//   |       |          |       |
//   |     3 o----------|-------o 2
//   |      /           |      /
//   |     /            |     /
//   |    /             |    /
//   |   /              |   /
//   |  /               |  /
//   | /                | /
//   |/                 |/
//   o------------------o
//  0                    1
//
//
//   Quadratic 20-Node Hexahedron node locations:
//
//           7         18         6
//            o--------o---------o
//           /|                 /|
//          / |                / |
//         /  |               /  |
//      19o   |            17o   |
//       /  15o             /    o14
//      /     |            /     |
//   4 /      | 16        /      |
//    o---------o--------o 5     |
//    |       |       10 |       |
//    |     3 o-------o--|-------o 2
//    |      /           |      /
//    |     /            |     /
//  12o    /             o13  /
//    |   o11            |   o9
//    |  /               |  /
//    | /                | /
//    |/                 |/
//    o---------o--------o
//   0          8         1
//
//
//   Quadratic 27-Node Hexahedron additional node locations:
//
//
//            x--------x---------x
//           /|                 /|
//          / |                / |
//         /  |   22          /  |
//        x   |    o         x   |
//       /    x       o26   /    x     Node #20 is at centroid of element
//      /     |            /     |
//     /      |           /      |     "QUAD_9" beginning with nodes
//    x---------x--------x       |      0,1,5,4 has node 25 at center....
//    | 23o   |          |   o24 |
//    |       x-------x--|-------x
//    |      /           |      /
//    |     /  25        |     /
//    x    /    o        x    /
//    |   x        o21   |   x
//    |  /               |  /
//    | /                | /
//    |/                 |/
//    x---------x--------x
*/
//***************************************************************************
template <>
struct topology_data<topology::HEX_8>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static constexpr topology::topology_t value = topology::HEX_8;
  static constexpr topology::topology_t base = topology::HEX_8;

  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::FACE_RANK;
  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2,
                                                                  topology::LINE_2};

  static constexpr bool has_homogeneous_faces = true;
  static constexpr bool is_shell = false;
  static constexpr unsigned dimension = 3;
  static constexpr unsigned num_nodes = 8;
  static constexpr unsigned num_vertices = 8;
  static constexpr unsigned num_edges = 12;
  static constexpr unsigned num_faces = 6;
  static constexpr unsigned num_permutations = 24;
  static constexpr unsigned num_positive_permutations = 24;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       false,  // 2d
                                                       true};  // 3d

  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_4,
                                                                  topology::QUAD_4,
                                                                  topology::QUAD_4,
                                                                  topology::QUAD_4,
                                                                  topology::QUAD_4,
                                                                  topology::QUAD_4};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,
                                                           1, 2,
                                                           2, 3,
                                                           3, 0,
                                                           4, 5,
                                                           5, 6,
                                                           6, 7,
                                                           7, 4,
                                                           0, 4,
                                                           1, 5,
                                                           2, 6,
                                                           3, 7};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 4, 8, 12, 16, 20, 24};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 5, 4,
                                                           1, 2, 6, 5,
                                                           2, 3, 7, 6,
                                                           0, 4, 7, 3,
                                                           0, 3, 2, 1,
                                                           4, 5, 6, 7};

  static constexpr unsigned permutation_node_ordinals_vector[] = {0, 1, 2, 3, 4, 5, 6, 7,
                                                                  0, 1, 5, 4, 3, 2, 6, 7,
                                                                  0, 4, 7, 3, 1, 5, 6, 2,
                                                                  1, 2, 3, 0, 5, 6, 7, 4,
                                                                  1, 2, 6, 5, 0, 3, 7, 4,
                                                                  1, 5, 4, 0, 2, 6, 7, 3,
                                                                  2, 3, 0, 1, 6, 7, 4, 5,
                                                                  2, 3, 7, 6, 1, 0, 4, 5,
                                                                  2, 6, 5, 1, 3, 7, 4, 0,
                                                                  3, 0, 1, 2, 7, 4, 5, 6,
                                                                  3, 0, 4, 7, 2, 1, 5, 6,
                                                                  3, 7, 6, 2, 0, 4, 5, 1,
                                                                  4, 0, 1, 5, 7, 3, 2, 6,
                                                                  4, 7, 3, 0, 5, 6, 2, 1,
                                                                  4, 7, 6, 5, 0, 3, 2, 1,
                                                                  5, 1, 2, 6, 4, 0, 3, 7,
                                                                  5, 4, 0, 1, 6, 7, 3, 2,
                                                                  5, 4, 7, 6, 1, 0, 3, 2,
                                                                  6, 2, 3, 7, 5, 1, 0, 4,
                                                                  6, 5, 1, 2, 7, 4, 0, 3,
                                                                  6, 5, 4, 7, 2, 1, 0, 3,
                                                                  7, 3, 0, 4, 6, 2, 1, 5,
                                                                  7, 6, 2, 3, 4, 5, 1, 0,
                                                                  7, 6, 5, 4, 3, 2, 1, 0};
};

template <>
struct topology_data<topology::HEX_20>
  : public topology_data<topology::HEX_8>
{
  static constexpr topology::topology_t value = topology::HEX_20;
  static constexpr unsigned num_nodes = 20;
  static constexpr unsigned num_permutations = 24;
  static constexpr unsigned num_positive_permutations = 24;

  static constexpr topology::topology_t edge_topology_vector[] = {topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3,
                                                                  topology::LINE_3};

  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_8,
                                                                  topology::QUAD_8,
                                                                  topology::QUAD_8,
                                                                  topology::QUAD_8,
                                                                  topology::QUAD_8,
                                                                  topology::QUAD_8};

  static constexpr unsigned edge_node_ordinals_offsets[] = {0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36};
  static constexpr unsigned edge_node_ordinals_vector[] = {0, 1,  8,
                                                           1, 2,  9,
                                                           2, 3,  10,
                                                           3, 0,  11,
                                                           4, 5,  16,
                                                           5, 6,  17,
                                                           6, 7,  18,
                                                           7, 4,  19,
                                                           0, 4,  12,
                                                           1, 5,  13,
                                                           2, 6,  14,
                                                           3, 7,  15};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 8, 16, 24, 32, 40, 48};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 5, 4,   8, 13, 16, 12,
                                                           1, 2, 6, 5,   9, 14, 17, 13,
                                                           2, 3, 7, 6,  10, 15, 18, 14,
                                                           0, 4, 7, 3,  12, 19, 15, 11,
                                                           0, 3, 2, 1,  11, 10,  9,  8,
                                                           4, 5, 6, 7,  16, 17, 18, 19};

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3, 4, 5, 6, 7,   8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    0, 1, 5, 4, 3, 2, 6, 7,   8, 13, 16, 12, 11,  9, 17, 19, 10, 14, 18, 15,
    0, 4, 7, 3, 1, 5, 6, 2,  12, 19, 15, 11,  8, 16, 18, 10, 13, 17, 14,  9,
    1, 2, 3, 0, 5, 6, 7, 4,   9, 10, 11,  8, 13, 14, 15, 12, 17, 18, 19, 16,
    1, 2, 6, 5, 0, 3, 7, 4,   9, 14, 17, 13,  8, 10, 18, 16, 11, 15, 19, 12,
    1, 5, 4, 0, 2, 6, 7, 3,  13, 16, 12,  8,  9, 17, 19, 11, 14, 18, 15, 10,
    2, 3, 0, 1, 6, 7, 4, 5,  10, 11,  8,  9, 14, 15, 12, 13, 18, 19, 16, 17,
    2, 3, 7, 6, 1, 0, 4, 5,  10, 15, 18, 14,  9, 11, 19, 17,  8, 12, 16, 13,
    2, 6, 5, 1, 3, 7, 4, 0,  14, 17, 13,  9, 10, 18, 16,  8, 15, 19, 12, 11,
    3, 0, 1, 2, 7, 4, 5, 6,  11,  8,  9, 10, 15, 12, 13, 14, 19, 16, 17, 18,
    3, 0, 4, 7, 2, 1, 5, 6,  11, 12, 19, 15, 10,  8, 16, 18,  9, 13, 17, 14,
    3, 7, 6, 2, 0, 4, 5, 1,  15, 18, 14, 10, 11, 19, 17,  9, 12, 16, 13,  8,
    4, 0, 1, 5, 7, 3, 2, 6,  12,  8, 13, 16, 19, 11,  9, 17, 15, 10, 14, 18,
    4, 7, 3, 0, 5, 6, 2, 1,  19, 15, 11, 12, 16, 18, 10,  8, 17, 14,  9, 13,
    4, 7, 6, 5, 0, 3, 2, 1,  19, 18, 17, 16, 12, 15, 14, 13, 11, 10,  9,  8,
    5, 1, 2, 6, 4, 0, 3, 7,  13,  9, 14, 17, 16,  8, 10, 18, 12, 11, 15, 19,
    5, 4, 0, 1, 6, 7, 3, 2,  16, 12,  8, 13, 17, 19, 11,  9, 18, 15, 10, 14,
    5, 4, 7, 6, 1, 0, 3, 2,  16, 19, 18, 17, 13, 12, 15, 14,  8, 11, 10,  9,
    6, 2, 3, 7, 5, 1, 0, 4,  14, 10, 15, 18, 17,  9, 11, 19, 13,  8, 12, 16,
    6, 5, 1, 2, 7, 4, 0, 3,  17, 13,  9, 14, 18, 16,  8, 10, 19, 12, 11, 15,
    6, 5, 4, 7, 2, 1, 0, 3,  17, 16, 19, 18, 14, 13, 12, 15,  9,  8, 11, 10,
    7, 3, 0, 4, 6, 2, 1, 5,  15, 11, 12, 19, 18, 10,  8, 16, 14,  9, 13, 17,
    7, 6, 2, 3, 4, 5, 1, 0,  18, 14, 10, 15, 19, 17,  9, 11, 16, 13,  8, 12,
    7, 6, 5, 4, 3, 2, 1, 0,  18, 17, 16, 19, 15, 14, 13, 12, 10,  9,  8, 11
  };
};

template <>
struct topology_data<topology::HEX_27>
  : public topology_data<topology::HEX_20>
{
  static constexpr topology::topology_t value = topology::HEX_27;
  static constexpr unsigned num_nodes = 27;

  static constexpr topology::topology_t face_topology_vector[] = {topology::QUAD_9,
                                                                  topology::QUAD_9,
                                                                  topology::QUAD_9,
                                                                  topology::QUAD_9,
                                                                  topology::QUAD_9,
                                                                  topology::QUAD_9};

  static constexpr unsigned face_node_ordinals_offsets[] = {0, 9, 18, 27, 36, 45, 54};
  static constexpr unsigned face_node_ordinals_vector[] = {0, 1, 5, 4,   8, 13, 16, 12,  25,
                                                           1, 2, 6, 5,   9, 14, 17, 13,  24,
                                                           2, 3, 7, 6,  10, 15, 18, 14,  26,
                                                           0, 4, 7, 3,  12, 19, 15, 11,  23,
                                                           0, 3, 2, 1,  11, 10,  9,  8,  21,
                                                           4, 5, 6, 7,  16, 17, 18, 19,  22};

  static constexpr unsigned permutation_node_ordinals_vector[] = {
    0, 1, 2, 3, 4, 5, 6, 7,   8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  20, 21, 22, 23, 24, 25, 26,
    0, 1, 5, 4, 3, 2, 6, 7,   8, 13, 16, 12, 11,  9, 17, 19, 10, 14, 18, 15,  20, 25, 26, 23, 24, 21, 22,
    0, 4, 7, 3, 1, 5, 6, 2,  12, 19, 15, 11,  8, 16, 18, 10, 13, 17, 14,  9,  20, 23, 24, 21, 22, 25, 26,
    1, 2, 3, 0, 5, 6, 7, 4,   9, 10, 11,  8, 13, 14, 15, 12, 17, 18, 19, 16,  20, 21, 22, 25, 26, 24, 23,
    1, 2, 6, 5, 0, 3, 7, 4,   9, 14, 17, 13,  8, 10, 18, 16, 11, 15, 19, 12,  20, 24, 23, 25, 26, 21, 22,
    1, 5, 4, 0, 2, 6, 7, 3,  13, 16, 12,  8,  9, 17, 19, 11, 14, 18, 15, 10,  20, 25, 26, 21, 22, 24, 23,
    2, 3, 0, 1, 6, 7, 4, 5,  10, 11,  8,  9, 14, 15, 12, 13, 18, 19, 16, 17,  20, 21, 22, 24, 23, 26, 25,
    2, 3, 7, 6, 1, 0, 4, 5,  10, 15, 18, 14,  9, 11, 19, 17,  8, 12, 16, 13,  20, 26, 25, 24, 23, 21, 22,
    2, 6, 5, 1, 3, 7, 4, 0,  14, 17, 13,  9, 10, 18, 16,  8, 15, 19, 12, 11,  20, 24, 23, 21, 22, 26, 25,
    3, 0, 1, 2, 7, 4, 5, 6,  11,  8,  9, 10, 15, 12, 13, 14, 19, 16, 17, 18,  20, 21, 22, 26, 25, 23, 24,
    3, 0, 4, 7, 2, 1, 5, 6,  11, 12, 19, 15, 10,  8, 16, 18,  9, 13, 17, 14,  20, 23, 24, 26, 25, 21, 22,
    3, 7, 6, 2, 0, 4, 5, 1,  15, 18, 14, 10, 11, 19, 17,  9, 12, 16, 13,  8,  20, 26, 25, 21, 22, 23, 24,
    4, 0, 1, 5, 7, 3, 2, 6,  12,  8, 13, 16, 19, 11,  9, 17, 15, 10, 14, 18,  20, 25, 26, 22, 21, 23, 24,
    4, 7, 3, 0, 5, 6, 2, 1,  19, 15, 11, 12, 16, 18, 10,  8, 17, 14,  9, 13,  20, 23, 24, 25, 26, 22, 21,
    4, 7, 6, 5, 0, 3, 2, 1,  19, 18, 17, 16, 12, 15, 14, 13, 11, 10,  9,  8,  20, 22, 21, 25, 26, 23, 24,
    5, 1, 2, 6, 4, 0, 3, 7,  13,  9, 14, 17, 16,  8, 10, 18, 12, 11, 15, 19,  20, 24, 23, 22, 21, 25, 26,
    5, 4, 0, 1, 6, 7, 3, 2,  16, 12,  8, 13, 17, 19, 11,  9, 18, 15, 10, 14,  20, 25, 26, 24, 23, 22, 21,
    5, 4, 7, 6, 1, 0, 3, 2,  16, 19, 18, 17, 13, 12, 15, 14,  8, 11, 10,  9,  20, 22, 21, 24, 23, 25, 26,
    6, 2, 3, 7, 5, 1, 0, 4,  14, 10, 15, 18, 17,  9, 11, 19, 13,  8, 12, 16,  20, 26, 25, 22, 21, 24, 23,
    6, 5, 1, 2, 7, 4, 0, 3,  17, 13,  9, 14, 18, 16,  8, 10, 19, 12, 11, 15,  20, 24, 23, 26, 25, 22, 21,
    6, 5, 4, 7, 2, 1, 0, 3,  17, 16, 19, 18, 14, 13, 12, 15,  9,  8, 11, 10,  20, 22, 21, 26, 25, 24, 23,
    7, 3, 0, 4, 6, 2, 1, 5,  15, 11, 12, 19, 18, 10,  8, 16, 14,  9, 13, 17,  20, 23, 24, 22, 21, 26, 25,
    7, 6, 2, 3, 4, 5, 1, 0,  18, 14, 10, 15, 19, 17,  9, 11, 16, 13,  8, 12,  20, 26, 25, 23, 24, 22, 21,
    7, 6, 5, 4, 3, 2, 1, 0,  18, 17, 16, 19, 15, 14, 13, 12, 10,  9,  8, 11,  20, 22, 21, 23, 24, 26, 25
  };
};

//***************************************************************************
// topology::SUPERFACE -- topology::SUPERFACE
//***************************************************************************

template <topology::topology_t Topology>
struct topology_data<Topology, typename std::enable_if< (Topology > topology::SUPEREDGE_START && Topology < topology::SUPEREDGE_END) >::type >
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  static constexpr topology::topology_t value = Topology;
  static constexpr topology::topology_t base = value;
  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::EDGE_RANK;
  static constexpr topology::rank_t side_rank = topology::INVALID_RANK;
  static constexpr unsigned num_nodes = Topology - topology::SUPEREDGE_START;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       true,   // 2d
                                                       true};  // 3d
};

template <topology::topology_t Topology>
struct topology_data<Topology, typename std::enable_if< (Topology > topology::SUPERFACE_START && Topology < topology::SUPERFACE_END) >::type >
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  static constexpr topology::topology_t value = Topology;
  static constexpr topology::topology_t base = value;
  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::FACE_RANK;
  static constexpr topology::rank_t side_rank = topology::INVALID_RANK;
  static constexpr unsigned num_nodes = Topology - topology::SUPERFACE_START;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       false,  // 1d
                                                       false,  // 2d
                                                       true};  // 3d
};

//***************************************************************************
// topology::SUPERELEMENT -- topology::SUPERELEMENT
//***************************************************************************

template <topology::topology_t Topology>
struct topology_data<Topology, typename std::enable_if< (Topology > topology::SUPERELEMENT_START) >::type >
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  static constexpr topology::topology_t value = Topology;
  static constexpr topology::topology_t base = value;
  static constexpr bool is_valid = true;
  static constexpr topology::rank_t rank = topology::ELEMENT_RANK;
  static constexpr topology::rank_t side_rank = topology::INVALID_RANK;
  static constexpr unsigned num_nodes = Topology - topology::SUPERELEMENT_START;

  static constexpr bool spatial_dimension_vector[4] = {false,  // 0d
                                                       true,   // 1d
                                                       true,   // 2d
                                                       true};  // 3d
};

}} // namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_TOPOLOGY_DATA_HPP
