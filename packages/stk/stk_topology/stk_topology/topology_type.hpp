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
#ifndef STKTOPOLOGY_TOPOLOGY_TYPE_HPP
#define STKTOPOLOGY_TOPOLOGY_TYPE_HPP

// IWYU pragma: private, include "stk_topology/topology.hpp"

#include <stk_topology/topology.hpp>

#include <stk_topology/topology_detail/fill_container.hpp>
#include <stk_topology/topology_detail/topology_data.hpp>
#include <stk_topology/topology_detail/meta_functions.hpp>
#include <stk_topology/topology_detail/equivalent_helper.hpp>

namespace stk {

//******************************************************************************
// struct topology::topology_type<Topology> is the compile time topology
//
// provides simple runtime methods for edge, face, and permutation information
//******************************************************************************
template <topology::topology_t Topology>
struct topology::topology_type
{
  STK_FUNCTION topology_type() {}

  typedef topology_detail::topology_data<Topology> data;
  typedef topology_type<Topology>         type;
  typedef topology_t                      value_type;

  static const topology_t value               = Topology;
  static const topology_t base                = data::base;
  static const bool is_valid                  = data::is_valid;
  static const rank_t rank                    = data::rank;
  static const rank_t side_rank               = data::side_rank;
  static const bool has_homogeneous_faces     = data::has_homogeneous_faces;
  static const bool is_shell                  = data::is_shell;
  static const unsigned dimension                  = data::dimension;
  static const unsigned num_nodes                  = data::num_nodes;
  static const unsigned num_vertices               = data::num_vertices;
  static const unsigned num_edges                  = data::num_edges;
  static const unsigned num_faces                  = data::num_faces;
  static const unsigned num_permutations           = data::num_permutations;
  static const unsigned num_positive_permutations  = data::num_positive_permutations;


  //***************************************************************************
  //static member functions
  //***************************************************************************

  /// name of the current topology
  static std::string name() { return topology(Topology).name(); }

  /// is the current topology defined on the given spatial dimension
  STK_FUNCTION
  static bool defined_on_spatial_dimension(unsigned spatial_dimension)
  {
    switch(spatial_dimension)
    {
    case 1: return topology_detail::defined_on_spatial_dimension_<data, 1>();
    case 2: return topology_detail::defined_on_spatial_dimension_<data, 2>();
    case 3: return topology_detail::defined_on_spatial_dimension_<data, 3>();
    default: break;
    }
    return false;
  }

  /// the topology of the edge at the given ordinal
  STK_FUNCTION
  static topology edge_topology(unsigned edge_ordinal = 0)
  {
    switch (edge_ordinal)
    {
      case 0:   return topology_detail::edge_topology_<data, 0 >();
      case 1:   return topology_detail::edge_topology_<data, 1 >();
      case 2:   return topology_detail::edge_topology_<data, 2 >();
      case 3:   return topology_detail::edge_topology_<data, 3 >();
      case 4:   return topology_detail::edge_topology_<data, 4 >();
      case 5:   return topology_detail::edge_topology_<data, 5 >();
      case 6:   return topology_detail::edge_topology_<data, 6 >();
      case 7:   return topology_detail::edge_topology_<data, 7 >();
      case 8:   return topology_detail::edge_topology_<data, 8 >();
      case 9:   return topology_detail::edge_topology_<data, 9 >();
      case 10:  return topology_detail::edge_topology_<data, 10>();
      case 11:  return topology_detail::edge_topology_<data, 11>();
      default: break;
    }

    return INVALID_TOPOLOGY;
  }

  /// the topology of the face at the given ordinal
  STK_FUNCTION
  static topology face_topology(unsigned face_ordinal = 0)
  {
    switch (face_ordinal)
    {
      case 0:  return topology_detail::face_topology_<data, 0>();
      case 1:  return topology_detail::face_topology_<data, 1>();
      case 2:  return topology_detail::face_topology_<data, 2>();
      case 3:  return topology_detail::face_topology_<data, 3>();
      case 4:  return topology_detail::face_topology_<data, 4>();
      case 5:  return topology_detail::face_topology_<data, 5>();
      default: break;
    }

    return INVALID_TOPOLOGY;
  }

  /// node ordinals that make up the given edge
  template <typename OrdinalOutputIterator>
  STK_FUNCTION
  static void edge_node_ordinals(unsigned edge_ordinal, OrdinalOutputIterator output_ordinals)
  {
    topology_detail::fill_ordinal_container<OrdinalOutputIterator> f(output_ordinals);

    switch (edge_ordinal)
    {
      case 0:  topology_detail::edge_node_ordinals_<data, 0 >(f); break;
      case 1:  topology_detail::edge_node_ordinals_<data, 1 >(f); break;
      case 2:  topology_detail::edge_node_ordinals_<data, 2 >(f); break;
      case 3:  topology_detail::edge_node_ordinals_<data, 3 >(f); break;
      case 4:  topology_detail::edge_node_ordinals_<data, 4 >(f); break;
      case 5:  topology_detail::edge_node_ordinals_<data, 5 >(f); break;
      case 6:  topology_detail::edge_node_ordinals_<data, 6 >(f); break;
      case 7:  topology_detail::edge_node_ordinals_<data, 7 >(f); break;
      case 8:  topology_detail::edge_node_ordinals_<data, 8 >(f); break;
      case 9:  topology_detail::edge_node_ordinals_<data, 9 >(f); break;
      case 10: topology_detail::edge_node_ordinals_<data, 10>(f); break;
      case 11: topology_detail::edge_node_ordinals_<data, 11>(f); break;
      default: break;
    }

    return;
  }

  /// the node ordinals that make up the given face
  template <typename OrdinalOutputIterator>
  STK_FUNCTION
  static void face_node_ordinals(unsigned face_ordinal, OrdinalOutputIterator output_ordinals)
  {
    topology_detail::fill_ordinal_container<OrdinalOutputIterator> f(output_ordinals);

    switch (face_ordinal)
    {
      case 0:  topology_detail::face_node_ordinals_<data, 0>(f); break;
      case 1:  topology_detail::face_node_ordinals_<data, 1>(f); break;
      case 2:  topology_detail::face_node_ordinals_<data, 2>(f); break;
      case 3:  topology_detail::face_node_ordinals_<data, 3>(f); break;
      case 4:  topology_detail::face_node_ordinals_<data, 4>(f); break;
      case 5:  topology_detail::face_node_ordinals_<data, 5>(f); break;
      default: break;
    }

    return;
  }

  /// the node ordinals of the topology in the given permutation order
  template <typename OrdinalOutputIterator>
  STK_FUNCTION
  static void permutation_node_ordinals(unsigned permutation_ordinal, OrdinalOutputIterator output_ordinals)
  {
    topology_detail::fill_ordinal_container<OrdinalOutputIterator> f(output_ordinals);

    switch (permutation_ordinal)
    {
      case 0:  topology_detail::permutation_node_ordinals_<data, 0 >(f); break;
      case 1:  topology_detail::permutation_node_ordinals_<data, 1 >(f); break;
      case 2:  topology_detail::permutation_node_ordinals_<data, 2 >(f); break;
      case 3:  topology_detail::permutation_node_ordinals_<data, 3 >(f); break;
      case 4:  topology_detail::permutation_node_ordinals_<data, 4 >(f); break;
      case 5:  topology_detail::permutation_node_ordinals_<data, 5 >(f); break;
      case 6:  topology_detail::permutation_node_ordinals_<data, 6 >(f); break;
      case 7:  topology_detail::permutation_node_ordinals_<data, 7 >(f); break;
      case 8:  topology_detail::permutation_node_ordinals_<data, 8 >(f); break;
      case 9:  topology_detail::permutation_node_ordinals_<data, 9 >(f); break;
      case 10: topology_detail::permutation_node_ordinals_<data, 10>(f); break;
      case 11: topology_detail::permutation_node_ordinals_<data, 11>(f); break;
      case 12: topology_detail::permutation_node_ordinals_<data, 12>(f); break;
      case 13: topology_detail::permutation_node_ordinals_<data, 13>(f); break;
      case 14: topology_detail::permutation_node_ordinals_<data, 14>(f); break;
      case 15: topology_detail::permutation_node_ordinals_<data, 15>(f); break;
      case 16: topology_detail::permutation_node_ordinals_<data, 16>(f); break;
      case 17: topology_detail::permutation_node_ordinals_<data, 17>(f); break;
      case 18: topology_detail::permutation_node_ordinals_<data, 18>(f); break;
      case 19: topology_detail::permutation_node_ordinals_<data, 19>(f); break;
      case 20: topology_detail::permutation_node_ordinals_<data, 20>(f); break;
      case 21: topology_detail::permutation_node_ordinals_<data, 21>(f); break;
      case 22: topology_detail::permutation_node_ordinals_<data, 22>(f); break;
      case 23: topology_detail::permutation_node_ordinals_<data, 23>(f); break;
      default: break;
    }

    return;
  }

  /// node that make up the given edge
  template <typename NodeArray, typename NodeOutputIterator>
  STK_FUNCTION
  static void edge_nodes(const NodeArray & nodes, unsigned edge_ordinal, NodeOutputIterator output_nodes)
  {
    topology_detail::fill_node_container<NodeArray, NodeOutputIterator> f(nodes, output_nodes);

    switch (edge_ordinal)
    {
      case 0:  topology_detail::edge_node_ordinals_<data, 0 >(f); break;
      case 1:  topology_detail::edge_node_ordinals_<data, 1 >(f); break;
      case 2:  topology_detail::edge_node_ordinals_<data, 2 >(f); break;
      case 3:  topology_detail::edge_node_ordinals_<data, 3 >(f); break;
      case 4:  topology_detail::edge_node_ordinals_<data, 4 >(f); break;
      case 5:  topology_detail::edge_node_ordinals_<data, 5 >(f); break;
      case 6:  topology_detail::edge_node_ordinals_<data, 6 >(f); break;
      case 7:  topology_detail::edge_node_ordinals_<data, 7 >(f); break;
      case 8:  topology_detail::edge_node_ordinals_<data, 8 >(f); break;
      case 9:  topology_detail::edge_node_ordinals_<data, 9 >(f); break;
      case 10: topology_detail::edge_node_ordinals_<data, 10>(f); break;
      case 11: topology_detail::edge_node_ordinals_<data, 11>(f); break;
      default: break;
    }

    return;
  }

  /// node that make up the given face
  template <typename NodeArray, typename NodeOutputIterator>
  STK_FUNCTION
  static void face_nodes(const NodeArray & nodes, unsigned face_ordinal, NodeOutputIterator output_nodes)
  {
    topology_detail::fill_node_container<NodeArray, NodeOutputIterator> f(nodes, output_nodes);

    switch (face_ordinal)
    {
      case 0:  topology_detail::face_node_ordinals_<data, 0>(f); break;
      case 1:  topology_detail::face_node_ordinals_<data, 1>(f); break;
      case 2:  topology_detail::face_node_ordinals_<data, 2>(f); break;
      case 3:  topology_detail::face_node_ordinals_<data, 3>(f); break;
      case 4:  topology_detail::face_node_ordinals_<data, 4>(f); break;
      case 5:  topology_detail::face_node_ordinals_<data, 5>(f); break;
      default: break;
    }

    return;
  }

  /// node that make up the given permutation
  template <typename NodeArray, typename NodeOutputIterator>
  STK_FUNCTION
  static void permutation_nodes(const NodeArray & nodes, unsigned permutation_ordinal, NodeOutputIterator output_nodes)
  {
    topology_detail::fill_node_container<NodeArray,NodeOutputIterator> f(nodes,output_nodes);

    switch (permutation_ordinal)
    {
      case 0:  topology_detail::permutation_node_ordinals_<data, 0 >(f); break;
      case 1:  topology_detail::permutation_node_ordinals_<data, 1 >(f); break;
      case 2:  topology_detail::permutation_node_ordinals_<data, 2 >(f); break;
      case 3:  topology_detail::permutation_node_ordinals_<data, 3 >(f); break;
      case 4:  topology_detail::permutation_node_ordinals_<data, 4 >(f); break;
      case 5:  topology_detail::permutation_node_ordinals_<data, 5 >(f); break;
      case 6:  topology_detail::permutation_node_ordinals_<data, 6 >(f); break;
      case 7:  topology_detail::permutation_node_ordinals_<data, 7 >(f); break;
      case 8:  topology_detail::permutation_node_ordinals_<data, 8 >(f); break;
      case 9:  topology_detail::permutation_node_ordinals_<data, 9 >(f); break;
      case 10: topology_detail::permutation_node_ordinals_<data, 10>(f); break;
      case 11: topology_detail::permutation_node_ordinals_<data, 11>(f); break;
      case 12: topology_detail::permutation_node_ordinals_<data, 12>(f); break;
      case 13: topology_detail::permutation_node_ordinals_<data, 13>(f); break;
      case 14: topology_detail::permutation_node_ordinals_<data, 14>(f); break;
      case 15: topology_detail::permutation_node_ordinals_<data, 15>(f); break;
      case 16: topology_detail::permutation_node_ordinals_<data, 16>(f); break;
      case 17: topology_detail::permutation_node_ordinals_<data, 17>(f); break;
      case 18: topology_detail::permutation_node_ordinals_<data, 18>(f); break;
      case 19: topology_detail::permutation_node_ordinals_<data, 19>(f); break;
      case 20: topology_detail::permutation_node_ordinals_<data, 20>(f); break;
      case 21: topology_detail::permutation_node_ordinals_<data, 21>(f); break;
      case 22: topology_detail::permutation_node_ordinals_<data, 22>(f); break;
      case 23: topology_detail::permutation_node_ordinals_<data, 23>(f); break;
      default: break;
    }

    return;
  }

  /// fill the output ordinals with the ordinals that make up the given sub topology
  template <typename OrdinalOutputIterator>
  STK_FUNCTION
  static void sub_topology_node_ordinals(unsigned sub_rank, unsigned sub_ordinal, OrdinalOutputIterator output_ordinals)
  {
    switch(sub_rank)
    {
    case NODE_RANK: *output_ordinals = sub_ordinal;                    break;
    case EDGE_RANK: edge_node_ordinals(sub_ordinal, output_ordinals);  break;
    case FACE_RANK: face_node_ordinals(sub_ordinal, output_ordinals);  break;
    default: break;
    }
  }

  /// fill the output nodes with the nodes that make up the given sub topology
  template <typename NodeArray, typename NodeOutputIterator>
  STK_FUNCTION
  static void sub_topology_nodes(const NodeArray & nodes, unsigned sub_rank, unsigned sub_ordinal, NodeOutputIterator output_nodes)
  {
    switch(sub_rank)
    {
    case NODE_RANK: *output_nodes = nodes[sub_ordinal];            break;
    case EDGE_RANK: edge_nodes(nodes, sub_ordinal, output_nodes);  break;
    case FACE_RANK: face_nodes(nodes, sub_ordinal, output_nodes);  break;
    default: break;
    }
  }

  /// do the two arrays defined equivalent entities (same nodes, but maybe a different permutation)
  /// return a struct containing a bool and permutation number from a to b
  template <typename NodeArrayA, typename NodeArrayB>
  STK_FUNCTION
  static EquivalentPermutation is_equivalent(const NodeArrayA & a, const NodeArrayB & b)
  {
    return topology_detail::is_equivalent_helper(type(), a, b, a[0]);
  }

  template <typename NodeArray>
  STK_FUNCTION
  static unsigned lexicographical_smallest_permutation( const NodeArray & nodes, bool only_positive_permutations = false)
  {
    return topology_detail::lexicographical_smallest_permutation_helper( type(), nodes, only_positive_permutations, nodes[0]);
  }

  template <typename NodeArray>
  STK_FUNCTION
  static unsigned lexicographical_smallest_permutation_preserve_polarity( const NodeArray & nodes, const NodeArray & element_nodes)
  {
    return topology_detail::lexicographical_smallest_permutation_preserve_polarity_helper( type(), nodes, element_nodes, nodes[0]);
  }

  STK_FUNCTION
  operator topology_t() const
  { return Topology; }

};

} //namespace stk

#undef STKTOPOLOGY_META_FUNCTION_SWITCH_WITH_FOR_EACH_FUNCTOR

#endif //STKTOPOLOGY_TOPOLOGY_TYPE_TCC
