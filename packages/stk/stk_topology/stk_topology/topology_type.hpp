// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
#ifndef STKTOPOLOGY_TOPOLOGY_TYPE_HPP
#define STKTOPOLOGY_TOPOLOGY_TYPE_HPP

// IWYU pragma: private, include "stk_topology/topology.hpp"

#include <stk_topology/topology.hpp>

#include <stk_topology/topology_detail/fill_container.hpp>
#include <stk_topology/topology_detail/topology_data.hpp>
#include <stk_topology/topology_detail/meta_functions.hpp>
#include <stk_topology/topology_detail/equivalent_helper.hpp>

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/assert.hpp>

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
  static const topology_t edge_topology       = data::edge_topology;
  static const bool has_homogeneous_faces     = data::has_homogeneous_faces;
  static const bool is_shell                  = data::is_shell;
  static const unsigned dimension                  = data::dimension;
  static const unsigned num_nodes                  = data::num_nodes;
  static const unsigned num_vertices               = data::num_vertices;
  static const unsigned num_edges                  = data::num_edges;
  static const unsigned num_faces                  = data::num_faces;
  static const unsigned num_permutations           = data::num_permutations;
  static const unsigned num_positive_permutations  = data::num_positive_permutations;

  typedef typename data::spatial_dimension_vector                        spatial_dimension_vector;
  typedef typename data::face_topology_vector                            face_topology_vector;
  typedef typename data::edge_node_ordinals_vector                       edge_node_ordinals_vector;
  typedef typename data::face_node_ordinals_vector                       face_node_ordinals_vector;
  typedef typename data::permutation_node_ordinals_vector                permutation_node_ordinals_vector;



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
    case 1: return topology_detail::defined_on_spatial_dimension_< type, 1>::value;
    case 2: return topology_detail::defined_on_spatial_dimension_< type, 2>::value;
    case 3: return topology_detail::defined_on_spatial_dimension_< type, 3>::value;
    default: break;
    }
    return false;
  }

  /// the topology of the face at the given ordinal
  STK_FUNCTION
  static topology face_topology(unsigned face_ordinal = 0)
  {
    switch (face_ordinal)
    {
      case 0:  return topology_detail::face_topology_<type,0 >::value;
      case 1:  return topology_detail::face_topology_<type,1 >::value;
      case 2:  return topology_detail::face_topology_<type,2 >::value;
      case 3:  return topology_detail::face_topology_<type,3 >::value;
      case 4:  return topology_detail::face_topology_<type,4 >::value;
      case 5:  return topology_detail::face_topology_<type,5 >::value;
      default: break;
    }

    return INVALID_TOPOLOGY;
  }

  /// node ordinals that make up the given edge
  template <typename OrdinalOutputIterator>
  static void edge_node_ordinals(unsigned edge_ordinal, OrdinalOutputIterator output_ordinals)
  {
    topology_detail::fill_ordinal_container<OrdinalOutputIterator> f(output_ordinals);

    switch (edge_ordinal)
    {
      case 0:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,0 >::type >( f ); break;
      case 1:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,1 >::type >( f ); break;
      case 2:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,2 >::type >( f ); break;
      case 3:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,3 >::type >( f ); break;
      case 4:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,4 >::type >( f ); break;
      case 5:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,5 >::type >( f ); break;
      case 6:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,6 >::type >( f ); break;
      case 7:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,7 >::type >( f ); break;
      case 8:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,8 >::type >( f ); break;
      case 9:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,9 >::type >( f ); break;
      case 10: boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,10>::type >( f ); break;
      case 11: boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,11>::type >( f ); break;
      default: break;
    }

    return;
  }

  /// the node ordinals that make up the given face
  template <typename OrdinalOutputIterator>
  static void face_node_ordinals(unsigned face_ordinal, OrdinalOutputIterator output_ordinals)
  {
    topology_detail::fill_ordinal_container<OrdinalOutputIterator> f(output_ordinals);

    switch (face_ordinal)
    {
      case 0:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,0 >::type >( f ); break;
      case 1:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,1 >::type >( f ); break;
      case 2:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,2 >::type >( f ); break;
      case 3:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,3 >::type >( f ); break;
      case 4:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,4 >::type >( f ); break;
      case 5:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,5 >::type >( f ); break;
      default: break;
    }

    return;
  }

  /// the node ordinals of the topology in the given permutation order
  template <typename OrdinalOutputIterator>
  static void permutation_node_ordinals(unsigned permutation_ordinal, OrdinalOutputIterator output_ordinals)
  {
    topology_detail::fill_ordinal_container<OrdinalOutputIterator> f(output_ordinals);

    switch (permutation_ordinal)
    {
      case 0:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,0 >::type >( f ); break;
      case 1:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,1 >::type >( f ); break;
      case 2:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,2 >::type >( f ); break;
      case 3:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,3 >::type >( f ); break;
      case 4:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,4 >::type >( f ); break;
      case 5:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,5 >::type >( f ); break;
      case 6:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,6 >::type >( f ); break;
      case 7:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,7 >::type >( f ); break;
      case 8:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,8 >::type >( f ); break;
      case 9:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,9 >::type >( f ); break;
      case 10: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,10>::type >( f ); break;
      case 11: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,11>::type >( f ); break;
      case 12: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,12>::type >( f ); break;
      case 13: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,13>::type >( f ); break;
      case 14: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,14>::type >( f ); break;
      case 15: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,15>::type >( f ); break;
      case 16: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,16>::type >( f ); break;
      case 17: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,17>::type >( f ); break;
      case 18: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,18>::type >( f ); break;
      case 19: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,19>::type >( f ); break;
      case 20: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,20>::type >( f ); break;
      case 21: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,21>::type >( f ); break;
      case 22: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,22>::type >( f ); break;
      case 23: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,23>::type >( f ); break;
      default: break;
    }

    return;
  }

  /// node that make up the given edge
  template <typename NodeArray, typename NodeOutputIterator>
  static void edge_nodes(const NodeArray & nodes, unsigned edge_ordinal, NodeOutputIterator output_nodes)
  {
    topology_detail::fill_node_container<NodeArray,NodeOutputIterator> f(nodes,output_nodes);

    switch (edge_ordinal)
    {
      case 0:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,0 >::type >( f ); break;
      case 1:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,1 >::type >( f ); break;
      case 2:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,2 >::type >( f ); break;
      case 3:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,3 >::type >( f ); break;
      case 4:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,4 >::type >( f ); break;
      case 5:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,5 >::type >( f ); break;
      case 6:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,6 >::type >( f ); break;
      case 7:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,7 >::type >( f ); break;
      case 8:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,8 >::type >( f ); break;
      case 9:  boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,9 >::type >( f ); break;
      case 10: boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,10>::type >( f ); break;
      case 11: boost::mpl::for_each< typename topology_detail::edge_node_ordinals_<type,11>::type >( f ); break;
      default: break;
    }

    return;
  }

  /// node that make up the given face
  template <typename NodeArray, typename NodeOutputIterator>
  static void face_nodes(const NodeArray & nodes, unsigned face_ordinal, NodeOutputIterator output_nodes)
  {
    topology_detail::fill_node_container<NodeArray,NodeOutputIterator> f(nodes,output_nodes);

    switch (face_ordinal)
    {
      case 0:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,0 >::type >( f ); break;
      case 1:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,1 >::type >( f ); break;
      case 2:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,2 >::type >( f ); break;
      case 3:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,3 >::type >( f ); break;
      case 4:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,4 >::type >( f ); break;
      case 5:  boost::mpl::for_each< typename topology_detail::face_node_ordinals_<type,5 >::type >( f ); break;
      default: break;
    }

    return;
  }

  /// node that make up the given permutation
  template <typename NodeArray, typename NodeOutputIterator>
  static void permutation_nodes(const NodeArray & nodes, unsigned permutation_ordinal, NodeOutputIterator output_nodes)
  {
    topology_detail::fill_node_container<NodeArray,NodeOutputIterator> f(nodes,output_nodes);

    switch (permutation_ordinal)
    {
      case 0:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,0 >::type >( f ); break;
      case 1:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,1 >::type >( f ); break;
      case 2:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,2 >::type >( f ); break;
      case 3:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,3 >::type >( f ); break;
      case 4:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,4 >::type >( f ); break;
      case 5:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,5 >::type >( f ); break;
      case 6:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,6 >::type >( f ); break;
      case 7:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,7 >::type >( f ); break;
      case 8:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,8 >::type >( f ); break;
      case 9:  boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,9 >::type >( f ); break;
      case 10: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,10>::type >( f ); break;
      case 11: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,11>::type >( f ); break;
      case 12: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,12>::type >( f ); break;
      case 13: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,13>::type >( f ); break;
      case 14: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,14>::type >( f ); break;
      case 15: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,15>::type >( f ); break;
      case 16: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,16>::type >( f ); break;
      case 17: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,17>::type >( f ); break;
      case 18: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,18>::type >( f ); break;
      case 19: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,19>::type >( f ); break;
      case 20: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,20>::type >( f ); break;
      case 21: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,21>::type >( f ); break;
      case 22: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,22>::type >( f ); break;
      case 23: boost::mpl::for_each< typename topology_detail::permutation_node_ordinals_<type,23>::type >( f ); break;
      default: break;
    }

    return;
  }

  /// fill the output ordinals with the ordinals that make up the given sub topology
  template <typename OrdinalOutputIterator>
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
  static void sub_topology_nodes(const NodeArray & nodes, unsigned sub_rank, unsigned sub_ordinal, NodeOutputIterator output_nodes)
  {
    switch(sub_rank)
    {
    case NODE_RANK: *output_nodes = nodes[sub_ordinal];                    break;
    case EDGE_RANK: edge_node_ordinals(nodes, sub_ordinal, output_nodes);  break;
    case FACE_RANK: face_node_ordinals(nodes, sub_ordinal, output_nodes);  break;
    default: break;
    }
  }

  /// how many 'sub topologies' does this topology have
  static unsigned num_sub_topology(unsigned sub_rank)
  {
    switch(sub_rank)
    {
    case NODE_RANK: return num_vertices;
    case EDGE_RANK: return num_edges;
    case FACE_RANK: return num_faces;
    default: break;
    }
    return 0;

  }


  /// what is the topology of the given sub topology
  static topology sub_topology(unsigned sub_rank, unsigned sub_ordinal = 0)
  {
    switch(sub_rank)
    {
    case NODE_RANK: return NODE;
    case EDGE_RANK: return edge_topology;
    case FACE_RANK: return face_topology(sub_ordinal);
    default: break;
    }
    return INVALID_TOPOLOGY;
  }

  /// do the two arrays defined equivalent entities (same nodes, but maybe a different permutation)
  /// return a pair<bool, permutation_number> bool and permutation number from a to b
  template <typename NodeArrayA, typename NodeArrayB>
  static std::pair<bool,unsigned> equivalent(const NodeArrayA & a, const NodeArrayB & b)
  {
    return topology_detail::equivalent_helper(type(),a,b,a[0]);
  }

  template <typename NodeArray>
  static unsigned lexicographical_smallest_permutation( const NodeArray & nodes, bool only_positive_permutations = false)
  {
    return topology_detail::lexicographical_smallest_permutation_helper( type(), nodes, only_positive_permutations, nodes[0]);
  }

  template <typename NodeArray>
  static unsigned lexicographical_smallest_permutation_preserve_polarity( const NodeArray & nodes, const NodeArray & element_nodes)
  {
    return topology_detail::lexicographical_smallest_permutation_preserve_polarity_helper( type(), nodes, element_nodes, nodes[0]);
  }

  operator topology_t() const
  { return Topology; }

  /// fill the output ordinals with the ordinals that make up the given side topology
  template <typename OrdinalOutputIterator>
  static void side_node_ordinals(unsigned side_ordinal, OrdinalOutputIterator output_ordinals)
  {
    sub_topology_node_ordinals( side_rank, side_ordinal, output_ordinals);
  }

  /// fill the output nodes with the nodes that make up the given side topology
  /// input 'nodes' is expected to be of length num_nodes.
  template <typename NodeArray, typename NodeOutputIterator>
  static void side_nodes(const NodeArray & nodes, unsigned side_ordinal, NodeOutputIterator output_nodes)
  {
    sub_topology_nodes( nodes, side_rank, side_ordinal, output_nodes);
  }

  /// how many 'side topologies' does this topology have
  static unsigned num_sides()
  {
    return num_sub_topology(side_rank);
  }


  /// what is the topology of the given side topology
  static topology side_topology(unsigned side_ordinal = 0)
  {
    return sub_topology(side_rank, side_ordinal);
  }

};

} //namespace stk

#undef STKTOPOLOGY_META_FUNCTION_SWITCH_WITH_FOR_EACH_FUNCTOR

#endif //STKTOPOLOGY_TOPOLOGY_TYPE_TCC
