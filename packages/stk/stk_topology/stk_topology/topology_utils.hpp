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

#ifndef STKTOPOLOGY_TOPOLOGY_TCC
#define STKTOPOLOGY_TOPOLOGY_TCC

// IWYU pragma: private, include "stk_topology/topology.hpp"
#include "stk_topology/topology_decl.hpp"
#include "stk_util/stk_config.h"
#include <type_traits>

namespace stk { namespace topology_detail {

struct num_nodes_impl {
  using result_type = unsigned;
  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_nodes; }
};

struct rank_impl {
  using result_type = topology::rank_t;
  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::rank; }
};

template <typename NodeArrayA, typename NodeArrayB>
struct is_equivalent_impl {
  using result_type = EquivalentPermutation;

  KOKKOS_FUNCTION
  is_equivalent_impl( const NodeArrayA &a , const NodeArrayB &b )
    : m_a(a), m_b(b)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::is_equivalent(m_a, m_b); }

  const NodeArrayA & m_a;
  const NodeArrayB & m_b;
};

template <typename NodeArray>
struct lexicographical_smallest_permutation_impl {
  using result_type = unsigned;

  KOKKOS_FUNCTION
  lexicographical_smallest_permutation_impl( const NodeArray &nodes , bool only_positive_permutations )
    : m_nodes(nodes), m_only_positive_permutations(only_positive_permutations)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const
  { return Topology::lexicographical_smallest_permutation(m_nodes, m_only_positive_permutations); }

  const NodeArray & m_nodes;
  bool              m_only_positive_permutations;
};

template <typename NodeArray>
struct lexicographical_smallest_permutation_preserve_polarity_impl {
  using result_type = unsigned;

  KOKKOS_FUNCTION
  lexicographical_smallest_permutation_preserve_polarity_impl( const NodeArray &nodes, const NodeArray &element_nodes)
    : m_nodes(nodes), m_element_nodes(element_nodes)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const
  { return Topology::lexicographical_smallest_permutation_preserve_polarity(m_nodes, m_element_nodes); }

  const NodeArray & m_nodes;
  const NodeArray & m_element_nodes;
};

// temporary helpers until shell side conversions are fully complete
template <typename Topology, typename = void>
struct HasMixedRankSides : std::false_type {};

template <typename Topology>
struct HasMixedRankSides<Topology, decltype(void(Topology::has_mixed_rank_sides))> : std::true_type {};

struct has_mixed_rank_sides_impl {
  using result_type = bool;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const {
    if constexpr (HasMixedRankSides<Topology>::value) {
      return Topology::has_mixed_rank_sides;
    }
    return false;
  }
};

struct has_homogeneous_faces_impl {
  using result_type = bool;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::has_homogeneous_faces; }
};

struct is_shell_impl {
  using result_type = bool;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::is_shell; }
};

struct side_rank_impl {
  using result_type = stk::topology::rank_t;

  KOKKOS_INLINE_FUNCTION
  side_rank_impl(unsigned ordinal)
    : m_ordinal(ordinal)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { 
    if (!has_mixed_rank_sides_impl{}(Topology{})) {
      return Topology::side_rank;
    } else {
      return Topology::side_topology_rank(m_ordinal);
    }
  }

  unsigned m_ordinal;
};

struct num_side_ranks_impl {
  using result_type = unsigned;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_side_ranks(); }
};

struct dimension_impl {
  using result_type = unsigned;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::dimension; }
};

struct num_vertices_impl {
  using result_type = unsigned;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_vertices; }
};

struct num_edges_impl {
  using result_type = unsigned;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_edges; }
};

struct num_faces_impl {
  using result_type = unsigned;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_faces; }
};

struct num_permutations_impl {
  using result_type = unsigned;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_permutations; }
};

struct num_positive_permutations_impl {
  using result_type = unsigned;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_positive_permutations; }
};

struct base_impl {
  using result_type = stk::topology;

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::base; }
};

struct edge_topology_impl {
  using result_type = stk::topology;

  KOKKOS_INLINE_FUNCTION
  edge_topology_impl(unsigned ordinal = 0)
    : m_ordinal(ordinal)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::edge_topology(m_ordinal); }

  unsigned m_ordinal;
};

struct defined_on_spatial_dimension_impl {
  using result_type = bool;

  KOKKOS_INLINE_FUNCTION
  defined_on_spatial_dimension_impl(unsigned ordinal)
    : m_ordinal(ordinal)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::defined_on_spatial_dimension(m_ordinal); }

  unsigned m_ordinal;
};

struct face_topology_impl {
  using result_type = stk::topology;

  KOKKOS_INLINE_FUNCTION
  face_topology_impl(unsigned ordinal)
    : m_ordinal(ordinal)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::face_topology(m_ordinal); }

  unsigned m_ordinal;
};

template <typename SideRanksOutputIterator>
struct side_ranks_impl {
  using result_type = void;
  KOKKOS_FUNCTION
  side_ranks_impl(SideRanksOutputIterator output_side_ranks)
    : m_output_side_ranks(output_side_ranks)
  {}
  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  void operator()(Topology) const { Topology::side_ranks(m_output_side_ranks); }

  SideRanksOutputIterator m_output_side_ranks;
};

template <typename OrdinalOutputIterator>
struct edge_node_ordinals_impl {
  using result_type = void;
  KOKKOS_FUNCTION
  edge_node_ordinals_impl(unsigned ordinal, OrdinalOutputIterator output_ordinals)
    : m_ordinal(ordinal)
    , m_output_ordinals(output_ordinals)
  {}
  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  void operator()(Topology) const { Topology::edge_node_ordinals(m_ordinal,m_output_ordinals); }
  unsigned m_ordinal;
  OrdinalOutputIterator m_output_ordinals;
};

template <typename OrdinalOutputIterator>
struct face_node_ordinals_impl {
  using result_type = void;
  KOKKOS_FUNCTION
  face_node_ordinals_impl(unsigned ordinal, OrdinalOutputIterator output_ordinals)
    : m_ordinal(ordinal)
    , m_output_ordinals(output_ordinals)
  {}
  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  void operator()(Topology) const { Topology::face_node_ordinals(m_ordinal,m_output_ordinals); }
  unsigned m_ordinal;
  OrdinalOutputIterator m_output_ordinals;
};

template <typename OrdinalOutputIterator>
struct permutation_node_ordinals_impl {
  using result_type = void;
  KOKKOS_FUNCTION
  permutation_node_ordinals_impl(unsigned ordinal, OrdinalOutputIterator output_ordinals)
    : m_ordinal(ordinal)
    , m_output_ordinals(output_ordinals)
  {}
  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  void operator()(Topology) const { Topology::permutation_node_ordinals(m_ordinal,m_output_ordinals); }
  unsigned m_ordinal;
  OrdinalOutputIterator m_output_ordinals;
};

template <typename NodeArray, typename NodeOutputIterator>
struct edge_nodes_impl {
  using result_type = void;
  KOKKOS_FUNCTION
  edge_nodes_impl(const NodeArray &nodes,
                  unsigned ordinal,
                  NodeOutputIterator output_ordinals)
    : m_nodes(nodes),
      m_ordinal(ordinal),
      m_output_ordinals(output_ordinals)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  void operator()(Topology) const { Topology::edge_nodes(m_nodes,m_ordinal,m_output_ordinals); }

  const NodeArray & m_nodes;
  unsigned m_ordinal;
  NodeOutputIterator m_output_ordinals;
};

template <typename NodeArray, typename NodeOutputIterator>
struct face_nodes_impl {
  using result_type = void;
  KOKKOS_FUNCTION
  face_nodes_impl(const NodeArray &nodes,
                  unsigned ordinal,
                  NodeOutputIterator output_ordinals)
    : m_nodes(nodes),
      m_ordinal(ordinal),
      m_output_ordinals(output_ordinals)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  void operator()(Topology) const { Topology::face_nodes(m_nodes,m_ordinal,m_output_ordinals); }

  const NodeArray & m_nodes;
  unsigned m_ordinal;
  NodeOutputIterator m_output_ordinals;
};

template <typename NodeArray, typename NodeOutputIterator>
struct permutation_nodes_impl {
  using result_type = void;
  KOKKOS_FUNCTION
  permutation_nodes_impl(const NodeArray &nodes,
                  unsigned ordinal,
                  NodeOutputIterator output_ordinals)
    : m_nodes(nodes),
      m_ordinal(ordinal),
      m_output_ordinals(output_ordinals)
  {}

  template <typename Topology>
  KOKKOS_INLINE_FUNCTION
  void operator()(Topology) const { Topology::permutation_nodes(m_nodes,m_ordinal,m_output_ordinals); }

  const NodeArray & m_nodes;
  unsigned m_ordinal;
  NodeOutputIterator m_output_ordinals;
};

} /*namespace stk::topology_detail*/
} /*namespace stk*/

#endif //STKTOPOLOGY_TOPOLOGY_TCC

