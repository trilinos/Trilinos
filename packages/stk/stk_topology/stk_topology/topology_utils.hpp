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

#ifndef STKTOPOLOGY_TOPOLOGY_TCC
#define STKTOPOLOGY_TOPOLOGY_TCC

// IWYU pragma: private, include "stk_topology/topology.hpp"

namespace stk { namespace topology_detail {

struct num_nodes_impl {
  using result_type = unsigned;
  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const
  { return Topology::num_nodes; }
};

struct rank_impl {
  using result_type = topology::rank_t;
  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const
  { return Topology::rank; }
};

template <typename NodeArrayA, typename NodeArrayB>
struct equivalent_impl {
  using result_type = std::pair<bool,unsigned>;

  equivalent_impl( const NodeArrayA &a , const NodeArrayB &b )
    : m_a(a), m_b(b)
  {}

  template <typename Topology>
  result_type operator()(Topology) const
  { return Topology::equivalent(m_a, m_b); }

  const NodeArrayA & m_a;
  const NodeArrayB & m_b;
};

template <typename NodeArray>
struct lexicographical_smallest_permutation_impl {
  using result_type = unsigned;

  lexicographical_smallest_permutation_impl( const NodeArray &nodes , bool only_positive_permutations )
    : m_nodes(nodes), m_only_positive_permutations(only_positive_permutations)
  {}

  template <typename Topology>
  result_type operator()(Topology) const
  { return Topology::lexicographical_smallest_permutation(m_nodes, m_only_positive_permutations); }

  const NodeArray & m_nodes;
  bool              m_only_positive_permutations;
};

template <typename NodeArray>
struct lexicographical_smallest_permutation_preserve_polarity_impl {
  using result_type = unsigned;

  lexicographical_smallest_permutation_preserve_polarity_impl( const NodeArray &nodes, const NodeArray &element_nodes)
    : m_nodes(nodes), m_element_nodes(element_nodes)
  {}

  template <typename Topology>
  result_type operator()(Topology) const
  { return Topology::lexicographical_smallest_permutation_preserve_polarity(m_nodes, m_element_nodes); }

  const NodeArray & m_nodes;
  const NodeArray & m_element_nodes;
};

struct has_homogeneous_faces_impl {
  using result_type = bool;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::has_homogeneous_faces; }
};

struct is_shell_impl {
  using result_type = bool;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::is_shell; }
};

struct side_rank_impl {
  using result_type = stk::topology::rank_t;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::side_rank; }
};

struct dimension_impl {
  using result_type = unsigned;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::dimension; }
};

struct num_vertices_impl {
  using result_type = unsigned;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_vertices; }
};

struct num_edges_impl {
  using result_type = unsigned;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_edges; }
};

struct num_faces_impl {
  using result_type = unsigned;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_faces; }
};

struct num_permutations_impl {
  using result_type = unsigned;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_permutations; }
};

struct num_positive_permutations_impl {
  using result_type = unsigned;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::num_positive_permutations; }
};

struct base_impl {
  using result_type = stk::topology;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::base; }
};

struct edge_topology_impl {
  using result_type = stk::topology;

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::edge_topology; }
};

struct defined_on_spatial_dimension_impl {
  using result_type = bool;

  STK_INLINE_FUNCTION
  defined_on_spatial_dimension_impl(unsigned ordinal)
    : m_ordinal(ordinal)
  {}

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::defined_on_spatial_dimension(m_ordinal); }

  unsigned m_ordinal;
};

struct face_topology_impl {
  using result_type = stk::topology;

  STK_INLINE_FUNCTION
  face_topology_impl(unsigned ordinal)
    : m_ordinal(ordinal)
  {}

  template <typename Topology>
  STK_INLINE_FUNCTION
  result_type operator()(Topology) const { return Topology::face_topology(m_ordinal); }

  unsigned m_ordinal;
};

template <typename OrdinalOutputIterator>
struct edge_node_ordinals_impl {
  using result_type = void;
  edge_node_ordinals_impl(unsigned ordinal, OrdinalOutputIterator output_ordinals)
    : m_ordinal(ordinal)
    , m_output_ordinals(output_ordinals)
  {}
  template <typename Topology>
  void operator()(Topology) const { Topology::edge_node_ordinals(m_ordinal,m_output_ordinals); }
  unsigned m_ordinal;
  OrdinalOutputIterator m_output_ordinals;
};

template <typename OrdinalOutputIterator>
struct face_node_ordinals_impl {
  using result_type = void;
  face_node_ordinals_impl(unsigned ordinal, OrdinalOutputIterator output_ordinals)
    : m_ordinal(ordinal)
    , m_output_ordinals(output_ordinals)
  {}
  template <typename Topology>
  void operator()(Topology) const
  { Topology::face_node_ordinals(m_ordinal,m_output_ordinals); }
  unsigned m_ordinal;
  OrdinalOutputIterator m_output_ordinals;
};

template <typename OrdinalOutputIterator>
struct permutation_node_ordinals_impl {
  using result_type = void;
  permutation_node_ordinals_impl(unsigned ordinal, OrdinalOutputIterator output_ordinals)
    : m_ordinal(ordinal)
    , m_output_ordinals(output_ordinals)
  {}
  template <typename Topology>
  void operator()(Topology) const
  { Topology::permutation_node_ordinals(m_ordinal,m_output_ordinals); }
  unsigned m_ordinal;
  OrdinalOutputIterator m_output_ordinals;
};

template <typename NodeArray, typename NodeOutputIterator>
struct edge_nodes_impl {
  using result_type = void;
  edge_nodes_impl(const NodeArray &nodes,
                  unsigned ordinal,
                  NodeOutputIterator output_ordinals)
    : m_nodes(nodes),
      m_ordinal(ordinal),
      m_output_ordinals(output_ordinals)
  {}

  template <typename Topology>
  void operator()(Topology) const { Topology::edge_nodes(m_nodes,m_ordinal,m_output_ordinals); }

  const NodeArray & m_nodes;
  unsigned m_ordinal;
  NodeOutputIterator m_output_ordinals;
};

template <typename NodeArray, typename NodeOutputIterator>
struct face_nodes_impl {
  using result_type = void;
  face_nodes_impl(const NodeArray &nodes,
                  unsigned ordinal,
                  NodeOutputIterator output_ordinals)
    : m_nodes(nodes),
      m_ordinal(ordinal),
      m_output_ordinals(output_ordinals)
  {}

  template <typename Topology>
  void operator()(Topology) const { Topology::face_nodes(m_nodes,m_ordinal,m_output_ordinals); }

  const NodeArray & m_nodes;
  unsigned m_ordinal;
  NodeOutputIterator m_output_ordinals;
};

template <typename NodeArray, typename NodeOutputIterator>
struct permutation_nodes_impl {
  using result_type = void;
  permutation_nodes_impl(const NodeArray &nodes,
                  unsigned ordinal,
                  NodeOutputIterator output_ordinals)
    : m_nodes(nodes),
      m_ordinal(ordinal),
      m_output_ordinals(output_ordinals)
  {}

  template <typename Topology>
  void operator()(Topology) const { Topology::permutation_nodes(m_nodes,m_ordinal,m_output_ordinals); }

  const NodeArray & m_nodes;
  unsigned m_ordinal;
  NodeOutputIterator m_output_ordinals;
};

} /*namespace stk::topology_detail*/


STK_INLINE_FUNCTION
unsigned topology::num_nodes() const
{
  using functor = topology_detail::num_nodes_impl;
  topology::apply_functor< functor > apply;
  if (m_value < END_TOPOLOGY)
      return apply(m_value);
  else if (is_superedge())
      return m_value - SUPEREDGE_START;
  else if (is_superface())
      return m_value - SUPERFACE_START;
  else
      return m_value - SUPERELEMENT_START;
}

inline
topology::rank_t topology::rank() const
{
  using functor = topology_detail::rank_impl;
  topology::apply_functor< functor > apply;
  if (m_value < END_TOPOLOGY)
      return apply(m_value);
  else if (is_superedge())
      return topology::EDGE_RANK;
  else if (is_superface())
      return topology::FACE_RANK;
  else if (is_superelement())
      return topology::ELEMENT_RANK;
  return topology::INVALID_RANK;
}

template <typename NodeArrayA, typename NodeArrayB>
inline
std::pair<bool,unsigned> topology::equivalent( const NodeArrayA &a, const NodeArrayB &b) const {
  using functor = topology_detail::equivalent_impl<NodeArrayA,NodeArrayB>;
  functor f(a,b);
  topology::apply_host_functor< functor > apply( f );
  return apply(m_value);
}

template <typename NodeArray>
inline
unsigned topology::lexicographical_smallest_permutation( const NodeArray &nodes, bool only_positive_permutations) const {
  using functor = topology_detail::lexicographical_smallest_permutation_impl< NodeArray >;
  functor f(nodes, only_positive_permutations);
  topology::apply_host_functor< functor > apply( f );
  return apply(m_value);
}

template <typename NodeArray>
inline
unsigned topology::lexicographical_smallest_permutation_preserve_polarity( const NodeArray &nodes, const NodeArray &element_nodes) const {
  using functor = topology_detail::lexicographical_smallest_permutation_preserve_polarity_impl< NodeArray >;
  functor f(nodes, element_nodes);
  topology::apply_host_functor< functor > apply( f );
  return apply(m_value);
}

STK_INLINE_FUNCTION
bool topology::has_homogeneous_faces() const {
  using functor = topology_detail::has_homogeneous_faces_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
bool topology::is_shell() const {
  using functor = topology_detail::is_shell_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
stk::topology::rank_t topology::side_rank() const {
  using functor = topology_detail::side_rank_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
unsigned topology::dimension() const {
  using functor = topology_detail::dimension_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
unsigned topology::num_vertices() const {
  using functor = topology_detail::num_vertices_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
unsigned topology::num_edges() const {
  using functor = topology_detail::num_edges_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
unsigned topology::num_faces() const {
  using functor = topology_detail::num_faces_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
unsigned topology::num_permutations() const {
  using functor = topology_detail::num_permutations_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
unsigned topology::num_positive_permutations() const {
  using functor = topology_detail::num_positive_permutations_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
stk::topology topology::base() const {
  using functor = topology_detail::base_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
stk::topology topology::edge_topology() const {
  using functor = topology_detail::edge_topology_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

STK_INLINE_FUNCTION
bool topology::defined_on_spatial_dimension(unsigned ordinal) const {
  using functor = topology_detail::defined_on_spatial_dimension_impl;
  functor f(ordinal);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

STK_INLINE_FUNCTION
stk::topology topology::face_topology(unsigned ordinal) const {
  using functor = topology_detail::face_topology_impl;
  functor f(ordinal);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename OrdinalOutputIterator>
void topology::edge_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  using functor = topology_detail::edge_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_host_functor< functor > apply( f );
  apply(m_value);
}

template <typename OrdinalOutputIterator>
void topology::face_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  using functor = topology_detail::face_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_host_functor< functor > apply( f );
  apply(m_value);
}

template <typename OrdinalOutputIterator>
void topology::permutation_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  using functor = topology_detail::permutation_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_host_functor< functor > apply( f );
  apply(m_value);
}

template <typename NodeArray, typename NodeOutputIterator>
void topology::edge_nodes(const NodeArray & nodes,
                          unsigned ordinal,
                          NodeOutputIterator output_ordinals) const
{
  using functor = topology_detail::edge_nodes_impl<NodeArray,NodeOutputIterator>;
  functor f(nodes,ordinal,output_ordinals);
  topology::apply_host_functor< functor > apply( f );
  apply(m_value);
}

template <typename NodeArray, typename NodeOutputIterator>
void topology::face_nodes(const NodeArray & nodes,
                          unsigned ordinal,
                          NodeOutputIterator output_ordinals) const
{
  using functor = topology_detail::face_nodes_impl<NodeArray,NodeOutputIterator>;
  functor f(nodes,ordinal,output_ordinals);
  topology::apply_host_functor< functor > apply( f );
  apply(m_value);
}

template <typename NodeArray, typename NodeOutputIterator>
void topology::permutation_nodes(const NodeArray & nodes,
                          unsigned ordinal,
                          NodeOutputIterator output_ordinals) const
{
  using functor = topology_detail::permutation_nodes_impl<NodeArray,NodeOutputIterator>;
  functor f(nodes,ordinal,output_ordinals);
  topology::apply_host_functor< functor > apply( f );
  apply(m_value);
}

} /*namespace stk*/

#endif //STKTOPOLOGY_TOPOLOGY_TCC

