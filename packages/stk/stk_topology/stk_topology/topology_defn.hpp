#ifndef TOPOLOGY_DEFN_HPP
#define TOPOLOGY_DEFN_HPP

// IWYU pragma: private, include "stk_topology/topology.hpp"
#include "stk_topology/topology_decl.hpp"
#include "stk_topology/topology_utils.hpp"
#include "stk_topology/apply_functor.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/stk_config.h"

#include <climits>

namespace stk {

KOKKOS_INLINE_FUNCTION
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

KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
EquivalentPermutation topology::is_equivalent( const NodeArrayA &a, const NodeArrayB &b) const {
  using functor = topology_detail::is_equivalent_impl<NodeArrayA, NodeArrayB>;
  functor f(a,b);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename NodeArray>
KOKKOS_INLINE_FUNCTION
unsigned topology::lexicographical_smallest_permutation( const NodeArray &nodes, bool only_positive_permutations) const {
  using functor = topology_detail::lexicographical_smallest_permutation_impl< NodeArray >;
  functor f(nodes, only_positive_permutations);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename NodeArray>
KOKKOS_INLINE_FUNCTION
unsigned topology::lexicographical_smallest_permutation_preserve_polarity( const NodeArray &nodes, const NodeArray &element_nodes) const {
  using functor = topology_detail::lexicographical_smallest_permutation_preserve_polarity_impl< NodeArray >;
  functor f(nodes, element_nodes);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename OrdinalOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::sub_topology_node_ordinals(unsigned sub_rank, unsigned sub_ordinal, OrdinalOutputIterator output_ordinals) const
{
  switch(sub_rank)
  {
  case NODE_RANK: *output_ordinals = sub_ordinal;                    break;
  case EDGE_RANK: edge_node_ordinals(sub_ordinal, output_ordinals);  break;
  case FACE_RANK: face_node_ordinals(sub_ordinal, output_ordinals);  break;
  default: break;
  }
}

template <typename NodeArray, typename NodeOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::sub_topology_nodes(const NodeArray & nodes, unsigned sub_rank, unsigned sub_ordinal, NodeOutputIterator output_nodes) const
{
  switch(sub_rank)
  {
  case NODE_RANK: *output_nodes = nodes[sub_ordinal];            break;
  case EDGE_RANK: edge_nodes(nodes, sub_ordinal, output_nodes);  break;
  case FACE_RANK: face_nodes(nodes, sub_ordinal, output_nodes);  break;
  default: break;
  }
}

KOKKOS_INLINE_FUNCTION
unsigned topology::num_sub_topology(unsigned sub_rank) const
{
  switch(sub_rank)
  {
  case NODE_RANK: return num_nodes();
  case EDGE_RANK: return num_edges();
  case FACE_RANK: return num_faces();
  default: break;
  }
  return 0;
}

KOKKOS_INLINE_FUNCTION
topology topology::sub_topology(unsigned sub_rank, unsigned sub_ordinal) const
{
  switch(sub_rank)
  {
  case NODE_RANK: return NODE;
  case EDGE_RANK: return edge_topology(sub_ordinal);
  case FACE_RANK: return face_topology(sub_ordinal);
  default: break;
  }
  return INVALID_TOPOLOGY;
}

template <typename OrdinalOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::side_node_ordinals(unsigned side_ordinal, OrdinalOutputIterator output_ordinals) const
{
  auto fix_ordinal = has_mixed_rank_sides() && side_ordinal >= num_sub_topology(side_rank());
  auto adjusted_ordinal = (fix_ordinal) ? side_ordinal - num_sub_topology(side_rank()) : side_ordinal;

  sub_topology_node_ordinals(side_rank(side_ordinal), adjusted_ordinal, output_ordinals);
}

template <typename NodeArray, typename NodeOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::side_nodes(const NodeArray & nodes, unsigned side_ordinal, NodeOutputIterator output_nodes) const
{
  auto fix_ordinal = has_mixed_rank_sides() && side_ordinal >= num_sub_topology(side_rank());
  auto adjusted_ordinal = (fix_ordinal) ? side_ordinal - num_sub_topology(side_rank()) : side_ordinal;

  sub_topology_nodes(nodes, side_rank(side_ordinal), adjusted_ordinal, output_nodes);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::num_sides() const
{
  unsigned num_sides_out = 0u;
  if (side_rank() != INVALID_RANK) {
    num_sides_out = side_rank() > NODE_RANK ? num_sub_topology(side_rank()) : num_vertices();

    if (has_mixed_rank_sides() && side_rank() > EDGE_RANK) {
      num_sides_out += num_sub_topology(EDGE_RANK);
    }
  }
  return num_sides_out;
}

KOKKOS_INLINE_FUNCTION
topology topology::side_topology(unsigned side_ordinal) const
{
  auto fix_ordinal = has_mixed_rank_sides() && side_ordinal >= num_sub_topology(side_rank());
  auto adjusted_ordinal = (fix_ordinal) ? side_ordinal - num_sub_topology(side_rank()) : side_ordinal;

  return sub_topology(side_rank(side_ordinal), adjusted_ordinal);
}

KOKKOS_INLINE_FUNCTION
bool topology::is_superelement() const
{
  return m_value > SUPERELEMENT_START;
}

KOKKOS_INLINE_FUNCTION
bool topology::is_superface() const
{
  return m_value > SUPERFACE_START && m_value < SUPERFACE_END;
}

KOKKOS_INLINE_FUNCTION
bool topology::is_superedge() const
{
  return m_value > SUPEREDGE_START && m_value < SUPEREDGE_END;
}

KOKKOS_INLINE_FUNCTION
bool topology::is_super_topology() const
{
  return is_superelement() || is_superface() || is_superedge();
}

KOKKOS_INLINE_FUNCTION
bool topology::has_homogeneous_faces() const {
  using functor = topology_detail::has_homogeneous_faces_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
bool topology::is_shell() const {
  using functor = topology_detail::is_shell_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
bool topology::has_mixed_rank_sides() const {
  using functor = topology_detail::has_mixed_rank_sides_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
stk::topology::rank_t topology::side_rank(unsigned ord) const {
  using functor = topology_detail::side_rank_impl;
  functor f(ord);
  topology::apply_functor< functor > apply(f);
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::num_side_ranks() const {
  using functor = topology_detail::num_side_ranks_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::side_ordinal(unsigned ranked_side_ordinal, rank_t rank) const {
  auto invalid_ordinal = UINT_MAX;

  if(ranked_side_ordinal >= num_sub_topology(rank)) {
    return invalid_ordinal;
  }

  if (num_side_ranks() == 2) {
    if(rank != stk::topology::FACE_RANK && rank != stk::topology::EDGE_RANK) {
      return invalid_ordinal;
    }

    if(rank == stk::topology::FACE_RANK) {
      return ranked_side_ordinal;
    } else if(rank == stk::topology::EDGE_RANK) {
      return ranked_side_ordinal + num_faces();
    }
  } else if (num_side_ranks() == 1) {
    if(rank != side_rank()) {
      return invalid_ordinal;
    }

    return ranked_side_ordinal;
  }

  return invalid_ordinal;
}

KOKKOS_INLINE_FUNCTION
void topology::ranked_side_ordinal(unsigned side_ordinal, unsigned& ranked_side_ordinal, rank_t& rank) const {
  auto invalid_ordinal = UINT_MAX;

  if(side_ordinal >= num_sides()) {
    ranked_side_ordinal = invalid_ordinal;
    rank = stk::topology::INVALID_RANK;
    return;
  }

  auto offset = num_sub_topology(side_rank());

  rank = side_rank(side_ordinal);
  ranked_side_ordinal = side_ordinal;

  if(has_mixed_rank_sides() && (side_ordinal >= offset)) {
    ranked_side_ordinal -= offset;
  }
}

template <typename SideRanksOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::side_ranks( SideRanksOutputIterator output_ranks) const {
  using functor = topology_detail::side_ranks_impl<SideRanksOutputIterator>;
  functor f(output_ranks);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::dimension() const {
  using functor = topology_detail::dimension_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::num_vertices() const {
  using functor = topology_detail::num_vertices_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::num_edges() const {
  using functor = topology_detail::num_edges_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::num_faces() const {
  using functor = topology_detail::num_faces_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::num_permutations() const {
  using functor = topology_detail::num_permutations_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
unsigned topology::num_positive_permutations() const {
  using functor = topology_detail::num_positive_permutations_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
bool topology::is_positive_polarity(unsigned permutation_ordinal) const {
  return (permutation_ordinal < num_positive_permutations());
}

KOKKOS_INLINE_FUNCTION
stk::topology topology::base() const {
  using functor = topology_detail::base_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
stk::topology topology::edge_topology() const {
  using functor = topology_detail::edge_topology_impl;
  topology::apply_functor< functor > apply;
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
stk::topology topology::edge_topology(unsigned ordinal) const {
  using functor = topology_detail::edge_topology_impl;
  functor f(ordinal);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
bool topology::defined_on_spatial_dimension(unsigned ordinal) const {
  using functor = topology_detail::defined_on_spatial_dimension_impl;
  functor f(ordinal);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

KOKKOS_INLINE_FUNCTION
stk::topology topology::face_topology(unsigned ordinal) const {
  using functor = topology_detail::face_topology_impl;
  functor f(ordinal);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename OrdinalOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::edge_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  using functor = topology_detail::edge_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

template <typename OrdinalOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::face_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  using functor = topology_detail::face_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

template <typename OrdinalOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::permutation_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  STK_NGP_ThrowAssert(m_value != QUAD_6 && m_value != WEDGE_12);
  using functor = topology_detail::permutation_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

template <typename NodeArray, typename NodeOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::edge_nodes(const NodeArray & nodes,
                          unsigned ordinal,
                          NodeOutputIterator output_ordinals) const
{
  using functor = topology_detail::edge_nodes_impl<NodeArray,NodeOutputIterator>;
  functor f(nodes,ordinal,output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

template <typename NodeArray, typename NodeOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::face_nodes(const NodeArray & nodes,
                          unsigned ordinal,
                          NodeOutputIterator output_ordinals) const
{
  using functor = topology_detail::face_nodes_impl<NodeArray,NodeOutputIterator>;
  functor f(nodes,ordinal,output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

template <typename NodeArray, typename NodeOutputIterator>
KOKKOS_INLINE_FUNCTION
void topology::permutation_nodes(const NodeArray & nodes,
                          unsigned ordinal,
                          NodeOutputIterator output_ordinals) const
{
  using functor = topology_detail::permutation_nodes_impl<NodeArray,NodeOutputIterator>;
  functor f(nodes,ordinal,output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

}

#endif // TOPOLOGY_DEFN_HPP
