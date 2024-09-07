#ifndef TOPOLOGY_DEFN_HPP
#define TOPOLOGY_DEFN_HPP

// IWYU pragma: private, include "stk_topology/topology.hpp"
#include "stk_topology/topology_decl.hpp"
#include "stk_topology/topology_utils.hpp"
#include "stk_topology/apply_functor.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "stk_util/stk_config.h"

namespace stk {

namespace impl {

// Temporary function used to identify the new SHELL_[TRI|QUAD]_ALL_FACE_SIDES
// Will be removed once a proper conversion is available
STK_INLINE_FUNCTION
bool is_temporary_shell_with_all_face_sides(topology::topology_t m_value) {
  return (m_value == topology::topology_t::SHELL_QUAD_4_ALL_FACE_SIDES ||
          m_value == topology::topology_t::SHELL_QUAD_8_ALL_FACE_SIDES ||
          m_value == topology::topology_t::SHELL_QUAD_9_ALL_FACE_SIDES ||
          m_value == topology::topology_t::SHELL_TRI_3_ALL_FACE_SIDES ||
          m_value == topology::topology_t::SHELL_TRI_4_ALL_FACE_SIDES ||
          m_value == topology::topology_t::SHELL_TRI_6_ALL_FACE_SIDES);
}
}

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

STK_INLINE_FUNCTION
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
STK_INLINE_FUNCTION
EquivalentPermutation topology::is_equivalent( const NodeArrayA &a, const NodeArrayB &b) const {
  using functor = topology_detail::is_equivalent_impl<NodeArrayA, NodeArrayB>;
  functor f(a,b);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename NodeArray>
STK_INLINE_FUNCTION
unsigned topology::lexicographical_smallest_permutation( const NodeArray &nodes, bool only_positive_permutations) const {
  using functor = topology_detail::lexicographical_smallest_permutation_impl< NodeArray >;
  functor f(nodes, only_positive_permutations);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename NodeArray>
STK_INLINE_FUNCTION
unsigned topology::lexicographical_smallest_permutation_preserve_polarity( const NodeArray &nodes, const NodeArray &element_nodes) const {
  using functor = topology_detail::lexicographical_smallest_permutation_preserve_polarity_impl< NodeArray >;
  functor f(nodes, element_nodes);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename OrdinalOutputIterator>
STK_INLINE_FUNCTION
void topology::sub_topology_node_ordinals(unsigned sub_rank, unsigned sub_ordinal, OrdinalOutputIterator output_ordinals) const
{
  switch(sub_rank)
  {
  case NODE_RANK: *output_ordinals = sub_ordinal;                    break;
  case EDGE_RANK: edge_node_ordinals(sub_ordinal, output_ordinals);  break;
  case FACE_RANK:
    if (is_shell_with_face_sides() && sub_ordinal >= num_faces()) {
      edge_node_ordinals(sub_ordinal - num_faces(), output_ordinals);
    } else {
      face_node_ordinals(sub_ordinal, output_ordinals);
    }
    break;
  default: break;
  }
}

template <typename NodeArray, typename NodeOutputIterator>
STK_INLINE_FUNCTION
void topology::sub_topology_nodes(const NodeArray & nodes, unsigned sub_rank, unsigned sub_ordinal, NodeOutputIterator output_nodes) const
{
  switch(sub_rank)
  {
  case NODE_RANK: *output_nodes = nodes[sub_ordinal];            break;
  case EDGE_RANK: edge_nodes(nodes, sub_ordinal, output_nodes);  break;
  case FACE_RANK:
    if (is_shell_side_ordinal(sub_ordinal)) {
      edge_nodes(nodes, sub_ordinal - num_faces(), output_nodes);
    } else {
      face_nodes(nodes, sub_ordinal, output_nodes);
    }
    break;
  default: break;
  }
}

STK_INLINE_FUNCTION
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

STK_INLINE_FUNCTION
topology topology::sub_topology(unsigned sub_rank, unsigned sub_ordinal) const
{
  switch(sub_rank)
  {
  case NODE_RANK: return NODE;
  case EDGE_RANK: return edge_topology(sub_ordinal);
  case FACE_RANK:
    if (is_shell_side_ordinal(sub_ordinal)) {
      return edge_topology(sub_ordinal - num_faces());
    }
    return face_topology(sub_ordinal);
  default: break;
  }
  return INVALID_TOPOLOGY;
}

template <typename OrdinalOutputIterator>
STK_INLINE_FUNCTION
void topology::side_node_ordinals(unsigned side_ordinal, OrdinalOutputIterator output_ordinals) const
{
  if (is_shell_side_ordinal(side_ordinal)) {
    sub_topology_node_ordinals(EDGE_RANK, side_ordinal-num_faces(), output_ordinals);
  } else {
    sub_topology_node_ordinals( side_rank(), side_ordinal, output_ordinals);
  }
}

template <typename NodeArray, typename NodeOutputIterator>
STK_INLINE_FUNCTION
void topology::side_nodes(const NodeArray & nodes, unsigned side_ordinal, NodeOutputIterator output_nodes) const
{
  if (is_shell_side_ordinal(side_ordinal)) {
    sub_topology_nodes( nodes, EDGE_RANK, side_ordinal-num_faces(), output_nodes);
  } else {
    sub_topology_nodes( nodes, side_rank(), side_ordinal, output_nodes);
  }
}

STK_INLINE_FUNCTION
unsigned topology::num_sides() const
{
  unsigned num_sides_out = 0u;
  if (side_rank() != INVALID_RANK) {
    num_sides_out = side_rank() > NODE_RANK? num_sub_topology(side_rank()) : num_vertices();

    if (is_shell_with_face_sides() &&
        !impl::is_temporary_shell_with_all_face_sides(m_value))
      num_sides_out += num_sub_topology(EDGE_RANK);
  }
  return num_sides_out;
}

STK_INLINE_FUNCTION
topology topology::side_topology(unsigned side_ordinal) const
{
  if (is_shell_side_ordinal(side_ordinal) && !impl::is_temporary_shell_with_all_face_sides(m_value))
    return shell_side_topology(side_ordinal-num_faces());

  return sub_topology(side_rank(), side_ordinal);
}

STK_INLINE_FUNCTION
bool topology::is_superelement() const
{
  return m_value > SUPERELEMENT_START;
}

STK_INLINE_FUNCTION
bool topology::is_superface() const
{
  return m_value > SUPERFACE_START && m_value < SUPERFACE_END;
}

STK_INLINE_FUNCTION
bool topology::is_superedge() const
{
  return m_value > SUPEREDGE_START && m_value < SUPEREDGE_END;
}

STK_INLINE_FUNCTION
bool topology::is_super_topology() const
{
  return is_superelement() || is_superface() || is_superedge();
}

STK_INLINE_FUNCTION
bool topology::is_shell_side_ordinal(unsigned ord) const
{
  return is_shell_with_face_sides() && ord >= num_faces();
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
bool topology::is_shell_with_face_sides() const {
  using functor = topology_detail::is_shell_with_face_sides_impl;
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
bool topology::is_positive_polarity(unsigned permutation_ordinal) const {
  return (permutation_ordinal < num_positive_permutations());
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
stk::topology topology::edge_topology(unsigned ordinal) const {
  using functor = topology_detail::edge_topology_impl;
  functor f(ordinal);
  topology::apply_functor< functor > apply( f );
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

STK_INLINE_FUNCTION
topology topology::shell_side_topology(unsigned ordinal) const {
  using functor = topology_detail::shell_side_topology_impl;
  functor f(ordinal);
  topology::apply_functor< functor > apply( f );
  return apply(m_value);
}

template <typename OrdinalOutputIterator>
STK_INLINE_FUNCTION
void topology::edge_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  using functor = topology_detail::edge_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

template <typename OrdinalOutputIterator>
STK_INLINE_FUNCTION
void topology::face_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  using functor = topology_detail::face_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

template <typename OrdinalOutputIterator>
STK_INLINE_FUNCTION
void topology::permutation_node_ordinals( unsigned ordinal, OrdinalOutputIterator output_ordinals) const
{
  STK_NGP_ThrowAssert(m_value != QUAD_6 && m_value != WEDGE_12);
  using functor = topology_detail::permutation_node_ordinals_impl<OrdinalOutputIterator>;
  functor f(ordinal, output_ordinals);
  topology::apply_functor< functor > apply( f );
  apply(m_value);
}

template <typename NodeArray, typename NodeOutputIterator>
STK_INLINE_FUNCTION
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
STK_INLINE_FUNCTION
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
STK_INLINE_FUNCTION
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
