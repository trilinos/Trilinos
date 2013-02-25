#ifndef STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP
#define STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP

#include <algorithm>

#include <boost/utility.hpp>

namespace stk { namespace topology_detail {

template <typename Topology, typename NodeArrayA, typename NodeArrayB, typename Node>
BOOST_GPU_ENABLED inline
typename boost::enable_if_c< (Topology::num_permutations > 0u), std::pair<bool,unsigned> >::type
equivalent_helper(Topology, const NodeArrayA &a, const NodeArrayB &b, Node)
{
  Node permutation[Topology::num_nodes];

  for (unsigned i=0; i<Topology::num_permutations; ++i) {
    Topology::permutation_nodes(a,i,permutation);

    if ( std::equal(permutation, permutation + Topology::num_nodes, &b[0]) )
      return std::make_pair(true,i);
  }
  return std::make_pair(false, 0);
}

template <typename Topology, typename NodeArrayA, typename NodeArrayB, typename Node>
BOOST_GPU_ENABLED inline
typename boost::enable_if_c< (Topology::num_permutations == 0u), std::pair<bool,unsigned> >::type
equivalent_helper(Topology, const NodeArrayA &a, const NodeArrayB &b, Node)
{
  if ( std::equal(&a[0], &a[0] + Topology::num_nodes, &b[0]) )
      return std::make_pair(true,0);
  return std::make_pair(false, 0);
}

template <typename Topology, typename NodeArray, typename Node>
BOOST_GPU_ENABLED inline
typename boost::enable_if_c< (Topology::num_permutations > 0u), unsigned >::type
lexicographical_smallest_permutation_helper(Topology, const NodeArray &nodes, bool only_positive_permutations, Node)
{
  Node permutation[Topology::num_nodes];

  const Node * nbegin = &nodes[0];
  const Node * nend = nbegin + Topology::num_nodes;

  unsigned min_permutation_index = 0;
  Node min_permutation[Topology::num_nodes];

  std::copy(nbegin,nend,min_permutation);

  if (!only_positive_permutations) {
    for (unsigned i=1; i<Topology::num_permutations; ++i) {
      Topology::permutation_nodes(nodes,i,permutation);

      if ( std::lexicographical_compare( permutation,     permutation     + Topology::num_nodes,
                                         min_permutation, min_permutation + Topology::num_nodes ) )
      {
        std::copy(permutation, permutation + Topology::num_nodes, min_permutation);
        min_permutation_index = i;
      }
    }
  }
  else {
    for (unsigned i=1; i<Topology::num_positive_permutations; ++i) {
      Topology::permutation_nodes(nodes,i,permutation);

      if ( std::lexicographical_compare( permutation,     permutation     + Topology::num_nodes,
                                         min_permutation, min_permutation + Topology::num_nodes ) )
      {
        std::copy(permutation, permutation + Topology::num_nodes, min_permutation);
        min_permutation_index = i;
      }
    }
  }

  return min_permutation_index;

}

template <typename Topology, typename NodeArray, typename Node>
BOOST_GPU_ENABLED inline
typename boost::enable_if_c< (Topology::num_permutations == 0u), unsigned >::type
lexicographical_smallest_permutation_helper(Topology, const NodeArray &, bool , Node)
{ return 0; }

}} // namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP

