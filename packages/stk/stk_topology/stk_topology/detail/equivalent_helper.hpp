#ifndef STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP
#define STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP

#include <algorithm>

namespace stk { namespace detail {

template <typename Topology, typename NodeArrayA, typename NodeArrayB, typename Node>
STKTOPOLOGY_INLINE_FUNCTION
std::pair<bool,int> equivalent_helper(Topology, const NodeArrayA &a, const NodeArrayB &b, Node)
{
  // size to num_nodes+1 to work around gcc4.4 error
  // when NodeArrayB is a char[]
  Node permutation[Topology::num_nodes+1];

  const Node * b_ptr= &b[0];
  const Node * b_end= b_ptr + Topology::num_nodes;

  for (int i=0; i<Topology::num_permutations; ++i) {
    Topology::permutation_nodes(a,i,permutation);

    if ( std::equal(b_ptr, b_end, permutation) )
      return std::make_pair(true,i);
  }

  return std::make_pair(false, -1);
}



}} // namespace stk::detail

#endif //STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP

