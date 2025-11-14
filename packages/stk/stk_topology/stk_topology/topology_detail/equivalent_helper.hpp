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

#ifndef STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP
#define STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP

#include "stk_topology/topology_decl.hpp"
#include "stk_util/stk_config.h"
#include <algorithm>
#include <type_traits>

namespace stk { namespace topology_detail {

template <class InputIt1, class InputIt2>
KOKKOS_INLINE_FUNCTION
bool stk_equal(InputIt1 first1, InputIt1 last1, InputIt2 first2)
{
  for (; first1 != last1; ++first1, ++first2) {
    if (!(*first1 == *first2)) {
      return false;
    }
  }
  return true;
}

template <typename Topology, typename NodeArrayA, typename NodeArrayB, typename Node>
KOKKOS_INLINE_FUNCTION
EquivalentPermutation is_equivalent_helper(Topology, const NodeArrayA &a, const NodeArrayB &b, Node)
{
  if constexpr (Topology::num_permutations == 0u)
  {
    if ( stk_equal(&a[0], &a[0] + Topology::num_nodes, &b[0]) ) {
      return EquivalentPermutation(true, 0);
    }
  }
  else
  {
    Node permutedNodes[Topology::num_nodes];
    for (unsigned i = 0; i < Topology::num_permutations; ++i) {
      Topology::permutation_nodes(&a[0], i, permutedNodes);
      if ( stk_equal(permutedNodes, permutedNodes + Topology::num_nodes, &b[0]) ) {
        return EquivalentPermutation(true, i);
      }
    }
  }
  return EquivalentPermutation(false, 0);
}

template<class InputIt, class OutputIt>
KOKKOS_INLINE_FUNCTION
OutputIt stk_copy(InputIt first, InputIt last, OutputIt d_first)
{
    while (first != last) {
        *d_first++ = *first++;
    }
    return d_first;
}

template<class InputIt1, class InputIt2>
KOKKOS_INLINE_FUNCTION
bool stk_lexicographical_compare(InputIt1 first1, InputIt1 last1,
                                 InputIt2 first2, InputIt2 last2)
{
    for ( ; (first1 != last1) && (first2 != last2); ++first1, (void) ++first2 ) {
        if (*first1 < *first2) return true;
        if (*first2 < *first1) return false;
    }
    return (first1 == last1) && (first2 != last2);
}

template <typename Topology, typename NodeArray, typename Node>
KOKKOS_INLINE_FUNCTION
unsigned lexicographical_smallest_permutation_helper(Topology, const NodeArray &nodes, bool only_positive_permutations, Node)
{
  if constexpr (Topology::num_permutations > 0u)
  {
    Node permutation[Topology::num_nodes];

    const Node * nbegin = &nodes[0];
    const Node * nend = nbegin + Topology::num_nodes;

    unsigned min_permutation_index = 0;
    Node min_permutation[Topology::num_nodes];

    stk_copy(nbegin,nend,min_permutation);

    if (!only_positive_permutations) {
      for (unsigned i=1; i<Topology::num_permutations; ++i) {
        Topology::permutation_nodes(nbegin,i,permutation);

        if ( stk_lexicographical_compare( permutation,     permutation     + Topology::num_nodes,
                                          min_permutation, min_permutation + Topology::num_nodes ) )
        {
          stk_copy(permutation, permutation + Topology::num_nodes, min_permutation);
          min_permutation_index = i;
        }
      }
    }
    else {
      for (unsigned i=1; i<Topology::num_positive_permutations; ++i) {
        Topology::permutation_nodes(nbegin,i,permutation);

        if ( stk_lexicographical_compare( permutation,     permutation     + Topology::num_nodes,
                                          min_permutation, min_permutation + Topology::num_nodes ) )
        {
          stk_copy(permutation, permutation + Topology::num_nodes, min_permutation);
          min_permutation_index = i;
        }
      }
    }
    return min_permutation_index;
  }
  return 0u;
}

template <typename Topology, typename NodeArray, typename Node>
KOKKOS_INLINE_FUNCTION
unsigned lexicographical_smallest_permutation_preserve_polarity_helper(Topology, const NodeArray &nodes, const NodeArray &element_nodes, Node)
{
  if constexpr (Topology::num_permutations > 0u)
  {
    Node permutation[Topology::num_nodes];

    unsigned my_permutation_index = Topology::num_permutations;

    //am i positive or negative
    for (unsigned i=0; i<Topology::num_permutations; ++i) {
        Topology::permutation_nodes(&nodes[0], i, permutation);

        if ( stk_equal( permutation, permutation + Topology::num_nodes, &element_nodes[0] ) )
        {
            my_permutation_index = i;
            break;
        }
    }
    unsigned min_permutation_index;
    if (my_permutation_index != Topology::num_permutations) {
        bool positive_polarity = my_permutation_index < Topology::num_positive_permutations;
        //if positive, only match 0-num_positive
        unsigned low = 0;
        unsigned high = Topology::num_positive_permutations;
        //if negative, only match num_positive - num_perm
        if (!positive_polarity) {
            low = Topology::num_positive_permutations;
            high = Topology::num_permutations;
        }
        min_permutation_index = low;
        Node min_permutation[Topology::num_nodes];
        Topology::permutation_nodes(&element_nodes[0], min_permutation_index, min_permutation);
        for (unsigned i=low + 1; i<high; ++i) {
            Topology::permutation_nodes(&element_nodes[0], i, permutation);

            if ( stk_lexicographical_compare( permutation,     permutation     + Topology::num_nodes,
                                              min_permutation, min_permutation + Topology::num_nodes ) )
            {
                stk_copy(permutation, permutation + Topology::num_nodes, min_permutation);
                min_permutation_index = i;
            }
        }
    }
    else {
        min_permutation_index = Topology::num_permutations;
    }
    return min_permutation_index;
  }
  return 0u;
}

}} // namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_EQUIVALENT_HELPER_HPP

