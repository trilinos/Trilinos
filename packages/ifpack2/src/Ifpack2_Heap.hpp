// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_HEAP_HPP
#define IFPACK2_HEAP_HPP

#include <algorithm>
#include "Teuchos_Array.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Ifpack2 {

template<typename Scalar, typename Ordinal>
struct greater_indirect {
  greater_indirect(const Teuchos::Array<Scalar>& vals)
  : m_vals(vals) {}
  ~greater_indirect(){}

  bool operator()(const Ordinal& lhs, const Ordinal& rhs) const
  { return Teuchos::ScalarTraits<Scalar>::magnitude(m_vals[lhs]) >
           Teuchos::ScalarTraits<Scalar>::magnitude(m_vals[rhs]); }

private:
  const Teuchos::Array<Scalar>& m_vals;
};//struct greater_indirect


/** Add idx to heap, don't assume heap occupies entire vector.
*/
template<typename Ordinal, typename SizeType>
void add_to_heap(const Ordinal& idx, Teuchos::Array<Ordinal>& heap, SizeType& heap_len)
{
  if (heap.size() == heap_len) heap.push_back(idx);
  else heap[heap_len] = idx;
  ++heap_len;
  std::push_heap(heap.begin(), heap.begin()+heap_len, std::greater<Ordinal>());
}

/** Add idx to heap, don't assume heap occupies entire vector.
    Also take custom comparator.
*/
template<typename Ordinal, typename SizeType, class Compare>
void add_to_heap(const Ordinal& idx, Teuchos::Array<Ordinal>& heap, SizeType& heap_len, Compare comp)
{
  if (heap.size() == heap_len) heap.push_back(idx);
  else heap[heap_len] = idx;
  ++heap_len;
  std::push_heap(heap.begin(), heap.begin()+heap_len, comp);
}

/** Remove heap root, don't shorten vector but update a heap_len parameter. */
template<typename Ordinal, typename SizeType>
void rm_heap_root(Teuchos::Array<Ordinal>& heap, SizeType& heap_len)
{
  std::pop_heap(heap.begin(), heap.begin()+heap_len, std::greater<Ordinal>());
  --heap_len;
}

/** Remove heap root, with custom comparator, don't assume heap occupies
  entire vector.
*/
template<typename Ordinal, typename SizeType, class Compare>
void rm_heap_root(Teuchos::Array<Ordinal>& heap, SizeType& heap_len, Compare comp)
{
  std::pop_heap(heap.begin(), heap.begin()+heap_len, comp);
  --heap_len;
}

}//namespace Ifpack2

#endif

