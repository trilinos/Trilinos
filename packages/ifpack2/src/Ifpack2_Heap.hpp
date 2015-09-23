/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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

