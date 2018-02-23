// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER


#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include "typedefs.hpp"

namespace BlockCrsTest {
  
  template<typename T1, typename T2, typename CompareType>
  KOKKOS_INLINE_FUNCTION
  static T1* lower_bound(T1* first, T1* last, const T2& val,
                         CompareType compare) {
    T1 *it;
    local_ordinal_type step = 0, count = last - first;
    while (count > 0) {
      it = first; step = count/2; it += step;
      if (compare(*it,val)) {
        first = ++it;
        count -= step + 1;
      } else {
        count = step;
      }
    }
    return first;
  }
  
  template<typename T1, typename T2>
  KOKKOS_FORCEINLINE_FUNCTION
  static T1* lower_bound(T1* first, T1* last, const T2& val) {
    return lower_bound(first, last, val, [](T1 left, T2 right) { return left < right; });
  }

  template<typename T1, typename T2>
  KOKKOS_INLINE_FUNCTION
  static void heapify(T1 *v, T2 n, T2 i) {
    T2 largest = i;
    T2 l = 2*i + 1;
    T2 r = 2*i + 2;

    if (l < n && v[l] > v[largest]) largest = l;
    if (r < n && v[r] > v[largest]) largest = r;
    if (largest != i) {
      // swap
      T1 tmp = v[i]; v[i] = v[largest]; v[largest] = tmp;
      heapify(v, n, largest);
    }
  }

  template<typename T1, typename T2>
  KOKKOS_INLINE_FUNCTION
  static void heap_sort(T1 *v, T2 n) { 
    for (T2 i=n/2-1;i>=0;--i) heapify(v, n, i);
    for (T2 i=n-1;i>=0;--i) {
      T1 tmp = v[0]; v[0] = v[i]; v[i] = tmp;
      heapify(v, i, 0);
    }
  }
  
}

#endif
  

