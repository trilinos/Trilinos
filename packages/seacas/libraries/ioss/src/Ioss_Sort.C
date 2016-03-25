// Copyright(C) 2016
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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

#include <Ioss_Sort.h>
#include <cstdlib>
#include <vector>

// The following 'qsort' routine is modified from Sedgewicks algorithm
// It selects the pivot based on the median of the left, right, and
// center values to try to avoid degenerate cases ocurring when a
// single value is chosen.  It performs a quicksort on intervals down
// to the QSORT_CUTOFF size and then performs a final insertion sort
// on the almost sorted final array.  Based on data in Sedgewick, the
// QSORT_CUTOFF value should be between 5 and 20.
//
// See Sedgewick for further details

namespace {
  const int QSORT_CUTOFF=12;

  template <typename INT>
  inline void SWAP(INT *V, size_t I, size_t J)
  {
    std::swap(V[I], V[J]);
  }

  template <typename INT>
  size_t median3(INT v[], size_t left, size_t right)
  {
    size_t center;
    center = (left + right) / 2;

    if (v[left] > v[center])
      SWAP(v, left, center);
    if (v[left] > v[right])
      SWAP(v, left, right);
    if (v[center] > v[right])
      SWAP(v, center, right);

    SWAP(v, center, right-1);
    return right-1;
  }

  template <typename INT>
  void qsort_int(INT v[], size_t left, size_t right)
  {
    size_t pivot;
    size_t i, j;

    if (left + QSORT_CUTOFF <= right) {
      pivot = median3(v, left, right);
      i = left;
      j = right - 1;

      for ( ; ; ) {
	while (v[++i] < v[pivot]);
	while (v[--j] > v[pivot]);
	if (i < j) {
	  SWAP(v, i, j);
	} else {
	  break;
	}
      }

      SWAP(v, i, right-1);
      qsort_int(v, left, i-1);
      qsort_int(v, i+1, right);
    }
  }

  template <typename INT>
  void isort_int(INT v[], size_t N)
  {
    size_t i,j;
    size_t ndx = 0;
    INT small;
    INT tmp;

    if (N <= 1) return;
    small = v[0];
    for (i = 1; i < N; i++) {
      if (v[i] < small) {
	small = v[i];
	ndx = i;
      }
    }
    /* Put smallest value in slot 0 */
    SWAP(v, 0, ndx);

    for (i=1; i <N; i++) {
      tmp = v[i];
      for (j=i; tmp < v[j-1]; j--) {
	v[j] = v[j-1];
      }
      v[j] = tmp;
    }
  }
}

namespace Ioss {
  // Explicit Template Instantiation
  // -- "Know" these are only types that Ioss will use with this.
  //    Avoids having implementation in header file.
  template void qsort<>(std::vector<int> &v);
  template void qsort<>(std::vector<int64_t> &v);
  template void qsort<>(std::vector<std::pair<int,int>> &v);
  template void qsort<>(std::vector<std::pair<int,long>> &v);
  template void qsort<>(std::vector<std::pair<long,int>> &v);
  template void qsort<>(std::vector<std::pair<int64_t,int64_t>> &v);

  template <typename INT>
  void qsort(std::vector<INT> &v)
  {
    if (v.size() <= 1) return;
    qsort_int(v.data(), 0, v.size()-1);
    isort_int(v.data(), v.size());
  }
}
