// Copyright(C) 2010 Sandia Corporation.
// 
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
/*!
 * The following 'indexed qsort' routine is modified from Sedgewicks
 * algorithm It selects the pivot based on the median of the left,
 * right, and center values to try to avoid degenerate cases ocurring
 * when a single value is chosen.  It performs a quicksort on
 * intervals down to the EX_QSORT_CUTOFF size and then performs a final
 * insertion sort on the almost sorted final array.  Based on data in
 * Sedgewick, the EX_QSORT_CUTOFF value should be between 5 and 20.
 *
 * See Sedgewick for further details
 * Define DEBUG_QSORT at the top of this file and recompile to compile
 * in code that verifies that the array is sorted.
 */

#include <iostream>
#include "index_sort.h"

#define EX_QSORT_CUTOFF 12

namespace {
  /* swap - interchange v[i] and v[j] */
  template <typename T>
  void ex_swap (T *v, int i, int j)
{
  T temp;

  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}

  template <typename T>
  int ex_int_median3(T *v, int iv[], int left, int right)
  {
    int center;
    center = (left + right) / 2;

    if (v[iv[left]] > v[iv[center]])
      ex_swap(iv, left, center);
    if (v[iv[left]] > v[iv[right]])
      ex_swap(iv, left, right);
    if (v[iv[center]] > v[iv[right]])
      ex_swap(iv, center, right);

    ex_swap(iv, center, right-1);
    return iv[right-1];
  }

  template <typename T>
  void ex_int_iqsort(T *v, int iv[], int left, int right)
  {
    int pivot;
    int i, j;
  
    if (left + EX_QSORT_CUTOFF <= right) {
      pivot = ex_int_median3(v, iv, left, right);
      i = left;
      j = right - 1;

      for ( ; ; ) {
	while (v[iv[++i]] < v[pivot]);
	while (v[iv[--j]] > v[pivot]);
	if (i < j) {
	  ex_swap(iv, i, j);
	} else {
	  break;
	}
      }

      ex_swap(iv, i, right-1);
      ex_int_iqsort(v, iv, left, i-1);
      ex_int_iqsort(v, iv, i+1, right);
    }
  }

  template <typename T>
  void ex_int_iisort(T *v, int iv[], int N)
  {
    int ndx = 0;
    int j;
    
    if (N == 0)
      return;

    double small = v[iv[0]];
    for (int i = 1; i < N; i++) {
      if (v[iv[i]] < small) {
	small = v[iv[i]];
	ndx = i;
      }
    }
    /* Put smallest value in slot 0 */
    ex_swap(iv, 0, ndx);

    for (int i=1; i <N; i++) {
      int tmp = iv[i];
      for (j=i; v[tmp] < v[iv[j-1]]; j--) {
	iv[j] = iv[j-1];
      }
      iv[j] = tmp;
    }
  }

  template <typename T>
  void ex_iqsort(T *v, int iv[], int N)
  {
    ex_int_iqsort(v, iv, 0, N-1);
    ex_int_iisort(v, iv, N);

#if defined(DEBUG_QSORT)
    fprintf(stderr, "Checking sort of %d values\n", N+1);
    int i;
    for (i=1; i < N; i++) {
      assert(v[iv[i-1]] <= v[iv[i]]);
    }
#endif
  }
}

void index_coord_sort(const std::vector<double> &xyz, std::vector<int> &index, int axis)
{
  // For now, let's extract the component we want to sort on into a separate vector.
  std::vector<double> comp(xyz.size()/3);
  size_t j = 0;
  for (size_t i=axis; i < xyz.size(); i+=3)
    comp[j++] = xyz[i];
  ex_iqsort(&comp[0], &index[0], (int)index.size());
}

void index_sort(const std::vector<int> &ids, std::vector<int> &index)
{
  index.resize(ids.size());
  for (size_t i=0; i < index.size(); i++)
    index[i] = i;
  
  ex_iqsort(&ids[0], &index[0], (int)index.size());
}
