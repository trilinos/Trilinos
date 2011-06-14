/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#if defined(DEBUG_QSORT)
#include <assert.h>
#endif

#include <string.h>

/*
 * The following 'indexed qsort' routine is modified from Sedgewicks
 * algorithm It selects the pivot based on the median of the left,
 * right, and center values to try to avoid degenerate cases ocurring
 * when a single value is chosen.  It performs a quicksort on
 * intervals down to the GDS_QSORT_CUTOFF size and then performs a final
 * insertion sort on the almost sorted final array.  Based on data in
 * Sedgewick, the GDS_QSORT_CUTOFF value should be between 5 and 20.
 *
 * See Sedgewick for further details
 * Define DEBUG_QSORT at the top of this file and recompile to compile
 * in code that verifies that the array is sorted.
 */

#define GDS_QSORT_CUTOFF 12
#define GDS_SWAP(V,I,J) do{int _t = V[I]; V[I] = V[J]; V[J] = _t;} while (0)

int gds_int_imedian3(int v[], int iv[], int left, int right)
{
  int center;
  center = (left + right) / 2;

  if (v[iv[left]] > v[iv[center]])
    GDS_SWAP(iv, left, center);
  if (v[iv[left]] > v[iv[right]])
    GDS_SWAP(iv, left, right);
  if (v[iv[center]] > v[iv[right]])
    GDS_SWAP(iv, center, right);

  GDS_SWAP(iv, center, right-1);
  return iv[right-1];
}

void gds_int_iqsort(int v[], int iv[], int left, int right)
{
  int pivot;
  int i, j;
  
  if (left + GDS_QSORT_CUTOFF <= right) {
    pivot = gds_int_imedian3(v, iv, left, right);
    i = left;
    j = right - 1;

    for ( ; ; ) {
      while (v[iv[++i]] < v[pivot]);
      while (v[iv[--j]] > v[pivot]);
      if (i < j) {
	GDS_SWAP(iv, i, j);
      } else {
	break;
      }
    }

    GDS_SWAP(iv, i, right-1);
    gds_int_iqsort(v, iv, left, i-1);
    gds_int_iqsort(v, iv, i+1, right);
  }
}

void gds_int_iisort(int v[], int iv[], int N)
{
  int i,j;
  int ndx = 0;
  int small;
  int tmp;
  
  if (N == 0) return;
  small = v[iv[0]];
  for (i = 1; i < N; i++) {
    if (v[iv[i]] < small) {
      small = v[iv[i]];
      ndx = i;
    }
  }
  /* Put smallest value in slot 0 */
  GDS_SWAP(iv, 0, ndx);

  for (i=1; i <N; i++) {
    tmp = iv[i];
    for (j=i; v[tmp] < v[iv[j-1]]; j--) {
      iv[j] = iv[j-1];
    }
    iv[j] = tmp;
  }
}

/*
 * Sort the values in 'v' using the index array 'iv'.
 * The values in 'v' will be unchanged, but the smallest
 * value will be v[iv[0]] and the largest v[iv[N-1]]
 *
 * The 'iv' array should be initialized 0..N-1 on entry.
 */
   
void gds_iqsort(int v[], int iv[], int N)
{
  if (N == 0) return;
  gds_int_iqsort(v, iv, 0, N-1);
  gds_int_iisort(v, iv, N);

#if defined(DEBUG_QSORT)
  fprintf(stderr, "Checking sort of %d values\n", N+1);
  int i;
  for (i=1; i < N; i++) {
    assert(v[iv[i-1]] <= v[iv[i]]);
  }
#endif
}

int  gds_ibin_search(int value, int v[], int iv[], int N)
{
  int l, r, x;
  l = 0;
  r = N-1;
  do {
    x = (l+r)/2;
     if (value < v[iv[x]]) {
      r = x-1;
    } else {
      l = x+1;
    }
  } while (value != v[iv[x]] && l <= r);

  if (value == v[iv[x]])
    return x;
  else
    return -1; /* to maintain compatibility with in_list() */
}

int  gds_bin_search(int value, int v[], int N)
{
  int l, r, x;
  l = 0;
  r = N-1;
  do {
    x = (l+r)/2;
    if (value < v[x]) {
      r = x-1;
    } else {
      l = x+1;
    }
  } while (value != v[x] && l <= r);

  if (value == v[x])
    return x;
  else
    return -1; /* to maintain compatibility with in_list() */
}

int gds_int_median3(int v[], int left, int right)
{
  int center;
  center = (left + right) / 2;

  if (v[left] > v[center])
    GDS_SWAP(v, left, center);
  if (v[left] > v[right])
    GDS_SWAP(v, left, right);
  if (v[center] > v[right])
    GDS_SWAP(v, center, right);

  GDS_SWAP(v, center, right-1);
  return right-1;
}

void gds_int_qsort(int v[], int left, int right)
{
  int pivot;
  int i, j;
  
  if (left + GDS_QSORT_CUTOFF <= right) {
    pivot = gds_int_median3(v, left, right);
    i = left;
    j = right - 1;

    for ( ; ; ) {
      while (v[++i] < v[pivot]);
      while (v[--j] > v[pivot]);
      if (i < j) {
	GDS_SWAP(v, i, j);
      } else {
	break;
      }
    }

    GDS_SWAP(v, i, right-1);
    gds_int_qsort(v, left, i-1);
    gds_int_qsort(v, i+1, right);
  }
}

void gds_int_isort(int v[], int N)
{
  int i,j;
  int ndx = 0;
  int small;
  int tmp;
  
  if (N == 0) return;
  small = v[0];
  for (i = 1; i < N; i++) {
    if (v[i] < small) {
      small = v[i];
      ndx = i;
    }
  }
  /* Put smallest value in slot 0 */
  GDS_SWAP(v, 0, ndx);

  for (i=1; i <N; i++) {
    tmp = v[i];
    for (j=i; tmp < v[j-1]; j--) {
      v[j] = v[j-1];
    }
    v[j] = tmp;
  }
}

/*
 * Sort the values in 'v' 
 */
   
void gds_qsort(int v[], int N)
{
  if (N == 0) return;
  gds_int_qsort(v, 0, N-1);
  gds_int_isort(v, N);

#if defined(DEBUG_QSORT)
  fprintf(stderr, "Checking sort of %d values\n", N+1);
  int i;
  for (i=1; i < N; i++) {
    assert(v[i-1] <= v[i]);
  }
#endif
}

