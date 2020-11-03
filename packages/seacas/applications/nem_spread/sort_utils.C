/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#define DEBUG_SORT 1
#if DEBUG_SORT
#include <cassert>
#include <cstdio>
#include <fmt/ostream.h>
#endif

#include <cstddef> // for size_t
#include <cstdint> // for int64_t

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

namespace {
  template <typename INT> void GDS_SWAP(INT *V, size_t I, size_t J)
  {
    INT _t = V[I];
    V[I]   = V[J];
    V[J]   = _t;
  }

  template <typename INT> INT gds_imedian3(INT v[], INT iv[], size_t left, size_t right)
  {
    size_t center;
    center = (left + right) / 2;

    if (v[iv[left]] > v[iv[center]]) {
      GDS_SWAP(iv, left, center);
    }
    if (v[iv[left]] > v[iv[right]]) {
      GDS_SWAP(iv, left, right);
    }
    if (v[iv[center]] > v[iv[right]]) {
      GDS_SWAP(iv, center, right);
    }

    GDS_SWAP(iv, center, right - 1);
    return iv[right - 1];
  }

  template <typename INT> void gds_iqsort(INT v[], INT iv[], size_t left, size_t right)
  {
    size_t pivot;
    size_t i, j;

    if (left + GDS_QSORT_CUTOFF <= right) {
      pivot = gds_imedian3(v, iv, left, right);
      i     = left;
      j     = right - 1;

      for (;;) {
        while (v[iv[++i]] < v[pivot]) {
          ;
        }
        while (v[iv[--j]] > v[pivot]) {
          ;
        }
        if (i < j) {
          GDS_SWAP(iv, i, j);
        }
        else {
          break;
        }
      }

      GDS_SWAP(iv, i, right - 1);
      gds_iqsort(v, iv, left, i - 1);
      gds_iqsort(v, iv, i + 1, right);
    }
  }

  template <typename INT> void gds_iisort(INT v[], INT iv[], size_t N)
  {
    size_t i, j;
    size_t ndx = 0;
    INT    low;
    INT    tmp;

    if (N == 0) {
      return;
    }
    low = v[iv[0]];
    for (i = 1; i < N; i++) {
      if (v[iv[i]] < low) {
        low = v[iv[i]];
        ndx = i;
      }
    }
    /* Put lowest value in slot 0 */
    GDS_SWAP(iv, 0, ndx);

    for (i = 1; i < N; i++) {
      tmp = iv[i];
      for (j = i; v[tmp] < v[iv[j - 1]]; j--) {
        iv[j] = iv[j - 1];
      }
      iv[j] = tmp;
    }
  }

  /*
   * Sort the values in 'v' using the index array 'iv'.
   * The values in 'v' will be unchanged, but the lowest
   * value will be v[iv[0]] and the largest v[iv[N-1]]
   *
   * The 'iv' array should be initialized 0..N-1 on entry.
   */

  template <typename INT> size_t gds_median3(INT v[], size_t left, size_t right)
  {
    size_t center;
    center = (left + right) / 2;

    if (v[left] > v[center]) {
      GDS_SWAP(v, left, center);
    }
    if (v[left] > v[right]) {
      GDS_SWAP(v, left, right);
    }
    if (v[center] > v[right]) {
      GDS_SWAP(v, center, right);
    }

    GDS_SWAP(v, center, right - 1);
    return right - 1;
  }

  template <typename INT> void gds_qsort(INT v[], size_t left, size_t right)
  {
    if (left + GDS_QSORT_CUTOFF <= right) {
      size_t pivot = gds_median3(v, left, right);
      size_t i     = left;
      size_t j     = right - 1;

      for (;;) {
        while (v[++i] < v[pivot]) {
          ;
        }
        while (v[--j] > v[pivot]) {
          ;
        }
        if (i < j) {
          GDS_SWAP(v, i, j);
        }
        else {
          break;
        }
      }

      GDS_SWAP(v, i, right - 1);
      gds_qsort(v, left, i - 1);
      gds_qsort(v, i + 1, right);
    }
  }

  template <typename INT> void gds_isort(INT v[], size_t N)
  {
    size_t i, j;
    size_t ndx = 0;
    INT    low;
    INT    tmp;

    if (N <= 1) {
      return;
    }
    low = v[0];
    for (i = 1; i < N; i++) {
      if (v[i] < low) {
        low = v[i];
        ndx = i;
      }
    }
    /* Put lowest value in slot 0 */
    GDS_SWAP(v, 0, ndx);

    for (i = 1; i < N; i++) {
      tmp = v[i];
      for (j = i; tmp < v[j - 1]; j--) {
        v[j] = v[j - 1];
      }
      v[j] = tmp;
    }
  }

  template <typename INT> void siftDown(INT *v, INT *iv, size_t start, size_t end)
  {
    size_t root = start;

    while (root * 2 + 1 < end) {
      size_t child = 2 * root + 1;
      if ((child + 1 < end) && (v[iv[child]] < v[iv[child + 1]])) {
        child += 1;
      }
      if (v[iv[root]] < v[iv[child]]) {
        GDS_SWAP(iv, child, root);
        root = child;
      }
      else {
        return;
      }
    }
  }
} // namespace

/*
 * Sort the values in 'v'
 */

template <typename INT> void indexed_sort(INT v[], INT iv[], size_t N)
{
  int64_t start, end;
  int64_t count = N;

  if (N <= 1) {
    return;
  }

  /* heapify */
  for (start = (count - 2) / 2; start >= 0; start--) {
    siftDown(v, iv, start, count);
  }

  for (end = count - 1; end > 0; end--) {
    GDS_SWAP(iv, end, 0);
    siftDown(v, iv, 0, end);
  }

#if DEBUG_SORT
  fmt::print(stderr, "Checking sort of {:n} values\n", count + 1);
  for (size_t i = 1; i < N; i++) {
    assert(v[iv[i - 1]] <= v[iv[i]]);
  }
#endif
}

template <typename INT> void gds_iqsort(INT v[], INT iv[], size_t N)
{
  if (N <= 1) {
    return;
  }
  gds_iqsort(v, iv, 0, N - 1);
  gds_iisort(v, iv, N);

#if defined(DEBUG_QSORT)
  fmt::print(stderr, "Checking sort of {:n} values\n", N + 1);
  size_t i;
  for (i = 1; i < N; i++) {
    assert(v[iv[i - 1]] <= v[iv[i]]);
  }
#endif
}

template <typename INT> void gds_qsort(INT v[], size_t N)
{
  if (N <= 1) {
    return;
  }
  gds_qsort(v, 0, N - 1);
  gds_isort(v, N);

#if defined(DEBUG_QSORT)
  fmt::print(stderr, "Checking sort of {:n} values\n", N + 1);
  for (size_t i = 1; i < N; i++) {
    assert(v[i - 1] <= v[i]);
  }
#endif
}

template void gds_iqsort(int v[], int iv[], size_t N);
template void gds_qsort(int v[], size_t N);

template void gds_iqsort(int64_t v[], int64_t iv[], size_t N);
template void gds_qsort(int64_t v[], size_t N);

template void indexed_sort(int v[], int iv[], size_t N);
template void indexed_sort(int64_t v[], int64_t iv[], size_t N);
