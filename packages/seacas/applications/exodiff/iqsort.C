// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#include "iqsort.h"
#include <cstdint>
#include <fmt/ostream.h>

namespace {
  template <typename INT> void swap_(INT v[], size_t i, size_t j);

  template <typename T, typename INT> INT median3(const T v[], INT iv[], size_t left, size_t right);

  template <typename T, typename INT> void iqsort(const T v[], INT iv[], size_t left, size_t right);

  template <typename T, typename INT> void iisort(const T v[], INT iv[], size_t N);

  template <typename T, typename INT> void check(const T v[], INT iv[], size_t N);
} // namespace

template <typename T, typename INT> void index_qsort(const T v[], INT iv[], size_t N)
{
  if (N <= 1) {
    return;
  }
  iqsort(v, iv, 0, N - 1);
  iisort(v, iv, N);
#if defined(DEBUG_QSORT)
  check(v, iv, N);
#endif
}

// The following are not part of the public interFace...

namespace {
  /* The following 'indexed qsort' routine is modified from Sedgewicks
   * algorithm It selects the pivot based on the median of the left,
   * right, and center values to try to avoid degenerate cases ocurring
   * when a single value is chosen.  It performs a quicksort on
   * intervals down to the QSORT_CUTOFF size and then performs a final
   * insertion sort on the almost sorted final array.  Based on data in
   * Sedgewick, the QSORT_CUTOFF value should be between 5 and 20.
   *
   * See Sedgewick for further details Define DEBUG_QSORT at the top of
   * this file and recompile to compile in code that verifies that the
   * array is sorted.
   */

#define QSORT_CUTOFF 12

  /* swap - interchange v[i] and v[j] */
  template <typename INT> void swap_(INT v[], size_t i, size_t j)
  {
    INT temp;

    temp = v[i];
    v[i] = v[j];
    v[j] = temp;
  }

  template <typename T, typename INT> INT median3(const T v[], INT iv[], size_t left, size_t right)
  {
    size_t center = (left + right) / 2;

    if (v[iv[left]] > v[iv[center]]) {
      swap_(iv, left, center);
    }
    if (v[iv[left]] > v[iv[right]]) {
      swap_(iv, left, right);
    }
    if (v[iv[center]] > v[iv[right]]) {
      swap_(iv, center, right);
    }

    swap_(iv, center, right - 1);
    return iv[right - 1];
  }

  template <typename T, typename INT> void iqsort(const T v[], INT iv[], size_t left, size_t right)
  {
    if (left + QSORT_CUTOFF <= right) {
      size_t pivot = median3(v, iv, left, right);
      size_t i     = left;
      size_t j     = right - 1;

      for (;;) {
        while (v[iv[++i]] < v[pivot]) {
          ;
        }
        while (v[iv[--j]] > v[pivot]) {
          ;
        }
        if (i < j) {
          swap_(iv, i, j);
        }
        else {
          break;
        }
      }

      swap_(iv, i, right - 1);
      iqsort(v, iv, left, i - 1);
      iqsort(v, iv, i + 1, right);
    }
  }

  template <typename T, typename INT> void iisort(const T v[], INT iv[], size_t N)
  {
    size_t j;
    size_t ndx = 0;

    T low = v[iv[0]];
    for (size_t i = 1; i < N; i++) {
      if (v[iv[i]] < low) {
        low = v[iv[i]];
        ndx = i;
      }
    }
    /* Put lowest value in slot 0 */
    swap_(iv, 0, ndx);

    for (size_t i = 1; i < N; i++) {
      INT tmp = iv[i];
      for (j = i; v[tmp] < v[iv[j - 1]]; j--) {
        iv[j] = iv[j - 1];
      }
      iv[j] = tmp;
    }
  }

  template <typename T, typename INT> void check(const T v[], INT iv[], size_t N)
  {
    fmt::print(stderr, "Checking sort of {:n} values\n", N + 1);
    size_t i;
    for (i = 1; i < N; i++) {
      SMART_ASSERT(v[iv[i - 1]] <= v[iv[i]]);
    }
  }
} // namespace

template void index_qsort(const int v[], int iv[], size_t N);
template void index_qsort(const double v[], int iv[], size_t N);

template void index_qsort(const int64_t v[], int64_t iv[], size_t N);
template void index_qsort(const double v[], int64_t iv[], size_t N);
