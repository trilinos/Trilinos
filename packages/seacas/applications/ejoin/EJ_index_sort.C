// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
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

#include <EJ_CodeTypes.h>
#include <cstdint>
#include <unistd.h>

#include "EJ_index_sort.h"

#define EX_QSORT_CUTOFF 12

namespace {
  /* swap - interchange v[i] and v[j] */
  template <typename T, typename INT> void ex_swap(T *v, INT i, INT j)
  {
    T temp;

    temp = v[i];
    v[i] = v[j];
    v[j] = temp;
  }

  template <typename T, typename INT> int ex_int_median3(T *v, INT iv[], size_t left, size_t right)
  {
    size_t center;
    center = (static_cast<ssize_t>(left) + static_cast<ssize_t>(right)) / 2;

    if (v[iv[left]] > v[iv[center]]) {
      ex_swap(iv, left, center);
    }
    if (v[iv[left]] > v[iv[right]]) {
      ex_swap(iv, left, right);
    }
    if (v[iv[center]] > v[iv[right]]) {
      ex_swap(iv, center, right);
    }

    ex_swap(iv, center, right - 1);
    return iv[right - 1];
  }

  template <typename T, typename INT> void ex_int_iqsort(T *v, INT iv[], size_t left, size_t right)
  {
    size_t pivot;
    size_t i;
    size_t j;

    if (left + EX_QSORT_CUTOFF <= right) {
      pivot = ex_int_median3(v, iv, left, right);
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
          ex_swap(iv, i, j);
        }
        else {
          break;
        }
      }

      ex_swap(iv, i, right - 1);
      ex_int_iqsort(v, iv, left, i - 1);
      ex_int_iqsort(v, iv, i + 1, right);
    }
  }

  template <typename T, typename INT> void ex_int_iisort(T *v, INT iv[], size_t N)
  {
    size_t ndx = 0;
    size_t j;

    if (N == 0) {
      return;
    }

    double small = v[iv[0]];
    for (size_t i = 1; i < N; i++) {
      if (v[iv[i]] < small) {
        small = v[iv[i]];
        ndx   = i;
      }
    }
    /* Put smallest value in slot 0 */
    ex_swap(iv, static_cast<size_t>(0), ndx);

    for (size_t i = 1; i < N; i++) {
      INT tmp = iv[i];
      for (j = i; v[tmp] < v[iv[j - 1]]; j--) {
        iv[j] = iv[j - 1];
      }
      iv[j] = tmp;
    }
  }

  template <typename T, typename INT> void ex__iqsort(T *v, INT iv[], size_t N)
  {
    if (N <= 1) {
      return;
    }

    ex_int_iqsort(v, iv, 0, N - 1);
    ex_int_iisort(v, iv, N);

#if defined(DEBUG_QSORT)
    fprintf(stderr, "Checking sort of %lu values\n", N + 1);
    int i;
    for (i = 1; i < N; i++) {
      assert(v[iv[i - 1]] <= v[iv[i]]);
    }
#endif
  }
} // namespace

template <typename INT>
void index_coord_sort(const std::vector<double> &xyz, std::vector<INT> &index, int axis)
{
  // For now, let's extract the component we want to sort on into a separate vector.
  std::vector<double> comp(xyz.size() / 3);
  size_t              j = 0;
  for (size_t i = axis; i < xyz.size(); i += 3) {
    comp[j++] = xyz[i];
  }
  ex__iqsort(&comp[0], &index[0], index.size());
}

template <typename INT> void index_sort(const std::vector<INT> &ids, std::vector<INT> &index)
{
  index.resize(ids.size());
  for (size_t i = 0; i < index.size(); i++) {
    index[i] = i;
  }

  ex__iqsort(&ids[0], &index[0], index.size());
}

template void index_coord_sort(const std::vector<double> &xyz, std::vector<int> &index, int axis);
template void index_coord_sort(const std::vector<double> &xyz, std::vector<int64_t> &index,
                               int axis);

template void index_sort(const std::vector<int> &ids, std::vector<int> &index);
template void index_sort(const std::vector<int64_t> &ids, std::vector<int64_t> &index);
