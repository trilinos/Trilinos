// Copyright(C) 1999-2020, 2022, 2023, 2024, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include <gtest/gtest.h>
#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif
#include <string>
#include <unordered_map>

#include <memory>
#include <string>
#if defined(_WIN32) && !defined(__MINGW32__)
#include <string.h>
#define strcasecmp  _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif
#include <ostream>
#include <vector>

#include <cstdarg> // for va_end, va_arg, va_list, etc
#include <cstddef> // for size_t
#include <cstdio>  // for stderr
#include <cstdlib> // for exit, malloc
#include <fmt/format.h>
#include <fmt/ostream.h>

namespace utest_util {

  template <typename INT>
  inline int numbermatch(INT *sidenodes, size_t i, size_t j, size_t k, size_t value)
  {
    if ((size_t)sidenodes[(i + j) % k] == value) {
      return 1;
    }
    return 0;
  }

  template <typename INT> inline void vec_free(std::vector<INT> &V)
  {
    V.clear();
    V.shrink_to_fit();
  }

#ifndef NDEBUG
  template <typename INT> inline void assert_sorted(INT *vector, size_t vecsize)
  {
    for (size_t i = 1; i < vecsize; i++) {
      assert(vector[i - 1] <= vector[i]);
    }
  }
#endif

  /* Function find_intersection() begins:
   *----------------------------------------------------------------------------
   * This function finds the intersection between two lists of integer values,
   * and returns the number of values in the intersection.
   *****************************************************************************/
  template <typename INT>
  inline size_t find_intersection(const INT set1[],  /* the first set of integers */
                                  const INT set2[],  /* the second set of integers */
                                  size_t    length1, /* the length of the first set */
                                  size_t    length2, /* the length of the second set */
                                  INT       inter_ptr[])   /* the values in the intersection */
  /*
   *
   *      Function which finds the intersection of two integer lists.
   *      The points in set1 that belong in the intersection set are
   *      returned in the vector inter_pts, starting at position inter_pts[0].
   *      Enough space in inter_pts[] (min(length1, length2)) must
   *      have already been allocated in the calling program before this
   *      function is called.
   *
   *      Know that set1 and set2 are monotonically increasing
   *
   *      On return, find_inter returns 0 if there is no intersection.
   *      It returns the number of points in the intersection, if there
   *      is an intersection.
   */

  {
    size_t counter = 0;
    size_t i       = 0;
    size_t j       = 0;

    while (i < length1 && j < length2) {
      if (set1[i] < set2[j]) {
        ++i;
      }
      else if (set2[j] < set1[i]) {
        ++j;
      }
      else {
        inter_ptr[counter++] = i;
        ++i;
        ++j;
      }
    }
    return counter;
  }

  /* Function in_list() begins:
   *----------------------------------------------------------------------------
   * This function searches a vector for the input value. If the value is
   * found in the vector then it's index in that vector is returned, otherwise
   * the function returns -1;
   *****************************************************************************/
  template <typename INT, typename INT2>
  inline int64_t in_list(INT value, size_t count, const INT2 *vector)
  {
    for (size_t i = 0; i < count; i++) {
      if (vector[i] == value) {
        return i;
      }
    }
    return -1;
  }

  template <typename INT, typename INT2>
  inline int64_t in_list(INT value, const std::vector<INT2> &vector)
  {
    size_t count = vector.size();
    for (size_t i = 0; i < count; i++) {
      if (vector[i] == value) {
        return i;
      }
    }
    return -1;
  }

  template <typename INT> inline void SWAP(INT &r, INT &s)
  {
    INT t = r;
    r     = s;
    s     = t;
  }

  /////////////////////
  /*
   * The following 'qsort' routine is modified from Sedgewicks
   * algorithm It selects the pivot based on the median of the left,
   * right, and center values to try to avoid degenerate cases occurring
   * when a single value is chosen.  It performs a quicksort on
   * intervals down to the GDS_QSORT_CUTOFF size and then performs a final
   * insertion sort on the almost sorted final array.  Based on data in
   * Sedgewick, the GDS_QSORT_CUTOFF value should be between 5 and 20.
   *
   * See Sedgewick for further details
   */

#define GDS_QSORT_CUTOFF 12

  template <typename INT> inline void ISWAP(INT *V, size_t I, size_t J)
  {
    INT _t = V[I];
    V[I]   = V[J];
    V[J]   = _t;
  }

  template <typename INT> inline size_t gds_median3(INT v[], size_t left, size_t right)
  {
    size_t center;
    center = (left + right) / 2;

    if (v[left] > v[center]) {
      ISWAP(v, left, center);
    }
    if (v[left] > v[right]) {
      ISWAP(v, left, right);
    }
    if (v[center] > v[right]) {
      ISWAP(v, center, right);
    }

    ISWAP(v, center, right - 1);
    return right - 1;
  }

  template <typename INT> inline void gds_qsort(INT v[], size_t left, size_t right)
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
          ISWAP(v, i, j);
        }
        else {
          break;
        }
      }

      ISWAP(v, i, right - 1);
      gds_qsort(v, left, i - 1);
      gds_qsort(v, i + 1, right);
    }
  }

  template <typename INT> inline void gds_isort(INT v[], size_t N)
  {
    if (N <= 1) {
      return;
    }

    size_t ndx       = 0;
    INT    small_val = v[0];
    for (size_t i = 1; i < N; i++) {
      if (v[i] < small_val) {
        small_val = v[i];
        ndx       = i;
      }
    }
    /* Put smallest value in slot 0 */
    ISWAP(v, 0, ndx);

    size_t j;
    for (size_t i = 1; i < N; i++) {
      INT tmp = v[i];
      for (j = i; tmp < v[j - 1]; j--) {
        assert(j >= 1);
        v[j] = v[j - 1];
      }
      v[j] = tmp;
    }
  }

  template <typename INT> inline void siftDowniii(INT *a, INT *b, INT *c, size_t start, size_t end)
  {
    size_t root = start;
    while (root * 2 + 1 < end) {
      size_t child = 2 * root + 1;
      if ((child + 1 < end) &&
          (a[child] < a[child + 1] || (a[child] == a[child + 1] && b[child] < b[child + 1]))) {
        child += 1;
      }
      if (a[root] < a[child]) {
        SWAP(a[child], a[root]);
        SWAP(b[child], b[root]);
        SWAP(c[child], c[root]);
        root = child;
      }
      else {
        return;
      }
    }
  }

  template <typename INT> inline void siftDown(INT *a, INT *b, size_t start, size_t end)
  {
    size_t root = start;

    while (root * 2 + 1 < end) {
      size_t child = 2 * root + 1;
      if ((child + 1 < end) && (a[child] < a[child + 1])) {
        child += 1;
      }
      if (a[root] < a[child]) {
        SWAP(a[child], a[root]);
        SWAP(b[child], b[root]);
        root = child;
      }
      else {
        return;
      }
    }
  }

  template <typename INT> inline void sort2(int64_t count, INT ra[], INT rb[])
  {
    if (count <= 1) {
      return;
    }
    /* heapify */
    for (int64_t start = (count - 2) / 2; start >= 0; start--) {
      siftDown(ra, rb, start, count);
    }

    for (size_t end = count - 1; end > 0; end--) {
      SWAP(ra[end], ra[0]);
      SWAP(rb[end], rb[0]);
      siftDown(ra, rb, 0, end);
    }
  }

#define QSORT_CUTOFF 12

  template <typename INT>
  inline int is_less_than4(INT ra1, INT rb1, INT rc1, INT rd1, INT ra2, INT rb2, INT rc2, INT rd2)
  {
    if (ra1 < ra2) {
      return 1;
    }
    if (ra1 > ra2) {
      return 0;
    }
    assert(ra1 == ra2);

    if (rb1 < rb2) {
      return 1;
    }
    if (rb1 > rb2) {
      return 0;
    }
    assert(rb1 == rb2);

    if (rc1 < rc2) {
      return 1;
    }
    if (rc1 > rc2) {
      return 0;
    }
    assert(rc1 == rc2);

    if (rd1 < rd2) {
      return 1;
    }

    return 0;
  }

  template <typename INT>
  inline int is_less_than4v(INT *v1, INT *v2, INT *v3, INT *v4, size_t i, size_t j)
  {
    if (v1[i] < v1[j]) {
      return 1;
    }
    if (v1[i] > v1[j]) {
      return 0;
    }
    assert(v1[i] == v1[j]);

    if (v2[i] < v2[j]) {
      return 1;
    }
    if (v2[i] > v2[j]) {
      return 0;
    }
    assert(v2[i] == v2[j]);

    if (v3[i] < v3[j]) {
      return 1;
    }
    if (v3[i] > v3[j]) {
      return 0;
    }
    assert(v3[i] == v3[j]);

    if (v4[i] < v4[j]) {
      return 1;
    }

    return 0;
  }

  template <typename INT> inline void swap4(INT *v1, INT *v2, INT *v3, INT *v4, size_t i, size_t j)
  {
    ISWAP(v1, i, j);
    ISWAP(v2, i, j);
    ISWAP(v3, i, j);
    ISWAP(v4, i, j);
  }

  template <typename INT>
  inline size_t internal_median3_4(INT *v1, INT *v2, INT *v3, INT *v4, size_t left, size_t right)
  {
    size_t center;
    center = (left + right) / 2;

    if (is_less_than4v(v1, v2, v3, v4, center, left)) {
      swap4(v1, v2, v3, v4, left, center);
    }
    if (is_less_than4v(v1, v2, v3, v4, right, left)) {
      swap4(v1, v2, v3, v4, left, right);
    }
    if (is_less_than4v(v1, v2, v3, v4, right, center)) {
      swap4(v1, v2, v3, v4, center, right);
    }

    swap4(v1, v2, v3, v4, center, right - 1);
    return right - 1;
  }

  template <typename INT>
  inline void internal_qsort_4(INT *v1, INT *v2, INT *v3, INT *v4, size_t left, size_t right)
  {
    if (left + QSORT_CUTOFF <= right) {
      size_t pivot = internal_median3_4(v1, v2, v3, v4, left, right);
      size_t i     = left;
      size_t j     = right - 1;

      for (;;) {
        while (is_less_than4v(v1, v2, v3, v4, ++i, pivot)) {
          ;
        }
        while (is_less_than4v(v1, v2, v3, v4, pivot, --j)) {
          ;
        }
        if (i < j) {
          swap4(v1, v2, v3, v4, i, j);
        }
        else {
          break;
        }
      }

      swap4(v1, v2, v3, v4, i, right - 1);
      internal_qsort_4(v1, v2, v3, v4, left, i - 1);
      internal_qsort_4(v1, v2, v3, v4, i + 1, right);
    }
  }

  template <typename INT> inline void internal_isort_4(INT *v1, INT *v2, INT *v3, INT *v4, size_t N)
  {
    size_t ndx = 0;
    for (size_t i = 1; i < N; i++) {
      if (is_less_than4v(v1, v2, v3, v4, i, ndx)) {
        ndx = i;
      }
    }
    /* Put small_valest value in slot 0 */
    swap4(v1, v2, v3, v4, 0, ndx);

    for (size_t i = 1; i < N; i++) {
      INT    small_val1 = v1[i];
      INT    small_val2 = v2[i];
      INT    small_val3 = v3[i];
      INT    small_val4 = v4[i];
      size_t j;
      for (j = i; is_less_than4(small_val1, small_val2, small_val3, small_val4, v1[j - 1],
                                v2[j - 1], v3[j - 1], v4[j - 1]);
           j--) {
        assert(j >= 1);
        v1[j] = v1[j - 1];
        v2[j] = v2[j - 1];
        v3[j] = v3[j - 1];
        v4[j] = v4[j - 1];
      }
      v1[j] = small_val1;
      v2[j] = small_val2;
      v3[j] = small_val3;
      v4[j] = small_val4;
    }
  }

  template <typename INT> inline int is_less_than2(INT ra1, INT rb1, INT ra2, INT rb2)
  {
    if (ra1 < ra2) {
      return 1;
    }
    if (ra1 > ra2) {
      return 0;
    }
    assert(ra1 == ra2);

    if (rb1 < rb2) {
      return 1;
    }

    return 0;
  }

  template <typename INT> inline int is_less_than2v(INT *v1, INT *v2, size_t i, size_t j)
  {
    if (v1[i] < v1[j]) {
      return 1;
    }
    if (v1[i] > v1[j]) {
      return 0;
    }
    assert(v1[i] == v1[j]);

    if (v2[i] < v2[j]) {
      return 1;
    }

    return 0;
  }

  template <typename INT> inline void swap2(INT *v1, INT *v2, size_t i, size_t j)
  {
    ISWAP(v1, i, j);
    ISWAP(v2, i, j);
  }

  template <typename INT>
  inline size_t internal_median3_2(INT *v1, INT *v2, size_t left, size_t right)
  {
    size_t center = (left + right) / 2;

    if (is_less_than2v(v1, v2, center, left)) {
      swap2(v1, v2, left, center);
    }
    if (is_less_than2v(v1, v2, right, left)) {
      swap2(v1, v2, left, right);
    }
    if (is_less_than2v(v1, v2, right, center)) {
      swap2(v1, v2, center, right);
    }

    swap2(v1, v2, center, right - 1);
    return right - 1;
  }

  template <typename INT> inline void internal_qsort_2(INT *v1, INT *v2, size_t left, size_t right)
  {
    if (left + QSORT_CUTOFF <= right) {
      size_t pivot = internal_median3_2(v1, v2, left, right);
      size_t i     = left;
      size_t j     = right - 1;

      for (;;) {
        while (is_less_than2v(v1, v2, ++i, pivot)) {
          ;
        }
        while (is_less_than2v(v1, v2, pivot, --j)) {
          ;
        }
        if (i < j) {
          swap2(v1, v2, i, j);
        }
        else {
          break;
        }
      }

      swap2(v1, v2, i, right - 1);
      internal_qsort_2(v1, v2, left, i - 1);
      internal_qsort_2(v1, v2, i + 1, right);
    }
  }

  template <typename INT> inline void internal_isort_2(INT *v1, INT *v2, size_t N)
  {
    size_t ndx = 0;
    for (size_t i = 1; i < N; i++) {
      if (is_less_than2v(v1, v2, i, ndx)) {
        ndx = i;
      }
    }

    /* Put small_valest value in slot 0 */
    swap2(v1, v2, 0, ndx);

    for (size_t i = 1; i < N; i++) {
      INT    small_val1 = v1[i];
      INT    small_val2 = v2[i];
      size_t j;
      for (j = i; is_less_than2(small_val1, small_val2, v1[j - 1], v2[j - 1]); j--) {
        assert(j >= 1);
        v1[j] = v1[j - 1];
        v2[j] = v2[j - 1];
      }
      v1[j] = small_val1;
      v2[j] = small_val2;
    }
  }

  /*
   * Sort the values in 'v'
   */

  template <typename INT> inline void qsort4(INT *v1, INT *v2, INT *v3, INT *v4, size_t N)
  {
    if (N <= 1) {
      return;
    }
    internal_qsort_4(v1, v2, v3, v4, 0, N - 1);
    internal_isort_4(v1, v2, v3, v4, N);

#if defined(DEBUG_QSORT)
    fmt::print(stderr, "Checking sort of {} values\n", (size_t)N + 1);
    for (size_t i = 1; i < N; i++) {
      assert(is_less_than4v(v1, v2, v3, v4, i - 1, i));
    }
#endif
  }

  template <typename INT> inline void qsort2(INT *v1, INT *v2, size_t N)
  {
    if (N <= 1) {
      return;
    }
    internal_qsort_2(v1, v2, 0, N - 1);
    internal_isort_2(v1, v2, N);

#if defined(DEBUG_QSORT)
    fmt::print(stderr, "Checking sort of {} values\n", (size_t)N + 1);
    for (size_t i = 1; i < N; i++) {
      assert(is_less_than2v(v1, v2, i - 1, i));
    }
#endif
  }

  template <typename INT> inline void sort3(int64_t count, INT ra[], INT rb[], INT rc[])
  {
    if (count <= 1) {
      return;
    }
    /* heapify */
    for (int64_t start = (count - 2) / 2; start >= 0; start--) {
      siftDowniii(ra, rb, rc, start, count);
    }

    for (size_t end = count - 1; end > 0; end--) {
      SWAP(ra[end], ra[0]);
      SWAP(rb[end], rb[0]);
      SWAP(rc[end], rc[0]);
      siftDowniii(ra, rb, rc, 0, end);
    }
  }

  template <typename INT> inline int64_t bin_search2(INT value, size_t num, INT List[])
  {
    /*
     * Searches a monotonic list of values for the value, value.
     * It returns the index of the first position found, which matches value.
     * The list is assumed to be monotonic, and
     * consist of elements list[0], ..., list[n-1].
     * If no position in list matches value, it returns the value -1.
     *
     */
    size_t top    = num - 1;
    size_t bottom = 0;
    while (bottom <= top) {
      size_t middle = (bottom + top) >> 1;
      INT    g_mid  = List[middle];
      if (value < g_mid) {
        top = middle - 1;
      }
      else if (value > g_mid) {
        bottom = middle + 1;
      }
      else {
        return middle; /* found */
      }
    }
    return -1;
  }

  template <typename INT> inline void gds_qsort(INT v[], size_t N)
  {
    if (N <= 1) {
      return;
    }
    gds_qsort(v, 0, N - 1);
    gds_isort(v, N);
  }

  /*****************************************************************************
   * Function to find the first and last entries of a given value that are
   * consecutively present in an integer array.
   * ASSUMES that 'vector' is sorted....
   *****************************************************************************/

  template <typename INT>
  inline void find_first_last(INT val, size_t vecsize, INT *vector, INT *first, INT *last)
  {
#ifndef NDEBUG
    assert_sorted(vector, vecsize);
#endif

    *first = -1;
    *last  = -1;

    /* See if value is in the vector */
    int64_t i = bin_search2(val, vecsize, vector);
    *first    = i; /* Save this location */

    if (i != -1) {
      /* Value is in vector, find first occurrence */
      while (i >= 0 && vector[i] == val) {
        i--;
      }
      i++;

      *last  = *first; /* Use saved location */
      *first = i;

      size_t ii;
      for (ii = (*last); ii < vecsize; ii++) {
        if (vector[ii] != val) {
          *last = ii - 1;
          break;
        }
      }

      if (ii == vecsize) {
        *last = vecsize - 1;
      }
    }
  }

} // namespace utest_util
