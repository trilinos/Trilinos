/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *      token_compare()
 *      strip_string()
 *      string_to_lower()
 *      clean_string()
 *      sort2()
 *      sort3()
 *      find_first_last()
 *      find_int()
 *      in_list()
 *      roundfloat()
 *      find_max()
 *      find_min()
 *      find_inter()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "elb_util.h"
#include <cassert> // for assert
#include <cctype>  // for isupper, tolower
#include <cmath>   // for ceil, floor
#include <cstddef> // for size_t
#include <cstring> // for strlen
#include <fmt/ostream.h>

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int token_compare(char *token, const char *key)
{

  size_t kcnt = 0;

  size_t key_len = strlen(key);

  for (size_t i1 = 0; i1 < strlen(token); i1++) {
    if (isupper(token[i1]) != 0) {
      token[i1] = tolower(token[i1]);
    }

    if (token[i1] != ' ') {
      if (token[i1] == key[kcnt]) {
        kcnt++;
        if (kcnt > key_len) {
          return 0;
        }
      }
      else {
        return 0;
      }
    }
    if (key[kcnt] == ' ') {
      kcnt++;
    }
  }

  if (kcnt == strlen(key)) {
    return 1;
  }
  return 0;
} /*--------------End token_compare()-----------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void strip_string(char inp_str[], const char *tokens)
{
  int i;
  int j;
  int itok;
  int ntokes;
  int bval;

  i      = 0;
  ntokes = strlen(tokens);

  while (inp_str[i] != '\0') {
    bval = 0;
    for (itok = 0; itok < ntokes; itok++) {
      if (inp_str[i] == tokens[itok]) {
        i++;
        bval = 1;
        break; /* out of for loop */
      }
    }
    if (bval == 0) {
      break; /* out of while loop */
    }
  }

  /* Move real part of string to the front */
  j = 0;
  while (inp_str[j + i] != '\0') {
    inp_str[j] = inp_str[j + i];
    j++;
  }
  inp_str[j] = inp_str[j + i];
  j--;

  /* Remove trailing tokens */
  while (j != -1) {
    bval = 0;
    for (itok = 0; itok < ntokes; itok++) {
      if (inp_str[j] == tokens[itok]) {
        bval = 1;
        j--;
        break; /* out of for loop */
      }
    }
    if (bval == 0) {
      break; /* out of while loop */
    }
  }

  inp_str[j + 1] = '\0';
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void string_to_lower(char in_string[], char cval)
{
  int len;
  int cnt;

  len = strlen(in_string);
  for (cnt = 0; cnt < len; cnt++) {
    if (in_string[cnt] == cval) {
      return;
    }

    if (isupper(in_string[cnt]) != 0) {
      in_string[cnt] = tolower(in_string[cnt]);
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void clean_string(char inp_str[], const char *tokens)
{
  int i;
  int j;
  int itok;
  int ntokes;
  int bval;
  int inplen;

  ntokes = strlen(tokens);
  inplen = strlen(inp_str);

  i    = 0;
  bval = 0;
  while (inp_str[i] != '\0') {
    for (itok = 0; itok < ntokes; itok++) {
      if (i < 0) {
        i = 0;
      }
      if (inp_str[i] == tokens[itok]) {
        /* Find out if the next character is also a token */
        for (j = 0; j < ntokes; j++) {
          if (inp_str[i + 1] == tokens[j]) {
            bval = 1;
            break;
          }
        }

        if (bval == 1) {
          for (j = i + 1; j < inplen; j++) {
            inp_str[j] = inp_str[j + 1];
          }

          inplen--;
          bval = 0;
          i--;
        }
      }
    }

    i++;

  } /* End "while(inp_str[i] != '\0')" */

} /*---------------- End clean_string() -----------------*/

namespace {
  /*
   * The following 'qsort' routine is modified from Sedgewicks
   * algorithm It selects the pivot based on the median of the left,
   * right, and center values to try to avoid degenerate cases ocurring
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

  template <typename INT> size_t gds_median3(INT v[], size_t left, size_t right)
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

  template <typename INT> void gds_qsort(INT v[], size_t left, size_t right)
  {
    size_t pivot;
    size_t i;
    size_t j;

    if (left + GDS_QSORT_CUTOFF <= right) {
      pivot = gds_median3(v, left, right);
      i     = left;
      j     = right - 1;

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

  template <typename INT> void gds_isort(INT v[], size_t N)
  {
    size_t i;
    size_t j;
    size_t ndx = 0;
    INT    small_val;
    INT    tmp;

    if (N <= 1) {
      return;
    }
    small_val = v[0];
    for (i = 1; i < N; i++) {
      if (v[i] < small_val) {
        small_val = v[i];
        ndx       = i;
      }
    }
    /* Put smallest value in slot 0 */
    ISWAP(v, 0, ndx);

    for (i = 1; i < N; i++) {
      tmp = v[i];
      for (j = i; tmp < v[j - 1]; j--) {
        v[j] = v[j - 1];
      }
      v[j] = tmp;
    }
  }

  template <typename INT> void siftDowniii(INT *a, INT *b, INT *c, size_t start, size_t end)
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

  template <typename INT> void assert_sorted(INT *vector, size_t vecsize)
  {
    for (size_t i = 1; i < vecsize; i++) {
      assert(vector[i - 1] <= vector[i]);
    }
  }
} // namespace

/*****************************************************************************
 * Function to find the first and last entries of a given value that are
 * consecutively present in an integer array.
 * ASSUMES that 'vector' is sorted....
 *****************************************************************************/
template void find_first_last(int val, size_t vecsize, int *vector, int *first, int *last);
template void find_first_last(int64_t val, size_t vecsize, int64_t *vector, int64_t *first,
                              int64_t *last);

template <typename INT>
void find_first_last(INT val, size_t vecsize, INT *vector, INT *first, INT *last)
{
  /* assert_sorted(vector, vecsize); */

  *first = -1;
  *last  = -1;

  /* See if value is in the vector */
  ssize_t i = bin_search2(val, vecsize, vector);
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

/*****************************************************************************
 * Find the value in the integer array. Return -1 if not found.
 * New 1/21/97: change this so that it cross references a second
 *              array with a second value
 *****************************************************************************/
template <typename INT>
ssize_t find_int(INT value1, INT value2, size_t start, size_t stop, INT *vector1, INT *vector2)
{
  for (size_t i = start; i <= stop; i++) {
    if ((vector1[i] == value1) && (vector2[i] == value2)) {
      return i;
    }
  }
  return -1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function in_list() begins:
 *----------------------------------------------------------------------------
 * This function searches a vector for the input value. If the value is
 * found in the vector then it's index in that vector is returned, otherwise
 * the function returns -1;
 *****************************************************************************/
template ssize_t in_list(int value, size_t count, int *vector);
template ssize_t in_list(int64_t value, size_t count, int64_t *vector);

template <typename INT> ssize_t in_list(INT value, size_t count, INT *vector)
{
  for (size_t i = 0; i < count; i++) {
    if (vector[i] == value) {
      return i;
    }
  }
  return -1;
}

template ssize_t in_list(int value, std::vector<int> vector);
template ssize_t in_list(int64_t value, std::vector<int64_t> vector);

template <typename INT> ssize_t in_list(INT value, std::vector<INT> vector)
{
  size_t count = vector.size();
  for (size_t i = 0; i < count; i++) {
    if (vector[i] == value) {
      return i;
    }
  }
  return -1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function roundfloat() begins:
 *----------------------------------------------------------------------------
 * This function rounds off the float "value" to the nearest integer,
 * and returns that integer.
 *****************************************************************************/
int roundfloat(float value)
{
  float high;
  float low;
  int   ans;

  high = std::ceil(value);
  low  = std::floor(value);

  if ((value - low) < (high - value)) {
    ans = static_cast<int>(low);
  }
  else {
    ans = static_cast<int>(high);
  }

  return ans;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function find_inter() begins:
 *----------------------------------------------------------------------------
 * This function finds the intersection between two lists of integer values,
 * and returns the number of values in the intersection.
 *****************************************************************************/
template <typename INT>
size_t find_inter(const INT set1[],  /* the first set of integers */
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

template size_t find_inter(const int set1[], const int set2[], size_t length1, size_t length2,
                           int inter_ptr[]);
template size_t find_inter(const int64_t set1[], const int64_t set2[], size_t length1,
                           size_t length2, int64_t inter_ptr[]);

#define QSORT_CUTOFF 12

namespace {
  template <typename INT>
  int is_less_than4(INT ra1, INT rb1, INT rc1, INT rd1, INT ra2, INT rb2, INT rc2, INT rd2)
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

  template <typename INT> int is_less_than4v(INT *v1, INT *v2, INT *v3, INT *v4, size_t i, size_t j)
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

  template <typename INT> void swap4(INT *v1, INT *v2, INT *v3, INT *v4, size_t i, size_t j)
  {
    ISWAP(v1, i, j);
    ISWAP(v2, i, j);
    ISWAP(v3, i, j);
    ISWAP(v4, i, j);
  }

  template <typename INT>
  size_t internal_median3_4(INT *v1, INT *v2, INT *v3, INT *v4, size_t left, size_t right)
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
  void internal_qsort_4(INT *v1, INT *v2, INT *v3, INT *v4, size_t left, size_t right)
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

  template <typename INT> void internal_isort_4(INT *v1, INT *v2, INT *v3, INT *v4, size_t N)
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

  template <typename INT> int is_less_than2(INT ra1, INT rb1, INT ra2, INT rb2)
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

  template <typename INT> int is_less_than2v(INT *v1, INT *v2, size_t i, size_t j)
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

  template <typename INT> void swap2(INT *v1, INT *v2, size_t i, size_t j)
  {
    ISWAP(v1, i, j);
    ISWAP(v2, i, j);
  }

  template <typename INT> size_t internal_median3_2(INT *v1, INT *v2, size_t left, size_t right)
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

  template <typename INT> void internal_qsort_2(INT *v1, INT *v2, size_t left, size_t right)
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

  template <typename INT> void internal_isort_2(INT *v1, INT *v2, size_t N)
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
        v1[j] = v1[j - 1];
        v2[j] = v2[j - 1];
      }
      v1[j] = small_val1;
      v2[j] = small_val2;
    }
  }
} // namespace

/*
 * Sort the values in 'v'
 */

template void qsort4(int *v1, int *v2, int *v3, int *v4, size_t N);
template void qsort4(int64_t *v1, int64_t *v2, int64_t *v3, int64_t *v4, size_t N);

template <typename INT> void qsort4(INT *v1, INT *v2, INT *v3, INT *v4, size_t N)
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

template void                qsort2(int *v1, int *v2, size_t N);
template void                qsort2(int64_t *v1, int64_t *v2, size_t N);
template <typename INT> void qsort2(INT *v1, INT *v2, size_t N)
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

template void                sort3(ssize_t count, int ra[], int rb[], int rc[]);
template void                sort3(ssize_t count, int64_t ra[], int64_t rb[], int64_t rc[]);
template <typename INT> void sort3(ssize_t count, INT ra[], INT rb[], INT rc[])
{
  if (count <= 1) {
    return;
  }
  /* heapify */
  for (ssize_t start = (count - 2) / 2; start >= 0; start--) {
    siftDowniii(ra, rb, rc, start, count);
  }

  for (size_t end = count - 1; end > 0; end--) {
    SWAP(ra[end], ra[0]);
    SWAP(rb[end], rb[0]);
    SWAP(rc[end], rc[0]);
    siftDowniii(ra, rb, rc, 0, end);
  }
}

template ssize_t                bin_search2(int value, size_t num, int List[]);
template ssize_t                bin_search2(int64_t value, size_t num, int64_t List[]);
template <typename INT> ssize_t bin_search2(INT value, size_t num, INT List[])
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

template void gds_qsort(int v[], size_t N);
template void gds_qsort(int64_t v[], size_t N);

template <typename INT> void gds_qsort(INT v[], size_t N)
{
  if (N <= 1) {
    return;
  }
  gds_qsort(v, 0, N - 1);
  gds_isort(v, N);
}
