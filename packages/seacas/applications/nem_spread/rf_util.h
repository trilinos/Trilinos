/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef _RF_UTIL_H
#define _RF_UTIL_H
#include <climits>

template <typename T> inline void my_swap(T &r, T &s)
{
  T tmp = r;
  r     = s;
  s     = tmp;
}

template <typename T, typename U> void siftDown(T *a, U *b, int64_t start, int64_t end)
{
  int64_t root = start;

  while (root * 2 + 1 < end) {
    int64_t child = 2 * root + 1;
    if ((child + 1 < end) && (a[child] < a[child + 1])) {
      child += 1;
    }
    if (a[root] < a[child]) {
      my_swap(a[child], a[root]);
      my_swap(b[child], b[root]);
      root = child;
    }
    else {
      return;
    }
  }
}

template <typename T, typename U> void my_sort(int64_t count, T ra[], U rb[])
{
  int64_t start, end;

  /* heapify */
  for (start = (count - 2) / 2; start >= 0; start--) {
    siftDown(ra, rb, start, count);
  }

  for (end = count - 1; end > 0; end--) {
    my_swap(ra[end], ra[0]);
    my_swap(rb[end], rb[0]);
    siftDown(ra, rb, 0, end);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
extern void print_line(const char *charstr, int ntimes);

#endif /* #ifndef _RF_UTIL_H */
