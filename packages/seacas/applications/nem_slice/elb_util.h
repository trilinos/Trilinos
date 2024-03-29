/*
 * Copyright(C) 1999-2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#pragma once

#include <cstddef> // for size_t
#include <cstdint>
#include <vector>

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#if !defined(__MINGW32__)
#define strcasecmp  stricmp
#define strncasecmp strnicmp
#endif
#endif

/* Function prototypes */
extern int token_compare(char       *token, /* The input character string */
                         const char *key    /* The key to compare with token */
);

extern void strip_string(char        inp_str[], /* The string to strip */
                         const char *tokens     /* The tokens to strip from the beginning and
                                                 * end of the input string */
);

extern void clean_string(char        inp_str[], /* The string to clean */
                         const char *tokens     /* The tokens to strip multiple copies of */
);

extern void string_to_lower(char in_string[], /* The string to convert to lower case */
                            char cval         /* Character where to stop */
);

template <typename INT> void gds_qsort(INT v[], size_t N);

template <typename INT> void qsort4(INT *v1, INT *v2, INT *v3, INT *v4, size_t N);

template <typename INT> void qsort2(INT *v1, INT *v2, size_t N);

template <typename INT> inline void SWAP(INT &r, INT &s)
{
  INT t = r;
  r     = s;
  s     = t;
}

template <typename INT> void siftDown(INT *a, INT *b, size_t start, size_t end)
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

template <typename INT> void sort2(int64_t count, INT ra[], INT rb[])
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

template <typename INT> void sort3(int64_t count, INT ra[], INT rb[], INT rc[]);

template <typename INT>
void find_first_last(INT val, size_t vecsize, INT *vector, INT *first, INT *last);

template <typename INT>
int64_t find_int(INT value1, INT value2, size_t start, size_t stop, INT *vector1, INT *vector2);

template <typename INT, typename INT2> int64_t in_list(INT value, size_t count, const INT2 *vector);

template <typename INT, typename INT2> int64_t in_list(INT value, const std::vector<INT2> &vector);

extern int roundfloat(float value /* the value to be rounded */
);

template <typename INT>
size_t find_inter(const INT set1[],     /* the first set of integers */
                  const INT set2[],     /* the second set of integers */
                  size_t    length1,    /* the length of the first set */
                  size_t    length2,    /* the length of the second set */
                  INT       inter_ptr[] /* the values in the intersection */
);

template <typename INT> int64_t bin_search2(INT value, size_t num, INT List[]);
