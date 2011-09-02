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

#ifndef _ELB_UTIL_CONST_H_
#define _ELB_UTIL_CONST_H_

/* Function prototypes */
extern
void print_usage(void);

extern
int token_compare(
  char *token,		/* The input character string */
  const char *key	/* The key to compare with token */
  );

extern
void strip_string(
  char inp_str[],	/* The string to strip */
  const char *tokens	/* The tokens to strip from the beginning and
                         * end of the input string */
  );

extern
void clean_string(
  char inp_str[],	/* The string to clean */
  const char *tokens	/* The tokens to strip multiple copies of */
  );

extern
void string_to_lower(
  char inp_str[],	/* The string to convert to lower case */
  const char cstop	/* Character where to stop */
  );

extern
void qsort4(int *v1, int *v2, int *v3, int *v4, int N);

extern
void qsort2(int *v1, int *v2, int N);

extern
void sort2_int_int(
  int  count,
  int *array1,
  int *array2
  );

extern
void sort3_int_int_int(
  int  count,
  int *array1,
  int *array2,
  int *array3
  );

extern
void sort4_iiii(
  int  count,
  int *array1,
  int *array2,
  int *array3,
  int *array4
  );

extern
void qsort4(
  int *v1,
  int *v2,
  int *v3,
  int *v4,
  int N
  );

extern
void find_first_last(
  int  value,
  int  vecsize,
  int *vector,
  int *first,
  int *last
  );

extern
int find_int(
  int  value1,
  int  value2,
  int  start,
  int  stop,
  int *vector1,
  int *vector2
  );

extern
int in_list(
  const int  search,		/* The value to search for */
  const int  count,		/* Number of elements in vector to search */
  int       *vector		/* The vector to search */
);

extern
int roundfloat(
  const float value		/* the value to be rounded */
);

extern int find_max (
  const int list_length,
  const int list[]
);

extern int find_min (
  const int list_length,
  const int list[]
);


extern
int find_inter (
  const int set1[],		/* the first set of integers */
  const int set2[],		/* the second set of integers */
  const int length1,		/* the length of the first set */
  const int length2,		/* the length of the second set */
  int inter_ptr[]		/* the values in the intersection */
);

#endif /* _ELB_UTIL_CONST_H_ */
