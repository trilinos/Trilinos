// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef _DR_UTIL_CONST_H_
#define _DR_UTIL_CONST_H_

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Function prototypes */
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
void string_to_lower(
  char inp_str[],	/* The string to convert to lower case */
  const char cstop	/* Character where to stop */
);

extern
void clean_string(
  char inp_str[],	/* The string to clean */
  const char *tokens	/* The tokens to strip multiple copies of */
);

extern
int in_list(
  const ZOLTAN_ID_TYPE  search,	/* The value to search for */
  const int  count,	/* Number of elements in vector to search */
  ZOLTAN_ID_TYPE       *vector	/* The vector to search */
);

extern
int in_list2(
  const int  search,	/* The value to search for */
  const int  count,	/* Number of elements in vector to search */
  int       *vector	/* The vector to search */
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
  const int set1[],             /* the first set of integers */
  const int set2[],             /* the second set of integers */
  const int length1,            /* the length of the first set */
  const int length2,            /* the length of the second set */
  const int prob_type,          /* value indicating known info about lists */
  int inter_ptr[]               /* the values in the intersection */
);

extern void quicksort_pointer_inc_id_id(
  int *sorted,
  ZOLTAN_ID_TYPE * val1,
  ZOLTAN_ID_TYPE  *val2,
  int start,
  int end
);

extern void safe_free(void **ptr);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* _DR_UTIL_CONST_H_ */
