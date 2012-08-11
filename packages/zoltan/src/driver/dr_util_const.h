/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


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
