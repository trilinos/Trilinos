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
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif
                                /* It would be nice if the routines in this
				   file didn't have include dependencies!   */
#include "pe_common.h"
#include "el_geom_const.h"
#include "rf_fem_const.h"
#include "rf_allo.h"

int bin_search_min(int list[], int num, int value);

/************ R O U T I N E S   I N   T H I S   F I L E  **********************

       NAME            		        TYPE
----------------------------------------------------------------------
	sortN_int_int		        void
        sortN_int_float			void
        sortN_int_floatlist		void
        sortN_int_double		void
        sortN_int_doublelist		void
	sort_int_ptr		        void
	sort_int_int_ptr	        void
	sort_int_int_int	        void
	sort_int			void
        find_max			void
	find_min			int
	find_inter			int
	find_inter_pos			int
	exchange_pointers		void
        in_list_mono                    int
        iindexx				void
        sort3				void
        bin_search2                     int
        find_range                      int
        bin_search_min                  int
 	print_line			void
	int_cmp                         int
        break_message_up                int

******************************************************************************/

/******************* PROTOTYPES FOR STATIC FUNCTIONS *************************/

static void ICOPY (int *i1, int *i2, int m);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifdef __STDC__

void sortN_int_int(int nval, int ra[], int narr, ...)

#else

void
sortN_int_int (va_alist)
va_dcl

#endif
/*****************************************************************************/

/*
*	Numerical Recipies in C source code
*	modified to have first argument an integer array (JS)
*
*	Sorts the array ra[0,..,(nval-1)] in ascending numerical order using
*	heapsort algorityhm, while making the corresponding rearrangement
*	of narr arrays rb[0,..,(nval-1)], rc[0,..,(nval-1)], ...
*/

{
#ifndef __STDC__
  int     nval, *ra, narr;
#endif
  va_list va;          /* Current pointer in the argument list */

  int   **rb;
  int     l, j, ir, i, k;
  int     rra;
  int    *rrb;



#ifdef __STDC__
  va_start(va, narr);
#else
  va_start(va);
  nval = va_arg(va, int);
  ra   = va_arg(va, int *);
  narr = va_arg(va, int);
#endif

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  rb = (int **) array_alloc(__FILE__, __LINE__, 1, narr, sizeof(int *));
  rrb = (int *) array_alloc(__FILE__, __LINE__, 1, nval, sizeof(int));
  /* get the other arrays */
  for (i = 0; i < narr; i++)
    rb[i] = va_arg(va, int *);
  va_end(va);

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][l];
    } else {
      rra=ra[ir];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][ir];
      ra[ir]=ra[0];
      for (k=0; k<narr; k++)
        rb[k][ir]=rb[k][0];
      if (--ir == 0) {
        ra[0]=rra;
        for (k=0; k<narr; k++)
          rb[k][0]=rrb[k];
        safe_free((void **) &rb);
        safe_free((void **) &rrb);
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        for (k=0; k<narr; k++)
          rb[k][i]=rb[k][j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    for (k=0; k<narr; k++)
      rb[k][i]=rrb[k];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifdef __STDC__

void sortN_int_float(int nval, int ra[], int narr, ...)

#else

void
sortN_int_float (va_alist)
va_dcl

#endif
/*****************************************************************************/

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*
*       Sorts the array ra[0,..,(nval-1)] in ascending numerical order using
*       heapsort algorityhm, while making the corresponding rearrangement
*       of narr arrays rb[0,..,(nval-1)], rc[0,..,(nval-1)].
*/

{
#ifndef __STDC__
  int     nval, *ra, narr;
#endif
  va_list va;          /* Current pointer in the argument list */

  float **rb;
  int     l, j, ir, i, k;
  int     rra;
  float  *rrb;



#ifdef __STDC__
  va_start(va, narr);
#else
  va_start(va);
  nval = va_arg(va, int);
  ra   = va_arg(va, int *);
  narr = va_arg(va, int);
#endif

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  rb = (float **) array_alloc(__FILE__, __LINE__, 1, narr, sizeof(float *));
  rrb = (float *) array_alloc(__FILE__, __LINE__, 1, nval, sizeof(float));
  /* get the other arrays */
  for (i = 0; i < narr; i++)
    rb[i] = va_arg(va, float *);
  va_end(va);

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][l];
    } else {
      rra=ra[ir];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][ir];
      ra[ir]=ra[0];
      for (k=0; k<narr; k++)
        rb[k][ir]=rb[k][0];
      if (--ir == 0) {
        ra[0]=rra;
        for (k=0; k<narr; k++)
          rb[k][0]=rrb[k];
        safe_free((void **) &rb);
        safe_free((void **) &rrb);
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        for (k=0; k<narr; k++)
          rb[k][i]=rb[k][j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    for (k=0; k<narr; k++)
      rb[k][i]=rrb[k];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sortN_int_floatlist(int nval, int ra[], int narr, float *rb[])

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*
*       Sorts the array ra[0,..,(nval-1)] in ascending numerical order using
*       heapsort algorityhm, while making the corresponding rearrangement
*       of narr arrays rb[0][0,..,(nval-1)], rb[1][0,..,(nval-1)], ...
*/

{
  int     l, j, ir, i, k;
  int     rra;
  float  *rrb;

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  /* allocate temp space to hold other values */
  rrb = (float *) array_alloc(__FILE__, __LINE__, 1, nval, sizeof(float));

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][l];
    } else {
      rra=ra[ir];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][ir];
      ra[ir]=ra[0];
      for (k=0; k<narr; k++)
        rb[k][ir]=rb[k][0];
      if (--ir == 0) {
        ra[0]=rra;
        for (k=0; k<narr; k++)
          rb[k][0]=rrb[k];
        safe_free((void **) &rrb);
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        for (k=0; k<narr; k++)
          rb[k][i]=rb[k][j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    for (k=0; k<narr; k++)
      rb[k][i]=rrb[k];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifdef __STDC__

void sortN_int_double(int nval, int ra[], int narr, ...)

#else

void
sortN_int_double (va_alist)
va_dcl

#endif
/*****************************************************************************/

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*
*       Sorts the array ra[0,..,nval] in ascending numerical order using
*       heapsort algorityhm, while making the corresponding rearrangement
*       of narr arrays rb[0,..,(nval-1)], rc[0,..,(nval-1)].
*/

{
#ifndef __STDC__
  int      nval, *ra, narr;
#endif
  va_list  va;          /* Current pointer in the argument list */

  double **rb;
  int      l, j, ir, i, k;
  int      rra;
  double  *rrb;



#ifdef __STDC__
  va_start(va, narr);
#else
  va_start(va);
  nval = va_arg(va, int);
  ra   = va_arg(va, int *);
  narr = va_arg(va, int);
#endif

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  rb = (double **) array_alloc(__FILE__, __LINE__, 1, narr, sizeof(double *));
  rrb = (double *) array_alloc(__FILE__, __LINE__, 1, nval, sizeof(double));
  /* get the other arrays */
  for (i = 0; i < narr; i++)
    rb[i] = va_arg(va, double *);
  va_end(va);

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][l];
    } else {
      rra=ra[ir];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][ir];
      ra[ir]=ra[0];
      for (k=0; k<narr; k++)
        rb[k][ir]=rb[k][0];
      if (--ir == 0) {
        ra[0]=rra;
        for (k=0; k<narr; k++)
          rb[k][0]=rrb[k];
        safe_free((void **) &rb);
        safe_free((void **) &rrb);
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        for (k=0; k<narr; k++)
          rb[k][i]=rb[k][j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    for (k=0; k<narr; k++)
      rb[k][i]=rrb[k];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sortN_int_doublelist(int nval, int ra[], int narr, double *rb[])

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*
*       Sorts the array ra[0,..,(nval-1)] in ascending numerical order using
*       heapsort algorityhm, while making the corresponding rearrangement
*       of narr arrays rb[0][0,..,(nval-1)], rb[0][0,..,(nval-1)], ...
*/

{
  int      l, j, ir, i, k;
  int      rra;
  double  *rrb;

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  rrb = (double *) array_alloc(__FILE__, __LINE__, 1, nval, sizeof(double));

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][l];
    } else {
      rra=ra[ir];
      for (k=0; k<narr; k++)
        rrb[k]=rb[k][ir];
      ra[ir]=ra[0];
      for (k=0; k<narr; k++)
        rb[k][ir]=rb[k][0];
      if (--ir == 0) {
        ra[0]=rra;
        for (k=0; k<narr; k++)
          rb[k][0]=rrb[k];
        safe_free((void **) &rrb);
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        for (k=0; k<narr; k++)
          rb[k][i]=rb[k][j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    for (k=0; k<narr; k++)
      rb[k][i]=rrb[k];
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sort_int_ptr(int nval, int ra[], char *rb[])
/*****************************************************************************/

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*
*       Sorts the array ra[0,..,(nval-1)] in ascending numerical order using
*       heapsort algorityhm, while making the corresponding rearrangement
*       of array rb[0,..,(nval-1)].
*/

{
  int     l, j, ir, i;
  int     rra;
  char   *rrb;

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      rrb=rb[l];
    } else {
      rra=ra[ir];
      rrb=rb[ir];
      ra[ir]=ra[0];
      rb[ir]=rb[0];
      if (--ir == 0) {
        ra[0]=rra;
        rb[0]=rrb;
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        rb[i]=rb[j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
  }
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sort_int_int_ptr(int nval, int ra[], int rb[], char *rc[])
/*****************************************************************************/

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*
*       Sorts the array ra[0,..,(nval-1)] in ascending numerical order using
*       heapsort algorityhm, while making the corresponding rearrangement
*       of arrays rb[0,..,(nval-1)] & rc[0,..,(nval-1)].
*/

{
  int     l, j, ir, i;
  int     rra, rrb;
  char   *rrc;

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      rrb=rb[l];
      rrc=rc[l];
    } else {
      rra=ra[ir];
      rrb=rb[ir];
      rrc=rc[ir];
      ra[ir]=ra[0];
      rb[ir]=rb[0];
      rc[ir]=rc[0];
      if (--ir == 0) {
        ra[0]=rra;
        rb[0]=rrb;
        rc[0]=rrc;
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        rb[i]=rb[j];
        rc[i]=rc[j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
    rc[i]=rrc;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sort_int_int_int(int nval, int ra[], int rb[], int rc[])
/*****************************************************************************/

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*
*       Sorts the array ra[0,..,(nval-1)] in ascending numerical order using
*       heapsort algorityhm, while making the corresponding rearrangement
*       of arrays rb[0,..,(nval-1)] & rc[0,..,(nval-1)].
*/

{
  int     l, j, ir, i;
  int     rra, rrb, rrc;

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      rrb=rb[l];
      rrc=rc[l];
    } else {
      rra=ra[ir];
      rrb=rb[ir];
      rrc=rc[ir];
      ra[ir]=ra[0];
      rb[ir]=rb[0];
      rc[ir]=rc[0];
      if (--ir == 0) {
        ra[0]=rra;
        rb[0]=rrb;
        rc[0]=rrc;
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        rb[i]=rb[j];
        rc[i]=rc[j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
    rc[i]=rrc;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sort_int_int(int nval, int ra[], int rb[])
/*****************************************************************************/

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array (JS)
*
*       Sorts the array ra[0,..,(nval-1)] in ascending numerical order using
*       heapsort algorityhm, while making the corresponding rearrangement
*       of array rb[0,..,(nval-1)] 
*/

{
  int     l, j, ir, i;
  int     rra, rrb;

  /*
   *  No need to sort if one or fewer items.
   */
  if (nval <= 1) return;

  l=nval >> 1;
  ir=nval-1;
  for (;;) {
    if (l > 0) {
      rra=ra[--l];
      rrb=rb[l];
    } else {
      rra=ra[ir];
      rrb=rb[ir];
      ra[ir]=ra[0];
      rb[ir]=rb[0];
      if (--ir == 0) {
        ra[0]=rra;
        rb[0]=rrb;
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        rb[i]=rb[j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sort_int(int n, int ra[])

/*
*	Numerical Recipies in C source code
*	modified to have first argument an integer array
*
*	Sorts the array ra[0,..,(n-1)] in ascending numerical order using
*	heapsort algorithm.
*
*/

{
  int   l, j, ir, i;
  int   rra;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  l=n >> 1;
  ir=n-1;
  for (;;) {
    if (l > 0)
      rra=ra[--l];
    else {
      rra=ra[ir];
      ra[ir]=ra[0];
      if (--ir == 0) {
	ra[0]=rra;
	return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	j += (i=j)+1;
      }
      else j=ir+1;
    }
    ra[i]=rra;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int find_max(int list_length, int *list)

     /*
       Function which finds the largest integer from a vector of integers

       Author:          Scott Hutchinson (1421)
       Date:            8 December 1992

       */

/*****************************************************************************/

{

  /* local variables */

  register int i, max;

 /*************************** execution begins *******************************/

  if (list_length > 0) {
    max = list[0];
    for (i = 1; i < list_length; i++)
      if (list[i] > max) max = list[i];
    return (max);
  } else
    return (INT_MIN);

} /* find_max */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int find_min(int list_length, int list[])

/*
 *     Function which finds the smallest integer from a vector of integers
 *
 *      Author:          Scott Hutchinson (1421)
 *      Date:            8 December 1992
 *
 */

{

  /* local variables */

  register int i, min;

  /************************** execution begins *******************************/
  if (list_length > 0) {
    min = list[0];
    for (i = 1; i < list_length; i++)
      if (list[i] < min)  min = list[i];
    return min;
  } else
      return INT_MAX;

} /* find_min ****************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int find_inter(int inter_ptr[], int set1[], int set2[], int length1,
	       int length2, int prob_type)

/*
 *
 *      Function which finds the intersection of two integer lists.
 *      The points in set1 that belong in the intersection set are
 *      returned in the vector inter_pts, starting at position inter_pts[0].
 *      Enough space in inter_pts[] (min(length1, length2)) must
 *      have already been allocated in the calling program before this
 *      function is called.
 *
 *        prob_type defines the problem to be addressed:
 *            0 = don't know anything about set1 or set2.
 *            1 = Know that set2 is monotonically increasing.
 *            2 = Know that set1 and set2 are monotonically increasing
 * NOTE: Only prob_type == 2 supported now...
 *
 *      On return, find_inter returns 0 if there is no intersection.
 *      It returns the number of points in the intersection, if there
 *      is an intersection.
 */

{
  int i, j, counter = 0;
  int max_set1, min_set1, max_set2, min_set2;

  /* Error check the arguments */
  if ((length1 <= 0) || (length2 <= 0)) return (counter);

  /*
   *  Find the maximum and the minimum of the two sets
   */
  max_set1 = set1[length1-1];
  min_set1 = set1[0];
  max_set2 = set2[length2-1];
  min_set2 = set2[0];
  /*
   *    Check for a possible overlaps in node numbers;
   *    If there is an overlap, then do the search using a linearly
   *    scaled method
   *
   */
  if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) )   {
    i = 0;
    j = 0;
    while (i < length1 && j < length2 &&
	   set1[i] <= max_set2 && set2[j] <= max_set1) {
      if (set1[i] < set2[j])
	++i;
      else if (set2[j] < set1[i])
	++j;
      else {
	inter_ptr[counter++] = i;
	++i;
	++j;
      }
    }
  }
  return counter;
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int find_inter_pos(int intersect[], int length1, int set1[], int length2,
                   int set2[], int prob_type)

/*
 *
 *	Function which finds the intersection of two integer lists.
 *	inter_ptr[] must be allocated as at least the same length as
 *	set1[]. The values in inter_ptr[] correspond to the location
 *	in set2[] of the value in the same position in set2[]. If the
 *	value in set1[] is not in set2[], then inter_ptr is set to -1
 *	for that position.
 *
 *	    prob_type defines the problem to be addressed:
 *		0 = don't know anything about set1 or set2.
 *		1 = Know that set2 is monotonically increasing.
 *		2 = Know that set1 and set2 are monotonically increasing
 *
 *      On return, find_inter_pos returns 0 if there is no intersection.
 *      It returns the number of points in the intersection, if there
 *      is an intersection.
 */

{

  /* Local variables */

  register int    i, j, count = 0;
  int             max_set1, min_set1, max_set2, min_set2;

/****************************** execution begins *****************************/

  /* Error check the arguments */
  if ((length1 <= 0) || (length2 <= 0)) return (count);

  /* Initialize the intersect vector */
  for(i=0; i < length1; i++)
    intersect[i] = -1;

  if (prob_type == 0 ) {

    /* find the maximum and the minimum of the two sets */
    max_set1 = find_max (length1, set1);
    min_set1 = find_min (length1, set1);
    max_set2 = find_max (length2, set2);
    min_set2 = find_min (length2, set2);

    /*  check for a possible overlaps in node numbers;
     *  If there is an overlap, then do the search
     */

    if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) )   {

      for (i = 0; i < length1; i++)
        for (j = 0; j < length2; j++)
          if (set1[i] == set2[j]) {
            intersect[i] = j;
            count++;
          }

    }
  } else if (prob_type == 1) {

    fprintf (stderr, "prob_type = 1 is unimplemented\n");
    exit(1);

  } else if (prob_type == 2) {

    /* Check for the possibility of an intersection */
    min_set1 = set1[0];
    max_set1 = set1[length1-1];
    min_set2 = set2[0];
    max_set2 = set2[length2-1];

    /* Search through the set */

    if( (max_set2 >= min_set1) && (min_set2 <= max_set1) ) {
      for(i=0, j=0; i < length1; i++) {
        while( (j < (length2-1)) && (set1[i] > set2[j]) ) j++;
        if (set1[i] == set2[j]) {
          intersect[i] = j;
          count++;
        }
      }
    }
  } else {

    fprintf (stderr, "prob_type = %d is unknown\n", prob_type);
    exit(1);

  }

  return count;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void exchange_pointers(void **pointer1, void **pointer2)

/*
 *	This function can be used to exchange the addresses that two pointers
 *       point to in the calling routine.  It should be called in the
 *       following way:
 *
 *	    double *pt1, *pt2;
 *	    exchange_pointers ( (void **) &pt1,  (void **) &pt2 );
 *
 *
 */

{
  void *temp_pointer;

  temp_pointer = *pointer1;
  *pointer1    = *pointer2;
  *pointer2    = temp_pointer;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int in_list_mono(int ivalue, int *ibegin, int iend, int ivector[])

/*
*	 This function searches an integer vector, ivector(ibegin:iend-1),
*     for the presence of a number, ivalue.  It returns the index of the
*     value, or -1, if the number, ivalue, is not found in the list.
*     ivector is assumed to be monotonic.  Therefore, the search is
*     stopped when ivector[i] > ibegin.  It is also assumed that this
*     function will be called repeatibly with monotonically increasing
*     values of ivalue.  Therefore, *ibegin is incremented at the end of
*     the routine to the last value in ivector that was checked.
*/

{
  register int i = *ibegin;

  if ((ivalue >= ivector[i]) && (ivalue <= ivector[iend-1])) {
    while (ivalue > ivector[i]) {
      if (i < iend-1) i++;
    }
    *ibegin = i;
    if (ivalue == ivector[i])
      return i;
    else
      return -1;
  } else
    return -1;

}

/*****************************************************************************/
/*****************************************************************************/

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void iindexx (unsigned int n, int arr[], unsigned int indx[])
/*
 *  Numerical Recipes routine to create an index vector for
 *  sorting
 *  - works on arr[1], ..., arr[n]
 *  - returns values in indx[1], ..., indx[n]
 *  - Values returned in indx[] range from 1 to n.
 *
 */
{
	unsigned int i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0, *istack;
	int a;

	istack= (int *) array_alloc (__FILE__, __LINE__, 1, (1 + NSTACK),
                                     sizeof(int));
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) {
                          fprintf(stderr, "NSTACK too small in iindexx.\n");
                          exit (1);
                        }
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	safe_free ((void **) &istack);
}
#undef M
#undef NSTACK
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software #,V4+!5,. */

/*****************************************************************************/
/*****************************************************************************/

void sort3(int n, int ra[], int rb[], int m)

/*
*	Numerical Recipies in C source code
*	modified to have first argument an integer array (JS)
*
*	Sorts the array ra[1,..,n] in ascending numerical order using quicksort
*	algorithm, while making the corresponding rearrangement of the
*	matrix rb[1,..,n*m].
*       rb is structured so that the first m entries correspond to ra[1].
*       These entries are not sorted amongst themselves.
*
*         ra[1]    -    rb[1],           rb[2],         ..., rb[m]
*         ra[2]    -    rb[m*(2-1)+1],   rb[m*(2-1)+1], ..., rb[m*(2-1)+m]
*         ...                              ...
*         ra[n]    -    rb[m*(n-1)+1],   rb[m*(n-1)+1], ..., rb[m*(n-1)+m]
*
*/

{
  extern void   iindexx (unsigned int n, int arr[],
			 unsigned int indx[]);
  unsigned int           j, *iwksp;
  int                   *wksp;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  iwksp = (unsigned int *) array_alloc (__FILE__, __LINE__, 1, n+1,
                                        sizeof(unsigned int));
  wksp  = (int *)           array_alloc (__FILE__, __LINE__, 1, (n*m),
                                         sizeof(int));
  iindexx (n, (int *) &ra[-1], iwksp-1);

  for (j=0; j<n; j++) wksp[j] = ra[j];
  for (j=0; j<n; j++) ra[j]   = wksp[iwksp[j]-1];

  for (j=0; j<(n*m); j++) wksp[j]=rb[j];

  switch (m) {
  case 2:
    for (j=0; j<n; j++) {
      rb[j*2  ] = wksp[(iwksp[j])*2 - 2];
      rb[j*2+1] = wksp[(iwksp[j])*2 - 1];
    }
    break;
  case 4:
    for (j=0; j<n; j++) {
      rb[j*4  ] = wksp[(iwksp[j])*4 - 4];
      rb[j*4+1] = wksp[(iwksp[j])*4 - 3];
      rb[j*4+2] = wksp[(iwksp[j])*4 - 2];
      rb[j*4+3] = wksp[(iwksp[j])*4 - 1];
    }
    break;
  default:
    for (j=0; j<n; j++) ICOPY (&rb[j*m], &wksp[(iwksp[j]-1)*m], m);
  }
  safe_free ((void **) &iwksp);
  safe_free ((void **)  &wksp);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void ICOPY (int *i1, int *i2, int m)

{
  int i; for (i = 0; i < m; i++) *i1++ = *i2++;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int bin_search2 (int value, int num, int List[])

/*
 * Searches a monotonic list of values for the value, value.
 * It returns the index of the first position found, which matches value.
 * The list is assumed to be monotonic, and
 * consist of elements list[0], ..., list[n-1].
 * If no position in list matches value, it returns the value -1.
 *
 */

{

 register int top, bottom = 0, middle, g_mid;

 /***** execution begins *****/

 top = num - 1;
 while (bottom <= top) {
   middle = (bottom + top) >> 1;
   g_mid = List[middle];
   if (value < g_mid)
     top = middle - 1;
   else if (value > g_mid)
     bottom = middle + 1;
   else
     return middle;     /* found */
 }

 return -1;

} /* bin_search2 */

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

int find_range (int start_value, int end_value, int list[], int num,
                int *start_pos, int *end_pos)

/*
 * find_range:
 *
 *      This routine searches a monotonically increasing list for the range
 * of values which fall between two given values, i.e.,
 *
 *  start_value <= list[i] <= end_value, for all i in the range
 *
 * The number of values which fall in the given range is returned.
 * The lowest value of i such that the above is true is returned in
 * start_pos, while the highest value of i is returned in end_pos
 *   If there is no overlap, find_range is set to 0, and *start_pos and
 * *end_pos will have unspecified values.
 *
 *   Input
 * ----------
 *  start_value = Lower value of the range sought
 *  end_value   = Upper value of the range sought
 *  list        = Vector of integer values to be searched
 *                (must be a monotonically increasing function)
 *                i.e., {list[0], ..., list[num-1]}
 *  num         = Length of the vector, list
 *
 *   Output
 * ----------
 *  (return)    =  Number of values in the vector list, which are in the
 *                 designated range.
 *  start_pos   =  Address of the value of the starting position within list[]
 *                 which is in the designated range.
 *   end_pos    =  Address of the value of the ending position within list[]
 *                 which is in the designated range.
 *
 */

{

  /* Local Variables: */

  int        start, end;

  /* Error check the argument list */

  if (end_value < start_value) return (0);
  if (num <= 0)                return (0);

  /* Check for Obvious Limits -saves doing a lot of special cases below */

  if (start_value > list[num - 1]) return (0);
  if (  end_value < list[0])       return (0);

  /* Find starting position in list */

  start = bin_search_min (list, num, start_value);

  /* Find ending position in list */

  end   = bin_search_min (&list[start], num - start, end_value);
  end  += start;

  if (list[start] < start_value) start++;
  if (start > end) return (0);

  *start_pos = start;
  *end_pos   = end;
  return (end - start + 1);

}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

int bin_search_min(int list[], int num, int value)

/*
 *     Searches a monotonic list of values for the value, value.
 * It returns the index of the first position found, which matches value.
 * The list is assumed to be monotonic, and
 * consist of elements list[0], ..., list[n-1].
 *   If no position in list matches value, it the next lowest position
 * in the list, or 0 if the value is lower than all values in list []
 * Error conditions are specified with a return value of -1.
 *
 *  example:
 * -------------
 *               0  1  2  3
 *    list [] = {3, 5, 7, 9}, num = 4
 *
 *       value       return_value
 *      ---------   --------------
 *          1             0 -- less than first element
 *          3             0 -- matches first element
 *          6             1 -- > list[1] and < list[2]
 *          13            3 -- > list[3]
 * ---------------------------------------------
 */

{
  int top, bottom, diff, middle, g_mid, g_bot, g_top;

  if (num <= 0) return (-1);
  bottom = 0;
  g_bot = list[0];
  if (value <= g_bot)
    return bottom;
  else {
    top    = num - 1;
    g_top  = list[top];
    if (value >=  g_top)
      return top;
    else {
      while ( (diff = top - bottom) > 1) {
        middle = bottom + diff/2;
        g_mid = list[middle];
        if (value > g_mid) {
          bottom = middle;
        } else if (value < g_mid) {
          top    = middle;
        } else {
          return middle;
        }
      }
    }
  }

  return bottom;

}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void print_line (char *charstr, int ntimes)

{
  int i;
  for (i = 0; i < ntimes; i++) (void) printf("%c", *charstr);
  (void) printf("\n");
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

 int int_cmp(int *i1, int *i2)

/*

  Function which compares two integers and returns negative if its first
  argument is less than the second, zero if equal and positive if greater.
  Designed to be used by the stdlib functions bsearch and qsort.

*/

{

  if (*i1 < *i2) return -1;
  else if (*i1 == *i2) return 0;
  else return 1;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int break_message_up (size_t unit_size, size_t num_units, size_t max_bytes,
		      int *start_pos[])

/*
 *   break_message_up:
 *
 *         This routine will break a long message up into chunks.  It returns
 * the number of chunks and the starting and ending positions, in terms of
 * unit numbers, of each chunk.
 * It assumes that the first unit of the message has a positional value of 0.
 *
 * NOTE:
 *   This routine does malloc memory.  This memory must be freed in the
 * calling routine, when it is no longer needed.  Use the following statement
 * to free the memory:
 *
 *       safe_free ((void **) &start_pos);
 *
 * Input
 * -------
 *
 * unit_size      = Size in bytes of each "unit" in the message.
 * num_units      = Number of units in the message.
 * max_bytes      = Max_bytes allowed in each message.
 *
 * Output
 * -------
 *
 * return_value   = Number of messages
 * start_pos []   = Starting positions of each message
 *                      (note: start_pos[0] = 0, by definition).
 *                    It is defined to have a length of num_mesg+1
 *                    start_pos[num_mesg] = num_units
 *
 *
 * Usage
 *----------
 *       To get the length of message i:
 *
 *                   mesg_length = start_pos[i+1] - start_pos[i]
 *
 */

{
  size_t  i, num_units_per_message, remainder, num_mesg;

 /*------------------------- begin execution --------------------------------*/

  if (num_units <= 0) {
    num_mesg = 0;
    *start_pos = NULL;
    return (num_mesg);
  }

  num_units_per_message = max_bytes / unit_size ;
  if (num_units < num_units_per_message) num_units_per_message = num_units;
  num_mesg   = num_units / num_units_per_message;
  remainder  = num_units % num_units_per_message;
  if (remainder != 0) num_mesg++;

  *start_pos = (int *) array_alloc (__FILE__, __LINE__, 1, (num_mesg + 1),
                                    sizeof (int));

  for (i = 0; i < num_mesg; i++) {
    (*start_pos)[i] =  i * num_units_per_message;
  }
  (*start_pos) [num_mesg] = num_units;

  return num_mesg;
}

/*****************************************************************************/
/*                END OF FILE rf_util.c					     */
/*****************************************************************************/
