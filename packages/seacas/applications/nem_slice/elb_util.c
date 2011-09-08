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

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *	token_compare()
 *	strip_string()
 *	string_to_lower()
 *	clean_string()
 *	sort2_int_int()
 *	sort3_int_int_int()
 *	sort4_iiii()
 *	find_first_last()
 *	find_int()
 *	in_list()
 *	roundfloat()
 *      find_max()
 *      find_min()
 *      find_inter()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>

#include "elb_const.h"
#include "elb_util_const.h"

int is_less_than4(int ra1,int rb1,int rc1,int rd1,
		  int ra2,int rb2,int rc2,int rd2);
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int token_compare(char *token, const char *key)
{

  int i1, key_len, kcnt=0;

  key_len = strlen(key);

  for(i1=0; i1 < strlen(token); i1++)
  {
    if(isupper(token[i1]))
      token[i1] = tolower(token[i1]);

    if(token[i1] != ' ')
    {
      if(token[i1] == key[kcnt])
      {
        kcnt++;
        if(kcnt > key_len)
          return 0;
      }
      else
        return 0;
    }
    if(key[kcnt] == ' ')
      kcnt++;
  }

  if(kcnt == strlen(key))
    return 1;
  else
    return 0;

} /*--------------End token_compare()-----------*/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void strip_string(char inp_str[], const char *tokens)
{
  int  i, j, itok, ntokes, bval;

  i = 0;
  ntokes = strlen(tokens);

  while(inp_str[i] != '\0')
  {
    bval = 0;
    for(itok=0; itok < ntokes; itok++)
    {
      if(inp_str[i] == tokens[itok])
      {
        i++;
        bval = 1;
        break; /* out of for loop */
      }
    }
    if(bval == 0)
      break; /* out of while loop */
  }

  /* Move real part of string to the front */
  j = 0;
  while(inp_str[j+i] != '\0')
  {
    inp_str[j] = inp_str[j+i];
    j++;
  }
  inp_str[j] = inp_str[j+i];
  j--;

  /* Remove trailing tokens */
  while(j != -1)
  {
    bval = 0;
    for(itok=0; itok < ntokes; itok++)
    {
      if(inp_str[j] == tokens[itok])
      {
        bval = 1;
        j--;
        break; /* out of for loop */
      }
    }
    if(bval == 0)
      break; /* out of while loop */
  }

  inp_str[j+1] = '\0';

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void string_to_lower(char in_string[], const char cval)
{
  int len, cnt;

  len = strlen(in_string);
  for(cnt=0; cnt < len; cnt++)
  {
    if(in_string[cnt] == cval)
      return;

    if(isupper(in_string[cnt]))
      in_string[cnt] = tolower(in_string[cnt]);
  }

  return;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void clean_string(char inp_str[], const char *tokens)
{
  int  i, j, itok, ntokes, bval, inplen;

  ntokes = strlen(tokens);
  inplen = strlen(inp_str);

  i = 0;
  bval = 0;
  while(inp_str[i] != '\0')
  {
    for(itok=0; itok < ntokes; itok++)
    {
      if (i < 0) i = 0;
      if(inp_str[i] == tokens[itok])
      {
        /* Find out if the next character is also a token */
        for(j=0; j < ntokes; j++)
        {
          if(inp_str[i+1] == tokens[j])
          {
            bval = 1;
            break;
          }
        }

        if(bval == 1)
        {
          for(j=i+1; j < inplen; j++)
            inp_str[j] = inp_str[j+1];

          inplen--;
          bval = 0;
          i--;
        }
      }
    }

    i++;

  } /* End "while(inp_str[i] != '\0')" */

  return;

} /*---------------- End clean_string() -----------------*/

/*
 * Function from "Numerical Methods in C"
 */
void sort2_int_int(int n, int *ra, int *rb)
{
  /* NOTE: This currently assumes ra[1] is first entry in array */
  int l,j,ir,i;
  int rra, rrb;

  if(n <= 1)return;

  l  = (n >> 1)+1;
  ir = n;
  for(;;)
  {
    if(l > 1)
    {
      rra = ra[--l];
      rrb = rb[l];
    }
    else
    {
      rra    = ra[ir];
      rrb    = rb[ir];
      ra[ir] = ra[1];
      rb[ir] = rb[1];

      if(--ir == 1)
      {
        ra[1] = rra;
        rb[1] = rrb;
        return;
      }
    }

    i = l;
    j = l << 1;

    while(j <= ir)
    {
      if(j < ir && ra[j] < ra[j+1]) ++j;
      if(rra < ra[j])
      {
        ra[i]  = ra[j];
        rb[i]  = rb[j];
        j     += (i=j);
      }
      else
        j = ir + 1;
    }
    ra[i] = rra;
    rb[i] = rrb;
  }

}

/*
 * Function from "Numerical Methods in C"
 * Note that this sort uses 2nd arg as primary key and 3rd arg as secondary key
 */
void sort3_int_int_int(int n, int *ra, int *rb, int *rc)
{
  int l,j,ir,i;
  int rra, rrb, rrc;

  if(n <= 1)return;

  l  = n >> 1;
  ir = n-1;
  for(;;)
  {
    if(l > 0)
    {
      rra = ra[--l];
      rrb = rb[l];
      rrc = rc[l];
    }
    else
    {
      rra    = ra[ir];
      rrb    = rb[ir];
      rrc    = rc[ir];
      ra[ir] = ra[0];
      rb[ir] = rb[0];
      rc[ir] = rc[0];

      if(--ir == 0)
      {
        ra[0] = rra;
        rb[0] = rrb;
        rc[0] = rrc;
        return;
      }
    }

    i = l;
    j = (l << 1)+1;

    while(j <= ir)
    {
      if(j < ir &&
	 (ra[j] < ra[j+1] || (ra[j]==ra[j+1] && rb[j] < rb[j+1]))) ++j;

      if(rra < ra[j] || (rra == ra[j] && rrb < rb[j]))
      {
        ra[i]  = ra[j];
        rb[i]  = rb[j];
        rc[i]  = rc[j];
        j     += (i=j)+1;
      }
      else
        j = ir + 1;
    }
    ra[i] = rra;
    rb[i] = rrb;
    rc[i] = rrc;
  }

}

/*****************************************************************************/
/* Function from "Numerical Methods in C"
 *****************************************************************************/
#if 1
void sort4_iiii(int n, int *ra, int *rb, int *rc, int *rd)
{
  int l,j,ir,i;
  int rra, rrb, rrc, rrd;

  if (n <= 1) return;

  l  = n >> 1;
  ir = n-1;
  for(;;)
    {
      if(l > 0)
	{
	  rra = ra[--l];
	  rrb = rb[l];
	  rrc = rc[l];
	  rrd = rd[l];
	}
      else
	{
	  rra    = ra[ir];
	  rrb    = rb[ir];
	  rrc    = rc[ir];
	  rrd    = rd[ir];
	  ra[ir] = ra[0];
	  rb[ir] = rb[0];
	  rc[ir] = rc[0];
	  rd[ir] = rd[0];

	  if(--ir == 0)
	    {
	      ra[0] = rra;
	      rb[0] = rrb;
	      rc[0] = rrc;
	      rd[0] = rrd;
	      return;
	    }
	}

      i = l;
      j = (l << 1) + 1;

      while (j <= ir)
	{
	  if (j < ir && is_less_than4(ra[j],   rb[j],   rc[j],   rd[j],
				     ra[j+1], rb[j+1], rc[j+1], rd[j+1]))
	    ++j;

	  if (is_less_than4(rra,   rrb,   rrc,   rrd,
			   ra[j], rb[j], rc[j], rd[j])) {
	    ra[i]  = ra[j];
	    rb[i]  = rb[j];
	    rc[i]  = rc[j];
	    rd[i]  = rd[j];
	    j     += (i=j)+1;
	  } else {
	    j = ir + 1;
	  }
	}
      ra[i] = rra;
      rb[i] = rrb;
      rc[i] = rrc;
      rd[i] = rrd;
    }
}
#else
void sort4_iiii(int n, int *ra, int *rb, int *rc, int *rd)
{
  int l,j,ir,i;
  int rra, rrb, rrc, rrd;

  if(n <= 1)return;

  l  = (n >> 1)+1;
  ir = n;
  for(;;)
  {
    if(l > 1)
    {
      rra = ra[--l];
      rrb = rb[l];
      rrc = rc[l];
      rrd = rd[l];
    }
    else
    {
      rra    = ra[ir];
      rrb    = rb[ir];
      rrc    = rc[ir];
      rrd    = rd[ir];
      ra[ir] = ra[1];
      rb[ir] = rb[1];
      rc[ir] = rc[1];
      rd[ir] = rd[1];

      if(--ir == 1)
      {
        ra[1] = rra;
        rb[1] = rrb;
        rc[1] = rrc;
        rd[1] = rrd;
        return;
      }
    }

    i = l;
    j = l << 1;

    while(j <= ir)
    {
      if(j < ir && ra[j] < ra[j+1]) ++j;
      if(rra < ra[j])
      {
        ra[i]  = ra[j];
        rb[i]  = rb[j];
        rc[i]  = rc[j];
        rd[i]  = rd[j];
        j     += (i=j);
      }
      else
        j = ir + 1;
    }
    ra[i] = rra;
    rb[i] = rrb;
    rc[i] = rrc;
    rd[i] = rrd;
  }

}
#endif

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

void assert_sorted(int *vector, int vecsize)
{
  int i;
  for (i=1; i < vecsize; i++)
    assert(vector[i-1] <= vector[i]);
}

/*****************************************************************************
 * Function to find the first and last entries of a given value that are
 * consecutively present in an integer array.
 * ASSUMES that 'vector' is sorted....
 *****************************************************************************/
void find_first_last(int val, int vecsize, int *vector, int *first, int *last)
{
  int i;

  /* assert_sorted(vector, vecsize); */
  
  *first = -1;
  *last  = -1;

  /* See if value is in the vector */
  i = bin_search2(val, vecsize, vector);
  *first = i; /* Save this location */
  
  if (i != -1) {
    /* Value is in vector, find first occurance */
    while (i >=0 && vector[i] == val) {
      i--;
    }
    i++;

    *last = *first; /* Use saved location */
    *first = i;

    for (i=(*last); i < vecsize; i++) {
      if (vector[i] != val) {
        *last = i-1;
        break;
      }
    }

    if (i == vecsize)
      *last = vecsize - 1;
  }
  return;
}

/*****************************************************************************
 * Find the value in the integer array. Return -1 if not found.
 * New 1/21/97: change this so that it cross references a second
 *              array with a second value
 *****************************************************************************/
int find_int(int value1, int value2, int start, int stop, int *vector1,
             int *vector2)
{
  int i;

  for(i=start; i <= stop; i++)
  {
    if((*(vector1+i) == value1) && (*(vector2+i) == value2))
      return i;
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
int in_list(const int value, const int count, int *vector)
{
  int i;

  for(i=0; i < count; i++)
  {
    if(*vector == value)
      return i;

    vector++;
  }

  return -1;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function roundfloat() begins:
 *----------------------------------------------------------------------------
 * This function rounds off the float "value" to the nearest integer,
 * and returns that interger.
 *****************************************************************************/
int roundfloat(const float value)
{
  float high, low;
  int ans;

  high = (float) ceil(value);
  low = (float) floor(value);

  if ((value - low) < (high - value))
    ans = (int) low;
  else
    ans = (int) high;

  return ans;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int find_max(const int list_length, const int list[])

     /*
       Function which finds the largest integer from a vector of integers

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

int find_min(const int list_length, const int list[])

/*
 *     Function which finds the smallest integer from a vector of integers
 *
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
/* Function find_inter() begins:
 *----------------------------------------------------------------------------
 * This function finds the intersection between two lists of interger values,
 * and returns the number of values in the intersection.
 *****************************************************************************/
int find_inter(const int set1[], const int set2[], const int length1,
               const int length2, int inter_ptr[])

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
  int counter = 0;
  int i = 0;
  int j = 0;

  while (i < length1 && j < length2) {
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

#if 0
  fprintf(stderr, "%d %d -- %d %d -- %d %d -- %d\n", length1, length2,
	  min_set1, max_set1, min_set2, max_set2, counter);
#endif
  return counter;

} /* find_inter **************************************************************/

#define QSORT_CUTOFF 12
#define SWAP(V, I,J) do{int _t = V[I]; V[I] = V[J]; V[J] = _t;} while (0)


int is_less_than4(int ra1,int rb1,int rc1,int rd1,
		 int ra2,int rb2,int rc2,int rd2)
{
  if (ra1 < ra2)
    return 1;
  else if (ra1 > ra2)
    return 0;
  assert(ra1 == ra2);

  if (rb1 < rb2)
    return 1;
  else if (rb1 > rb2)
    return 0;
  assert(rb1 == rb2);

  if (rc1 < rc2)
    return 1;
  else 
    return 0;
}
 
int is_less_than4v(int *v1, int *v2, int *v3, int *v4, int i, int j)
{
  if (v1[i] < v1[j])
    return 1;
  else if (v1[i] > v1[j])
    return 0;
  assert(v1[i] == v1[j]);

  if (v2[i] < v2[j])
    return 1;
  else if (v2[i] > v2[j])
    return 0;
  assert(v2[i] == v2[j]);

  if (v3[i] < v3[j])
    return 1;
  else 
    return 0;
}

void swap4(int *v1, int *v2, int *v3, int *v4, int i, int j)
{
  SWAP(v1, i, j);
  SWAP(v2, i, j);
  SWAP(v3, i, j);
  SWAP(v4, i, j);
}

int internal_median3_4(int *v1, int *v2, int *v3, int *v4, int left, int right)
{
  int center;
  center = (left + right) / 2;

  if (is_less_than4v(v1, v2, v3, v4, center, left))
    swap4(v1, v2, v3, v4, left, center);
  if (is_less_than4v(v1, v2, v3, v4, right, left))
    swap4(v1, v2, v3, v4, left, right);
  if (is_less_than4v(v1, v2, v3, v4, right, center))
    swap4(v1, v2, v3, v4, center, right);

  swap4(v1, v2, v3, v4, center, right-1);
  return right-1;
}

void internal_qsort_4(int *v1, int *v2, int *v3, int *v4, int left, int right)
{
  int pivot;
  int i, j;
  
  if (left + QSORT_CUTOFF <= right) {
    pivot = internal_median3_4(v1, v2, v3, v4, left, right);
    i = left;
    j = right - 1;

    for ( ; ; ) {
      while (is_less_than4v(v1, v2, v3, v4, ++i, pivot));
      while (is_less_than4v(v1, v2, v3, v4, pivot, --j));
      if (i < j) {
	swap4(v1, v2, v3, v4, i, j);
      } else {
	break;
      }
    }

    swap4(v1, v2, v3, v4, i, right-1);
    internal_qsort_4(v1, v2, v3, v4, left, i-1);
    internal_qsort_4(v1, v2, v3, v4, i+1, right);
  }
}

void internal_isort_4(int *v1, int *v2, int *v3, int *v4, int N)
{
  int i,j;
  int ndx = 0;
  int small1, small2, small3, small4;
  
  ndx = 0;
  for (i = 1; i < N; i++) {
    if (is_less_than4v(v1, v2, v3, v4, i, ndx)) {
      ndx = i;
    }
  }
  /* Put smallest value in slot 0 */
  swap4(v1, v2, v3, v4, 0, ndx);

  for (i=1; i <N; i++) {
    small1 = v1[i];
    small2 = v2[i];
    small3 = v3[i];
    small4 = v4[i];
    for (j=i; is_less_than4(small1, small2, small3, small4,
			    v1[j-1], v2[j-1], v3[j-1], v4[j-1]); j--) {
      v1[j] = v1[j-1];
      v2[j] = v2[j-1];
      v3[j] = v3[j-1];
      v4[j] = v4[j-1];
    }
    v1[j] = small1;
    v2[j] = small2;
    v3[j] = small3;
    v4[j] = small4;
  }
}

/*
 * Sort the values in 'v' 
 */
   
void qsort4(int *v1, int *v2, int *v3, int *v4, int N)
{
  internal_qsort_4(v1, v2, v3, v4, 0, N-1);
  internal_isort_4(v1, v2, v3, v4, N);

#if defined(DEBUG_QSORT)
  fprintf(stderr, "Checking sort of %d values\n", N+1);
  int i;
  for (i=1; i < N; i++) {
    assert(is_less_than4v(v1, v2, v3, v4, i-1, i));
  }
#endif
}

int is_less_than2(int ra1,int rb1, int ra2,int rb2)
{
  if (ra1 < ra2)
    return 1;
  else if (ra1 > ra2)
    return 0;
  assert(ra1 == ra2);

  if (rb1 < rb2)
    return 1;
  else 
    return 0;
}
 
int is_less_than2v(int *v1, int *v2, int i, int j)
{
  if (v1[i] < v1[j])
    return 1;
  else if (v1[i] > v1[j])
    return 0;
  assert(v1[i] == v1[j]);

  if (v2[i] < v2[j])
    return 1;
  else 
    return 0;
}

void swap2(int *v1, int *v2, int i, int j)
{
  SWAP(v1, i, j);
  SWAP(v2, i, j);
}

int internal_median3_2(int *v1, int *v2, int left, int right)
{
  int center;
  center = (left + right) / 2;

  if (is_less_than2v(v1, v2, center, left))
    swap2(v1, v2, left, center);
  if (is_less_than2v(v1, v2, right, left))
    swap2(v1, v2, left, right);
  if (is_less_than2v(v1, v2, right, center))
    swap2(v1, v2, center, right);

  swap2(v1, v2, center, right-1);
  return right-1;
}

void internal_qsort_2(int *v1, int *v2, int left, int right)
{
  int pivot;
  int i, j;
  
  if (left + QSORT_CUTOFF <= right) {
    pivot = internal_median3_2(v1, v2, left, right);
    i = left;
    j = right - 1;

    for ( ; ; ) {
      while (is_less_than2v(v1, v2, ++i, pivot));
      while (is_less_than2v(v1, v2, pivot, --j));
      if (i < j) {
	swap2(v1, v2, i, j);
      } else {
	break;
      }
    }

    swap2(v1, v2, i, right-1);
    internal_qsort_2(v1, v2, left, i-1);
    internal_qsort_2(v1, v2, i+1, right);
  }
}

void internal_isort_2(int *v1, int *v2, int N)
{
  int i,j;
  int ndx = 0;
  int small1, small2;
  
  ndx = 0;
  for (i = 1; i < N; i++) {
    if (is_less_than2v(v1, v2, i, ndx)) {
      ndx = i;
    }
  }

  /* Put smallest value in slot 0 */
  swap2(v1, v2, 0, ndx);

  for (i=1; i <N; i++) {
    small1 = v1[i];
    small2 = v2[i];
    for (j=i; is_less_than2(small1, small2, v1[j-1], v2[j-1]); j--) {
      v1[j] = v1[j-1];
      v2[j] = v2[j-1];
    }
    v1[j] = small1;
    v2[j] = small2;
  }
}

/*
 * Sort the values in 'v' 
 */
   
void qsort2(int *v1, int *v2, int N)
{
  internal_qsort_2(v1, v2, 0, N-1);
  internal_isort_2(v1, v2, N);

#if defined(DEBUG_QSORT)
  fprintf(stderr, "Checking sort of %d values\n", N+1);
  int i;
  for (i=1; i < N; i++) {
    assert(is_less_than2v(v1, v2, i-1, i));
  }
#endif
}

