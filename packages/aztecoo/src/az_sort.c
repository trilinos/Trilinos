/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"
#include <string.h>

int divider;
int type_size;

/* ----------------- External Definitions ---------------------------------*/

extern void move_dble(double first[], double second[], unsigned int 
		      first_length, unsigned int second_length);

extern void move_ints(int first[], int second[], unsigned int first_length, 
               unsigned int second_length);

/* ------------------------------------------------------------------------*/

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_quick_find(int key, int list[], int length, int shift, int bins[])

/*******************************************************************************

  Find 'key' in 'list' and return the indices number. On exit, AZ_find_index()
  returns:
            -1 ==> key not found
             i ==> list[i] = key

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  key:             Element to be search for in list.

  list:            List (assumed to be in ascending order) to be searched.

  length:          Length of list.

  shift:           Integer used to compute an indices into bins[]. Computed by
                   AZ_init_quick_find().

  bins:            Used to obtain a sublist where should contain 'key'. In
                   particular, list[bins[k] ... bins[k+1]-1] should contain
                   'key'. Computed by AZ_init_quick_find().

*******************************************************************************/

{

  /* local variables */

  int i, loc, oldkey;

  /**************************** execution begins ******************************/

  if (length == 0)            return -1;
  if (key > list[length - 1]) return -1;

  oldkey = key;
  key   -= list[0];

  if (key < 0) return -1;

  loc = key >> shift;

  i = AZ_find_index(oldkey, &(list[bins[loc]]), bins[loc + 1] - bins[loc]);

  if (i == -1) return -1;

  return (i + bins[loc]);

} /* AZ_quick_find */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_init_quick_find(int list[], int length, int *shift, int *bins)

/*******************************************************************************

  The value 'shift' and the array 'bins' are initialized so that subsequent
  calls to AZ_quick_find() will work properly. In particular, the array 'bins'
  and which should be 1/4 the size of list is set so that

    1)  range>>shift     > length/4
              and
        range>>(shift+1) < length/4

    where range = list[length-1] - list[0].

    2)  list[j] = value  ==> bins[k] <= j < bins[k+1]
                             where k = (value-list[0]) >> shift

 Note: list[] is sorted in ascending order. This routine is used in conjunction
 with AZ_quick_find(). The idea is to use bins[] to get a good initial guess as
 to the location of 'value' in list[].

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  list:            List (assumed to be in ascending order) to be searched.

  length:          Length of list.

  shift:           Integer used to compute an indices into bins[]. Computed by
                   AZ_init_quick_find().

  bins:            Used to obtain a sublist where should contain 'key'. In
                   particular, list[bins[k] ... bins[k+1]-1] should contain
                   'key'. Computed by AZ_init_quick_find().

*******************************************************************************/

{

  /* local variables */

  register int i, j = 0;
  int          range, temp;

  /**************************** execution begins ******************************/

  if (length == 0) return;

  range  = list[length - 1] - list[0];
  *shift = 0;

  while ((range >> (*shift)) > length / 4)
    (*shift)++;

  bins[j++] = 0;

  for (i = 0; i < length; i++) {
    temp = list[i] - list[0];

    while ((temp >> (*shift)) >= j)
      bins[j++] = i;
  }

  bins[j] = length;

} /* AZ_init_quick_find */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sort(int list[], int N, int list2[], double list3[])

/*******************************************************************************

  This routine was taken from Knuth: Sorting and Searching. It puts the input
  data list into a heap and then sorts it.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     double, maximum value in vector 'vec'
  ============

  Parameter list:
  ===============

  list:            On input, values to be sorted. On output, sorted values
                   (i.e., list[i] <= list[i+1]).

  N:               length of vector 'vec'.

  list2:           If on input,
                   a) list2 = NULL: it is unchanged on output,
                   b) list2 is a list associated with 'list':
                   on output, if list[k] on input is now element 'j' on output,
                   list2[j] on output is list2[k].

  list3:           If on input,
                   a) list3 = NULL: it is unchanged on output,
                   b) list3 is a list associated with 'list':
                   on output, list3[j] is assigned the input value of list3[k],
                   if list[j] has been assigned the input value of list[k].

*******************************************************************************/

{

  /* local variables */

  int    l, r, RR, K, j, i, flag;
  int    RR2;
  double RR3;

  /**************************** execution begins ******************************/

  if (N <= 1) return;

  l   = N / 2 + 1;
  r   = N - 1;
  l   = l - 1;
  RR  = list[l - 1];
  K   = list[l - 1];

  if ((list2 != NULL) && (list3 != NULL)) {
    RR2 = list2[l - 1];
    RR3 = list3[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list2[i - 1] = list2[j - 1];
            list3[i - 1] = list3[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list2[i - 1] = RR2;
      list3[i - 1] = RR3;

      if (l == 1) {
        RR  = list [r];
        RR2 = list2[r];
        RR3 = list3[r];

        K = list[r];
        list[r ] = list[0];
        list2[r] = list2[0];
        list3[r] = list3[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR2 = list2[l - 1];
        RR3 = list3[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list2[0] = RR2;
    list3[0] = RR3;
  }
  else if (list2 != NULL) {
    RR2 = list2[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list2[i - 1] = list2[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list2[i - 1] = RR2;

      if (l == 1) {
        RR  = list [r];
        RR2 = list2[r];

        K = list[r];
        list[r ] = list[0];
        list2[r] = list2[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR2 = list2[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list2[0] = RR2;
  }
  else if (list3 != NULL) {
    RR3 = list3[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list3[i - 1] = list3[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list3[i - 1] = RR3;

      if (l == 1) {
        RR  = list [r];
        RR3 = list3[r];

        K = list[r];
        list[r ] = list[0];
        list3[r] = list3[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR3 = list3[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list3[0] = RR3;

  }
  else {
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;

      if (l == 1) {
        RR  = list [r];

        K = list[r];
        list[r ] = list[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
  }

} /* AZ_sort */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sortqlists(char a[], int b[], int indices[], int length, int type_length,
                int ind_length)

/*******************************************************************************

  Sort the set of sublists in a[]. In particular a[] contains a sequence of k
  sublists (the length of sublist j is given by b[j] unless b=0 then the list is
  of size 1). The array indices[] contains the numbers 0 to k-1. On output, the
  sublists in a[] are ordered such that the jth list on input becomes the
  indices[j]th list on output.

  IMPORTANT: This routine assumes that indices[] contains 2 sequencies of numbers
  that are ordered but intertwined. For example,

         indices:   4 5 0 6 1 2 3 7
         seq 1 =>     0   1 2 3
         seq 2 => 4 5   6       7

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  a:               Set of lists to be sorted.
                   List 0 consists of a[0-->b[0]-1], list 1 consists of
                   a[b[0]-->b[0]+b[1]-1], list 2 consists of
                   a[b[0]+b[1]-->b[0]+b[1]+b[2]-1], etc.
                   On output, a[] is reordered such that list j on input
                   becomes list indices[j].

  b:               b[i] is the length of list i in a[]. If b = 0 then all the
                   sublists are of length 1.

  indices:           Contains the numbers 0 --> number of lists. It is assumed
                   that indices[] contains 2 sequencies of numbers that are
                   ordered but intertwined as discussed in the comments above.
                   On output, list j on input becomes list indices[j].

  length:          The length of a[].

  type_length:     This is the size of each element in a[]. Normally, if we are
                   sorting integers ... type_length = 4.

  ind_length:      Number of sublists to be sorted.

*******************************************************************************/

{

  /* local variables */

  char *abuffer = 0;
  int   next_ind = 0;
  int   i, midpoint, real_lists, the_first, buf_len;

  /**************************** execution begins ******************************/

  type_size = type_length;

  /*
   * Divider is the indices of the first sublist in sequence 2.  In the example
   * above, divider = 4.
   */

  divider = -1;
  for (i = 0; i < ind_length; i++) {
    if (indices[i] != next_ind) {
      divider = indices[i];
      break;
    }

    next_ind++;
  }

  if (divider != -1) {

    /* allocate a temporary buffer ... try and make it as big as possible */

    buf_len = type_length * length/2;
    while (abuffer == 0) {
      buf_len = buf_len / 2;
      abuffer = (char *) AZ_allocate(buf_len*sizeof(char));
    }

    /*
     * Effectively merge together lists that are already in the correct order.
     * That is, two consecutive lists which are also to appear consecutively on
     * output are grouped together. This grouping is accomplished by changing
     * the arrays: indices and b.
     */

    AZ_change_it(indices, ind_length, &the_first, &real_lists, b);

    /* now sort the merged sublists */

    if (type_length == 4)
      AZ_sort_ints(a, indices, 0, type_length*length-1, b, &midpoint, real_lists,
                   abuffer, buf_len, the_first, 0);
    else if (type_length == 8)
      AZ_sort_dble(a, indices, 0, type_length*length-1, b, &midpoint, real_lists,
                   abuffer, buf_len, the_first, 0);
    else
      (void) AZ_printf_err("ERROR: unknown type size in AZ_sortqlists\n");

    AZ_free(abuffer);

    /*
     * Reverse the changes made by AZ_change_it() so that the arrays: indices and
     * b have the same values that they had before AZ_change_it()
     */

    AZ_reverse_it(indices, ind_length, the_first, real_lists, b);
  }

} /* AZ_sortqlists */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void move_dble(double first[], double second[], unsigned int first_length,
               unsigned int second_length)

/*******************************************************************************

  Swap the contents of the arrays 'first' and 'second'. NOTE: THE ARRAY 'second'
  FOLLOWS IMMEDIATELY AFTER THE ARRAY 'first' IN STORAGE. For example,

          ----------------------------------------------
          |first                             |second   |
          ----------------------------------------------

  On output, the two arrays are swapped.

          ----------------------------------------------
          |second   |first                             |
          ----------------------------------------------

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  first:           Pointer to first array (as illustrated above).

  second:          Pointer to second array (as illustrated above).

  first_length:    Length of array 'first'.

  second_length:   Length of array 'second'.

*******************************************************************************/

{

  /* local variables */

  int    working = 1, tlength, i;
  double temp;

  /**************************** execution begins ******************************/

  if ((first_length == 0) || (second_length == 0)) working = 0;
  while (working) {
    if (first_length < second_length) tlength = first_length;
    else                              tlength = second_length;

    for (i = 0; i < tlength; i++) {
      temp      = first[i];
      first[i]  = second[i];
      second[i] = temp;
    }

    if (first_length  > second_length) {
      first         = &(first[second_length]);
      first_length -= second_length;
    }
    else if (first_length < second_length) {
      first          = &(first[first_length]);
      second         = &(second[first_length]);
      second_length -= first_length;
    }
    else
      working = 0;
  }

} /* move_dble */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void move_ints(int first[], int second[], unsigned int first_length, 
               unsigned int second_length)

/*******************************************************************************

  Swap the contents of the arrays 'first' and 'second'. NOTE: THE ARRAY 'second'
  FOLLOWS IMMEDIATELY AFTER THE ARRAY 'first' IN STORAGE. For example,

          ----------------------------------------------
          |first                             |second   |
          ----------------------------------------------

  On output, the two arrays are swapped.

          ----------------------------------------------
          |second   |first                             |
          ----------------------------------------------

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  first:           Pointer to first array (as illustrated above).

  second:          Pointer to second array (as illustrated above).

  first_length:    Length of array 'first'.

  second_length:   Length of array 'second'.

*******************************************************************************/

{

  /* local variables */

  int i, tlength, temp, working = 1;

  /**************************** execution begins ******************************/

  if ((first_length == 0) || (second_length == 0)) working = 0;
  while (working) {
    if (first_length < second_length) tlength = first_length;
    else                              tlength = second_length;

    for (i = 0; i < tlength; i++) {
      temp      = first[i];
      first[i]  = second[i];
      second[i] = temp;
    }

    if (first_length  > second_length) {
      first         = &(first[second_length]);
      first_length -= second_length;
    }
    else if (first_length < second_length) {
      first          = &(first[first_length]);
      second         = &(second[first_length]);
      second_length -= first_length;
    }
    else
      working = 0;
  }

} /* move_ints */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_change_it(int indx[], int length, int *first, int *total, int b[])

/*******************************************************************************

  Modify indx[] and b[] to effectively merge the lists found in indx[]. That is,
  indx[] and b[] will correspond to a set of merged lists where all consecutive
  elements belonging to the same subsequence (see AZ_sortqlists() comments) are
  grouped together.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  indx:            On input, indx[] contains the numbers 0 --> 'length'-1. It is
                   assumed that indices[] contains 2 sequencies of numbers that
                   are ordered but intertwined as discussed in the comments for
                   AZ_sortqlists(). After AZ_sortqlists() reorders list j on
                   input will become list indices[j].
                   On output to AZ_change_it(),
                     indx[i] - indx[i-1] i > 0
                     indx[i]             i = 0
                   gives the length of the ith merged list.

  length:          The length of indx[].

  first:           Indicates if the first merged list belongs to the first
                   subsequence or the second subsequence. In particular ...
                   On output,
                     first = 0   ==> the first merged list will be the first
                                     list after sorting.
                     first = 1   ==> the first merged list will appear
                                     immediately after the first subsequence
                                     when the sorting is finished.
                   NOTE: IT IS ASSUMED THAT THE MERGED LISTS ALTERNATE. That is,
                   if the first merged lists belongs to the first subsequence,
                   then all the even merged lists belong to the second
                   subsequence and all the odd merged lists belong to the first
                   subsequence.

  total:           On output, number of merged lists.

  b:               On input, b[i] indicates the length of list i.
                   On output,
                       b[indx[i-1]]     i > 0
                       b[i]             i = 0
                  gives the length of merged list i.

*******************************************************************************/

{

  /* local variables */

  int i, ii = 0, count = 0;

  /**************************** execution begins ******************************/

  if (indx[0] == 0) *first = 0;
  else              *first = 1;

  for (i = 1; i < length; i++) {
    if (indx[i-1] < divider) {
      if (indx[i] >= divider)
        indx[ii++] = i;
    }
    else if (indx[i] < divider)
      indx[ii++] = i;
  }

  *total         = ii + 1;
  indx[*total-1] = length;

  if (b != 0) {
    for (i = 1; i < *total; i++) {
      count = 0;

      for (ii = indx[i-1]; ii < indx[i]; ii++)
        count += b[ii];

      count       *= type_size;
      b[indx[i-1]] = count;
    }

    count = 0;
    for (ii = 0; ii < indx[0]; ii++)
      count += b[ii];

    count *= type_size;
    b[0]   = count;
  }

} /* AZ_change_it */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_reverse_it(int indx[], int length, int first, int total, int b[])

/*******************************************************************************

  Reverse the changes performed in AZ_change_it(). That is, the arrays 'indx'
  and 'b' should be identical before and after the following code sequence:

        AZ_change_it(indx, b)
        AZ_reverse_it(indx, b)

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  indx:            On input, indx[] contains the numbers 0 --> 'length'-1. It is
                   assumed that indices[] contains 2 sequencies of numbers that
                   are ordered but intertwined as discussed in the comments for
                   AZ_sortqlists(). After AZ_sortqlists() reorders list j on
                   input will become list indices[j].
                   On output to AZ_change_it(),
                     indx[i] - indx[i-1] i > 0
                     indx[i]             i = 0
                   gives the length of the ith merged list.

  length:          The length of indx[].

  first:           Indicates if the first merged list belongs to the first
                   subsequence or the second subsequence. In particular ...
                   On output,
                     first = 0   ==> the first merged list will be the first
                                     list after sorting.
                     first = 1   ==> the first merged list will appear
                                     immediately after the first subsequence
                                     when the sorting is finished.
                   NOTE: IT IS ASSUMED THAT THE MERGED LISTS ALTERNATE. That is,
                   if the first merged lists belongs to the first subsequence,
                   then all the even merged lists belong to the second
                   subsequence and all the odd merged lists belong to the first
                   subsequence.

  total:           On output, number of merged lists.

  b:               On input, b[i] indicates the length of list i.
                   On output,
                       b[indx[i-1]]     i > 0
                       b[i]             i = 0
                  gives the length of merged list i.

*******************************************************************************/

{

  /* local variables */

  int i, j, ii, count;
  int seq1, seq2;
  int cur, toggle, sub_length;

  /**************************** execution begins ******************************/

  if (b != 0) {
    count = 0;
    for (ii = 1; ii < indx[0]; ii++)
      count += b[ii];
    count *= type_size;

    b[0] = (b[0] - count) / type_size;

    for (i = 1; i < total; i++) {
      count = 0;
      for (ii = indx[i-1]+1; ii < indx[i]; ii++)
        count += b[ii];
      count *= type_size;

      b[indx[i-1]] = (b[indx[i-1]] - count) / type_size;
    }
  }

  seq1 = divider - 1;
  seq2 = cur = length - 1;

  if (first == 0) {
    if (total%2 == 0) toggle = 1;
    else              toggle = 0;
  }
  else {
    if (total%2 == 0) toggle = 0;
    else              toggle = 1;
  }

  for (i = total-1; i > 0; i--) {
    sub_length = indx[i] - indx[i-1];

    if (toggle)
      for (j = 0; j < sub_length; j++)
        indx[cur--] = seq2--;
    else
      for (j = 0; j < sub_length; j++)
        indx[cur--] = seq1--;

    toggle = 1 - toggle;
  }

  sub_length = indx[0];

  if (toggle)
    for (j = 0; j < sub_length; j++)
      indx[cur--] = seq2--;
  else
    for (j = 0; j < sub_length; j++ )
      indx[cur--] = seq1--;

} /* AZ_reverse_it */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sort_dble(char a[], int indx[], int start, int end, int b[], int *mid,
               int real_lists, char buffer[], int buf_len, int the_first,
               int ind_indices)

/*******************************************************************************

  Sort the set of merged lists in a[]. In particular a[] contains a sequence
  of 'real_lists' merged lists of the form
        list 0   list 1  list 2   list 3
  Each merged list is in fact comprised of a bunch of sublists which consist
  of double precision numbers.  On output, AZ_sort_dble will put all the even
  lists first followed by all the odd lists (e.g. list 0, list 2, list4, ...
  list1, list3, list 5, ...) if 'the_first' == 0. Otherwise, the odd lists are
  first followed by the even lists.

  The algorithm proceeds as follows:
      1) sort as many merged lists as we can using a buffer (of size buf_len)
         and simply copy what does not belong in the first subsequence
         (either even or odd depending on 'the_first') into the buffer.
      2) split the remaining merged lists into 2 sublists.
      3) use AZ_sort_dble() recursively for each of the 2 sublists.
      4) merge the results of the 2 recursive calls so that the entire
         list is sorted.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  a:               Set of merged lists to be sorted.
                   List 0 consists of a[0-->b[0]-1], list 1 consists of
                   a[b[0]-->b[0]+b[1]-1], list 2 consists of
                   a[b[0]+b[1]-->b[0]+b[1]+b[2]-1], etc.
                   On output, a[] is reordered such that all the even (or
                   all the odd if the_first == 1) come first as described above.

  indx:            indx[i] - indx[i-1] i > 0 (and indx[i] for i = 0)
                   gives the number of sublists comprising the ith merged list.

  start:           a[start] is the first element of the first list to be sorted.

  end:             a[end] is the last element of the last list to be sorted.

  b:               The length of the ith merged list in a[] is given by
                       if (b==0) && (i == 0) indx[0]
                       if (b==0) && (i >  0) indx[i] - indx[i-1]
                       if (b!=0) && (i == 0) b[0]
                       if (b!=0) && (i >  0) b[indx[i-1]]

  mid:             On output, a[mid] is the first element of the second
                   subsequence after sorting (i.e. the first element in the
                   first even merged list if 'the_first' == 0).

  real_lists:      Number of merged lists to be sorted.

  buffer:          A workspace array used to copy values during the sorting.

  buf_len:         The length of the workspace array 'buffer'.

  the_first:       Indicates whether we put the evens or the odd lists first.
                   That is, a[start ... end] contains a series of merged lists
                   to be sorted (l0 l1 l2 l3 l4 l5 l6 ... )
                   the_first = 0 on input, indicates that on output the list
                   will be sorted as (l0 l2 l4 l6 ... l1 l3 l5 ...)
                   the_first = 1 on input, indicates that on output the list
                   will be sorted as (l1 l3 l5 l7 ... l0 l2 l4 ...)

  ind_indices:       Index into  indx[] pointing to the first sublist within
                   the first merged list that will be sorted.

*******************************************************************************/

{

  /* local variables */

  int pre_mid, mid1, mid2, real1, real2;
  int count, rest, i;
  int first2;

  /**************************** execution begins ******************************/

  /* sort as many lists as we can directly using the buffer */

  AZ_direct_sort(b, indx, buffer, a, &start, buf_len, &ind_indices, &the_first,
                 &real_lists, &pre_mid);

  if (real_lists > 2 ) {

    /* split the rest of the list */

    real1 = real_lists >> 1;
    real2 = real_lists - real1;

    if ((real1%2) == 0) first2 = the_first;
    else                first2 = 1 - the_first;

    if (b == 0) {
      count = indx[ind_indices+real1-1];

      if (ind_indices != 0)
        count -= indx[ind_indices-1];
      count *= type_size;
    }
    else {
      count = 0;
      for (i = ind_indices; i < ind_indices+real1-1; i++)
        count += b[indx[i]];

      if (ind_indices == 0) count += b[0];
      else                count += b[indx[ind_indices-1]];
    }

    AZ_sort_dble(a, indx, start, start+count-1, b, &mid1, real1, buffer,
                 buf_len, the_first, ind_indices);

    AZ_sort_dble(a, indx, start+count, end, b, &mid2, real2, buffer, buf_len,
                 first2, ind_indices + real1);

    /* merge the lists */

    if (start+count-1-mid1 < 0 )
      *mid = mid2;
    else if (mid2-1-(start+count) < 0)
      *mid = mid1;
    else {
      move_dble((double *) &(a[mid1]), (double *) &(a[start+count]),
      (start+count-mid1) / sizeof(double),
      (mid2-(start+count)) / sizeof(double));

      *mid = mid1 + (mid2 - (start + count));
    }
  }

  else {
    *mid = start;
    if (real_lists == 2) {
      if (ind_indices ==0) {
        if (b == 0) count = type_size * indx[0];
        else        count = b[0];
      }
      else {
        if (b == 0) count = type_size * (indx[ind_indices] - indx[ind_indices-1]);
        else        count = b[indx[ind_indices-1]];
      }

      rest = end - (start + count) + 1;

      if (the_first == 0)
        *mid = start + count;
      else {
        *mid = start + rest;

        move_dble((double *) &(a[start]), (double *) &(a[start+count]),
        count / sizeof(double), rest / sizeof(double));
      }
    }
    else if (real_lists == 1) {
      if (the_first == 1) *mid = start;
      else                *mid = end + 1;
    }
  }

  if (start > pre_mid) {
    if (*mid == start) {
      *mid = pre_mid;
    }
    else {
      move_dble((double *) &(a[pre_mid]), (double *) &(a[start]),
      (start-pre_mid) / sizeof(double),
      (*mid-start) / sizeof(double));

      *mid = pre_mid + *mid - start;
    }
  }

} /* AZ_sort_dble */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sort_ints(char a[], int indx[], int start, int end, int b[], int *mid,
                  int real_lists, char buffer[], int buf_len, int the_first,
                  int ind_indices)

/*******************************************************************************

  Sort the set of merged lists in a[]. In particular a[] contains a sequence of
  'real_lists' merged lists of the form
        list 0   list 1  list 2   list 3
  Each merged list is in fact comprised of a bunch of sublists which consist of
  integer numbers.  On output, AZ_sort_dble will put all the even lists first
  followed by all the odd lists (e.g. list 0, list 2, list4, ...  list1, list3,
  list 5, ...) if 'the_first' == 0. Otherwise, the odd lists are first followed
  by the even lists.

  The algorithm proceeds as follows:
      1) sort as many merged lists as we can using a buffer (of size buf_len)
         and simply copy what does not belong in the first subsequence
         (either even or odd depending on 'the_first') into the buffer.
      2) split the remaining merged lists into 2 sublists.
      3) use AZ_sort_dble() recursively for each of the 2 sublists.
      4) merge the results of the 2 recursive calls so that the entire
         list is sorted.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  a:               Set of merged lists to be sorted.
                   List 0 consists of a[0-->b[0]-1], list 1 consists of
                   a[b[0]-->b[0]+b[1]-1], list 2 consists of
                   a[b[0]+b[1]-->b[0]+b[1]+b[2]-1], etc.
                   On output, a[] is reordered such that all the even (or
                   all the odd if the_first == 1) come first as described above.

  indx:            indx[i] - indx[i-1] i > 0 (and indx[i] for i = 0)
                   gives the number of sublists comprising the ith merged list.

  start:           a[start] is the first element of the first list to be sorted.

  end:             a[end] is the last element of the last list to be sorted.

  b:               The length of the ith merged list in a[] is given by
                       if (b==0) && (i == 0) indx[0]
                       if (b==0) && (i >  0) indx[i] - indx[i-1]
                       if (b!=0) && (i == 0) b[0]
                       if (b!=0) && (i >  0) b[indx[i-1]]

  mid:             On output, a[mid] is the first element of the second
                   subsequence after sorting (i.e. the first element in the
                   first even merged list if 'the_first' == 0).

  real_lists:      Number of merged lists to be sorted.

  buffer:          A workspace array used to copy values during the sorting.

  buf_len:         The length of the workspace array 'buffer'.

  the_first:       Indicates whether we put the evens or the odd lists first.
                   That is, a[start ... end] contains a series of merged lists
                   to be sorted (l0 l1 l2 l3 l4 l5 l6 ... )
                   the_first = 0 on input, indicates that on output the list
                   will be sorted as (l0 l2 l4 l6 ... l1 l3 l5 ...)
                   the_first = 1 on input, indicates that on output the list
                   will be sorted as (l1 l3 l5 l7 ... l0 l2 l4 ...)

  ind_indices:       Index into  indx[] pointing to the first sublist within
                   the first merged list that will be sorted.

*******************************************************************************/

{

  /* local variables */

  int pre_mid, mid1, mid2, real1, real2;
  int count, rest, i;
  int first2;

  /**************************** execution begins ******************************/

  /* sort as many lists as we can directly using the buffer */

  AZ_direct_sort(b, indx, buffer, a, &start, buf_len, &ind_indices, &the_first,
                 &real_lists, &pre_mid);

  if (real_lists > 2 ) {

    /* split the rest of the list */

    real1 = real_lists >> 1;
    real2 = real_lists - real1;

    if ( (real1%2) == 0) first2 = the_first;
    else                 first2 = 1 - the_first;

    if (b == 0) {
      count = indx[ind_indices+real1-1];

      if (ind_indices != 0)
        count -= indx[ind_indices-1];

      count *= type_size;
    }
    else {
      count = 0;
      for (i = ind_indices; i < ind_indices+real1-1; i++)
        count += b[indx[i]];

      if (ind_indices == 0) count += b[0];
      else                count += b[indx[ind_indices-1]];
    }

    AZ_sort_ints(a, indx, start, start+count-1, b, &mid1, real1, buffer,
                 buf_len, the_first, ind_indices);

    AZ_sort_ints(a, indx, start+count, end, b, &mid2, real2, buffer, buf_len,
                 first2, ind_indices+real1);

    /* merge the lists */

    if(start+count-1-mid1 < 0 )
      *mid = mid2;
    else if (mid2-1-(start+count) < 0)
      *mid = mid1;
    else {
      move_ints((int *) &(a[mid1]), (int *) &(a[start+count]),
      (start+count-mid1) / sizeof(int),
      (mid2-(start+count)) / sizeof(int));

      *mid = mid1 + (mid2 - (start + count));
    }
  }
  else {
    *mid = start;
    if (real_lists == 2) {
      if (ind_indices == 0) {
        if (b == 0) count = type_size * indx[0];
        else        count = b[0];
      }
      else {
        if (b == 0) count =type_size * (indx[ind_indices] - indx[ind_indices-1]);
        else        count = b[indx[ind_indices-1]];
      }

      rest = end - (start + count) + 1;

      if (the_first == 0)
        *mid = start+count;
      else {
        *mid = start + rest;
        move_ints((int *) &(a[start]), (int *) &(a[start+count]),
        count / sizeof(int), rest / sizeof(int));
      }
    }

    else if (real_lists == 1) {
      if (the_first == 1) *mid = start;
      else                *mid = end + 1;
    }
  }

  if (start > pre_mid) {
    if (*mid == start) {
      *mid = pre_mid;
    }
    else {
      move_ints((int *) &(a[pre_mid]), (int *) &(a[start]),
      (start-pre_mid) / sizeof(int), (*mid-start) / sizeof(int));
      *mid = pre_mid + *mid - start;
    }
  }

} /* AZ_sort_ints */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_direct_sort(int b[], int indx[], char buffer[], char a[], int *start,
                    int buf_len, int *ind_indices, int *the_first,
                    int *real_lists, int *firstmid)

/*******************************************************************************

  Sort as many lists as we can using the buffer. In particular, we go through
  the lists (starting with indx[0]). If a list belongs to the first subsequence
  then we move it to its proper location. If a list belongs to the second
  subsequence we copy it to the buffer. We continue doing this until there is
  not enough room in the buffer. At this point, we copy the buffer contents back
  to the data array starting immediately after the last element put in its
  proper location of the first subsequence.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  b:               Contains lengths of merged lists (see AZ_change_it()).

  indx:            Indicates the number of sublists in this merged list.
                   indx[i] - indx[i-1]    i > 0
                   indx[i]                i = 0

  buffer:          Workspace array.

  a:               Data array containing lists to be sorted.

  start:           On input, starting location in a[] of information to be
                   sorted. On output, start indicates the location of the first
                   list in a[] which was not sorted or copied into the buffer.

  buf_len:         Length of workspace array (buffer[]).

  ind_indices:       Index of the first merged list not to be sorted.
                   On input, first list to be sorted with buffer.
                   On output, first list that was not sorted with buffer.

  the_first:       On input, the_first = 0   ==> first merged list belongs
                                                 to first subsequence
                             the_first = 1   ==> first merged list belongs
                                                 to second subsequence

  real_lists:      On input, number of merged lists to be sorted.
                   On output, number of merged lists remaining to be sorted.

  firstmid:        On output, pointer to first element in a[] that was stored
                   in the buffer (hence belongs to the second subsequence).

*******************************************************************************/

{

  /* local variables */

  int   buffer_full = 0, cur, i, si, flag;
  char *ptr1, *ptr2;
  unsigned int buf_end = 0, thelength;

  /**************************** execution begins ******************************/

  cur = *start;
  i   = *ind_indices;
  si  = *start;

  if (*the_first == 0) flag = 0;
  else                 flag = 1;

  while (buffer_full == 0) {
    if (flag == 1) {

      /* must move this to the buffer */

      if (i ==0) {
        if (b == 0) thelength = type_size * indx[0];
        else        thelength = b[0];
      }
      else {
        if (b == 0) thelength = type_size * (indx[i] - indx[i-1]);
        else        thelength = b[indx[i-1]];
      }

      if  ( (int) (buf_end+thelength) > buf_len )
        buffer_full = 1;
      else {
        ptr1 = (char *) &(buffer[buf_end]);
        ptr2 = (char *) &(a[si]);

        memcpy(ptr1, ptr2, thelength);
        buf_end += thelength;
        si      += thelength;
      }
      flag = 0;
    }
    else {

      /* move this stuff to the front */

      if (i ==0) {
        if (b == 0) thelength = type_size * indx[0];
        else        thelength = b[0];
      }
      else {
        if (b == 0) thelength = type_size * (indx[i] - indx[i-1]);
        else        thelength = b[indx[i-1]];
      }

      if (cur != si) {
        ptr1 = (char *) &(a[cur]);
        ptr2 = (char *) &(a[si]);
        memmove(ptr1, ptr2, thelength);
      }

      cur += thelength;
      si  += thelength;

      flag = 1;
    }

    i++;
    if (buffer_full == 0) {
      if ((i- (*ind_indices)) == (*real_lists) ) {
        buffer_full = 1;
        i++;
      }
    }
  }

  i--;

  /* empty the buffer back into the list */

  *firstmid = cur;
  ptr1      = (char *) &(a[cur]);
  memcpy(ptr1, buffer, buf_end);
  cur += buf_end;

  *real_lists = *real_lists - (i - *ind_indices);
  *start      = cur;
  *the_first  = 1;
  *ind_indices  = i;

} /* AZ_direct_sort */
