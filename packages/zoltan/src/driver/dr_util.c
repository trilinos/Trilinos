/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "dr_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Author(s):   Gary L. Hennigan (SNL 9221)
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *      token_compare()
 *      strip_string()
 *      string_to_lower()
 *      clean_string()
 *      in_list()
 *      find_min()
 *      find_max()
 *      find_inter()
 *      sort_int()
 *      sort2_index()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int token_compare(char *token, const char *key)
{

  unsigned int i1, key_len, kcnt=0;

  key_len = strlen(key);

  for(i1=0; i1 < strlen(token); i1++)
  {
    if(isupper((int)(token[i1])))
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

    if(isupper((int)(in_string[cnt])))
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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function in_list() begins:
 *----------------------------------------------------------------------------
 * This function searches a vector for the input value. If the value is
 * found in the vector then it's index in that vector is returned, otherwise
 * the function returns -1;
 *****************************************************************************/
int in_list(const ZOLTAN_ID_TYPE value, const int count, ZOLTAN_ID_TYPE *vector)
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

int in_list2(const int value, const int count, int *vector)
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
ZOLTAN_ID_TYPE find_max(const ZOLTAN_ID_TYPE list_length, const ZOLTAN_ID_TYPE list[])

     /*
       Function which finds the largest integer from a vector of integers

       Author:          Scott Hutchinson (1421)
       Date:            8 December 1992

       */

/*****************************************************************************/

{

  /* local variables */

  register ZOLTAN_ID_TYPE i, max;

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

ZOLTAN_ID_TYPE find_min(const ZOLTAN_ID_TYPE list_length, const ZOLTAN_ID_TYPE list[])

/*
 *     Function which finds the smallest integer from a vector of integers
 *
 *      Author:          Scott Hutchinson (1421)
 *      Date:            8 December 1992
 *
 */

{

  /* local variables */

  register ZOLTAN_ID_TYPE i, min;

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
 * This function finds the intersection between two lists of integer values,
 * and returns the number of values in the intersection.
 *****************************************************************************/
int find_inter(const ZOLTAN_ID_TYPE set1[], const ZOLTAN_ID_TYPE set2[], const int length1,
               const int length2, const int prob_type, int inter_ptr[])

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
 *
 *      On return, find_inter returns 0 if there is no intersection.
 *      It returns the number of points in the intersection, if there
 *      is an intersection.
 */

{

  /* Local variables */

  register int    i, j, counter = 0;
  ZOLTAN_ID_TYPE        max_set1, min_set1, max_set2, min_set2;

/****************************** execution begins *****************************/

  /* Error check the arguments */
  if ((length1 <= 0) || (length2 <= 0)) return (counter);

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
          if (set1[i] == set2[j])  inter_ptr[counter++] = i;

    }
  } else if (prob_type == 1) {

    fprintf (stderr, "prob_type = 1 is unimplemented\n");
    return -1;

  } else if (prob_type == 2) {
    /*
     *    Find the maximum and the minimum of the two sets
     */
    max_set1 =  set1[length1-1];
    min_set1 =  set1[0];
    max_set2 =  set2[length2-1];
    min_set2 =  set2[0];
    /*
     *    Check for a possible overlaps in node numbers;
     *    If there is an overlap, then do the search using a linearly
     *    scaled method
     *
     */
    if ( (max_set2 >= min_set1) && (min_set2 <= max_set1) )   {
      for (i = 0, j = 0; i < length1; i++) {
        while ( (j < (length2-1)) && (set1[i] > set2[j]) ) j++;
        if (set1[i] == set2[j]) inter_ptr[counter++] = i;
      }
    }
  } else {

    fprintf (stderr, "prob_type = %d is unknown\n", prob_type);
    return -1;

  }

  return counter;

} /* find_inter **************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void safe_free(void **ptr) {
  if (*ptr != NULL) {
    free(*ptr);
    *ptr = NULL;
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void sort2_index(int n, ZOLTAN_ID_TYPE ra[], ZOLTAN_ID_TYPE sa[], int indx[])

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array
*
*       Sorts the array ra[0,..,(n-1)] in ascending numerical order using
*       heapsort algorithm.  Use array sa as secondary sort key; that is,
*       if (ra[i] == ra[j]), then compare sa[i], sa[j] to determine order.
*       Array ra is not reorganized.  An index array indx is built that
*       gives the new order.
*
*/

{
  int   l, j, ir, i, irra;
  ZOLTAN_ID_TYPE   rra;
  ZOLTAN_ID_TYPE   ssa;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  l=n >> 1;
  ir=n-1;
  for (;;) {
    if (l > 0) {
      rra=ra[indx[--l]];
      ssa=sa[indx[l]];
      irra = indx[l];
    }
    else {
      rra=ra[indx[ir]];
      ssa=sa[indx[ir]];
      irra=indx[ir];
     
      indx[ir]=indx[0];
      if (--ir == 0) {
        indx[0]=irra;
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && 
          ((ra[indx[j]] <  ra[indx[j+1]]) || 
           (ra[indx[j]] == ra[indx[j+1]] && sa[indx[j]] < sa[indx[j+1]])))
        ++j;
      if ((rra <  ra[indx[j]]) ||
          (rra == ra[indx[j]] && ssa < sa[indx[j]])) {
        indx[i] = indx[j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    indx[i]=irra;
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void sort_index(int n, int ra[], int indx[])

/*
*       Numerical Recipies in C source code
*       modified to have first argument an integer array
*
*       Sorts the array ra[0,..,(n-1)] in ascending numerical order using
*       heapsort algorithm. 
*       Array ra is not reorganized.  An index array indx is built that
*       gives the new order.
*
*/

{
  int   l, j, ir, i;
  int   rra, irra;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  l=n >> 1;
  ir=n-1;
  for (;;) {
    if (l > 0) {
      rra=ra[indx[--l]];
      irra = indx[l];
    }
    else {
      rra=ra[indx[ir]];
      irra=indx[ir];
     
      indx[ir]=indx[0];
      if (--ir == 0) {
        indx[0]=irra;
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && 
          (ra[indx[j]] <  ra[indx[j+1]]))
        ++j;
      if (rra <  ra[indx[j]]) {
        indx[i] = indx[j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    indx[i]=irra;
  }
}

void sort_id_type_index(int n, ZOLTAN_ID_TYPE ra[], ZOLTAN_ID_TYPE indx[])
{
  ZOLTAN_ID_TYPE   l, j, ir, i;
  ZOLTAN_ID_TYPE   irra;
  ZOLTAN_ID_TYPE  rra;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  l=n >> 1;
  ir=n-1;
  for (;;) {
    if (l > 0) {
      rra=ra[indx[--l]];
      irra = indx[l];
    }
    else {
      rra=ra[indx[ir]];
      irra=indx[ir];
     
      indx[ir]=indx[0];
      if (--ir == 0) {
        indx[0]=irra;
        return;
      }
    }
    i=l;
    j=(l << 1)+1;
    while (j <= ir) {
      if (j < ir && 
          (ra[indx[j]] <  ra[indx[j+1]]))
        ++j;
      if (rra <  ra[indx[j]]) {
        indx[i] = indx[j];
        j += (i=j)+1;
      }
      else j=ir+1;
    }
    indx[i]=irra;
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
