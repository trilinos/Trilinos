// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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

/* Sorting pointers in increasing order. Sort key is ZOLTAN_ID_TYPE. Sub key is
ZOLTAN_ID_TYPE. */
static void quickpart_pointer_inc_id_id (
  int *sorted, ZOLTAN_ID_TYPE *val1, ZOLTAN_ID_TYPE *val2,
  int start, int end, int *equal, int *larger)
{
int i, next; 
ZOLTAN_ID_TYPE key1, key1_next, key2, key2_next;

  i = (end + start) / 2;
  key1 = val1 ? val1[sorted[i]] : 1;
  key2 = val2 ? val2[sorted[i]] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++) {
     next = sorted[i];
     key1_next = val1 ? val1[next] : 1;
     key2_next = val2 ? val2[next] : 1;
     if (key1_next < key1 || (key1_next == key1 && key2_next < key2)) {
        sorted[i]           = sorted[*larger];
        sorted[(*larger)++] = sorted[*equal];
        sorted[(*equal)++]  = next;
     }
     else if (key1_next == key1  &&  key2_next == key2) {
        sorted[i]           = sorted[*larger];
        sorted[(*larger)++] = next;
     }
  }
}


/* Sorts in increasing order with primary key val1 and secondary key val2.
   The arrays val1 and val2 are not rearranged; rather the index array
   sorted is rearranged based on values in val1 and val2. */
void quicksort_pointer_inc_id_id(
  int *sorted,   /* index array that is rearranged; should be initialized
                    so that sorted[i] == i. */
  ZOLTAN_ID_TYPE * val1,     /* array of primary key values. */
  ZOLTAN_ID_TYPE  *val2,     /* array of secondary key values. */
  int start,     /* first array position to be considered for sorting. */
  int end        /* last array position to be considered for sorting. */
)
{
int  equal, larger;

  if (start < end) {
     quickpart_pointer_inc_id_id (sorted,val1,val2,start,end,&equal,&larger);
     quicksort_pointer_inc_id_id (sorted, val1, val2, start, equal-1);
     quicksort_pointer_inc_id_id (sorted, val1, val2, larger, end);
  }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
