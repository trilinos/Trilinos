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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"



/****************************************************************************/



/* Sorting pointers in decreasing order. Criteria (key) is float. */
static void quickpart_pointer_dec_float (
  int *sorted, float *val, int start, int end, int* equal, int* smaller
)
{
  int   i, next;
  float key = (val ? val[sorted[(end+start)/2]] : 1.0);

  *equal = *smaller = start;
  for (i = start; i <= end; i++) {
    next = sorted[i];
    if ((val ? val[next] : 1.0) > key) {
      sorted[i]            = sorted[*smaller];
      sorted[(*smaller)++] = sorted[*equal];
      sorted[(*equal)++]   = next;
    }
    else if ((val ? val[next] : 1.0) == key) {
      sorted[i]            = sorted[*smaller];
      sorted[(*smaller)++] = next;
    }
  }
}



void Zoltan_quicksort_pointer_dec_float (
  int *sorted, float* val, int start, int end
)
{
  int  equal, smaller;

  if (start < end) {
    quickpart_pointer_dec_float (sorted, val, start, end, &equal, &smaller);
    Zoltan_quicksort_pointer_dec_float (sorted, val, start,   equal-1);
    Zoltan_quicksort_pointer_dec_float (sorted, val, smaller, end);
  }
}



/****************************************************************************/
/* Sorting pointers in decreasing order. Sort key is float. Sub key is int. */
static void quickpart_pointer_dec_float_int (
  int *sorted, float *val1, int *val2, int start, int end, int *equal,
  int *smaller
)
{
  int i, next, key2, key2_next;
  float key1, key1_next;

  i = (end + start) / 2;
  key1 = val1 ? val1[sorted[i]] : 1.0;
  key2 = val2 ? val2[sorted[i]] : 1;

  *equal = *smaller = start;
  for (i = start; i <= end; i++) {
    next = sorted[i];
    key1_next = val1 ? val1[next] : 1.0;
    key2_next = val2 ? val2[next] : 1;
    if (key1_next > key1 || (key1_next == key1 && key2_next > key2)) {
      sorted[i]            = sorted[*smaller];
      sorted[(*smaller)++] = sorted[*equal];
      sorted[(*equal)++]   = next;
    }
    else if (key1_next == key1 && key2_next == key2) {
      sorted[i]            = sorted[*smaller];
      sorted[(*smaller)++] = next;
    }
  }
}



void Zoltan_quicksort_pointer_dec_float_int (
  int *sorted, float* val1, int *val2, int start, int end
)
{
  int  equal, smaller;

  if (start < end) {
    quickpart_pointer_dec_float_int(sorted,val1,val2,start,end,&equal,&smaller);
    Zoltan_quicksort_pointer_dec_float_int (sorted,val1,val2,start,equal-1);
    Zoltan_quicksort_pointer_dec_float_int (sorted,val1,val2,smaller,end);
  }
}



/****************************************************************************/
/* Sorting pointers in increasing order. Sort key is int. Sub key is int. */
static void quickpart_pointer_inc_int_int (
  int *sorted, int *val1, int *val2, int start, int end, int *equal, int *larger)
{
  int i, next, key1, key1_next, key2, key2_next;

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
void Zoltan_quicksort_pointer_inc_int_int (
  int *sorted,   /* index array that is rearranged; should be initialized
                    so that sorted[i] == i. */
  int* val1,     /* array of primary key values. */
  int *val2,     /* array of secondary key values. */
  int start,     /* first array position to be considered for sorting. */
  int end        /* last array position to be considered for sorting. */
)
{
  int  equal, larger;

  if (start < end) {
    quickpart_pointer_inc_int_int (sorted,val1,val2,start,end,&equal,&larger);
    Zoltan_quicksort_pointer_inc_int_int (sorted, val1, val2, start, equal-1);
    Zoltan_quicksort_pointer_inc_int_int (sorted, val1, val2, larger, end);
  }
}



/****************************************************************************/
/* Sorting values in increasing order. Criteria is int. */
static void quickpart_list_inc_int (
  int *list, int start, int end, int *equal, int *larger)
{
  int i, key, change;

  key = list ? list[(end+start)/2] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++)
    if (list[i] < key) {
      change            = list[i];
      list[i]           = list[*larger];
      list[(*larger)++] = list[*equal];
      list[(*equal)++]  = change;
    }
    else if (list[i] == key) {
      list[i]           = list[*larger];
      list[(*larger)++] = key;
      }
}



void Zoltan_quicksort_list_inc_int (int* list, int start, int end)
{
  int  equal, larger;

  if (start < end) {
    quickpart_list_inc_int (list, start, end, &equal, &larger);
    Zoltan_quicksort_list_inc_int (list, start, equal-1);
    Zoltan_quicksort_list_inc_int (list, larger,end);
  }
}



/****************************************************************************/
/* Sorting pointers in increasing order. Criteria are multiple ints. */
static void quickpart_pointer_inc_int_mult (
  int *sorted, int start, int end, int *equal, int *larger, int *index,
  int *data)
{
  int i, j_new, j_old, new, old;

  *equal = *larger = start;
  for (i = start; i <= end; i++) {
    new   = sorted[i];
    old   = sorted[*equal];
    j_old = index[old];
    j_new = index[new];
    while (j_old < index[old+1] && j_new < index[new+1]
     && data[j_old] == data[j_new]) {
      j_new++;
      j_old++;
    }
    if (j_old == index[old+1] && j_new == index[new+1]) {
      sorted[i]           = sorted[*larger];
      sorted[(*larger)++] = new;
    }
    else if (j_new == index[new+1] || (j_old < index[old+1]
     && data[j_new] < data[j_old])) {
      sorted[i]           = sorted[*larger];
       sorted[(*larger)++] = old;
       sorted[(*equal)++]  = new;
    }
  }
}



void Zoltan_quicksort_pointer_inc_int_mult (int *sorted, int start, int end,
 int *index, int *data)
{
  int  equal, larger;

  if (start < end) {
    quickpart_pointer_inc_int_mult(sorted,start,end,&equal,&larger,index,data);
    Zoltan_quicksort_pointer_inc_int_mult (sorted, start, equal-1,index, data);
    Zoltan_quicksort_pointer_inc_int_mult (sorted, larger,end,    index, data);
    }
}

/****************************************************************************/



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
