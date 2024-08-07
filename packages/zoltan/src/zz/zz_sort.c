// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_sort.h"

#include "zz_const.h"


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

/****************************************************************************/
/* Sorting pointers in decreasing order. Criteria (key) is double. */
static void quickpart_pointer_dec_double (
  int *sorted, double *val, int start, int end, int* equal, int* smaller
)
{
int   i, next;
double key = (val ? val[sorted[(end+start)/2]] : 1.0);

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


void Zoltan_quicksort_pointer_dec_double (
  int *sorted, double* val, int start, int end
)
{
int  equal, smaller;

  if (start < end) {
     quickpart_pointer_dec_double (sorted, val, start, end, &equal, &smaller);
     Zoltan_quicksort_pointer_dec_double (sorted, val, start,   equal-1);
     Zoltan_quicksort_pointer_dec_double (sorted, val, smaller, end);
  }
}
/****************************************************************************/
/****************************************************************************/
/* Sort in increasing order by first calling the decreasing sort,
   then reverse the order in linear time. */
void Zoltan_quicksort_pointer_inc_double (
  int *sorted, double* val, int start, int end
)
{
  int i, j;
  double temp;

  /* sort in decreasing order */
  Zoltan_quicksort_pointer_dec_double(sorted, val, start, end);
  /* reverse order */
  for (i=start, j=end; i<j; i++, j--){
    temp = sorted[i];
    sorted[i] = sorted[j];
    sorted[j] = temp;
  }
}

/****************************************************************************/
/* Sort in increasing order by first calling the decreasing sort,
   then reverse the order in linear time. */
void Zoltan_quicksort_pointer_inc_float (
  int *sorted, float* val, int start, int end
)
{
  int i, j;
  float temp;

  /* sort in decreasing order */
  Zoltan_quicksort_pointer_dec_float(sorted, val, start, end);
  /* reverse order */
  for (i=start, j=end; i<j; i++, j--){
    temp = sorted[i];
    sorted[i] = sorted[j];
    sorted[j] = temp;
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
void Zoltan_quicksort_pointer_inc_int_int(
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

/* Same code as above except that the primary key is ZOLTAN_GNO_TYPE */

static void quickpart_pointer_inc_gno_int(
  int *sorted, ZOLTAN_GNO_TYPE *val1, int *val2, int start, int end, int *equal, int *larger)
{
int i, next;
int key2, key2_next;
ZOLTAN_GNO_TYPE key1, key1_next; 

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

void Zoltan_quicksort_pointer_inc_gno_int(
  int *sorted, ZOLTAN_GNO_TYPE *val1, int *val2, int start, int end)
{
int  equal, larger;

  if (start < end) {
     quickpart_pointer_inc_gno_int (sorted,val1,val2,start,end,&equal,&larger);
     Zoltan_quicksort_pointer_inc_gno_int (sorted, val1, val2, start, equal-1);
     Zoltan_quicksort_pointer_inc_gno_int (sorted, val1, val2, larger, end);
  }
}

/****************************************************************************/



/* Sorting values in array list in increasing order. Criteria is int. */
static void quickpart_list_inc_one_int (
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

void Zoltan_quicksort_list_inc_one_int (int* list, int start, int end)
{
int  equal, larger;

  if (start < end) {
    quickpart_list_inc_one_int (list, start, end, &equal, &larger);
    Zoltan_quicksort_list_inc_one_int (list, start,  equal-1);
    Zoltan_quicksort_list_inc_one_int (list, larger, end);
  }
}

/* Also rearrange values in array parlist to match the new order of list. */
static void quickpart_list_inc_int (
  int *list, int *parlist, int start, int end, int *equal, int *larger)
{
int i, key, change, parchange;

  key = list ? list[(end+start)/2] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++)
    if (list[i] < key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parlist[*equal];
      parlist[(*equal)] = parchange;
      change            = list[i];
      list[i]           = list[*larger];
      list[(*larger)++] = list[*equal];
      list[(*equal)++]  = change;
    }
    else if (list[i] == key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parchange;
      list[i]           = list[*larger];
      list[(*larger)++] = key;
    }
}

void Zoltan_quicksort_list_inc_int (int* list, int *parlist, int start, int end)
{
int  equal, larger;

  if (start < end) {
    quickpart_list_inc_int (list, parlist, start, end, &equal, &larger);
    Zoltan_quicksort_list_inc_int (list, parlist, start,  equal-1);
    Zoltan_quicksort_list_inc_int (list, parlist, larger, end);
  }
}

/* Exact same code except the list to be sorted is ZOLTAN_GNO_TYPEs */

static void quickpart_list_inc_gno (
  ZOLTAN_GNO_TYPE *list, int *parlist, int start, int end, int *equal, int *larger)
{
int i, parchange;
ZOLTAN_GNO_TYPE key, change;

  key = list ? list[(end+start)/2] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++)
    if (list[i] < key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parlist[*equal];
      parlist[(*equal)] = parchange;
      change            = list[i];
      list[i]           = list[*larger];
      list[(*larger)++] = list[*equal];
      list[(*equal)++]  = change;
    }
    else if (list[i] == key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parchange;
      list[i]           = list[*larger];
      list[(*larger)++] = key;
    }
}

void Zoltan_quicksort_list_inc_gno(ZOLTAN_GNO_TYPE * list, int *parlist, int start, int end)
{
int  equal, larger;

  if (start < end) {
    quickpart_list_inc_gno(list, parlist, start, end, &equal, &larger);
    Zoltan_quicksort_list_inc_gno(list, parlist, start,  equal-1);
    Zoltan_quicksort_list_inc_gno(list, parlist, larger, end);
  }
}

/* Exact same code except the list to be sorted is short*/

static void quickpart_list_inc_short(
  short *list, int *parlist, int start, int end, int *equal, int *larger)
{
int i, parchange;
short key, change;

  key = list ? list[(end+start)/2] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++)
    if (list[i] < key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parlist[*equal];
      parlist[(*equal)] = parchange;
      change            = list[i];
      list[i]           = list[*larger];
      list[(*larger)++] = list[*equal];
      list[(*equal)++]  = change;
    }
    else if (list[i] == key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parchange;
      list[i]           = list[*larger];
      list[(*larger)++] = key;
    }
}

void Zoltan_quicksort_list_inc_short(short *list, int *parlist, int start, int end)
{
int  equal, larger;

  if (start < end) {
    quickpart_list_inc_short(list, parlist, start, end, &equal, &larger);
    Zoltan_quicksort_list_inc_short(list, parlist, start,  equal-1);
    Zoltan_quicksort_list_inc_short(list, parlist, larger, end);
  }
}

/* Exact same code except the list to be sorted is long*/

static void quickpart_list_inc_long (
  long *list, int *parlist, int start, int end, int *equal, int *larger)
{
int i, parchange;
long key, change;

  key = list ? list[(end+start)/2] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++)
    if (list[i] < key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parlist[*equal];
      parlist[(*equal)] = parchange;
      change            = list[i];
      list[i]           = list[*larger];
      list[(*larger)++] = list[*equal];
      list[(*equal)++]  = change;
    }
    else if (list[i] == key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parchange;
      list[i]           = list[*larger];
      list[(*larger)++] = key;
    }
}

void Zoltan_quicksort_list_inc_long(long *list, int *parlist, int start, int end)
{
int  equal, larger;

  if (start < end) {
    quickpart_list_inc_long(list, parlist, start, end, &equal, &larger);
    Zoltan_quicksort_list_inc_long(list, parlist, start,  equal-1);
    Zoltan_quicksort_list_inc_long(list, parlist, larger, end);
  }
}

/* Exact same code except the list to be sorted is long long*/

static void quickpart_list_inc_long_long (
  int64_t *list, int *parlist, int start, int end, int *equal, int *larger)
{
int i, parchange;
long long key, change;

  key = list ? list[(end+start)/2] : 1;

  *equal = *larger = start;
  for (i = start; i <= end; i++)
    if (list[i] < key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parlist[*equal];
      parlist[(*equal)] = parchange;
      change            = list[i];
      list[i]           = list[*larger];
      list[(*larger)++] = list[*equal];
      list[(*equal)++]  = change;
    }
    else if (list[i] == key) {
      parchange         = parlist[i];
      parlist[i]        = parlist[*larger];
      parlist[(*larger)]= parchange;
      list[i]           = list[*larger];
      list[(*larger)++] = key;
    }
}

void Zoltan_quicksort_list_inc_long_long(int64_t*list, int *parlist, int start, int end)
{
int  equal, larger;

  if (start < end) {
    quickpart_list_inc_long_long(list, parlist, start, end, &equal, &larger);
    Zoltan_quicksort_list_inc_long_long(list, parlist, start,  equal-1);
    Zoltan_quicksort_list_inc_long_long(list, parlist, larger, end);
  }
}



/****************************************************************************/



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
