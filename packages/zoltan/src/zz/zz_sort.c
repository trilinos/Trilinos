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

#include "hypergraph.h"

/****************************************************************************/

static void quickpart_dec_float (
int *sorted, float *val, int start, int end, int* equal, int* smaller
)
{ int   i, next;
  float key=(val?val[sorted[(end+start)/2]]:1.0);

  (*equal) = (*smaller) = start;
  for (i=start; i<=end; i++)
  { next = sorted[i];
    if ((val?val[next]:1.0) > key)
    { sorted[i] = sorted[(*smaller)];
      sorted[(*smaller)++] = sorted[(*equal)];
      sorted[(*equal)++] = next;
    }
    else if ((val?val[next]:1.0) == key)
    { sorted[i] = sorted[(*smaller)];
      sorted[(*smaller)++] = next;
  } }
}

void quicksort_dec_float (
int *sorted, float* val, int start, int end
)
{ int  equal, smaller;

  if (start < end)
  { quickpart_dec_float (sorted,val,start,end,&equal,&smaller);
    quicksort_dec_float (sorted,val,start,equal-1);
    quicksort_dec_float (sorted,val,smaller,end);
  }
}

/****************************************************************************/

static void quickpart_dec_float_int (
int *sorted, float *val1, int *val2, int start, int end, int *equal, int *smaller
)
{ int	i, next, key2, key2_next;
  float key1, key1_next;

  i = (end+start)/2;
  key1 = val1?val1[sorted[i]]:1.0;
  key2 = val2?val2[sorted[i]]:1;

  (*equal) = (*smaller) = start;
  for (i=start; i<=end; i++)
  { next = sorted[i];
    key1_next = val1?val1[next]:1.0;
    key2_next = val2?val2[next]:1;
    if (key1_next>key1 || (key1_next==key1 && key2_next>key2))
    { sorted[i] = sorted[(*smaller)];
      sorted[(*smaller)++] = sorted[(*equal)];
      sorted[(*equal)++] = next;
    }
    else if (key1_next==key1 && key2_next==key2)
    { sorted[i] = sorted[(*smaller)];
      sorted[(*smaller)++] = next;
  } }
}

void quicksort_dec_float_int (
int *sorted, float* val1, int *val2, int start, int end
)
{ int  equal, smaller;

  if (start < end)
  { quickpart_dec_float_int (sorted,val1,val2,start,end,&equal,&smaller);
    quicksort_dec_float_int (sorted,val1,val2,start,equal-1);
    quicksort_dec_float_int (sorted,val1,val2,smaller,end);
  }
}

/****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
