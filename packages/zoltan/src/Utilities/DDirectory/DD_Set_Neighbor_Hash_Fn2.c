// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include <stdio.h>
#include <stdlib.h>

#include "zoltan_dd_const.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif



typedef struct
   {
   int high;
   int low;
   int proc;
   } Range_Info;

struct dd_nh2_struct {
  Range_Info *ptr;
  int nproc;
  int low_limit;
  int high_limit;
  int debug_level;
  int count;
};

static unsigned int dd_nh2 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc, struct dd_nh2_struct* hashdata, ZOLTAN_HASH_FN *fn);

static int compare_sort   (const void *a, const void *b);
static int compare_search (const void *a, const void *b);

static void dd_nh2_cleanup (struct dd_nh2_struct* hashdata);



/*************  Zoltan_DD_Set_Hash_Fn2()  ***********************/
/*
**  These routines allow the user to specify a beginning and ending GID
**  (number) per processor to associate an arbitrary GID to its orginal
**  owner processor.  It assumes the GIDs are consecutive numbers.
**  It assumes that GIDs primarily stay near their original owner. It
**  requires that the number of high, low, & proc entries is 'n'. It
**  assumes the GID length is 1.
*/



int Zoltan_DD_Set_Neighbor_Hash_Fn2 (
 Zoltan_DD_Directory *dd,     /* directory state information */
 int *proc,                   /* list of processors for following info */
 int *low,                    /* lowest GID for corresponding processor */
 int *high,                   /* highest GID for corresponding processor */
 int n)                       /* number of processors in above lists */
   {
   int i;
   char *yo = "Zoltan_DD_Set_Hash_Fn2";
   struct dd_nh2_struct *hashdata;


   if (dd == NULL || proc == NULL || low == NULL || high == NULL)  {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument");
      return ZOLTAN_FATAL;
   }

  hashdata = (struct dd_nh2_struct*) ZOLTAN_MALLOC(sizeof(struct dd_nh2_struct));
  if (hashdata == NULL) {
    ZOLTAN_PRINT_ERROR (0, yo, "Memory error");
    return ZOLTAN_FATAL;
  }


   /* register functions for automatic invocation */
   dd->hash    = (DD_Hash_fn*) &dd_nh2;
   dd->hashdata = hashdata;
   dd->hashfn    = NULL;
   dd->cleanup = (DD_Cleanup_fn*)&dd_nh2_cleanup;

   /* malloc and initialize storage for range information structures */
   hashdata->ptr = (Range_Info*)  ZOLTAN_MALLOC (n * sizeof (Range_Info));
   if (hashdata->ptr == NULL)  {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to Malloc range info");
      return ZOLTAN_MEMERR;
   }
   for (i = 0;  i < n; i++)  {
      hashdata->ptr[i].high = high[i] ;
      hashdata->ptr[i].low  = low [i] ;
      hashdata->ptr[i].proc = (proc[i] < n) ? proc[i] : 0;
   }

   /* do not assume user lists were ordered */
   qsort (hashdata->ptr, n, sizeof (Range_Info), compare_sort);

   hashdata->low_limit   = hashdata->ptr[0].low;
   hashdata->high_limit  = hashdata->ptr[n-1].high;
   hashdata->debug_level = dd->debug_level;
   hashdata->count       = n;
   hashdata->nproc       = dd->nproc;

   return ZOLTAN_OK;
   }



static unsigned int dd_nh2 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int junk, struct dd_nh2_struct* hashdata, ZOLTAN_HASH_FN *fn)
   {
   Range_Info *p;
   char *yo = "dd_ny2";
   int id = (signed) *gid;

   /* check if gid is out of range */
   if (id > hashdata->high_limit || id < hashdata->low_limit)
      return id % hashdata->nproc;

   p = (Range_Info*) bsearch (gid, hashdata->ptr, hashdata->count, sizeof (Range_Info),
    compare_search);

   if (p == NULL) {             /* shouldn't happen */
      if (hashdata->debug_level > 1)
         ZOLTAN_PRINT_ERROR (0, yo, "C function bsearch returned NULL.") ;
      return id % hashdata->nproc;
   }

   return p->proc;
   }







static int compare_sort (const void *a, const void *b)
    {
    if (((Range_Info*) a)->low < ((Range_Info *) b)->low) return -1;
    if (((Range_Info*) a)->low > ((Range_Info *) b)->low) return  1;
    else return 0 ;
    }


static int compare_search (const void *a, const void *b)
    {
    int temp = (signed) *((ZOLTAN_ID_TYPE *) a);
    if (temp < ((Range_Info*) b)->low)  return -1;
    if (temp > ((Range_Info*) b)->high) return  1;
    else return 0 ;
    }





static void dd_nh2_cleanup (struct dd_nh2_struct *hashdata)
{
  if (hashdata == NULL) return;
  ZOLTAN_FREE (&hashdata->ptr);
  ZOLTAN_FREE (&hashdata);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
