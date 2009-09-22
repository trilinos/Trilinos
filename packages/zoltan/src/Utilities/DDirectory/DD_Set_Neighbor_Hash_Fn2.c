/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"


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
 unsigned int nproc, struct dd_nh2_struct* data);

static int compare_sort   (const void *a, const void *b);
static int compare_search (const void *a, const void *b);

static void dd_nh2_cleanup (struct dd_nh2_struct* data);



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
   struct dd_nh2_struct *data;


   if (dd == NULL || proc == NULL || low == NULL || high == NULL)  {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument");
      return ZOLTAN_FATAL;
   }

  data = (struct dd_nh2_struct*) ZOLTAN_MALLOC(sizeof(struct dd_nh2_struct));
  if (data == NULL) {
    ZOLTAN_PRINT_ERROR (0, yo, "Memory error");
    return ZOLTAN_FATAL;
  }


   /* register functions for automatic invocation */
   dd->hash    = (DD_Hash_fn*) &dd_nh2;
   dd->cleanup = (DD_Cleanup_fn*)&dd_nh2_cleanup;
   dd->data = data;

   /* malloc and initialize storage for range information structures */
   data->ptr = (Range_Info*)  ZOLTAN_MALLOC (n * sizeof (Range_Info));
   if (data->ptr == NULL)  {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to Malloc range info");
      return ZOLTAN_MEMERR;
   }
   for (i = 0;  i < n; i++)  {
      data->ptr[i].high = high[i] ;
      data->ptr[i].low  = low [i] ;
      data->ptr[i].proc = (proc[i] < n) ? proc[i] : 0;
   }

   /* do not assume user lists were ordered */
   qsort (data->ptr, n, sizeof (Range_Info), compare_sort);

   data->low_limit   = data->ptr[0].low;
   data->high_limit  = data->ptr[n-1].high;
   data->debug_level = dd->debug_level;
   data->count       = n;
   data->nproc       = dd->nproc;

   return ZOLTAN_OK;
   }



static unsigned int dd_nh2 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int junk, struct dd_nh2_struct* data)
   {
   Range_Info *p;
   char *yo = "dd_ny2";
   int id = (signed) *gid;

   /* check if gid is out of range */
   if (id > data->high_limit || id < data->low_limit)
      return id % data->nproc;

   p = (Range_Info*) bsearch (gid, data->ptr, data->count, sizeof (Range_Info),
    compare_search);

   if (p == NULL) {             /* shouldn't happen */
      if (data->debug_level > 1)
         ZOLTAN_PRINT_ERROR (0, yo, "C function bsearch returned NULL.") ;
      return id % data->nproc;
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





static void dd_nh2_cleanup (struct dd_nh2_struct *data)
{
  if (data == NULL) return;
  ZOLTAN_FREE (&data->ptr);
  ZOLTAN_FREE (&data);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
