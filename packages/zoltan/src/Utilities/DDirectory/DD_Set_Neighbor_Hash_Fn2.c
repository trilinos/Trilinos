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

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>

#include "DD.h"



/*  NOTE: See file, README, for associated documentation. (RTH)  */


typedef struct
   {
   int high ;
   int low ;
   int proc ;
   } Range_Info ;


static unsigned int dd_nh2 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc) ;

static int compare_sort   (const void *a, const void *b) ;
static int compare_search (const void *a, const void *b) ;

static void dd_nh2_cleanup (void) ;

static Range_Info *ptr ;
static int         nproc ;
static int         low_limit ;
static int         high_limit ;
static int         debug_level ;
static int         count ;


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
   int i ;
   char *yo = "Zoltan_DD_Set_Hash_Fn2" ;

   if (dd == NULL || proc == NULL || low == NULL || high == NULL)
      {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument") ;
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   /* register functions for automatic invocation */
   dd->hash    = dd_nh2 ;
   dd->cleanup = dd_nh2_cleanup ;

   /* malloc and initialize storage for range information structures */
   ptr = (Range_Info *)  ZOLTAN_MALLOC (n * sizeof (Range_Info)) ;
   if (ptr == NULL)
      {
      ZOLTAN_PRINT_ERROR (dd->my_proc, yo, "Unable to Malloc range info") ;
      return ZOLTAN_DD_MEMORY_ERROR ;
      }
   for (i = 0 ;  i < n ; i++)
      {
      ptr[i].high = high[i] ;
      ptr[i].low  = low [i] ;
      ptr[i].proc = (proc[i] < n) ? proc[i] : 0 ;
      }

   /* do not assume user lists were ordered */
   qsort (ptr, n, sizeof (Range_Info), compare_sort) ;

   low_limit   = ptr[0].low ;
   high_limit  = ptr[n-1].high ;
   debug_level = dd->debug_level ;
   count       = n ;
   nproc       = dd->nproc ;

   return ZOLTAN_DD_NORMAL_RETURN ;
   }



static unsigned int dd_nh2 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int junk)
   {
   Range_Info *p ;
   char *yo = "dd_ny2" ;
   int id = (signed) *gid ;

   /* check if gid is out of range */
   if (id > high_limit || id < low_limit)
      {
      return id % nproc ;
      }

   p = (Range_Info *) bsearch (gid, ptr, count, sizeof (Range_Info),
    compare_search) ;

   if (p == NULL)              /* shouldn't happen */
      {
      if (debug_level > 1)
         ZOLTAN_PRINT_ERROR (0, yo, "C function bsearch returned NULL.") ;
      return id % nproc ;
      }

   return p->proc ;
   }







static int compare_sort (const void *a, const void *b)
    {
    if (((Range_Info *) a)->low < ((Range_Info *) b)->low) return -1 ;
    if (((Range_Info *) a)->low > ((Range_Info *) b)->low) return  1 ;
    else return 0 ;
    }


static int compare_search (const void *a, const void *b)
    {
    int temp = (signed) *((ZOLTAN_ID_TYPE *) a) ;
    if (temp < ((Range_Info *) b)->low)  return -1 ;
    if (temp > ((Range_Info *) b)->high) return  1 ;
    else return 0 ;
    }





static void dd_nh2_cleanup (void)
   {
   ZOLTAN_FREE (&ptr) ;
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
