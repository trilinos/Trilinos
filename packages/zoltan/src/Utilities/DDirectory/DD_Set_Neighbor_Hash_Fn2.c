/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
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

#include "DD_Const.h"



/*  NOTE: See file, README, for associated documentation. (RTH)  */


typedef struct
   {
   int high ;
   int low ;
   int proc ;
   } Range_Info ;


static unsigned int dd_nh2 (LB_ID_PTR gid, int gid_length,
 unsigned int nproc) ;

static int compare_sort   (const void *a, const void *b) ;
static int compare_search (const void *a, const void *b) ;

static void dd_nh2_cleanup (void) ;




/*************  Zoltan_DD_Set_Hash_Fn2()  ***********************/
/*
**  These routines allow the user to specify a beginning and ending
**  number per processor to associate an arbitrary GID to its orginal
**  owner processor.  It assumes the GIDs are consecutive numbers.
**  It assumes that GIDs primarily stay near their original owner.
*/

static Range_Info *ptr ;
static int nproc ;

int Zoltan_DD_Set_Neighbor_Hash2_Fn (Zoltan_DD_Directory *dd, int *proc,
 int *low, int *high)
   {
   int i ;

   if (dd == NULL || proc == NULL || low == NULL || high == NULL)
      return ZOLTAN_DD_INPUT_ERROR ;

   dd->hash = dd_nh2 ;
   dd->cleanup = dd_nh2_cleanup ;

   nproc = dd->nproc ;

   ptr = (Range_Info *)  LB_MALLOC (dd->nproc * sizeof (Range_Info)) ;

   for (i = 0 ;  i < dd->nproc ; i++)
       {
       ptr[i].high = high[i] ;
       ptr[i].low  = low [i] ;
       ptr[i].proc = proc[i] ;
       }

   qsort (ptr, dd->nproc, sizeof (Range_Info), compare_sort) ;

   return ZOLTAN_DD_NORMAL_RETURN ;
   }

static unsigned int dd_nh2 (LB_ID_PTR gid, int gid_length,
 unsigned int junk)
   {
   Range_Info *temp ;
 
   /* check if gid is out of range */
   if (*gid > ptr[nproc-1].high || *gid < ptr[0].low)
      return nproc-1 ;

   temp = (Range_Info *) bsearch (gid, ptr, nproc,
    sizeof (Range_Info), compare_search) ;

   if (temp == NULL)              /* shouldn't happen */
      return ZOLTAN_DD_NO_PROC ;

   return temp->proc ;
   }


static int compare_sort (const void *a, const void *b)
    {
    if (((Range_Info *) a)->low < ((Range_Info *) b)->low) return -1 ;
    if (((Range_Info *) a)->low > ((Range_Info *) b)->low) return  1 ;
    else return 0 ;
    }


static int compare_search (const void *a, const void *b)
    {
    if (*((LB_ID_TYPE *) a) < ((Range_Info *) b)->low)  return -1 ;
    if (*((LB_ID_TYPE *) a) > ((Range_Info *) b)->high) return  1 ;
    else return 0 ;
    }


static void dd_nh2_cleanup (void)
   {
   LB_FREE (ptr) ;
   }
