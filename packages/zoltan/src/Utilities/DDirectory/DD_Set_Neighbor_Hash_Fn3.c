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



/*  NOTE: See file, README, for associated documentation. (RTH) */


static unsigned int dd_nh3 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc) ;

static int remainder ;
static int average ;
static int breakpt ;
static int total_ ;


/*************  Zoltan_DD_Set_Hash_Fn3() ***********************/
/*
**  These routines associate the first n=groupsize GIDs to proc 0, the
**  next n to proc 1, etc.  It assumes the GIDs are consecutive numbers.
**  It assumes that GIDs primarily stay near their original owner. The
**  GID length is assumed to be 1. GIDs outside of range are evenly
**  distributed among the processors via modulo(nproc).  This method
**  is designed for Trilinos/Epetra linear map.
*/


int Zoltan_DD_Set_Neighbor_Hash_Fn3 (
 Zoltan_DD_Directory *dd,          /* directory state information */
 int total)                        /* total number of GIDS */
   {
   char *yo = "Zoltan_DD_Set_Hash_Fn3" ;

   if (dd == NULL || total < 1)
      {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument") ;
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   total_    = total ;
   average   = total / dd->nproc ;
   remainder = total % dd->nproc ;
   breakpt   = (average+1) * remainder ;

   dd->hash    = dd_nh3 ;
   dd->cleanup = NULL ;                 /* no need to free anything */

   return ZOLTAN_DD_NORMAL_RETURN ;
   }



static unsigned int dd_nh3 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc)
   {
   int id = (signed) *gid ;
   if (id < breakpt)
      return  id/(average+1) ;
   if (id < total_)
      return remainder + (id-breakpt)/average ;

   return 0 ;                    /* error, gid is out of range */
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
