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

#include "DD.h"



/*  NOTE: See file, README, for associated documentation. (RTH) */


static unsigned int dd_nh1 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc) ;

static int max_gid ;
static int groupsize ;


/*************  Zoltan_DD_Set_Hash_Fn1() ***********************/
/*
**  These routines associate the first n=groupsize GIDs to proc 0, the
**  next n to proc 1, etc.  It assumes the GIDs are consecutive numbers.
**  It assumes that GIDs primarily stay near their original owner. The
**  GID length is assumed to be 1. GIDs outside of range are evenly
**  distributed among the processors via modulo(nproc).
*/


int Zoltan_DD_Set_Neighbor_Hash_Fn1 (
 Zoltan_DD_Directory *dd,          /* directory state information */
 int size)                         /* number of reserved GIDs per CPU */
   {
   char *yo = "Zoltan_DD_Set_Hash_Fn1" ;

   if (dd == NULL || size < 1)
      {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument") ;
      return ZOLTAN_DD_INPUT_ERROR ;
      }

   groupsize   = size ;
   dd->hash    = dd_nh1 ;
   dd->cleanup = NULL ;                 /* no need to free anything */

   max_gid     = size * dd->nproc ;     /* larger GIDs out of range */

   return ZOLTAN_DD_NORMAL_RETURN ;
   }



static unsigned int dd_nh1 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc)
   {
   return (*gid < max_gid) ? (*gid / groupsize) : (*gid % nproc) ;
   }
