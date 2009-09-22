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


/*  NOTE: See file, README, for associated documentation. (RTH) */

static unsigned int dd_hash_user (ZOLTAN_ID_PTR, int, unsigned int,
				  unsigned int (*data) (ZOLTAN_ID_PTR, int, unsigned int));

/*************  Zoltan_DD_Set_Hash_Fn()  ***********************/


int Zoltan_DD_Set_Hash_Fn (
 Zoltan_DD_Directory *dd,              /* directory state information */
 unsigned int (*hash) (ZOLTAN_ID_PTR, int, unsigned int))
   {
   char *yo = "Zoltan_DD_Set_Hash_Fn";

   /* input sanity checking */
   if (dd == NULL || hash == NULL)  {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument");
      return ZOLTAN_FATAL ;
   }

   dd->hash = (DD_Hash_fn*)dd_hash_user;
   dd->data = (void*)hash;
   dd->cleanup = (DD_Cleanup_fn*) &Zoltan_DD_default_cleanup;

   if (dd->debug_level > 0)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "Successful");

   return ZOLTAN_OK;
   }

static unsigned int dd_hash_user (ZOLTAN_ID_PTR gid, int gid_length, unsigned int nproc,
				  unsigned int (*data) (ZOLTAN_ID_PTR, int, unsigned int)) {
  return (*data)(gid, gid_length, nproc);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
