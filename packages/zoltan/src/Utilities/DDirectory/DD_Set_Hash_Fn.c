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


/*  NOTE: See file, README, for associated documentation. (RTH) */


static unsigned int dd_hash_user (
  ZOLTAN_ID_PTR gid, int gid_length, 
  unsigned int nproc,
  void *data,
  ZOLTAN_HASH_FN *fn)
{
  return (*fn)(gid, gid_length, nproc);
}

/*************  Zoltan_DD_Set_Hash_Fn()  ***********************/


int Zoltan_DD_Set_Hash_Fn (
 Zoltan_DD_Directory *dd,              /* directory state information */
 ZOLTAN_HASH_FN *hash)
{
   char *yo = "Zoltan_DD_Set_Hash_Fn";

   /* input sanity checking */
   if (dd == NULL || hash == NULL)  {
      ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument");
      return ZOLTAN_FATAL ;
   }

   dd->hash = (DD_Hash_fn*)dd_hash_user;
   dd->hashdata = NULL;
   dd->hashfn = hash;
   dd->cleanup = (DD_Cleanup_fn*) NULL; 

   if (dd->debug_level > 0)
      ZOLTAN_PRINT_INFO (dd->my_proc, yo, "Successful");

   return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
