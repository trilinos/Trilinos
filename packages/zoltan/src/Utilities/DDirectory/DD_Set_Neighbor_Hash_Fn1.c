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

struct dd_nh1_struct {
 int max_gid ;
 int groupsize ;
};

static unsigned int dd_nh1 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc, struct dd_nh1_struct* hashdata, ZOLTAN_HASH_FN *fn) ;



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
   char *yo = "Zoltan_DD_Set_Neighbor_Hash_Fn1";
   struct dd_nh1_struct *hashdata;

   if (dd == NULL) {
     ZOLTAN_PRINT_ERROR (0, yo, "NULL DDirectory pointer");
     return ZOLTAN_FATAL;
   }
   if (size < 1) {
     ZOLTAN_PRINT_WARN (0, yo, "Invalid input argument; size < 1");
     return ZOLTAN_WARN;
   }

   hashdata = (struct dd_nh1_struct*) ZOLTAN_MALLOC(sizeof(struct dd_nh1_struct));
   if (hashdata == NULL) {
     ZOLTAN_PRINT_ERROR (0, yo, "Memory error");
     return ZOLTAN_FATAL;
   }

   hashdata->groupsize   = size;
   dd->hash    = (DD_Hash_fn*) &dd_nh1;
   dd->hashfn  = NULL;
   dd->hashdata = hashdata;
   dd->cleanup = (DD_Cleanup_fn*) &Zoltan_DD_default_cleanup;
   hashdata->max_gid = size * dd->nproc;     /* larger GIDs out of range */

   return ZOLTAN_OK;
   }



static unsigned int dd_nh1 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc, struct dd_nh1_struct* hashdata, ZOLTAN_HASH_FN *fn)
   {
   int id = (signed) *gid;
   return (unsigned int) ((id < hashdata->max_gid) ? (id / hashdata->groupsize)
                                                   : (id % (int)nproc));
   }

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
