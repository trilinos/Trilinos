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


struct dd_nh3_struct {
  int remainder;
  int average;
  int breakpt;
  int total_;};


static unsigned int dd_nh3 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc, struct dd_nh3_struct* hashdata, ZOLTAN_HASH_FN *);



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
  char *yo = "Zoltan_DD_Set_Hash_Fn3";
  struct dd_nh3_struct *hashdata;

  if (dd == NULL || total < 1) {
    ZOLTAN_PRINT_ERROR (0, yo, "Invalid input argument");
    return ZOLTAN_FATAL;
  }

  hashdata = (struct dd_nh3_struct*) ZOLTAN_MALLOC(sizeof(struct dd_nh3_struct));
  if (hashdata == NULL) {
    ZOLTAN_PRINT_ERROR (0, yo, "Memory error");
    return ZOLTAN_FATAL;
  }

  hashdata->total_    = total;
  hashdata->average   = total / dd->nproc;
  hashdata->remainder = total % dd->nproc;
  hashdata->breakpt   = (hashdata->average+1) * hashdata->remainder;

  dd->hash    = (DD_Hash_fn*) &dd_nh3;
  dd->hashdata    = hashdata;
  dd->hashfn  = NULL;
  dd->cleanup = (DD_Cleanup_fn*)&Zoltan_DD_default_cleanup;

  return ZOLTAN_OK;
}


static unsigned int dd_nh3 (ZOLTAN_ID_PTR gid, int gid_length,
 unsigned int nproc, struct dd_nh3_struct * hashdata, ZOLTAN_HASH_FN *fn)
{
  int id = (signed) *gid;
  if (id < hashdata->breakpt)
    return  id/(hashdata->average+1);
  if (id < hashdata->total_)
    return hashdata->remainder + (id-hashdata->breakpt)/hashdata->average;

  return 0;                    /* error, gid is out of range */
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
