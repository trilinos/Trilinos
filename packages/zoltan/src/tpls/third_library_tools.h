// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __THIRD_LIBRARY_TOOLS_H
#define __THIRD_LIBRARY_TOOLS_H

#include <limits.h>
#include "zoltan_comm.h"
#include "third_library_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Misc. local constants */
#define CHUNKSIZE 20  /* Number of nodes to allocate in initial chunk. */
#define REALLOC_FACTOR 1.5  /* Increase size by this factor if too small. */


/* Data structures used in ParMetis interface routines */
/* An array of this data structure works with a parallel array of
 * ZOLTAN_ID_PTR called proc_list_nbor containing the global IDs of the
 * neighboring object.
 * This separate array is needed to prevent individual mallocs of
 * neighboring global IDs.
 */
struct Edge_Info {
  ZOLTAN_ID_PTR my_gid;  /* Pointer to the Global id of local vtx */
  int my_gno;        /* Global number of local vtx */
  int nbor_proc;     /* Proc id for the neighboring proc */
  int *adj;          /* Pointer to adjcny array */
};

struct Hash_Node {
  ZOLTAN_ID_PTR gid;     /* Pointer to a Global id */
  int gno;           /* Global number */
  struct Hash_Node * next;
};

#ifdef __cplusplus
}
#endif

#endif /* __THIRD_LIBRARY_TOOLS_H */
