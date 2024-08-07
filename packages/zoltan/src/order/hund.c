// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "key_params.h"
#include "params_const.h"
#include "ha_const.h"
#include "zoltan_dd.h"
#include "phg_tree.h"
#include "order_const.h"

#define SIZE_PART 1000

/* #define SIMPLE_HUND */

#define SIMPLE_HUND

#ifndef SIMPLE_HUND
/* Removed CCOLAMD from BSD version. 
int Zoltan_CColAMD(
  ZZ *zz,               
  struct Zoltan_DD_Struct *dd_constraint,
  int nPart,
  int *num_obj,
  ZOLTAN_ID_PTR *gid,
  ZOLTAN_ID_PTR *rank
  );
*/
#else /* SIMPLE_HUND */
static int
HUND_Order_simple(ZZ* zz, int num_local_gid, ZOLTAN_ID_PTR part, int numParts, int* sizeParts, ZOLTAN_ID_PTR dperm); /* result is stored into a DD */
#endif

int Zoltan_HUND(
  ZZ *zz,               /* Zoltan structure */
  int num_gid_entries, /* # of entries for a global id */
  int num_obj,		/* Number of objects to order */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
                        /* The application must allocate enough space */
  ZOLTAN_ID_PTR rank  /* rank[i] is the rank of gids[i] */
  /* int *iperm           /\* iperm[rank[i]]=i, only for sequential ordering *\/ */
)
{
  char *yo = "Zoltan_HUND";
  int ierr = ZOLTAN_OK;
  int changes;
  int numGidEntries, numLidEntries;
  int numImport;
  int numExport;
  ZOLTAN_ID_TYPE *part = NULL;
  ZOLTAN_ID_TYPE *dperm = NULL;
  Zoltan_PHG_LB_Data *data = NULL;
  ZOLTAN_ID_PTR local_gid;
  int num_local_gid;
  ZOLTAN_ID_PTR gidImport, lidImport, gidExport, lidExport;
  int *procImport, *partImport, *procExport, *partExport;
  int numPart = 1;
  int numLocObj, numGlbObj;
  char partArg[MAX_PARAM_STRING_LEN];

  ZOLTAN_TRACE_ENTER(zz, yo);

  numLocObj = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  CHECK_IERR;

  MPI_Allreduce(&numLocObj, &numGlbObj, 1, MPI_INT, MPI_SUM, zz->Communicator);
  numPart = numGlbObj/SIZE_PART;

  /* Round up to the next highest power of 2 */
  /* From http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2 */
  numPart--;
  numPart |= numPart >> 1;
  numPart |= numPart >> 2;
  numPart |= numPart >> 4;
  numPart |= numPart >> 8;
  numPart |= numPart >> 16;

  numPart++; /* aprun -n 4 ./pddrive -r 1 -c 2 */

  /* First, we have to compute an hypergraph partitioning */
  Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
  sprintf (partArg, "%d", numPart);
  Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", partArg);
  Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.2");
  Zoltan_Set_Param(zz, "PHG_KEEP_TREE", "1"); /* We need the separator tree */
  Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTITION_ASSIGNMENTS");
  /* We want to minimize the number of columns in the separator,
     not the number of nnz */
  Zoltan_Set_Param(zz, "PHG_CUT_OBJECTIVE", "HYPEREDGES");

  ierr = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
			     &changes,        /* 1 if partitioning was changed, 0 otherwise */
			     &numGidEntries,  /* Number of integers used for a global ID */
			     &numLidEntries,  /* Number of integers used for a local ID */
			     &numImport,      /* Number of vertices to be sent to me */
        &gidImport,  /* Global IDs of vertices to be sent to me */
        &lidImport,   /* Local IDs of vertices to be sent to me */
        &procImport,    /* Process rank for source of each incoming vertex */
        &partImport,   /* New partition for each incoming vertex */
        &numExport,      /* Number of vertices I must send to other processes*/
        &gidExport,  /* Global IDs of the vertices I must send */
        &lidExport,   /* Local IDs of the vertices I must send */
        &procExport,    /* Process to which I send each of the vertices */
        &partExport);  /* Partition to which each vertex will belong */

  ZOLTAN_FREE(&gidExport);
  ZOLTAN_FREE(&lidExport);
  ZOLTAN_FREE(&procExport);

  /* Second compute ordering from partitioning */
  data = (Zoltan_PHG_LB_Data*)zz->LB.Data_Structure;
  if (zz->Proc == 0) {
    int i;
    fprintf (stderr, "Block sizes (%d - %d)\n", data->numParts, numPart);
    for (i = 1 ; (i <= data->numParts) && (i < 10) ; ++i) {
      fprintf (stderr, "%d ", data->sizeParts[data->numParts-i]);
    }
    fprintf (stderr, "\n");
  }

  ierr = Zoltan_DD_GetLocalKeys(data->ddHedge, &local_gid, &num_local_gid);
  CHECK_IERR;

#ifdef SIMPLE_HUND
  part = (ZOLTAN_ID_TYPE*) ZOLTAN_MALLOC(num_local_gid*sizeof(ZOLTAN_ID_TYPE));
  if (num_local_gid && part == NULL) MEMORY_ERROR;

  Zoltan_DD_Find (data->ddHedge, local_gid, part, NULL, NULL, num_local_gid, NULL);

  dperm = (ZOLTAN_ID_TYPE*) ZOLTAN_MALLOC(num_local_gid*sizeof(ZOLTAN_ID_TYPE));
  if (num_local_gid && dperm == NULL) MEMORY_ERROR;

  ierr = HUND_Order_simple(zz,  num_local_gid, part, data->numParts, data->sizeParts, dperm);
  CHECK_IERR;
#else /* SIMPLE_HUND */
  /* Removed CCOLAMD from BSD version.
  ierr = Zoltan_CColAMD(zz, data->ddHedge, data->numParts, &num_local_gid, &local_gid, &dperm);
  CHECK_IERR;
  */
#endif /* SIMPLE_HUND */

  Zoltan_DD_Update (data->ddHedge, local_gid, (ZOLTAN_ID_PTR)dperm, NULL,  NULL /* part */, num_local_gid);
  Zoltan_DD_Find (data->ddHedge, gids, rank, NULL, NULL, num_obj, NULL);

  Zoltan_DD_Destroy(&data->ddHedge);

 End:
  ZOLTAN_FREE(&dperm);
  ZOLTAN_FREE(&part);
  ZOLTAN_FREE(&local_gid);
  ZOLTAN_FREE(&data->sizeParts);
  ZOLTAN_FREE(&partExport);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

#ifdef SIMPLE_HUND
static int
HUND_Order_simple(ZZ* zz, int num_local_gid, ZOLTAN_ID_TYPE * part, int numParts, int* sizeParts, ZOLTAN_ID_TYPE* dperm) /* result is stored into a DD */
{
  char *yo = "HUND_Order_simple";
  int ierr = ZOLTAN_OK;
  int *local_size = NULL;
  int *offset = NULL;
  int i;

  ZOLTAN_TRACE_ENTER(zz, yo);

  local_size = (int*) ZOLTAN_CALLOC(numParts,sizeof(int));
  if (numParts && local_size == NULL) MEMORY_ERROR;

  for (i = 0 ; i < num_local_gid ; ++i)
    local_size[part[i]]++;

  offset = (int*) ZOLTAN_CALLOC(numParts,sizeof(int));
  if (numParts && offset == NULL) MEMORY_ERROR;

  /* Compute offset for numbering */
  MPI_Scan (local_size, offset, numParts, MPI_INT, MPI_SUM, zz->Communicator);

  for (i = 0 ; i < numParts ; ++i) {
    offset[i] -= local_size[i];
  }

  local_size[0] = 0;
  for (i=1 ; i < numParts ; ++i) {
    local_size[i] = sizeParts[i-1] + local_size[i-1];
  }

  for (i = 1 ; i < numParts ; ++i) {
    offset[i] += local_size[i];
  }

  ZOLTAN_FREE(&local_size);

  for (i = 0 ; i < num_local_gid ; ++i)
    dperm[i] = (offset[part[i]]++);
  ZOLTAN_FREE(&offset);

  End:
  ZOLTAN_FREE(&offset);
  ZOLTAN_FREE(&local_size);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}
#endif

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
