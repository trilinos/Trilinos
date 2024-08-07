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


#include <stdio.h>
#include <memory.h>
#include "zz_util_const.h"
#include "zz_rand.h"
#include "params_const.h"
#include "all_allo_const.h"

/* local function prototypes */
static void cyclic_part(ZZ *zz, int num_obj, int wtflag, float *wgts, 
            float *part_sizes, int *newparts);

/* Cyclic (round-robin) partitioning method.
   Based on the Block code; could be consolidated.
 */

int Zoltan_Cyclic(
  ZZ *zz,                       /* The Zoltan structure.                     */
  float *part_sizes,            /* Input:  Array of size zz->LB.Num_Global_Parts
                                   containing the percentage of work to be
                                   assigned to each partition.               */
  int *num_import,              /* Return -1. We use only export lists. */
  ZOLTAN_ID_PTR *import_global_ids, /* Not used. */
  ZOLTAN_ID_PTR *import_local_ids,  /* Not used. */
  int **import_procs,           /* Not used. */
  int **import_to_part,         /* Not used. */
  int *num_export,              /* Output: Number of objects to export. */
  ZOLTAN_ID_PTR *export_global_ids, /* Output: GIDs to export. */
  ZOLTAN_ID_PTR *export_local_ids,  /* Output: LIDs to export. */
  int **export_procs,           /* Output: Processsors to export to. */
  int **export_to_part          /* Output: Partitions to export to. */
)
{
  int ierr = ZOLTAN_OK;
  int i, count, num_obj;
  int wtflag = 0;
  ZOLTAN_ID_PTR global_ids = NULL;
  ZOLTAN_ID_PTR local_ids = NULL; 
  int *parts = NULL;
  int *newparts = NULL;
  float *wgts = NULL;
  static char *yo = "Zoltan_Cyclic";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* No import lists computed. */
  *num_import = -1;
  *export_global_ids = *export_local_ids = NULL;
  *export_procs = *export_to_part = NULL;

  /* Get list of local objects. */
  if (zz->Obj_Weight_Dim > 1) {
    ierr = ZOLTAN_FATAL;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                      "OBJ_WEIGHT_DIM > 1 not supported by LB_METHOD BLOCK.");
    goto End;
  }
  wtflag = (zz->Obj_Weight_Dim>0 ? 1 : 0);
  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, wtflag,
                             &wgts, &parts);

  /* Compute the new partition numbers. */
  newparts = (int *) ZOLTAN_MALLOC(num_obj * sizeof(int));
  if (num_obj && (!newparts)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  cyclic_part(zz, num_obj, wtflag, wgts, part_sizes, newparts);

  /* Check how many partition numbers changed. */
  count=0;
  for (i=0; i<num_obj; i++){
    if (newparts[i] != parts[i])
      ++count;
  }
  (*num_export) = count;

  /* Allocate export lists. */
  if ((*num_export) > 0) {
    if (!Zoltan_Special_Malloc(zz, (void **)export_global_ids, (*num_export),
                               ZOLTAN_SPECIAL_MALLOC_GID)
     || !Zoltan_Special_Malloc(zz, (void **)export_local_ids, (*num_export),
                               ZOLTAN_SPECIAL_MALLOC_LID)
     || !Zoltan_Special_Malloc(zz, (void **)export_procs, (*num_export),
                               ZOLTAN_SPECIAL_MALLOC_INT)
     || !Zoltan_Special_Malloc(zz, (void **)export_to_part, (*num_export),
                               ZOLTAN_SPECIAL_MALLOC_INT)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  /* Loop over objects and fill export lists. */
  count=0;
  for (i=0; i<num_obj; i++){
    if (newparts[i] != parts[i]){
      /* export_global_ids[count] = global_ids[i]; */
      ZOLTAN_SET_GID(zz, &((*export_global_ids)[count*zz->Num_GID]),
                     &global_ids[i*zz->Num_GID]);
      if (local_ids)
        /* export_local_ids[count] = local_ids[i]; */
        ZOLTAN_SET_LID(zz, &((*export_local_ids)[count*zz->Num_LID]),
                       &local_ids[i*zz->Num_LID]);
      /* Set new partition number. */
      (*export_to_part)[count] = newparts[i];
      /* Processor is derived from partition number. */
      (*export_procs)[count] = Zoltan_LB_Part_To_Proc(zz, 
                     (*export_to_part)[count], &global_ids[i*zz->Num_GID]);

      ++count;
    }
  }

End:
  /* Free local memory, but not export lists. */
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&parts);
  ZOLTAN_FREE(&newparts);
  if (wtflag) ZOLTAN_FREE(&wgts);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}


/* Core function to compute partition numbers. */
/* Pure cyclic. Ignore weights!                */
/* Output: newparts contains the new partition numbers. */

static void cyclic_part(ZZ *zz, int num_obj, int wtflag, float *wgts, 
            float *part_sizes, int *newparts)
{
  ZOLTAN_GNO_TYPE n, scan_sum;
  MPI_Datatype gno_mpi_type;
  int i, part=0;
  const int k = zz->LB.Num_Global_Parts;

  gno_mpi_type = Zoltan_mpi_gno_type();

  /* Compute offset for my proc */
  n = (ZOLTAN_GNO_TYPE)num_obj;
  MPI_Scan(&n, &scan_sum, 1, gno_mpi_type, MPI_SUM, zz->Communicator);
  part = (int)(scan_sum - num_obj); 
  part = (part % k);

  /* Loop over objects and assign parts. */
  for (i=0; i<num_obj; i++){
    newparts[i] = part;
    if (++part >= k)
      part= 0;
  }
 
}
