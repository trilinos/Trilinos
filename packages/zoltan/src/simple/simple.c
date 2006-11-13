/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
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
#include <memory.h>
#include "zz_const.h"
#include "zz_rand.h"
#include "params_const.h"
#include "all_allo_const.h"

/* Simple partitioning method.
 * Consider all objects as a linear sequence and do
 * a simple parallel sequence (1d) partitioning stategy.
 */


int Zoltan_Simple(
  ZZ *zz,                       /* The Zoltan structure.                     */
  float *part_sizes,            /* Input:  Array of size zz->LB.Num_Global_Parts
                                   containing the percentage of work to be
                                   assigned to each partition.               */
  int *num_import,              /* Return -1. Simple uses only export lists. */
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
  int part, max_export;
  int wtflag;
  double wtsum, scansum[zz->Num_Proc+1];
  ZOLTAN_ID_PTR global_ids = NULL;
  ZOLTAN_ID_PTR local_ids = NULL; 
  int *parts = NULL;
  float *wgts = NULL;
  static char *yo = "Zoltan_Simple";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* No import lists computed. */
  *num_import = -1;

  /* Get list of local objects. */
  wtflag = (zz->Obj_Weight_Dim>0 ? 1 : 0);
  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, 0,
                             &wgts, &parts);

  /* TODO: Estimate number of objects to export. */
  max_export = num_obj;

  /* Allocate export lists. */
  *export_global_ids = *export_local_ids = NULL;
  *export_procs = *export_to_part = NULL;
  if (max_export > 0) {
    if (!Zoltan_Special_Malloc(zz, (void **)export_global_ids, max_export,
                               ZOLTAN_SPECIAL_MALLOC_GID)
     || !Zoltan_Special_Malloc(zz, (void **)export_local_ids, max_export,
                               ZOLTAN_SPECIAL_MALLOC_LID)
     || !Zoltan_Special_Malloc(zz, (void **)export_procs, max_export,
                               ZOLTAN_SPECIAL_MALLOC_INT)
     || !Zoltan_Special_Malloc(zz, (void **)export_to_part, max_export,
                               ZOLTAN_SPECIAL_MALLOC_INT)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  /* Sum up local object weights. */
  if (wtflag){
    wtsum = 0.0;
    for (i=0; i<num_obj; i++)
      wtsum += wgts[i];
  }
  else
    wtsum = num_obj;

  /* Cumulative global wtsum */
  MPI_Allgather(&wtsum, 1, MPI_DOUBLE, &scansum[1], 1, MPI_DOUBLE, 
                zz->Communicator);
  /* scansum = sum of weights on lower processors, excluding self. */
  scansum[0] = 0.;
  for (i=1; i<=zz->Num_Proc; i++)
    scansum[i] += scansum[i-1]; 

  /* Overwrite part_sizes with cumulative sum (inclusive) part_sizes. */
  /* A cleaner way is to make a copy, but this works. */
  for (i=1; i<zz->LB.Num_Global_Parts; i++)
    part_sizes[i] += part_sizes[i-1]; 

  /* Loop over objects and assign partition. */
  count = 0;
  part = 0;
  wtsum = scansum[zz->Proc];
  for (i=0; i<num_obj; i++){
    /* wtsum is now sum of all lower-ordered object */
    /* determine new partition number for this object, 
       using the "center of gravity" */
    while (part<zz->LB.Num_Global_Parts-1 && (wtsum+0.5*(wtflag? wgts[i]: 1.0)) 
           > part_sizes[part]*scansum[zz->Num_Proc])
      part++;
    if (part != parts[i]){
      /* export_global_ids[count] = global_ids[i]; */
      ZOLTAN_SET_GID(zz, &((*export_global_ids)[count*zz->Num_GID]),
                     &global_ids[i*zz->Num_GID]);
      if (local_ids)
        /* export_local_ids[count] = local_ids[i]; */
        ZOLTAN_SET_LID(zz, &((*export_local_ids)[count*zz->Num_LID]),
                       &local_ids[i*zz->Num_LID]);
      /* Set new partition number. */
      (*export_to_part)[count] = part;
      /* Processor is derived from partition number. */
      (*export_procs)[count] = Zoltan_LB_Part_To_Proc(zz, 
                     (*export_to_part)[count], &global_ids[i*zz->Num_GID]);

/*
      printf("Debug: export gid %u to part %d and proc %d.\n", (*export_global_ids)[count], (*export_to_part)[count], (*export_procs)[count]);
*/
      ++count;
    }
    wtsum += (wtflag? wgts[i] : 1.0);
  }
  (*num_export) = count;

End:
  /* Free local memory, but not export lists. */
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&parts);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

