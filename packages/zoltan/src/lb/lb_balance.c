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


#include "zz_const.h"
#include "key_params.h"
#include "timer_const.h"
#include "ha_const.h"
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for performing load balancing with Zoltan.
 *  These functions are all callable by the application.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_LB(ZZ *, int, int *, int *, int *,
  int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **,
  int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **);

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_LB_Partition(
  ZZ *zz, 
  int *changes,
  int *num_gid_entries,
  int *num_lid_entries,
  int *num_import_objs,
  ZOLTAN_ID_PTR *import_global_ids,
  ZOLTAN_ID_PTR *import_local_ids,
  int **import_procs,
  int **import_to_part,
  int *num_export_objs,
  ZOLTAN_ID_PTR *export_global_ids,
  ZOLTAN_ID_PTR *export_local_ids,
  int **export_procs,
  int **export_to_part
)
{
/*
 * Wrapper around Zoltan_LB to generate partition information.
 * Arguments correspond directly with arguments of Zoltan_LB.
 */

char *yo = "Zoltan_LB_Partition";
int ierr = ZOLTAN_OK;    /* Error code */

  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = Zoltan_LB(zz, 1, changes, num_gid_entries, num_lid_entries,
           num_import_objs, import_global_ids, import_local_ids,
           import_procs, import_to_part, 
           num_export_objs, export_global_ids, 
           export_local_ids, export_procs, export_to_part);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Balance(
  ZZ *zz, 
  int *changes,
  int *num_gid_entries,
  int *num_lid_entries,
  int *num_import_objs,
  ZOLTAN_ID_PTR *import_global_ids,
  ZOLTAN_ID_PTR *import_local_ids,
  int **import_procs,
  int *num_export_objs,
  ZOLTAN_ID_PTR *export_global_ids,
  ZOLTAN_ID_PTR *export_local_ids,
  int **export_procs
)
{
/*
 * Wrapper around Zoltan_LB for backward compatibility with
 * previous Zoltan versions.  
 * Appropriate only when (# requested partitions == # processors), uniformly
 * distributed.
 * Arguments correspond directly with arguments of Zoltan_LB.
 */

char *yo = "Zoltan_LB_Balance";
int ierr = ZOLTAN_OK;    /* Error code */
int *import_to_part = NULL;    /* Array used as dummy arg in partitioning. */
int *export_to_part = NULL;    /* Array used as dummy arg in partitioning. */

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Determine whether partition parameters were set.  Report error if
   * values are unreasonable. */
  if ((zz->LB.Num_Global_Parts_Param != -1 && 
       zz->LB.Num_Global_Parts_Param != zz->Num_Proc) ||
      (zz->LB.Num_Local_Parts_Param != -1 &&
       zz->LB.Num_Local_Parts_Param != 1)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "Non-uniform distribution of partitions over processors is specified; "
      "use Zoltan_LB_Partition.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }
  
  ierr = Zoltan_LB(zz, 0, changes, num_gid_entries, num_lid_entries,
           num_import_objs, import_global_ids, import_local_ids,
           import_procs, &import_to_part, 
           num_export_objs, export_global_ids, 
           export_local_ids, export_procs, &export_to_part);


End:
  /* Not returning import/export partition information; free it if allocated. */
  if (import_to_part != NULL) 
    Zoltan_Special_Free(zz, (void **) &import_to_part, 
                        ZOLTAN_SPECIAL_MALLOC_INT);
  if (export_to_part != NULL) 
    Zoltan_Special_Free(zz, (void **) &export_to_part, 
                        ZOLTAN_SPECIAL_MALLOC_INT);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return(ierr);
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_LB(
  ZZ *zz, 
  int include_parts,             /* Flag indicating whether to generate
                                    partition informtion;
                                    0 if called by Zoltan_LB_Balance,
                                    1 if called by Zoltan_LB_Partition.       */
  int *changes,                  /* Set to zero or one depending on if 
                                    Zoltan determines a new
                                    decomposition or not:
                                    zero - No changes to the decomposition
                                           were made by the load-balancing
                                           algorithm; migration is not needed.
                                    one  - A new decomposition is suggested
                                           by the load-balancer; migration is
                                           needed to establish the new
                                           decomposition.                     */
  int *num_gid_entries,          /* The number of array entries in a global ID;
                                    set to be the max over all processors in
                                    zz->Communicator of the parameter
                                    Num_Global_ID_Entries.                    */
  int *num_lid_entries,          /* The number of array entries in a local ID;
                                    set to be the max over all processors in
                                    zz->Communicator of the parameter
                                    Num_Local_ID_Entries.                     */
  int *num_import_objs,          /* The number of non-local objects in the
                                    processor's new decomposition.            */
  ZOLTAN_ID_PTR *import_global_ids,/* Array of global IDs for non-local objects
                                    (i.e., objs to be imported) in
                                    the processor's new decomposition.        */
  ZOLTAN_ID_PTR *import_local_ids,   /* Array of local IDs for non-local objects
                                    (i.e., objs to be imported) in
                                    the processor's new decomposition.        */
  int **import_procs,            /* Array of processor IDs for processors 
                                    currently owning non-local objects (i.e.,
                                    objs to be imported) in this processor's
                                    new decomposition.                        */
  int **import_to_part,          /* Partition to which the objects should be
                                    imported.                                 */
  int *num_export_objs,          /* The number of local objects that need to
                                    be exported from the processor to establish
                                    the new decomposition.                    */
  ZOLTAN_ID_PTR *export_global_ids,/* Array of global IDs for objects that need
                                    to be exported (assigned and sent to other
                                    processors) to establish the new 
                                    decomposition.                            */
  ZOLTAN_ID_PTR *export_local_ids,   /* Array of local IDs for objects that need
                                    to be exported (assigned and sent to other
                                    processors) to establish the new 
                                    decomposition.                            */
  int **export_procs,            /* Array of destination processor IDs for
                                    objects that need to be exported 
                                    to establish the new decomposition.       */
  int **export_to_part           /* Partition to which objects should be 
                                    exported.                                 */
)
{
/*
 * Main load-balancing routine.
 * Input:  a Zoltan structure with appropriate function pointers set.
 * Output: 
 *   changes
 *   num_import_objs
 *   import_global_ids
 *   import_local_ids
 *   import_procs
 *   import_to_part
 *   num_export_objs
 *   export_global_ids
 *   export_local_ids
 *   export_procs
 *   export_to_part
 * Return values:
 *   Zoltan error code.
 */

char *yo = "Zoltan_LB";
int gmax;    /* Maximum number of imported/exported objects 
                over all processors.                       */
int error;    /* Error code */
double start_time, end_time;
double lb_time[2] = {0.0,0.0};
char msg[256];
int comm[3],gcomm[3]; 
float *part_sizes = NULL;
int part_dim;

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (zz->Proc == zz->Debug_Proc && zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS)
    Zoltan_Print_Key_Params(zz);

  start_time = Zoltan_Time(zz->Timer);

  /* 
   * Compute Max number of array entries per ID over all processors.
   * Compute Max number of return arguments for Zoltan_LB_Balance.
   * This is a sanity-maintaining step; we don't want different
   * processors to have different values for these numbers.
   */
  comm[0] = zz->Num_GID;
  comm[1] = zz->Num_LID;
  comm[2] = zz->LB.Return_Lists;
  MPI_Allreduce(comm, gcomm, 3, MPI_INT, MPI_MAX, zz->Communicator);
  zz->Num_GID = *num_gid_entries = gcomm[0];
  zz->Num_LID = *num_lid_entries = gcomm[1];
  zz->LB.Return_Lists = gcomm[2];

  /* assume no changes */
  *changes = 0;

  *num_import_objs = *num_export_objs = 0;
  *import_global_ids = NULL;
  *import_local_ids = NULL;
  *import_procs = NULL;
  *import_to_part = NULL;
  *export_global_ids = NULL;
  *export_local_ids = NULL;
  *export_procs = NULL;
  *export_to_part = NULL;

  /*
   *  Return if this processor is not in the Zoltan structure's
   *  communicator.
   */

  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_OK);
  }

  if (zz->LB.Method == NONE) {
    if (zz->Proc == zz->Debug_Proc && zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS)
      printf("%s Balancing method selected == NONE; no balancing performed\n",
              yo);

    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_WARN);
  }

  /*
   *  Construct the heterogenous machine description.
   */

  error = Zoltan_Build_Machine_Desc(zz);

  if (error == ZOLTAN_FATAL){
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (error);
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done machine description");

  /* Since generating a new partition, need to free old mapping vector */
  ZOLTAN_FREE(&zz->LB.Remap);

  error = Zoltan_LB_Build_PartDist(zz);
  if (error != ZOLTAN_OK && error != ZOLTAN_WARN) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (error);
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    int i, np, fp;
    for (i = 0; i < zz->Num_Proc; i++) {
      Zoltan_LB_Proc_To_Part(zz, i, &np, &fp);
      printf("%d Proc_To_Part Proc %d NParts %d FPart %d\n", 
             zz->Proc, i, np, fp);
    }
  }

  /*
   * Generate partitions sizes.
   */
  part_dim = ((zz->Obj_Weight_Dim > 0) ? zz->Obj_Weight_Dim : 1);

  part_sizes = (float *) ZOLTAN_MALLOC(sizeof(float) * part_dim 
                                     * zz->LB.Num_Global_Parts);
  if (part_sizes == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_MEMERR);
  }

  /* Get partition sizes. */
  Zoltan_LB_Get_Part_Sizes(zz, zz->LB.Num_Global_Parts, part_dim,
    part_sizes);

  /*
   * Call the actual load-balancing function.
   */

  error = zz->LB.LB_Fn(zz, part_sizes,
                       num_import_objs, import_global_ids, import_local_ids,
                       import_procs, import_to_part, 
                       num_export_objs, export_global_ids, export_local_ids, 
                       export_procs, export_to_part);

  ZOLTAN_FREE(&part_sizes);

  if (error == ZOLTAN_FATAL || error == ZOLTAN_MEMERR){
    sprintf(msg, "Partitioning routine returned code %d.", error);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (error);
  }
  else if (error){
    if (zz->Debug_Level >ZOLTAN_DEBUG_NONE) {
      sprintf(msg, "Partitioning routine returned code %d.", error);
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done partitioning");

  if (*num_import_objs >= 0)
    MPI_Allreduce(num_import_objs, &gmax, 1, MPI_INT, MPI_MAX, 
                zz->Communicator);
  else /* use export data */
    MPI_Allreduce(num_export_objs, &gmax, 1, MPI_INT, MPI_MAX, 
                zz->Communicator);

  if (gmax == 0) {

    /*
     *  Decomposition was not changed by the load balancing; no migration
     *  is needed.
     */

    if (zz->Proc == zz->Debug_Proc && zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS)
      printf("%s No changes to the decomposition due to partitioning; "
             "no migration is needed.\n", yo);

    /*
     *  Reset num_import_objs and num_export_objs; don't want to return
     *  -1 for the arrays that weren't returned by ZOLTAN_LB_FN.
     */

    *num_import_objs = *num_export_objs = 0;

    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_OK);
  }

  /*
   *  Check whether we know the import data, export data, or both.
   *
   *  If we were given the import data,
   *  we know what the new decomposition should look like on the
   *  processor, but we don't know which of our local objects we have
   *  to export to other processors to establish the new decomposition.
   *  Reverse the argument if we were given the export data.
   *
   *  Unless we were given both maps, compute the inverse map.
   */

  if (*num_import_objs >= 0){
    if (*num_export_objs >= 0) {
      /* Both maps already available; nothing to do. */;

      if (zz->LB.Return_Lists == ZOLTAN_LB_NO_LISTS) {
        /* This condition should never happen!! */
        /* Methods should not return arrays if no lists are requested. */
        *num_import_objs = *num_export_objs = -1;
        Zoltan_LB_Free_Part(import_global_ids, import_local_ids, import_procs,
                            import_to_part);
        Zoltan_LB_Free_Part(export_global_ids, export_local_ids, export_procs,
                            export_to_part);
        ZOLTAN_PRINT_WARN(zz->Proc, yo, 
                      "Method returned lists, but no lists requested.");
      }
    }
    else if (zz->LB.Return_Lists == ZOLTAN_LB_ALL_LISTS || 
             zz->LB.Return_Lists == ZOLTAN_LB_EXPORT_LISTS) {
      /* Export lists are requested; compute export map */
      error = Zoltan_Invert_Lists(zz, *num_import_objs, *import_global_ids, 
                                      *import_local_ids, *import_procs,
                                      *import_to_part,
                                      num_export_objs, export_global_ids,
                                      export_local_ids, export_procs,
                                      export_to_part);
      if (error != ZOLTAN_OK && error != ZOLTAN_WARN) {
        sprintf(msg, "Error building return arguments; "
                     "%d returned by Zoltan_Compute_Destinations\n", error);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return error;
      }
      if (zz->LB.Return_Lists == ZOLTAN_LB_EXPORT_LISTS) {
        /* Method returned import lists, but only export lists were desired. */
        /* Import lists not needed; free them. */
        *num_import_objs = -1;
        Zoltan_LB_Free_Part(import_global_ids, import_local_ids, import_procs,
                            import_to_part);
      }
    }
  }
  else { /* (*num_import_objs < 0) */
    if (*num_export_objs >= 0) {
      /* Only export lists have been returned. */
      if (zz->LB.Return_Lists == ZOLTAN_LB_ALL_LISTS || 
          zz->LB.Return_Lists == ZOLTAN_LB_IMPORT_LISTS) {
        /* Compute import map */
        error = Zoltan_Invert_Lists(zz, *num_export_objs, *export_global_ids, 
                                        *export_local_ids, *export_procs,
                                        *export_to_part,
                                        num_import_objs, import_global_ids,
                                        import_local_ids, import_procs, 
                                        import_to_part);

        if (error != ZOLTAN_OK && error != ZOLTAN_WARN) {
          sprintf(msg, "Error building return arguments; "
                       "%d returned by Zoltan_Compute_Destinations\n", error);
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
          ZOLTAN_TRACE_EXIT(zz, yo);
          return error;
        }
        if (zz->LB.Return_Lists == ZOLTAN_LB_IMPORT_LISTS) {
          /* Method returned export lists, but only import lists are desired. */
          /* Export lists not needed; free them. */
          *num_export_objs = -1;
          Zoltan_LB_Free_Part(export_global_ids, export_local_ids, export_procs,
                             export_to_part);
        }
      }
    }
    else {  /* *num_export_objs < 0 && *num_import_objs < 0) */
      if (zz->LB.Return_Lists) {
        /* No map at all available */
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Load-balancing function returned "
               "neither import nor export data.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_WARN;
      }
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done building return arguments");

  end_time = Zoltan_Time(zz->Timer);
  lb_time[0] = end_time - start_time;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST) {
    int i;
    Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
    printf("ZOLTAN: Objects to be imported to Proc %d\n", zz->Proc);
    for (i = 0; i < *num_import_objs; i++) {
      printf("    Obj: ");
      ZOLTAN_PRINT_GID(zz, &((*import_global_ids)[i*zz->Num_GID]));
      printf("  To partition: %4d", 
             (*import_to_part != NULL ? (*import_to_part)[i] 
                                      : zz->Proc));
      printf("  From processor: %4d\n", (*import_procs)[i]);
    }
    printf("\n");
    printf("ZOLTAN: Objects to be exported from Proc %d\n", zz->Proc);
    for (i = 0; i < *num_export_objs; i++) {
      printf("    Obj: ");
      ZOLTAN_PRINT_GID(zz, &((*export_global_ids)[i*zz->Num_GID]));
      printf("  To partition: %4d",
             (*export_to_part != NULL ? (*export_to_part)[i] 
                                      : (*export_procs)[i]));
      printf("  To processor: %4d\n", (*export_procs)[i]);
    }
    Zoltan_Print_Sync_End(zz->Communicator, TRUE);
  }

  /*
   *  If the Help_Migrate flag is set, perform migration for the application.
   */

  if (zz->Migrate.Auto_Migrate) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Begin auto-migration");

    start_time = Zoltan_Time(zz->Timer);
    error = Zoltan_Migrate(zz,
                            *num_import_objs, *import_global_ids,
                            *import_local_ids, *import_procs, *import_to_part,
                            *num_export_objs, *export_global_ids,
                            *export_local_ids, *export_procs, *export_to_part);
    if (error != ZOLTAN_OK && error != ZOLTAN_WARN) {
      sprintf(msg, "Error in auto-migration; %d returned from "
                    "Zoltan_Help_Migrate\n", error);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return error;
    }
    end_time = Zoltan_Time(zz->Timer);
    lb_time[1] = end_time - start_time;

    ZOLTAN_TRACE_DETAIL(zz, yo, "Done auto-migration");
  }
  
  /* Print timing info */
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ZTIME) {
    if (zz->Proc == zz->Debug_Proc) {
      printf("ZOLTAN Times:  \n");
    }
    Zoltan_Print_Stats (zz->Communicator, zz->Debug_Proc, lb_time[0], 
                   "ZOLTAN     Partition:     ");
    if (zz->Migrate.Auto_Migrate)
      Zoltan_Print_Stats (zz->Communicator, zz->Debug_Proc, lb_time[1], 
                      "ZOLTAN     Migrate: ");
  }

  *changes = 1;

  ZOLTAN_TRACE_EXIT(zz, yo);
  if (error)
    return (error);
  else
    return (ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Build_PartDist(ZZ *zz)
{
char *yo = "Zoltan_LB_Build_PartDist";
int ierr = ZOLTAN_OK;
int inflag[6], outflag[6] = {0,0,-1,0,0,0};
int global_parts_set = 0;   /* number of procs on which NUM_GLOBAL_PARTITIONS 
                               parameter was set. */
int local_parts_set = 0;    /* number of procs on which NUM_LOCAL_PARTITIONS
                               parameter was set. */
int max_global_parts = 0;   /* Max value of Num_Global_Parts_Param on all 
                               procs. */
int sum_local_parts = 0;    /* Sum of Num_Local_Parts over all procs.
                               Procs on which NUM_LOCAL_PARTITIONS was not
                               set assume zero parts on them.  Thus,
                               sum_local_parts may be < max_global_parts. */
int remaining_procs;        /* Num of procs not setting NUM_LOCAL_PARTITIONS */
int avail_local_parts;      /* max_global_parts - sum_local_parts */
int num_proc = zz->Num_Proc;
int *pdist;
int local_parts = 0;
int *local_parts_params = NULL;
int i, j, cnt, pcnt;
int frac = 0, mod = 0;
MPI_Op op;
MPI_User_function Zoltan_PartDist_MPIOp;

  MPI_Op_create(&Zoltan_PartDist_MPIOp,1,&op);

  /* Check whether global parts or local parts parameters were used. */
  inflag[0] = (zz->LB.Num_Global_Parts_Param != -1); 
  inflag[1] = (zz->LB.Num_Local_Parts_Param != -1); 
  inflag[2] =  zz->LB.Num_Global_Parts_Param;
  inflag[3] = ((zz->LB.Num_Local_Parts_Param == -1) 
                    ? 0 : zz->LB.Num_Local_Parts_Param);
  inflag[4] = (zz->LB.Num_Global_Parts_Param != zz->LB.Prev_Global_Parts_Param);
  inflag[5] = (zz->LB.Num_Local_Parts_Param != zz->LB.Prev_Local_Parts_Param);

  MPI_Allreduce(inflag, outflag, 6, MPI_INT, op, zz->Communicator);
  MPI_Op_free(&op);

  if (!outflag[4] && !outflag[5]) {
    /* Parameter values have not changed since last invocation of Zoltan. */
    /* Do not have to change PartDist or Num_Global_Parts. */
    goto End;
  }

  /* Since PartDist is changing, can't reuse old partitions.
   * Free LB.Data_Structure to prevent reuse. 
   * Also free LB.PartDist and LB.ProcDist.
   */
  if (zz->LB.Free_Structure != NULL)
    zz->LB.Free_Structure(zz);
  ZOLTAN_FREE(&(zz->LB.PartDist));
  ZOLTAN_FREE(&(zz->LB.ProcDist));

  zz->LB.Prev_Global_Parts_Param = zz->LB.Num_Global_Parts_Param;
  zz->LB.Prev_Local_Parts_Param = zz->LB.Num_Local_Parts_Param;

  global_parts_set = outflag[0];  /* Sum of inflag[0] */
  local_parts_set = outflag[1];   /* Sum of inflag[1] */
  max_global_parts = outflag[2];  /* Max of inflag[2] */
  sum_local_parts = outflag[3];   /* Sum of inflag[3] */

  /* Check whether any parameters were set;
   * No need to build the PartDist array if not. 
   */
  if ((!global_parts_set || (max_global_parts==num_proc)) && !local_parts_set) {
    /* Number of parts == number of procs, uniformly distributed; */
    zz->LB.Num_Global_Parts = num_proc;
    zz->LB.Single_Proc_Per_Part = 1;
  }

  else {
    /* Either NUM_GLOBAL_PARTITIONS is set != num_proc or NUM_LOCAL_PARTITIONS
     * is set.  Build PartDist, distributing partitions to processors as 
     * specified. 
     */

    /* error checking. */
    if (local_parts_set) {
      if (!global_parts_set) 
        max_global_parts = sum_local_parts;
      else if (sum_local_parts > max_global_parts) {
        char emsg[256];
        sprintf(emsg, 
                "Sum of NUM_LOCAL_PARTITIONS %d > NUM_GLOBAL_PARTITIONS %d", 
                sum_local_parts, max_global_parts);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, emsg);
        ierr = ZOLTAN_FATAL;
        goto End;
      }
      else if (sum_local_parts < max_global_parts && 
               local_parts_set == num_proc) {
        char emsg[256];
        sprintf(emsg, 
                "Sum of NUM_LOCAL_PARTITIONS %d < NUM_GLOBAL_PARTITIONS %d", 
                sum_local_parts, max_global_parts);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, emsg);
        ierr = ZOLTAN_FATAL;
        goto End;
      }
    }

    if (max_global_parts == 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zero partitions requested");
      ierr = ZOLTAN_FATAL;
      goto End;
    }

    /* Allocate space for PartDist. */
    zz->LB.PartDist = (int *) ZOLTAN_MALLOC((max_global_parts+1)*sizeof(int));
    if (zz->LB.PartDist == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      goto End;
    }

    pdist = zz->LB.PartDist;
 
    /* Compute the PartDist array. */

    if (!local_parts_set) {
      if (max_global_parts > num_proc) {
        /* NUM_LOCAL_PARTITIONS is not set; NUM_GLOBAL_PARTITIONS > num_proc. */
        /* Even distribution of partitions to processors. */
        zz->LB.Single_Proc_Per_Part = 1;
        frac = max_global_parts / num_proc;
        mod  = max_global_parts % num_proc;

        for (cnt = 0, i = 0; i < num_proc; i++) {
          local_parts = frac + ((num_proc - i) <= mod);
          for (j = 0; j < local_parts; j++)
            pdist[cnt++] = i;
        }
        pdist[cnt] = num_proc;
      }
      else { /* num_proc < max_global_parts */
        /* NUM_LOCAL_PARTITIONS is not set; NUM_GLOBAL_PARTITIONS < num_proc. */
        /* Even distribution of processors to partitions. */
        zz->LB.Single_Proc_Per_Part = 0;  /* Parts are spread across procs */
        pdist[0] = 0;
        frac = num_proc / max_global_parts;
        mod  = num_proc % max_global_parts;
        for (i = 1; i < max_global_parts; i++)
          pdist[i] = pdist[i-1] + frac + ((max_global_parts - i) <= mod);
        pdist[max_global_parts] = num_proc;
      }
    }
    else /* local_parts_set */ {

      /* NUM_LOCAL_PARTITIONS is set on at least some processors. */
      /* Distribute partitions to processors to match NUM_LOCAL_PARTITIONS
         where specified; distribute remaining partitions 
         to processors that didn't specify NUM_LOCAL_PARTITIONS */

      zz->LB.Single_Proc_Per_Part = 1;

      /* Gather the parameter values from all processors. */
      local_parts_params = (int *) ZOLTAN_MALLOC((num_proc+1)* sizeof(int));
      MPI_Allgather(&(zz->LB.Num_Local_Parts_Param), 1, MPI_INT, 
                    local_parts_params, 1, MPI_INT, zz->Communicator);

      /* Compute number of parts not specified by NUM_LOCAL_PARTITIONS */
      /* In MPI_Allreduce above, processors not specifying NUM_LOCAL_PARTITIONS
       * specified contributed one partition to sum_local_parts.  */

      remaining_procs = num_proc - local_parts_set;
      avail_local_parts = max_global_parts - sum_local_parts;
      if (remaining_procs > 0) {
        frac = avail_local_parts / remaining_procs;
        mod  = avail_local_parts % remaining_procs;
      }

      for (cnt = 0, pcnt = 0, i = 0; i < num_proc; i++)
        if (local_parts_params[i] != -1) {
          /* Fill in processor for its NUM_LOCAL_PARTITIONS partitions. */
          for (j = 0; j < local_parts_params[i]; j++)
            pdist[cnt++] = i;
        }
        else {
          /* Equally distribute avail_local_parts among remaining_procs. */
          local_parts = frac + ((remaining_procs - pcnt) <= mod);
          for (j = 0; j < local_parts; j++)
            pdist[cnt++] = i;
          pcnt++;
        }
  
      pdist[cnt] = num_proc;
      ZOLTAN_FREE(&local_parts_params);
    }

    /* Reset Num_Global_Parts.  */
    zz->LB.Num_Global_Parts = max_global_parts;

    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL && zz->LB.PartDist != NULL) {
      printf("[%1d] Debug: LB.PartDist = ", zz->Proc);
      for (i=0; i<=zz->LB.Num_Global_Parts; i++)
        printf("%d ", zz->LB.PartDist[i]);
      printf("\n");
    }
  }

End:
  return ierr;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void Zoltan_PartDist_MPIOp(
  void *in, 
  void *inout, 
  int *len, 
  MPI_Datatype *dptr)
{
int *int_in = (int *) in;
int *int_inout = (int *) inout;

  int_inout[0] += int_in[0];
  int_inout[1] += int_in[1];
  int_inout[2] = ((int_in[2] > int_inout[2]) ? int_in[2] : int_inout[2]);
  int_inout[3] += int_in[3];
  int_inout[4] += int_in[4];
  int_inout[5] += int_in[5];
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
