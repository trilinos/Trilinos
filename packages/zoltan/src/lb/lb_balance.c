/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "lb_util_const.h"
#include "params_const.h"
#include "timer_const.h"
#include "ha_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for performing partitioning with Zoltan.
 *  These functions are all callable by the application.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Balance(
  LB *lb, 
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
                                    lb->Communicator of the parameter
                                    Num_Global_ID_Entries.                    */
  int *num_lid_entries,          /* The number of array entries in a local ID;
                                    set to be the max over all processors in
                                    lb->Communicator of the parameter
                                    Num_Local_ID_Entries.                     */
  int *num_import_objs,          /* The number of non-local objects in the
                                    processor's new decomposition.            */
  ZOLTAN_ID_PTR *import_global_ids,  /* Array of global IDs for non-local objects
                                    (i.e., objs to be imported) in
                                    the processor's new decomposition.        */
  ZOLTAN_ID_PTR *import_local_ids,   /* Array of local IDs for non-local objects
                                    (i.e., objs to be imported) in
                                    the processor's new decomposition.        */
  int **import_procs,            /* Array of processor IDs for processors 
                                    currently owning non-local objects (i.e.,
                                    objs to be imported) in this processor's
                                    new decomposition.                        */
  int *num_export_objs,          /* The number of local objects that need to
                                    be exported from the processor to establish
                                    the new decomposition.                    */
  ZOLTAN_ID_PTR *export_global_ids,  /* Array of global IDs for objects that need
                                    to be exported (assigned and sent to other
                                    processors) to establish the new 
                                    decomposition.                            */
  ZOLTAN_ID_PTR *export_local_ids,   /* Array of local IDs for objects that need
                                    to be exported (assigned and sent to other
                                    processors) to establish the new 
                                    decomposition.                            */
  int **export_procs             /* Array of destination processor IDs for
                                    objects that need to be exported 
                                    to establish the new decomposition.       */
)
{
/*
 * Main user-call for load-balancing.
 * Input:  a load-balancing structure with appropriate function pointers set.
 * Output: 
 *   num_import_objs
 *   import_global_ids
 *   import_local_ids
 *   import_procs
 *   num_export_objs
 *   export_global_ids
 *   export_local_ids
 *   export_procs
 * Return values:
 *   Return zero if the decomposition is unchanged by load-balancing;
 *   Return one otherwise.
 */

char *yo = "LB_Balance";
int gmax;    /* Maximum number of imported/exported objects 
                over all processors.                       */
int error;    /* Error code */
double start_time, end_time;
double lb_time[2] = {0.0,0.0};
char msg[256];
int comm[3],gcomm[3]; 

  LB_TRACE_ENTER(lb, yo);

  if (lb->Proc == lb->Debug_Proc && lb->Debug_Level >= LB_DEBUG_PARAMS)
    LB_Print_Key_Params(lb);

  start_time = LB_Time(lb->Timer);

  /* 
   * Compute Max number of array entries per ID over all processors.
   * Compute Max number of return arguments for LB_Balance.
   * This is a sanity-maintaining step; we don't want different
   * processors to have different values for these numbers.
   */
  comm[0] = lb->Num_GID;
  comm[1] = lb->Num_LID;
  comm[2] = lb->Return_Lists;
  MPI_Allreduce(comm, gcomm, 3, MPI_INT, MPI_MAX, lb->Communicator);
  lb->Num_GID = *num_gid_entries = gcomm[0];
  lb->Num_LID = *num_lid_entries = gcomm[1];
  lb->Return_Lists = gcomm[2];

  /* assume no changes */
  *changes = 0;

  *num_import_objs = *num_export_objs = 0;
  *import_global_ids = NULL;
  *import_local_ids = NULL;
  *import_procs = NULL;
  *export_global_ids = NULL;
  *export_local_ids = NULL;
  *export_procs = NULL;

  /*
   *  Return if this processor is not in the load-balancing structure's
   *  communicator.
   */

  if (LB_PROC_NOT_IN_COMMUNICATOR(lb)) {
    LB_TRACE_EXIT(lb, yo);
    return (ZOLTAN_OK);
  }

  if (lb->Method == NONE) {
    if (lb->Proc == lb->Debug_Proc && lb->Debug_Level >= LB_DEBUG_PARAMS)
      printf("%s Balancing method selected == NONE; no balancing performed\n",
              yo);

    LB_TRACE_EXIT(lb, yo);
    return (ZOLTAN_WARN);
  }

  /*
   *  Construct the heterogenous machine description.
   */

  error = LB_Build_Machine_Desc(lb);

  if (error == ZOLTAN_FATAL){
    LB_TRACE_EXIT(lb, yo);
    return (error);
  }

  LB_TRACE_DETAIL(lb, yo, "Done machine description");

  /*
   * Call the actual load-balancing function.
   */

  error = lb->LB_Fn(lb, num_import_objs, import_global_ids, import_local_ids,
          import_procs, num_export_objs, export_global_ids, 
          export_local_ids, export_procs);

  if (error == ZOLTAN_FATAL){
    sprintf(msg, "Balancing routine returned error code %d.", error);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_TRACE_EXIT(lb, yo);
    return (error);
  }
  else if (error){
    sprintf(msg, "Balancing routine returned error code %d.", error);
    LB_PRINT_WARN(lb->Proc, yo, msg);
  }

  LB_TRACE_DETAIL(lb, yo, "Done load balancing");

  if (*num_import_objs >= 0)
    MPI_Allreduce(num_import_objs, &gmax, 1, MPI_INT, MPI_MAX, 
                lb->Communicator);
  else /* use export data */
    MPI_Allreduce(num_export_objs, &gmax, 1, MPI_INT, MPI_MAX, 
                lb->Communicator);

  if (gmax == 0) {

    /*
     *  Decomposition was not changed by the load balancing; no migration
     *  is needed.
     */

    if (lb->Proc == lb->Debug_Proc && lb->Debug_Level >= LB_DEBUG_PARAMS)
      printf("%s No changes to the decomposition due to load-balancing; "
             "no migration is needed.\n", yo);

    /*
     *  Reset num_import_objs and num_export_objs; don't want to return
     *  -1 for the arrays that weren't returned by LB_FN.
     */

    *num_import_objs = *num_export_objs = 0;

    LB_TRACE_EXIT(lb, yo);
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

      if (lb->Return_Lists == LB_NO_LISTS) {
        /* This condition should never happen!! */
        /* Methods should not return arrays if no lists are requested. */
        *num_import_objs = *num_export_objs = -1;
        LB_Free_Data(import_global_ids, import_local_ids, import_procs,
                     export_global_ids, export_local_ids, export_procs);
        LB_PRINT_WARN(lb->Proc, yo, 
                      "Method returned lists, but no lists requested.");
      }
    }
    else if (lb->Return_Lists == LB_ALL_LISTS || 
             lb->Return_Lists == LB_EXPORT_LISTS) {
      /* Export lists are requested; compute export map */
      error = LB_Compute_Destinations(lb,
                                      *num_import_objs, *import_global_ids, 
                                      *import_local_ids, *import_procs,
                                      num_export_objs, export_global_ids,
                                      export_local_ids, export_procs);
      if (error != ZOLTAN_OK && error != ZOLTAN_WARN) {
        sprintf(msg, "Error building return arguments; "
                     "%d returned by LB_Compute_Destinations\n", error);
        LB_PRINT_ERROR(lb->Proc, yo, msg);
        LB_TRACE_EXIT(lb, yo);
        return error;
      }
      if (lb->Return_Lists == LB_EXPORT_LISTS) {
        /* Method returned import lists, but only export lists were desired. */
        /* Import lists not needed; free them. */
        *num_import_objs = -1;
        LB_Free_Data(import_global_ids, import_local_ids, import_procs,
                     NULL, NULL, NULL);
      }
    }
  }
  else { /* (*num_import_objs < 0) */
    if (*num_export_objs >= 0) {
      /* Only export lists have been returned. */
      if (lb->Return_Lists == LB_ALL_LISTS || 
          lb->Return_Lists == LB_IMPORT_LISTS) {
        /* Compute import map */
        error = LB_Compute_Destinations(lb, 
                                        *num_export_objs, *export_global_ids, 
                                        *export_local_ids, *export_procs,
                                        num_import_objs, import_global_ids,
                                        import_local_ids, import_procs);

        if (error != ZOLTAN_OK && error != ZOLTAN_WARN) {
          sprintf(msg, "Error building return arguments; "
                       "%d returned by LB_Compute_Destinations\n", error);
          LB_PRINT_ERROR(lb->Proc, yo, msg);
          LB_TRACE_EXIT(lb, yo);
          return error;
        }
        if (lb->Return_Lists == LB_IMPORT_LISTS) {
          /* Method returned export lists, but only import lists are desired. */
          /* Export lists not needed; free them. */
          *num_export_objs = -1;
          LB_Free_Data(NULL, NULL, NULL, 
                       export_global_ids, export_local_ids, export_procs);
        }
      }
    }
    else {  /* *num_export_objs < 0 && *num_import_objs < 0) */
      if (lb->Return_Lists) {
        /* No map at all available */
        LB_PRINT_ERROR(lb->Proc, yo, "Load-balancing function returned "
               "neither import nor export data.");
        LB_TRACE_EXIT(lb, yo);
        return ZOLTAN_WARN;
      }
    }
  }

  LB_TRACE_DETAIL(lb, yo, "Done building return arguments");

  end_time = LB_Time(lb->Timer);
  lb_time[0] = end_time - start_time;

  if (lb->Debug_Level >= LB_DEBUG_LIST) {
    int i;
    LB_Print_Sync_Start(lb->Communicator, TRUE);
    printf("ZOLTAN: Objects to be imported to Proc %d\n", lb->Proc);
    for (i = 0; i < *num_import_objs; i++) {
      printf("    Obj: ");
      LB_PRINT_GID(lb, &((*import_global_ids)[i*lb->Num_GID]));
      printf("  From processor: %4d\n", (*import_procs)[i]);
    }
    printf("\n");
    printf("ZOLTAN: Objects to be exported from Proc %d\n", lb->Proc);
    for (i = 0; i < *num_export_objs; i++) {
      printf("    Obj: ");
      LB_PRINT_GID(lb, &((*export_global_ids)[i*lb->Num_GID]));
      printf("  Destination: %4d\n", (*export_procs)[i]);
    }
    LB_Print_Sync_End(lb->Communicator, TRUE);
  }

  /*
   *  If the Help_Migrate flag is set, perform migration for the application.
   */

  if (lb->Migrate.Auto_Migrate) {
    LB_TRACE_DETAIL(lb, yo, "Begin auto-migration");

    start_time = LB_Time(lb->Timer);
    error = LB_Help_Migrate(lb,
                            *num_import_objs, *import_global_ids,
                            *import_local_ids, *import_procs,
                            *num_export_objs, *export_global_ids,
                            *export_local_ids, *export_procs);
    if (error != ZOLTAN_OK && error != ZOLTAN_WARN) {
      sprintf(msg, "Error in auto-migration; %d returned from "
                    "LB_Help_Migrate\n", error);
      LB_PRINT_ERROR(lb->Proc, yo, msg);
      LB_TRACE_EXIT(lb, yo);
      return error;
    }
    end_time = LB_Time(lb->Timer);
    lb_time[1] = end_time - start_time;

    LB_TRACE_DETAIL(lb, yo, "Done auto-migration");
  }
  
  /* Print timing info */
  if (lb->Debug_Level >= LB_DEBUG_ZTIME) {
    if (lb->Proc == lb->Debug_Proc) {
      printf("ZOLTAN Times:  \n");
    }
    LB_Print_Stats (lb->Communicator, lb->Debug_Proc, lb_time[0], 
                   "ZOLTAN     Balance:     ");
    if (lb->Migrate.Auto_Migrate)
      LB_Print_Stats (lb->Communicator, lb->Debug_Proc, lb_time[1], 
                      "ZOLTAN     HelpMigrate: ");
  }

  *changes = 1;

  LB_TRACE_EXIT(lb, yo);
  if (error)
    return (error);
  else
    return (ZOLTAN_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Free_Data(
  ZOLTAN_ID_PTR *import_global_ids, /* Array of global IDs for non-local objects 
                                    assigned to this processor in the new
                                    decomposition.                           */
  ZOLTAN_ID_PTR *import_local_ids,  /* Array of local IDs for non-local objects
                                    assigned to the processor in the new
                                    decomposition.                           */
  int **import_procs,           /* Array of processor IDs of processors owning
                                   the non-local objects that are assigned to
                                   this processor in the new decomposition.  */
  ZOLTAN_ID_PTR *export_global_ids, /* Array of global IDs of
                                   objects to be exported to other processors
                                   to establish the new decomposition.       */
  ZOLTAN_ID_PTR *export_local_ids,  /* Array of local IDs of
                                   objects to be exported to other processors
                                   to establish the new decomposition.       */
  int **export_procs            /* Array of processor IDs
                                   to which objects will be exported 
                                   to establish the new decomposition.       */
)
{
/*
 *  Routine to free the arrays returning the results of the load balancing.
 */

  LB_FREE(import_global_ids);
  LB_FREE(import_local_ids);
  LB_FREE(import_procs);
  LB_FREE(export_global_ids);
  LB_FREE(export_local_ids);
  LB_FREE(export_procs);

  return (ZOLTAN_OK);

}
