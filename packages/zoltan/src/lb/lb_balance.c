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

int Zoltan_LB_Balance(
  ZZ *zz, 
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
 * Input:  a Zoltan structure with appropriate function pointers set.
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

char *yo = "Zoltan_LB_Balance";
int gmax;    /* Maximum number of imported/exported objects 
                over all processors.                       */
int error;    /* Error code */
double start_time, end_time;
double lb_time[2] = {0.0,0.0};
char msg[256];
int comm[3],gcomm[3]; 

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
  *export_global_ids = NULL;
  *export_local_ids = NULL;
  *export_procs = NULL;

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

  /*
   * Call the actual load-balancing function.
   */

  error = zz->LB.LB_Fn(zz, num_import_objs, import_global_ids, import_local_ids,
          import_procs, num_export_objs, export_global_ids, 
          export_local_ids, export_procs);

  if (error == ZOLTAN_FATAL){
    sprintf(msg, "Balancing routine returned error code %d.", error);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (error);
  }
  else if (error){
    sprintf(msg, "Balancing routine returned error code %d.", error);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done load balancing");

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
      printf("%s No changes to the decomposition due to load-balancing; "
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
        Zoltan_LB_Free_Data(import_global_ids, import_local_ids, import_procs,
                            export_global_ids, export_local_ids, export_procs);
        ZOLTAN_PRINT_WARN(zz->Proc, yo, 
                      "Method returned lists, but no lists requested.");
      }
    }
    else if (zz->LB.Return_Lists == ZOLTAN_LB_ALL_LISTS || 
             zz->LB.Return_Lists == ZOLTAN_LB_EXPORT_LISTS) {
      /* Export lists are requested; compute export map */
      error = Zoltan_Compute_Destinations(zz,
                                      *num_import_objs, *import_global_ids, 
                                      *import_local_ids, *import_procs,
                                      num_export_objs, export_global_ids,
                                      export_local_ids, export_procs);
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
        Zoltan_LB_Free_Data(import_global_ids, import_local_ids, import_procs,
                            NULL, NULL, NULL);
      }
    }
  }
  else { /* (*num_import_objs < 0) */
    if (*num_export_objs >= 0) {
      /* Only export lists have been returned. */
      if (zz->LB.Return_Lists == ZOLTAN_LB_ALL_LISTS || 
          zz->LB.Return_Lists == ZOLTAN_LB_IMPORT_LISTS) {
        /* Compute import map */
        error = Zoltan_Compute_Destinations(zz, 
                                        *num_export_objs, *export_global_ids, 
                                        *export_local_ids, *export_procs,
                                        num_import_objs, import_global_ids,
                                        import_local_ids, import_procs);

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
          Zoltan_LB_Free_Data(NULL, NULL, NULL, 
                             export_global_ids, export_local_ids, export_procs);
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
      printf("  From processor: %4d\n", (*import_procs)[i]);
    }
    printf("\n");
    printf("ZOLTAN: Objects to be exported from Proc %d\n", zz->Proc);
    for (i = 0; i < *num_export_objs; i++) {
      printf("    Obj: ");
      ZOLTAN_PRINT_GID(zz, &((*export_global_ids)[i*zz->Num_GID]));
      printf("  Destination: %4d\n", (*export_procs)[i]);
    }
    Zoltan_Print_Sync_End(zz->Communicator, TRUE);
  }

  /*
   *  If the Help_Migrate flag is set, perform migration for the application.
   */

  if (zz->Migrate.Auto_Migrate) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Begin auto-migration");

    start_time = Zoltan_Time(zz->Timer);
    error = Zoltan_Help_Migrate(zz,
                            *num_import_objs, *import_global_ids,
                            *import_local_ids, *import_procs,
                            *num_export_objs, *export_global_ids,
                            *export_local_ids, *export_procs);
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
                   "ZOLTAN     Balance:     ");
    if (zz->Migrate.Auto_Migrate)
      Zoltan_Print_Stats (zz->Communicator, zz->Debug_Proc, lb_time[1], 
                      "ZOLTAN     HelpMigrate: ");
  }

  *changes = 1;

  ZOLTAN_TRACE_EXIT(zz, yo);
  if (error)
    return (error);
  else
    return (ZOLTAN_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_LB_Free_Data(
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

  ZOLTAN_FREE(import_global_ids);
  ZOLTAN_FREE(import_local_ids);
  ZOLTAN_FREE(import_procs);
  ZOLTAN_FREE(export_global_ids);
  ZOLTAN_FREE(export_local_ids);
  ZOLTAN_FREE(export_procs);

  return (ZOLTAN_OK);

}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
