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
#include "comm_const.h"
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines implementing the migration-help tools.
 *  These functions are all callable by the application.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Compute_Destinations(
  LB *lb,                    /* Load balancing structure for current balance.*/
  int num_import,            /* Number of non-local objects assigned to the 
                                processor in the new decomposition.          */
  LB_GID *import_global_ids, /* Array of global IDs for non-local objects 
                                assigned to this processor in the new
                                decomposition.                               */
  LB_LID *import_local_ids,  /* Array of local IDs for non-local objects
                                assigned to the processor in the new
                                decomposition.                               */
  int *import_procs,         /* Array of processor IDs of processors owning
                                the non-local objects that are assigned to
                                this processor in the new decomposition.     */
  int *num_export,           /* Returned value:  Number of objs to be exported
                                to other processors to establish the new
                                decomposition.                               */
  LB_GID **export_global_ids,/* Returned value:  Array of global IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  LB_LID **export_local_ids, /* Returned value:  Array of local IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  int **export_procs         /* Returned value:  Array of processor IDs
                                to which objects will be exported 
                                to establish the new decomposition.          */
)
{
/*
 *  Routine to compute the inverse map:  Given, for each processor, a list 
 *  of non-local objects assigned to the processor, compute the list of objects
 *  that processor needs to export to other processors to establish the new
 *  decomposition.
 */

char *yo = "LB_Compute_Destinations";
char msg[256];
int *proc_list = NULL;      /* List of processors from which objs are to be 
                               imported.                                    */
COMM_OBJ *comm_plan;        /* Object returned communication routines  */
LB_TAG *import_objs = NULL; /* Array of import objects used to request objs
                               from other processors.                       */
LB_TAG *export_objs = NULL; /* Array of export objects describing which objs
                               must be sent to other processors.            */
int msgtag, msgtag2;        /* Message tags for communication routines */
int i;
int ierr = LB_OK;


  LB_TRACE_ENTER(lb, yo);
  /*
   *  Return if this processor is not in the load-balancing structure's
   *  communicator.
   */

  if (LB_PROC_NOT_IN_COMMUNICATOR(lb)) {
    LB_TRACE_EXIT(lb, yo);
    return (LB_OK);
  }

  /*
   *  Build processor's list of requests for non-local objs.
   */

  if (num_import > 0) {
    proc_list = (int *) LB_MALLOC(num_import*sizeof(int));
    if (!proc_list) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }
    import_objs = (LB_TAG *) LB_MALLOC(num_import*sizeof(LB_TAG));
    if (!import_objs) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }

    for (i = 0; i < num_import; i++) {
      proc_list[i] = import_procs[i];

      LB_SET_GID(import_objs[i].Global_ID, import_global_ids[i]);
      LB_SET_LID(import_objs[i].Local_ID, import_local_ids[i]);
      import_objs[i].Proc = lb->Proc;
    }
  }

  /*
   *  Compute communication map and num_export, the number of objs this
   *  processor has to export to establish the new decomposition.
   */

  msgtag = 32767;
  ierr = LB_Comm_Create(&comm_plan, num_import, proc_list, lb->Communicator, 
                        msgtag, lb->Deterministic, num_export);
  if (ierr != LB_COMM_OK && ierr != LB_COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Create.",
            (ierr == LB_COMM_MEMERR ? "LB_COMM_MEMERR" : "LB_COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&proc_list);
    LB_FREE(&import_objs);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == LB_COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }
  

  LB_TRACE_DETAIL(lb, yo, "Done comm create");

  /*
   *  Allocate space for the object tags that need to be exported.  Communicate
   *  to get the list of objects to be exported.
   */

  if (*num_export > 0) {
    export_objs = (LB_TAG *) LB_MALLOC((*num_export)*sizeof(LB_TAG));
    if (!export_objs) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_FREE(&import_objs);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }
    if (!LB_Special_Malloc(lb,(void **)export_global_ids,*num_export,
                           LB_SPECIAL_MALLOC_GID)) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_FREE(&import_objs);
      LB_FREE(&export_objs);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }
    if (!LB_Special_Malloc(lb,(void **)export_local_ids,*num_export,
                           LB_SPECIAL_MALLOC_LID)) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_FREE(&import_objs);
      LB_FREE(&export_objs);
      LB_Special_Free(lb,(void **)export_global_ids,LB_SPECIAL_MALLOC_GID);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }
    if (!LB_Special_Malloc(lb,(void **)export_procs,*num_export,
                           LB_SPECIAL_MALLOC_INT)) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_FREE(&import_objs);
      LB_FREE(&export_objs);
      LB_Special_Free(lb,(void **)export_global_ids,LB_SPECIAL_MALLOC_GID);
      LB_Special_Free(lb,(void **)export_local_ids,LB_SPECIAL_MALLOC_LID);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }

  }
  else {
    export_objs = NULL;
    *export_global_ids = NULL;
    *export_local_ids = NULL;
    *export_procs = NULL;
  }

  msgtag2 = 32766;
  ierr = LB_Comm_Do(comm_plan, msgtag2, (char *) import_objs, 
                    (int) sizeof(LB_TAG), (char *) export_objs);
  if (ierr != LB_COMM_OK && ierr != LB_COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Do.", 
            (ierr == LB_COMM_MEMERR ? "LB_COMM_MEMERR" : "LB_COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&proc_list);
    LB_FREE(&import_objs);
    LB_FREE(&export_objs);
    LB_Comm_Destroy(&comm_plan);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == LB_COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }

  LB_TRACE_DETAIL(lb, yo, "Done comm_do");

  /*
   *  Put the exported LB_TAGs into the output format.
   */

  for (i = 0; i < *num_export; i++) {
    LB_SET_GID((*export_global_ids)[i], export_objs[i].Global_ID);
    LB_SET_LID((*export_local_ids)[i], export_objs[i].Local_ID);
    (*export_procs)[i]      = export_objs[i].Proc;
  }

  LB_FREE(&proc_list);
  LB_FREE(&import_objs);
  LB_FREE(&export_objs);
  
  LB_Comm_Destroy(&comm_plan);

  LB_TRACE_DETAIL(lb, yo, "Done comm destroy");

  LB_TRACE_EXIT(lb, yo);
  return (ierr);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Help_Migrate(
  LB *lb,                    /* Load balancing structure for current balance.*/
  int num_import,            /* Number of non-local objects assigned to the 
                                processor in the new decomposition.          */
  LB_GID *import_global_ids, /* Array of global IDs for non-local objects 
                                assigned to this processor in the new
                                decomposition; this field can be NULL if 
                                the application does not provide import IDs. */
  LB_LID *import_local_ids,  /* Array of local IDs for non-local objects
                                assigned to the processor in the new
                                decomposition; this field can be NULL if the 
                                application does not provide import IDs.     */
  int *import_procs,         /* Array of processor IDs of processors owning
                                the non-local objects that are assigned to
                                this processor in the new decomposition; this
                                field can be NULL if the application does
                                not provide import IDs.                      */
  int num_export,            /* Number of objs to be exported
                                to other processors to establish the new
                                decomposition.                               */
  LB_GID *export_global_ids, /* Array of global IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  LB_LID *export_local_ids,  /* Array of local IDs of
                                objects to be exported to other processors
                                to establish the new decomposition.          */
  int *export_procs          /* Array of processor IDs
                                to which objects will be exported 
                                to establish the new decomposition.          */
)
{
/*
 *  Routine to help perform migration.  If migration pre-processing routine
 *  (LB_PRE_MIGRATE_FN) is specified, this routine first calls that function.
 *  It then calls a function to obtain the size of the migrating objects
 *  (LB_OBJ_SIZE_FN).  The routine next calls an application-specified
 *  object packing routine (LB_PACK_OBJ_FN) for each object
 *  to be exported.  It develops the needed communication map to move the
 *  objects to other processors.  It performs the communication according
 *  to the map, and then calls an application-specified object unpacking
 *  routine (LB_UNPACK_OBJ_FN) for each object imported.
 */

char *yo = "LB_Help_Migrate";
char msg[256];
int size;                /* size (in bytes) of the object data for export.  */
int id_size;             /* size (in bytes) of LB_GID + padding for 
                            alignment                                       */
char *export_buf = NULL; /* buffer for packing export data.                 */
char *import_buf = NULL; /* buffer for receiving imported data.             */
char *tmp;               /* temporary pointer into buffers.                 */
int i;                   /* loop counter.                                   */
int tmp_import;          /* number of objects to be imported.               */
int *proc_list = NULL;   /* list of processors to which this proc exports.  */
LB_GID global_id;        /* tmp global ID for unpacking objects.            */
LB_GID *tmp_id;          /* pointer to storage for an LB_GID in comm buf    */
COMM_OBJ *comm_plan;     /* Object returned by communication routines       */
int msgtag, msgtag2;     /* Tags for communication routines                 */
int ierr = 0;

  LB_TRACE_ENTER(lb, yo);

  /*
   *  Return if this processor is not in the load-balancing structure's
   *  communicator.
   */

  if (LB_PROC_NOT_IN_COMMUNICATOR(lb)) {
    LB_TRACE_EXIT(lb, yo);
    return (LB_OK);
  }


  if (lb->Migrate.Get_Obj_Size == NULL) {
    LB_PRINT_ERROR(lb->Proc, yo, "Must register an "
           "LB_OBJ_SIZE_FN_TYPE function to use the migration-help tools.");
    LB_TRACE_EXIT(lb, yo);
    return (LB_FATAL);
  }

  if (lb->Migrate.Pack_Obj == NULL) {
    LB_PRINT_ERROR(lb->Proc, yo, "Must register an "
           "LB_PACK_OBJ_FN_TYPE function to use the migration-help tools.");
    LB_TRACE_EXIT(lb, yo);
    return (LB_FATAL);
  }

  if (lb->Migrate.Unpack_Obj == NULL) {
    LB_PRINT_ERROR(lb->Proc, yo, "Must register an "
         "LB_UNPACK_OBJ_FN_TYPE function to use the migration-help tools.");
    LB_TRACE_EXIT(lb, yo);
    return (LB_FATAL);
  }

  if (lb->Migrate.Pre_Migrate != NULL) {
    lb->Migrate.Pre_Migrate(lb->Migrate.Pre_Migrate_Data,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, &ierr);
    if (ierr) {
      LB_PRINT_ERROR(lb->Proc, yo, "Error returned from user defined "
                      "Migrate.Pre_Migrate function.");
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }
  }

  LB_TRACE_DETAIL(lb, yo, "Done pre-migration processing");

  /*
   * For each object, allow space for its LB_GID and its data.
   * Zoltan will pack the LB_GIDs; the application must pack the data
   * through the pack routine.  Zoltan needs the LB_GIDs for unpacking,
   * as the order of the data received during communication is not 
   * necessarily the same order as import_global_ids[].
   */
  id_size = sizeof(LB_GID) + LB_pad_for_alignment(sizeof(LB_GID));
  size = lb->Migrate.Get_Obj_Size(lb->Migrate.Get_Obj_Size_Data, &ierr)
       + id_size;
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, "Error returned from user defined "
                    "Migrate.Get_Obj_Size function.");
    LB_TRACE_EXIT(lb, yo);
    return (LB_FATAL);
  }


  if (num_export > 0) {
    export_buf = (char *) LB_MALLOC(num_export*size);
    if (!export_buf) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }

    proc_list = (int *) LB_MALLOC(num_export*sizeof(int));
    if (!proc_list) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&export_buf);
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }


    /*
     *  Pack the proc_list (to create the map) and the objects for export.
     */
  
    tmp = export_buf;
    for (i = 0; i < num_export; i++) {
      proc_list[i] = export_procs[i];

      /* Pack the object's global ID */
      tmp_id = (LB_GID *) tmp;
      LB_SET_GID(*tmp_id, export_global_ids[i]);
    
      /* Pack the object's data */
      lb->Migrate.Pack_Obj(lb->Migrate.Pack_Obj_Data, export_global_ids[i],
                           export_local_ids[i], export_procs[i], size,
                           tmp+id_size, &ierr);
      if (ierr) {
        LB_PRINT_ERROR(lb->Proc, yo, "Error returned from user defined "
                        "Migrate.Pack_Obj function.");
        LB_FREE(&export_buf);
        LB_FREE(&proc_list);
        LB_TRACE_EXIT(lb, yo);
        return (LB_FATAL);
      }
      tmp += size;
    }
  }

  LB_TRACE_DETAIL(lb, yo, "Done packing objects");

  /*
   *  Compute communication map and tmp_import, the number of objs this
   *  processor has to import to establish the new decomposition.
   */

  msgtag = 32767;
  ierr = LB_Comm_Create(&comm_plan, num_export, proc_list, lb->Communicator, 
                        msgtag, lb->Deterministic, &tmp_import);
  if (ierr != LB_COMM_OK && ierr != LB_COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Create.", 
            (ierr == LB_COMM_MEMERR ? "LB_COMM_MEMERR" : "LB_COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&export_buf);
    LB_FREE(&proc_list);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == LB_COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }
  if (tmp_import != num_import) {
    sprintf(msg, "tmp_import %d != num_import %d.", tmp_import, num_import);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
  }

  if (num_import > 0) {
    import_buf = (char *) LB_MALLOC(num_import*size);
    if (!import_buf) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&export_buf);
      LB_FREE(&proc_list);
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }

  }

  /*
   *  Send the export data using the communication plan.
   */

  msgtag2 = 32766;
  ierr = LB_Comm_Do(comm_plan, msgtag2, export_buf, size, import_buf);
  if (ierr != LB_COMM_OK && ierr != LB_COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Do.", 
            (ierr == LB_COMM_MEMERR ? "LB_COMM_MEMERR" : "LB_COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&proc_list);
    LB_FREE(&export_buf);
    LB_FREE(&import_buf);
    LB_Comm_Destroy(&comm_plan);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == LB_COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }

  /*
   *  Free whatever memory we can.
   */

  LB_Comm_Destroy(&comm_plan);
  LB_FREE(&proc_list);
  LB_FREE(&export_buf);

  LB_TRACE_DETAIL(lb, yo, "Done communication");

  /* 
   *  Perform application-specified processing before unpacking the data.
   */
  if (lb->Migrate.Mid_Migrate != NULL) {
    lb->Migrate.Mid_Migrate(lb->Migrate.Mid_Migrate_Data,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, &ierr);
    if (ierr) {
      LB_PRINT_ERROR(lb->Proc, yo, "Error returned from user defined "
                      "Migrate.Mid_Migrate function.");
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }

    LB_TRACE_DETAIL(lb, yo, "Done mid-migration processing");
  }

  /*
   *  Unpack the object data.
   */

  tmp = import_buf;
  for (i = 0; i < num_import; i++) {

    /* Unpack the object's global ID */
    tmp_id = (LB_GID *) tmp;
    LB_SET_GID(global_id, *tmp_id);

    /* Unpack the object's data */
    lb->Migrate.Unpack_Obj(lb->Migrate.Unpack_Obj_Data, global_id, size,
                           tmp+id_size, &ierr);
    if (ierr) {
      LB_PRINT_ERROR(lb->Proc, yo, "Error returned from user defined "
                      "Migrate.Unpack_Obj function.");
      LB_FREE(&import_buf);
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }
    tmp += size;
  }

  LB_FREE(&import_buf);

  LB_TRACE_DETAIL(lb, yo, "Done unpacking objects");

  if (lb->Migrate.Post_Migrate != NULL) {
    lb->Migrate.Post_Migrate(lb->Migrate.Post_Migrate_Data,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, &ierr);
    if (ierr) {
      LB_PRINT_ERROR(lb->Proc, yo, "Error returned from user defined "
                      "Migrate.Post_Migrate function.");
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }

    LB_TRACE_DETAIL(lb, yo, "Done post-migration processing");
  }

  LB_TRACE_EXIT(lb, yo);
  return (LB_OK);
}
