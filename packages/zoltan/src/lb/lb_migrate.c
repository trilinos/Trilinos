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
  LB *lb,                      /* Load balancing structure.                  */
  int num_import,              /* Number of non-local objects assigned to the 
                                  processor in the new decomposition.        */
  LB_ID_PTR import_global_ids, /* Array of global IDs for non-local objects 
                                  assigned to this processor in the new
                                  decomposition.                             */
  LB_ID_PTR import_local_ids,  /* Array of local IDs for non-local objects
                                  assigned to the processor in the new
                                  decomposition.                             */
  int *import_procs,           /* Array of processor IDs of processors owning
                                  the non-local objects that are assigned to
                                  this processor in the new decomposition.   */
  int *num_export,             /* Returned value:  Number of objs to be exported
                                  to other processors to establish the new
                                  decomposition.                             */
  LB_ID_PTR *export_global_ids,/* Returned value:  Array of global IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  LB_ID_PTR *export_local_ids, /* Returned value:  Array of local IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  int **export_procs           /* Returned value:  Array of processor IDs
                                  to which objects will be exported 
                                  to establish the new decomposition.        */
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
int *import_proc_list = NULL;
                            /* Array containing owning processor IDs of import
                               objects; used to request objs from other procs.*/
int msgtag, msgtag2;        /* Message tags for communication routines */
int num_gid_entries, num_lid_entries;  /* Length of global and local ids */
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
   *  Check that all procs use the same id types.
   */

  MPI_Allreduce(&(lb->Num_GID), &num_gid_entries, 1,
                MPI_INT, MPI_MAX, lb->Communicator);
  if (lb->Num_GID != num_gid_entries){
    sprintf(msg, "Inconsistent global id sizes: Num_GID = %d "
      "but global max is %d\n", lb->Num_GID, num_gid_entries);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_TRACE_EXIT(lb, yo);
    return LB_FATAL;
  }

  MPI_Allreduce(&(lb->Num_LID), &num_lid_entries, 1,
                MPI_INT, MPI_MAX, lb->Communicator);
  if (lb->Num_LID != num_lid_entries){
    sprintf(msg, "Inconsistent local id sizes: Num_LID = %d "
      "but global max is %d\n", lb->Num_LID, num_lid_entries);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_TRACE_EXIT(lb, yo);
    return LB_FATAL;
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
    import_proc_list = (int *) LB_MALLOC(num_import * sizeof(int));
    if (!import_proc_list) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }

    for (i = 0; i < num_import; i++) {
      proc_list[i] = import_procs[i];
      import_proc_list[i] = lb->Proc;
    }
  }

  /*
   *  Compute communication map and num_export, the number of objs this
   *  processor has to export to establish the new decomposition.
   */

  msgtag = 32767;
  ierr = LB_Comm_Create(&comm_plan, num_import, proc_list, lb->Communicator, 
                        msgtag, num_export);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Create.",
            (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&proc_list);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }
  

  LB_TRACE_DETAIL(lb, yo, "Done comm create");

  /*
   *  Allocate space for the object tags that need to be exported.  Communicate
   *  to get the list of objects to be exported.
   */

  if (*num_export > 0) {
    if (!LB_Special_Malloc(lb,(void **)export_global_ids,*num_export,
                           LB_SPECIAL_MALLOC_GID)) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }
    if (!LB_Special_Malloc(lb,(void **)export_local_ids,*num_export,
                           LB_SPECIAL_MALLOC_LID)) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_Special_Free(lb,(void **)export_global_ids,LB_SPECIAL_MALLOC_GID);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }
    if (!LB_Special_Malloc(lb,(void **)export_procs,*num_export,
                           LB_SPECIAL_MALLOC_INT)) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&proc_list);
      LB_Special_Free(lb,(void **)export_global_ids,LB_SPECIAL_MALLOC_GID);
      LB_Special_Free(lb,(void **)export_local_ids,LB_SPECIAL_MALLOC_LID);
      LB_TRACE_EXIT(lb, yo);
      return (LB_MEMERR);
    }
  }
  else {
    *export_global_ids = NULL;
    *export_local_ids = NULL;
    *export_procs = NULL;
  }

  /*
   *  Use the communication plan to send global IDs, local IDs, and processor
   *  numbers.  Do in separate communications to avoid a memory copy and to
   *  simplify implementation when a data type is added to the comm. package
   *  (to support heterogeneous computing).
   */

  msgtag2 = 32766;
  ierr = LB_Comm_Do(comm_plan, msgtag2, (char *) import_global_ids, 
                    (int) (sizeof(LB_ID_TYPE)*(num_gid_entries)), 
                    (char *) *export_global_ids);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Do.", 
            (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&proc_list);
    LB_Comm_Destroy(&comm_plan);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }

  if (num_lid_entries) {
    msgtag2--;
    ierr = LB_Comm_Do(comm_plan, msgtag2, (char *) import_local_ids, 
                      (int) (sizeof(LB_ID_TYPE)*num_lid_entries), 
                      (char *) *export_local_ids);
    if (ierr != COMM_OK && ierr != COMM_WARN) {
      sprintf(msg, "Error %s returned from LB_Comm_Do.", 
              (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
      LB_PRINT_ERROR(lb->Proc, yo, msg);
      LB_FREE(&proc_list);
      LB_Comm_Destroy(&comm_plan);
      LB_TRACE_EXIT(lb, yo);
      return (ierr == COMM_MEMERR ? LB_MEMERR : LB_FATAL);
    }
  }

  msgtag2--;
  ierr = LB_Comm_Do(comm_plan, msgtag2, (char *) import_proc_list, 
                    (int) sizeof(int), (char *) *export_procs);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Do.", 
            (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&proc_list);
    LB_Comm_Destroy(&comm_plan);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }
  LB_TRACE_DETAIL(lb, yo, "Done comm_do");

  LB_FREE(&proc_list);
  LB_FREE(&import_proc_list);
  
  LB_Comm_Destroy(&comm_plan);

  LB_TRACE_DETAIL(lb, yo, "Done comm destroy");

  LB_TRACE_EXIT(lb, yo);
  return (ierr);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int LB_Help_Migrate(
  LB *lb,                      /* Load balancing structure.                  */
  int num_import,              /* Number of non-local objects assigned to the 
                                  processor in the new decomposition.        */
  LB_ID_PTR import_global_ids, /* Array of global IDs for non-local objects 
                                  assigned to this processor in the new
                                  decomposition; this field can be NULL if 
                                  the application doesn't provide import IDs.*/
  LB_ID_PTR import_local_ids,  /* Array of local IDs for non-local objects
                                  assigned to the processor in the new
                                  decomposition; this field can be NULL if the 
                                  application does not provide import IDs.   */
  int *import_procs,           /* Array of processor IDs of processors owning
                                  the non-local objects that are assigned to
                                  this processor in the new decomposition; this
                                  field can be NULL if the application does
                                  not provide import IDs.                    */
  int num_export,              /* Number of objs to be exported
                                  to other processors to establish the new
                                  decomposition.                             */
  LB_ID_PTR export_global_ids, /* Array of global IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  LB_ID_PTR export_local_ids,  /* Array of local IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  int *export_procs            /* Array of processor IDs
                                  to which objects will be exported 
                                  to establish the new decomposition.        */
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
int num_gid_entries, num_lid_entries;  /* lengths of global & local ids */
int *sizes = NULL;       /* sizes (in bytes) of the object data for export. */
int id_size;             /* size (in bytes) of LB_GID + padding for 
                            alignment                                       */
int tag_size;            /* size (in bytes) of LB_GID + one int (for message size) */
char *export_buf = NULL; /* buffer for packing export data.                 */
char *import_buf = NULL; /* buffer for receiving imported data.             */
char *tmp;               /* temporary pointer into buffers.                 */
int i;                   /* loop counter.                                   */
int tmp_import;          /* number of objects to be imported.               */
int *proc_list = NULL;   /* list of processors to which this proc exports.  */
LB_ID_PTR tmp_id = NULL; /* pointer to storage for a global ID in comm buf  */
LB_ID_PTR lid;           /* temporary pointer to a local ID; used to pass
                            NULL to query functions when NUM_LID_ENTRIES=0. */
COMM_OBJ *comm_plan;     /* Object returned by communication routines       */
int msgtag, msgtag2;     /* Tags for communication routines                 */
int total_send_size;     /* Total size of outcoming message (in #items)     */
int total_recv_size;     /* Total size of incoming message (in #items)      */
int size;                /* size (in bytes) of an object                    */
int aligned_int;         /* size of an int padded for alignment             */
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

  /*
   *  Check that all procs use the same id types.
   */

  MPI_Allreduce(&(lb->Num_GID), &num_gid_entries, 1,
                MPI_INT, MPI_MAX, lb->Communicator);
  if (lb->Num_GID != num_gid_entries){
    sprintf(msg, "Inconsistent global id sizes: Num_GID = %d "
      "but global max is %d\n", lb->Num_GID, num_gid_entries);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_TRACE_EXIT(lb, yo);
    return LB_FATAL;
  }

  MPI_Allreduce(&(lb->Num_LID), &num_lid_entries, 1,
                MPI_INT, MPI_MAX, lb->Communicator);
  if (lb->Num_LID != num_lid_entries){
    sprintf(msg, "Inconsistent local id sizes: Num_LID = %d "
      "but global max is %d\n", lb->Num_LID, num_lid_entries);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_TRACE_EXIT(lb, yo);
    return LB_FATAL;
  }

  /*
   *  Check that all necessary query functions are available.
   */

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
                            num_gid_entries, num_lid_entries,
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
   * For each object, allow space for its global ID and its data plus 
   * one int (for the object data size).
   * Zoltan will pack the global IDs; the application must pack the data
   * through the pack routine.  Zoltan needs the global IDs for unpacking,
   * as the order of the data received during communication is not 
   * necessarily the same order as import_global_ids[].
   * Zoltan also needs to communicate the sizes of the objects because
   * only the sender knows the size of each object.
   */
  if (num_export > 0) {
    sizes = (int *) LB_MALLOC(num_export * sizeof(int));
    if (!sizes) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }
  }
  id_size = num_gid_entries * sizeof(LB_ID_TYPE) 
          + LB_pad_for_alignment(num_gid_entries * sizeof(LB_ID_TYPE));
  /* Note that the padding for alignment is not strictly necessary 
     when LB_ID_TYPE is int or unsigned int. */
  aligned_int = sizeof(int) + LB_pad_for_alignment(sizeof(int));
  tag_size = id_size + aligned_int;
  total_send_size = 0;

  for (i=0; i<num_export; i++){
    lid = (num_lid_entries ? &(export_local_ids[i*num_lid_entries]) : NULL);
    size = lb->Migrate.Get_Obj_Size(lb->Migrate.Get_Obj_Size_Data, 
                 num_gid_entries, num_lid_entries,
                 &(export_global_ids[i*num_gid_entries]), 
                 lid, &ierr);
    if (ierr) {
      LB_PRINT_ERROR(lb->Proc, yo, "Error returned from user defined "
                      "Migrate.Get_Obj_Size function.");
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }
    sizes[i] = size + LB_pad_for_alignment(size);
    total_send_size += sizes[i] + tag_size;
  }


  if (num_export > 0) {
    export_buf = (char *) LB_MALLOC(total_send_size);
    if (!export_buf) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&sizes);
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }

    proc_list = (int *) LB_MALLOC(num_export*sizeof(int));
    if (!proc_list) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&export_buf);
      LB_FREE(&sizes);
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
      tmp_id = (LB_ID_PTR) tmp;
      LB_SET_GID(lb, tmp_id, &(export_global_ids[i*num_gid_entries]));
      tmp += id_size;
    
      /* Pack the object's size */
      *((int *)tmp) = sizes[i];
      tmp += aligned_int;

      if (lb->Debug_Level >= LB_DEBUG_ALL){
        printf("[%1d] DEBUG in %s: Packing object with gid ", lb->Proc, yo);
        LB_PRINT_GID(lb, tmp_id);
        printf("size = %d bytes\n", sizes[i]); 
      }

      /* Pack the object's data */
      lid = (num_lid_entries ? &(export_local_ids[i*num_lid_entries]) : NULL);
      lb->Migrate.Pack_Obj(lb->Migrate.Pack_Obj_Data, 
                           num_gid_entries, num_lid_entries,
                           &(export_global_ids[i*num_gid_entries]),
                           lid, export_procs[i], sizes[i], tmp, &ierr);
      if (ierr) {
        LB_PRINT_ERROR(lb->Proc, yo, "Error returned from user defined "
                        "Migrate.Pack_Obj function.");
        LB_FREE(&sizes);
        LB_FREE(&export_buf);
        LB_FREE(&proc_list);
        LB_TRACE_EXIT(lb, yo);
        return (LB_FATAL);
      }
      tmp += sizes[i];
    }
  }

  LB_TRACE_DETAIL(lb, yo, "Done packing objects");

  /*
   *  Compute communication map and tmp_import, the number of objs this
   *  processor has to import to establish the new decomposition.
   */

  msgtag = 32767;
  ierr = LB_Comm_Create(&comm_plan, num_export, proc_list, lb->Communicator, 
                        msgtag, &tmp_import);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Create.", 
            (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&sizes);
    LB_FREE(&export_buf);
    LB_FREE(&proc_list);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }
  if (tmp_import != num_import) {
    sprintf(msg, "tmp_import %d != num_import %d.", tmp_import, num_import);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
  }

  /* Modify sizes[] to contain message sizes, not object sizes */
  for (i=0; i<num_export; i++)
    sizes[i] += tag_size;
  msgtag--;
  ierr = LB_Comm_Resize(comm_plan, sizes, msgtag, &total_recv_size);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Create.", 
            (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&sizes);
    LB_FREE(&export_buf);
    LB_FREE(&proc_list);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }

  if (num_import > 0) {
    import_buf = (char *) LB_MALLOC(total_recv_size);
    if (!import_buf) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&sizes);
      LB_FREE(&export_buf);
      LB_FREE(&proc_list);
      LB_TRACE_EXIT(lb, yo);
      return (LB_FATAL);
    }

  }

  /*
   *  Send the export data using the communication plan.
   */

  msgtag2 = 32765;
  ierr = LB_Comm_Do(comm_plan, msgtag2, export_buf, 1, import_buf);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    sprintf(msg, "Error %s returned from LB_Comm_Do.", 
            (ierr == COMM_MEMERR ? "COMM_MEMERR" : "COMM_FATAL"));
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    LB_FREE(&sizes);
    LB_FREE(&proc_list);
    LB_FREE(&export_buf);
    LB_FREE(&import_buf);
    LB_Comm_Destroy(&comm_plan);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == COMM_MEMERR ? LB_MEMERR : LB_FATAL);
  }

  /*
   *  Free whatever memory we can.
   */

  LB_Comm_Destroy(&comm_plan);
  LB_FREE(&proc_list);
  LB_FREE(&export_buf);
  LB_FREE(&sizes);

  LB_TRACE_DETAIL(lb, yo, "Done communication");

  /* 
   *  Perform application-specified processing before unpacking the data.
   */
  if (lb->Migrate.Mid_Migrate != NULL) {
    lb->Migrate.Mid_Migrate(lb->Migrate.Mid_Migrate_Data,
                            num_gid_entries, num_lid_entries,
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
    tmp_id = (LB_ID_PTR) tmp;
    tmp += id_size;

    /* Unpack the object's size */
    size = *((int *)tmp);
    tmp += aligned_int;

    if (lb->Debug_Level >= LB_DEBUG_ALL){
      printf("[%1d] DEBUG in %s: Unpacking object with gid ", lb->Proc, yo);
      LB_PRINT_GID(lb, tmp_id);
      printf("size = %d bytes\n", size);
    }

    /* Unpack the object's data */
    lb->Migrate.Unpack_Obj(lb->Migrate.Unpack_Obj_Data, num_gid_entries,
                           tmp_id, size, tmp, &ierr);
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
                             num_gid_entries, num_lid_entries,
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
