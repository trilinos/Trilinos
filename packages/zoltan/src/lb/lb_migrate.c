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

#include "zz_const.h"
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

int Zoltan_Compute_Destinations(
  ZZ *zz,                      /* Zoltan structure.                  */
  int num_import,              /* Number of non-local objects assigned to the 
                                  processor in the new decomposition.        */
  ZOLTAN_ID_PTR import_global_ids, /* Array of global IDs for non-local objects 
                                  assigned to this processor in the new
                                  decomposition.                             */
  ZOLTAN_ID_PTR import_local_ids,  /* Array of local IDs for non-local objects
                                  assigned to the processor in the new
                                  decomposition.                             */
  int *import_procs,           /* Array of processor IDs of processors owning
                                  the non-local objects that are assigned to
                                  this processor in the new decomposition.   */
  int *num_export,             /* Returned value:  Number of objs to be exported
                                  to other processors to establish the new
                                  decomposition.                             */
  ZOLTAN_ID_PTR *export_global_ids,/* Returned value:  Array of global IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  ZOLTAN_ID_PTR *export_local_ids, /* Returned value:  Array of local IDs of
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

char *yo = "Zoltan_Compute_Destinations";
char msg[256];
int *proc_list = NULL;      /* List of processors from which objs are to be 
                               imported.                                    */
ZOLTAN_COMM_OBJ *comm_plan;        /* Object returned communication routines  */
int *import_proc_list = NULL;
                            /* Array containing owning processor IDs of import
                               objects; used to request objs from other procs.*/
int msgtag, msgtag2;        /* Message tags for communication routines */
int num_gid_entries, num_lid_entries;  /* Length of global and local ids */
int i;
int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);
  /*
   *  Return if this processor is not in the Zoltan structure's
   *  communicator.
   */

  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_OK);
  }

  /*
   *  Check that all procs use the same id types.
   */

  MPI_Allreduce(&(zz->Num_GID), &num_gid_entries, 1,
                MPI_INT, MPI_MAX, zz->Communicator);
  if (zz->Num_GID != num_gid_entries){
    sprintf(msg, "Inconsistent global id sizes: Num_GID = %d "
      "but global max is %d\n", zz->Num_GID, num_gid_entries);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }

  MPI_Allreduce(&(zz->Num_LID), &num_lid_entries, 1,
                MPI_INT, MPI_MAX, zz->Communicator);
  if (zz->Num_LID != num_lid_entries){
    sprintf(msg, "Inconsistent local id sizes: Num_LID = %d "
      "but global max is %d\n", zz->Num_LID, num_lid_entries);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }

  /*
   *  Build processor's list of requests for non-local objs.
   */

  if (num_import > 0) {
    proc_list = (int *) ZOLTAN_MALLOC(num_import*sizeof(int));
    if (!proc_list) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_MEMERR);
    }
    import_proc_list = (int *) ZOLTAN_MALLOC(num_import * sizeof(int));
    if (!import_proc_list) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&proc_list);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_MEMERR);
    }

    for (i = 0; i < num_import; i++) {
      proc_list[i] = import_procs[i];
      import_proc_list[i] = zz->Proc;
    }
  }

  /*
   *  Compute communication map and num_export, the number of objs this
   *  processor has to export to establish the new decomposition.
   */

  msgtag = 32767;
  ierr = Zoltan_Comm_Create(&comm_plan, num_import, proc_list, zz->Communicator, 
                        msgtag, num_export);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Create.",
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&proc_list);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }
  

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done comm create");

  /*
   *  Allocate space for the object tags that need to be exported.  Communicate
   *  to get the list of objects to be exported.
   */

  if (*num_export > 0) {
    if (!Zoltan_Special_Malloc(zz,(void **)export_global_ids,*num_export,
                           ZOLTAN_SPECIAL_MALLOC_GID)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&proc_list);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_MEMERR);
    }
    if (!Zoltan_Special_Malloc(zz,(void **)export_local_ids,*num_export,
                           ZOLTAN_SPECIAL_MALLOC_LID)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&proc_list);
      Zoltan_Special_Free(zz,(void **)export_global_ids,ZOLTAN_SPECIAL_MALLOC_GID);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_MEMERR);
    }
    if (!Zoltan_Special_Malloc(zz,(void **)export_procs,*num_export,
                           ZOLTAN_SPECIAL_MALLOC_INT)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&proc_list);
      Zoltan_Special_Free(zz,(void **)export_global_ids,ZOLTAN_SPECIAL_MALLOC_GID);
      Zoltan_Special_Free(zz,(void **)export_local_ids,ZOLTAN_SPECIAL_MALLOC_LID);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_MEMERR);
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
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) import_global_ids, 
                    (int) (sizeof(ZOLTAN_ID_TYPE)*(num_gid_entries)), 
                    (char *) *export_global_ids);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Do.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&proc_list);
    Zoltan_Comm_Destroy(&comm_plan);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  if (num_lid_entries) {
    msgtag2--;
    ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) import_local_ids, 
                      (int) (sizeof(ZOLTAN_ID_TYPE)*num_lid_entries), 
                      (char *) *export_local_ids);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      sprintf(msg, "Error %s returned from Zoltan_Comm_Do.", 
              (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      ZOLTAN_FREE(&proc_list);
      Zoltan_Comm_Destroy(&comm_plan);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }
  }

  msgtag2--;
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) import_proc_list, 
                    (int) sizeof(int), (char *) *export_procs);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Do.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&proc_list);
    Zoltan_Comm_Destroy(&comm_plan);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }
  ZOLTAN_TRACE_DETAIL(zz, yo, "Done comm_do");

  ZOLTAN_FREE(&proc_list);
  ZOLTAN_FREE(&import_proc_list);
  
  Zoltan_Comm_Destroy(&comm_plan);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done comm destroy");

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_Help_Migrate(
  ZZ *zz,                      /* Zoltan structure.                  */
  int num_import,              /* Number of non-local objects assigned to the 
                                  processor in the new decomposition.        */
  ZOLTAN_ID_PTR import_global_ids, /* Array of global IDs for non-local objects 
                                  assigned to this processor in the new
                                  decomposition; this field can be NULL if 
                                  the application doesn't provide import IDs.*/
  ZOLTAN_ID_PTR import_local_ids,  /* Array of local IDs for non-local objects
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
  ZOLTAN_ID_PTR export_global_ids, /* Array of global IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  ZOLTAN_ID_PTR export_local_ids,  /* Array of local IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  int *export_procs            /* Array of processor IDs
                                  to which objects will be exported 
                                  to establish the new decomposition.        */
)
{
/*
 *  Routine to help perform migration.  If migration pre-processing routine
 *  (ZOLTAN_PRE_MIGRATE_FN) is specified, this routine first calls that fn.
 *  It then calls a function to obtain the size of the migrating objects
 *  (ZOLTAN_OBJ_SIZE_FN).  The routine next calls an application-specified
 *  object packing routine (ZOLTAN_PACK_OBJ_FN) for each object
 *  to be exported.  It develops the needed communication map to move the
 *  objects to other processors.  It performs the communication according
 *  to the map, and then calls an application-specified object unpacking
 *  routine (ZOLTAN_UNPACK_OBJ_FN) for each object imported.
 */

char *yo = "Zoltan_Help_Migrate";
char msg[256];
int num_gid_entries, num_lid_entries;  /* lengths of global & local ids */
int *sizes = NULL;       /* sizes (in bytes) of the object data for export. */
int id_size;             /* size (in bytes) of ZOLTAN_GID + padding for 
                            alignment                                       */
int tag_size;            /* size (in bytes) of ZOLTAN_GID + one int 
                            (for message size) */
char *export_buf = NULL; /* buffer for packing export data.                 */
char *import_buf = NULL; /* buffer for receiving imported data.             */
char *tmp;               /* temporary pointer into buffers.                 */
int i;                   /* loop counter.                                   */
int tmp_import;          /* number of objects to be imported.               */
int *proc_list = NULL;   /* list of processors to which this proc exports.  */
ZOLTAN_ID_PTR tmp_id = NULL; /* pointer to storage for a global ID in comm  
                                buf  */
ZOLTAN_ID_PTR lid;       /* temporary pointer to a local ID; used to pass
                            NULL to query functions when NUM_LID_ENTRIES=0. */
ZOLTAN_COMM_OBJ *comm_plan;     /* Object returned by communication routines*/
int msgtag, msgtag2;     /* Tags for communication routines                 */
int total_send_size;     /* Total size of outcoming message (in #items)     */
int total_recv_size;     /* Total size of incoming message (in #items)      */
int size;                /* size (in bytes) of an object                    */
int aligned_int;         /* size of an int padded for alignment             */
int ierr = 0;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /*
   *  Return if this processor is not in the Zoltan structure's
   *  communicator.
   */

  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_OK);
  }

  /*
   *  Return if export and import lists are not provided (through faulty
   *  value of RETURN_LISTS parameter).
   */

  if (num_export == -1) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Export lists must be "
                                 "provided; change RETURN_LISTS parameter.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_FATAL);
  }

  /*
   *  Check that all procs use the same id types.
   */

  MPI_Allreduce(&(zz->Num_GID), &num_gid_entries, 1,
                MPI_INT, MPI_MAX, zz->Communicator);
  if (zz->Num_GID != num_gid_entries){
    sprintf(msg, "Inconsistent global id sizes: Num_GID = %d "
      "but global max is %d\n", zz->Num_GID, num_gid_entries);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }

  MPI_Allreduce(&(zz->Num_LID), &num_lid_entries, 1,
                MPI_INT, MPI_MAX, zz->Communicator);
  if (zz->Num_LID != num_lid_entries){
    sprintf(msg, "Inconsistent local id sizes: Num_LID = %d "
      "but global max is %d\n", zz->Num_LID, num_lid_entries);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }

  /*
   *  Check that all necessary query functions are available.
   */

  if (zz->Get_Obj_Size == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register an "
           "ZOLTAN_OBJ_SIZE_FN_TYPE function to use the migration-help tools.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  if (zz->Pack_Obj == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register an "
           "ZOLTAN_PACK_OBJ_FN_TYPE function to use the migration-help tools.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  if (zz->Unpack_Obj == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register an "
         "ZOLTAN_UNPACK_OBJ_FN_TYPE function to use the migration-help tools.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  if (zz->Migrate.Pre_Migrate != NULL) {
    zz->Migrate.Pre_Migrate(zz->Migrate.Pre_Migrate_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                      "Migrate.Pre_Migrate function.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done pre-migration processing");

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
    sizes = (int *) ZOLTAN_MALLOC(num_export * sizeof(int));
    if (!sizes) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }
  }
  id_size = Zoltan_Align(num_gid_entries * sizeof(ZOLTAN_ID_TYPE));
  /* Note that alignment is not strictly necessary 
     when ZOLTAN_ID_TYPE is int or unsigned int. */
  aligned_int = Zoltan_Align(sizeof(int));
  tag_size = id_size + aligned_int;
  total_send_size = 0;

  for (i=0; i<num_export; i++){
    lid = (num_lid_entries ? &(export_local_ids[i*num_lid_entries]) : NULL);
    size = zz->Get_Obj_Size(zz->Get_Obj_Size_Data, 
                 num_gid_entries, num_lid_entries,
                 &(export_global_ids[i*num_gid_entries]), 
                 lid, &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                      "Get_Obj_Size function.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }
    sizes[i] = Zoltan_Align(size);
    total_send_size += sizes[i] + tag_size;
  }


  if (num_export > 0) {
    export_buf = (char *) ZOLTAN_MALLOC(total_send_size);
    if (!export_buf) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&sizes);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }

    proc_list = (int *) ZOLTAN_MALLOC(num_export*sizeof(int));
    if (!proc_list) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&export_buf);
      ZOLTAN_FREE(&sizes);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }


    /*
     *  Pack the proc_list (to create the map) and the objects for export.
     */
  
    tmp = export_buf;
    for (i = 0; i < num_export; i++) {
      proc_list[i] = export_procs[i];

      /* Pack the object's global ID */
      tmp_id = (ZOLTAN_ID_PTR) tmp;
      ZOLTAN_SET_GID(zz, tmp_id, &(export_global_ids[i*num_gid_entries]));
      tmp += id_size;
    
      /* Pack the object's size */
      *((int *)tmp) = sizes[i];
      tmp += aligned_int;

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] DEBUG in %s: Packing object with gid ", zz->Proc, yo);
        ZOLTAN_PRINT_GID(zz, tmp_id);
        printf("size = %d bytes\n", sizes[i]); 
      }

      /* Pack the object's data */
      lid = (num_lid_entries ? &(export_local_ids[i*num_lid_entries]) : NULL);
      zz->Pack_Obj(zz->Pack_Obj_Data, 
                           num_gid_entries, num_lid_entries,
                           &(export_global_ids[i*num_gid_entries]),
                           lid, export_procs[i], sizes[i], tmp, &ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                        "Pack_Obj function.");
        ZOLTAN_FREE(&sizes);
        ZOLTAN_FREE(&export_buf);
        ZOLTAN_FREE(&proc_list);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ZOLTAN_FATAL);
      }
      tmp += sizes[i];
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done packing objects");

  /*
   *  Compute communication map and tmp_import, the number of objs this
   *  processor has to import to establish the new decomposition.
   */

  msgtag = 32767;
  ierr = Zoltan_Comm_Create(&comm_plan, num_export, proc_list, zz->Communicator, 
                        msgtag, &tmp_import);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Create.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&sizes);
    ZOLTAN_FREE(&export_buf);
    ZOLTAN_FREE(&proc_list);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }
  if ((num_import != -1) && (tmp_import != num_import)) {
    sprintf(msg, "tmp_import %d != num_import %d.", tmp_import, num_import);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
  }

  /* Modify sizes[] to contain message sizes, not object sizes */
  for (i=0; i<num_export; i++)
    sizes[i] += tag_size;
  msgtag--;
  ierr = Zoltan_Comm_Resize(comm_plan, sizes, msgtag, &total_recv_size);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Create.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&sizes);
    ZOLTAN_FREE(&export_buf);
    ZOLTAN_FREE(&proc_list);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  if (tmp_import > 0) {
    import_buf = (char *) ZOLTAN_MALLOC(total_recv_size);
    if (!import_buf) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&sizes);
      ZOLTAN_FREE(&export_buf);
      ZOLTAN_FREE(&proc_list);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }

  }

  /*
   *  Send the export data using the communication plan.
   */

  msgtag2 = 32765;
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, export_buf, 1, import_buf);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Do.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&sizes);
    ZOLTAN_FREE(&proc_list);
    ZOLTAN_FREE(&export_buf);
    ZOLTAN_FREE(&import_buf);
    Zoltan_Comm_Destroy(&comm_plan);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  /*
   *  Free whatever memory we can.
   */

  Zoltan_Comm_Destroy(&comm_plan);
  ZOLTAN_FREE(&proc_list);
  ZOLTAN_FREE(&export_buf);
  ZOLTAN_FREE(&sizes);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done communication");

  /* 
   *  Perform application-specified processing before unpacking the data.
   */
  if (zz->Migrate.Mid_Migrate != NULL) {
    zz->Migrate.Mid_Migrate(zz->Migrate.Mid_Migrate_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                      "Migrate.Mid_Migrate function.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }

    ZOLTAN_TRACE_DETAIL(zz, yo, "Done mid-migration processing");
  }

  /*
   *  Unpack the object data.
   */

  tmp = import_buf;
  for (i = 0; i < tmp_import; i++) {

    /* Unpack the object's global ID */
    tmp_id = (ZOLTAN_ID_PTR) tmp;
    tmp += id_size;

    /* Unpack the object's size */
    size = *((int *)tmp);
    tmp += aligned_int;

    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
      printf("[%1d] DEBUG in %s: Unpacking object with gid ", zz->Proc, yo);
      ZOLTAN_PRINT_GID(zz, tmp_id);
      printf("size = %d bytes\n", size);
    }

    /* Unpack the object's data */
    zz->Unpack_Obj(zz->Unpack_Obj_Data, num_gid_entries,
                           tmp_id, size, tmp, &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                      "Unpack_Obj function.");
      ZOLTAN_FREE(&import_buf);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }
    tmp += size;
  }

  ZOLTAN_FREE(&import_buf);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done unpacking objects");

  if (zz->Migrate.Post_Migrate != NULL) {
    zz->Migrate.Post_Migrate(zz->Migrate.Post_Migrate_Data,
                             num_gid_entries, num_lid_entries,
                             num_import, import_global_ids,
                             import_local_ids, import_procs,
                             num_export, export_global_ids,
                             export_local_ids, export_procs, &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                      "Migrate.Post_Migrate function.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }

    ZOLTAN_TRACE_DETAIL(zz, yo, "Done post-migration processing");
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ZOLTAN_OK);
}
