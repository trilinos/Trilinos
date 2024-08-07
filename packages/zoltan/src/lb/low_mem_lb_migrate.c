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
#include "zz_util_const.h"

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

static void printobj(int me, int nprocs, 
   int num_import, unsigned int *import_global_ids,
   unsigned int *import_local_ids, int *import_procs, int *import_to_part,
   int num_export, unsigned int *export_global_ids,
   unsigned int *export_local_ids, int *export_procs, int *export_to_part)
{
int i, j;

  for (i=0; i < nprocs; i++){

    if (i == me){
    fprintf(stderr,"(%d) %d exports %d imports\n",me,num_export,num_import);
    fprintf(stderr,"(%d) -> ",me);
    for (j=0 ; j < num_export; j++){
      fprintf(stderr,"%d/%d/%d/%d ",export_global_ids[j],export_local_ids[j],export_procs[j],export_to_part[j]);
    }
    fprintf(stderr,"\n");
    fprintf(stderr,"(%d) >- ",me);
    for (j=0 ; j < num_import; j++){
      fprintf(stderr,"%d/%d/%d/%d ",import_global_ids[j],import_local_ids[j],import_procs[j],import_to_part[j]);
    }
    fprintf(stderr,"\n");
    }

    MPI_Barrier(zoltan_get_global_comm());
    MPI_Barrier(zoltan_get_global_comm());
  }

}


static int check_input(ZZ *, int, int *);
static int actual_arrays(ZZ *, int, int, int, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
  int *, int *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **, int *);

static int get_int_array(ZZ *zz, int len, int **buf, int *oldlen);
static int get_char_array(ZZ *zz, int len, char **buf, int *oldlen);
static int get_gid_array(ZZ *zz, int len, ZOLTAN_ID_TYPE **buf, int *oldlen);
static void reset_arrays();

static int pack_message_for_proc(ZZ *zz, int dest_proc, int num,
             int num_export, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
             int *to_proc, int *to_part,
             int tag_size, int id_size, int aligned_int,
             char **send_buf, int *send_size,
             ZOLTAN_ID_PTR *egids, ZOLTAN_ID_PTR *elids, int **eprocs, int **eparts);

static int message_receive(MPI_Comm comm, int src, char *buf, int len);
static int message_send(MPI_Comm comm, int dest, char *buf, int len);
static int message_wait(int len);

static int create_import_lists(ZZ *zz, int msgLen, int src, int num,
                  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int *procs, int *parts,
                  int *objectCount,
                  ZOLTAN_ID_PTR *inGids, ZOLTAN_ID_PTR *inLids, int **inProcs, int **inParts);

/* message passing arrays */
static int *_sizes = NULL;
static int _sizes_len = 0;
static int *_idx = NULL;
static int _idx_len = 0;
static char *_export_buf= NULL;
static int _export_buf_len = 0;
static char *_import_buf= NULL;
static int _import_buf_len = 0;
static ZOLTAN_ID_PTR _gid_list = NULL;
static int _gid_list_len = 0;
static ZOLTAN_ID_PTR _lid_list = NULL;
static int _lid_list_len = 0;
static int *_proc_list = NULL;
static int _proc_list_len = 0;
static int *_part_list = NULL;
static int _part_list_len = 0;

/* mid-migrate arrays */
static ZOLTAN_ID_PTR _import_gids= NULL;
static int _import_gids_len = 0;
static ZOLTAN_ID_PTR _import_lids = NULL;
static int _import_lids_len = 0;
static int *_import_procs = NULL;
static int _import_procs_len = 0;
static int *_import_parts = NULL;
static int _import_parts_len = 0;

static int ackTag = 10000, dataTag = 20000, sizeTag=30000;
static MPI_Request _request;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Low_Mem_Migrate(
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
  int *import_to_part,         /* Array of partition numbers to which imported
                                  objects should be assigned.                */
  int num_export,              /* Number of objs to be exported
                                  to other processors to establish the new
                                  decomposition.                             */
  ZOLTAN_ID_PTR export_global_ids, /* Array of global IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  ZOLTAN_ID_PTR export_local_ids,  /* Array of local IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  int *export_procs,           /* Array of processor IDs
                                  to which objects will be exported 
                                  to establish the new decomposition.        */
  int *export_to_part          /* Array of partition numbers to which exported
                                  objects should be assigned.                */
)
{
/*
 *  Low memory version of Migrate.  Buffers are allocated to send to one process
 *  and to receive from one process.  The send of objects to another process is
 *  an MPI "ready" send.  The send does not happen until the receive is posted,
 *  so that MPI does not have to allocate a buffer on the send side.
 *
 *  Uses only the PP (proc/part) versions of pre-, mid-, post-migrate query functions.
 *
 *  Mid migrate query is called after each receive/send pair, not after all
 *     sends and receives are completed.  It contains the list only of those
 *     objects that were received and sent in that step.
 *
 *  zdrive's mid migrate function assumes that all objects have been sent out, so
 *     this migrate utility does not work with zdrive.  zdrive overwrites the
 *     sent out objects with the received objects, before the objects have been
 *     sent out.
 *
 *  Unpack is called after each receive/send pair, not after all
 *     sends and receives are completed.  It contains only the objects
 *     that were received in that step.
 */

char *yo = "Zoltan_Lo_Mem_Migrate";
int ierr = ZOLTAN_OK;

int num_gid_entries, num_lid_entries;  /* lengths of global & local ids */
int *sizes = NULL;       /* sizes (in bytes) of the object data for export. */
int id_size;             /* size (in bytes) of ZOLTAN_GID + padding for 
                            alignment                                       */
int tag_size;            /* size (in bytes) of ZOLTAN_GID + one int 
                            (for message size) */
int aligned_int;         /* size of an int padded for alignment             */
ZOLTAN_COMM_OBJ *imp_plan = NULL; /* Comm obj built from import lists. */
ZOLTAN_COMM_OBJ *exp_plan = NULL; /* Comm obj built from export lists. */
int include_parts = 0;   /* flag indicating whether partition info is
                            provided */
int actual_num_exp = 0;
int actual_exp_allocated = 0;
ZOLTAN_ID_PTR actual_exp_gids = NULL;    /* Arrays containing only objs to  */
ZOLTAN_ID_PTR actual_exp_lids = NULL;    /* actually be packed.  Objs that  */
int *actual_exp_procs = NULL;            /* are changing partition but not  */
int *actual_exp_to_part = NULL;          /* processor may not be included.  */
int actual_num_imp = 0;
int actual_imp_allocated = 0;
ZOLTAN_ID_PTR actual_imp_gids = NULL;    /* Arrays containing only objs to  */
ZOLTAN_ID_PTR actual_imp_lids = NULL;    /* actually be imported. Objs that  */
int *actual_imp_procs = NULL;            /* are changing partition but not  */
int *actual_imp_to_part = NULL;          /* processor may not be included.  */

int i, k, n;
int left, right, start_offset;
int *dummyPtr=NULL, *idx=NULL, *numGIDs=NULL;
int idx_cnt, tmp_size;
int destProcMap=-1;
int numDestProcs;
int nprocs = zz->Num_Proc;
int rank = zz->Proc;
MPI_Comm comm = zz->Communicator;
MPI_Request request;
MPI_Status status;
int need_import_lists = 0;
int need_export_lists = 0;
ZOLTAN_ID_PTR mid_migrate_export_gids = NULL; 
ZOLTAN_ID_PTR mid_migrate_export_lids = NULL;
int *mid_migrate_export_procs = NULL; 
int *mid_migrate_export_parts = NULL;
ZOLTAN_ID_PTR mid_migrate_import_gids = NULL; 
ZOLTAN_ID_PTR mid_migrate_import_lids = NULL;
int *mid_migrate_import_procs = NULL; 
int *mid_migrate_import_parts = NULL;
int inBytes, outBytes, inCount, outCount;
char *inBuf=NULL, *outBuf=NULL, *tmp;
int msgtag, msgtag2;
ZOLTAN_ID_PTR tmp_id;
intptr_t dummyVal;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /*
   *  Return if this processor is not in the Zoltan structure's
   *  communicator.
   */

  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz)) {
    goto End;
  }

  /*
   *  Check that all procs use the same id types.
   */

  ierr = check_input(zz, 
                    ((num_export >= 0 && export_to_part) || 
                     (num_import >= 0 && import_to_part)),
                     &include_parts);
  if (ierr != ZOLTAN_OK) 
    goto End;

  num_gid_entries = zz->Num_GID;
  num_lid_entries = zz->Num_LID;

  /*
   *  Check that all necessary query functions are available.
   */

  if (zz->Get_Obj_Size == NULL && zz->Get_Obj_Size_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
           "ZOLTAN_OBJ_SIZE_FN or ZOLTAN_OBJ_SIZE_MULTI_FN function "
           "to use the migration-help tools.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->Pack_Obj == NULL && zz->Pack_Obj_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
           "ZOLTAN_PACK_OBJ_FN or ZOLTAN_PACK_OBJ_MULTI_FN function "
           "to use the migration-help tools.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->Unpack_Obj == NULL && zz->Unpack_Obj_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
         "ZOLTAN_UNPACK_OBJ_FN or ZOLTAN_UNPACK_OBJ_MULTI_FN function "
         "to use the migration-help tools.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /* Create lists of objects to be exported */

  if (num_export >= 0) {

    if (zz->Migrate.Mid_Migrate_PP != NULL) {
      need_export_lists = 1;
    }

    /* Build the actual export arrays */
    ierr = actual_arrays(zz, num_gid_entries, num_lid_entries,
                         num_export, export_global_ids, export_local_ids, 
                         export_procs, export_to_part, 
                         &actual_num_exp, &actual_exp_gids, &actual_exp_lids,
                         &actual_exp_procs, &actual_exp_to_part,
                         &actual_exp_allocated);
    if (ierr < 0) 
      goto End;

  }

  else if (num_import >= 0) {

    if (zz->Migrate.Mid_Migrate_PP != NULL) {
      need_import_lists = 1;
    }

    /* Build the actual import arrays */
    ierr = actual_arrays(zz, num_gid_entries, num_lid_entries,
                         num_import, import_global_ids, import_local_ids, 
                         import_procs, import_to_part, 
                         &actual_num_imp, &actual_imp_gids, &actual_imp_lids,
                         &actual_imp_procs, &actual_imp_to_part,
                         &actual_imp_allocated);
    if (ierr < 0) 
      goto End;
    
    /* Compute communication map based on imports.  */
    msgtag = 32767;
    ierr = Zoltan_Comm_Create(&imp_plan, actual_num_imp, actual_imp_procs,
                              zz->Communicator, msgtag, &actual_num_exp);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc,yo,"Error returned from Zoltan_Comm_Create.");
      goto End;
    }

    /* Compute actual export lists for packing objects */
    if (actual_num_exp > 0) {
      actual_exp_allocated = 1;
      actual_exp_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, actual_num_exp);
      actual_exp_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, actual_num_exp);
      actual_exp_procs = (int *) ZOLTAN_MALLOC(sizeof(int) * actual_num_exp);
      if (include_parts)
        actual_exp_to_part = (int *) ZOLTAN_MALLOC(sizeof(int)*actual_num_exp);
      if (actual_exp_gids == NULL ||
          (num_lid_entries && actual_exp_lids == NULL) ||
          actual_exp_procs == NULL || 
          (import_to_part != NULL && actual_exp_to_part == NULL)) {
        Zoltan_Multifree(__FILE__, __LINE__, 4, 
                         &actual_exp_gids, &actual_exp_lids, 
                         &actual_exp_procs, &actual_exp_to_part);
        ierr = ZOLTAN_MEMERR;
        goto End;
      }

    }

    msgtag2 = 32766;
    ierr = Zoltan_Comm_Do(imp_plan, msgtag2, (char *) actual_imp_gids,
                    (int) (sizeof(ZOLTAN_ID_TYPE)*(num_gid_entries)),
                    (char *) actual_exp_gids);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
      goto End;
    }

    if (num_lid_entries) {
      msgtag2--;
      ierr = Zoltan_Comm_Do(imp_plan, msgtag2, (char *) actual_imp_lids,
                      (int) (sizeof(ZOLTAN_ID_TYPE)*num_lid_entries),
                      (char *) actual_exp_lids);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
        goto End;
      }
    }

    Zoltan_Comm_Info(imp_plan, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, actual_exp_procs, NULL);
 
    if (include_parts) {
      msgtag2--;
      ierr = Zoltan_Comm_Do(imp_plan, msgtag2, (char *) actual_imp_to_part,
                      (int) sizeof(int), (char *) actual_exp_to_part);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
        goto End;
      }
    }
   
    /* Create inverse plan (i.e., plan based on exports) so can set 
     * variable sizes. 
     * (Zoltan_Comm_Do_Reverse(imp_plan, ...) allows sending variable 
     * but does not tell how large to allocate receive buffer.
     */
    ierr = Zoltan_Comm_Invert_Plan(&imp_plan);
    if (ierr != ZOLTAN_OK){
      goto End;
    }
    exp_plan = imp_plan;
    imp_plan = NULL;
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Import or export lists needed.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }


  /* Create a searchable set of the process ranks we send to
   */

  destProcMap = Zoltan_Map_Create(zz, nprocs, sizeof(int), 1, 0);

  if (destProcMap < 0){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to create first map.");
    ierr = ZOLTAN_FATAL; 
    goto End;
  }

  for (i=0; i < actual_num_exp; i++){
     ierr = Zoltan_Map_Add(zz, destProcMap, (char *)(actual_exp_procs + i), dummyVal);
     if (ierr != ZOLTAN_OK){
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to add proc to first map.");
       goto End;
     }
  }

  numDestProcs = Zoltan_Map_Size(zz, destProcMap);

  ierr = Zoltan_Map_Destroy(zz, destProcMap);
  if (ierr != ZOLTAN_OK){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to destroy first map.");
    goto End;
  }
  destProcMap = -1;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_TRACE_DETAIL){
     printf("%s (%d) will send objects to %d processes\n", yo, zz->Proc, numDestProcs);
  }

  /* Now create searchable map from process rank to number of GIDs 
   */

  destProcMap = Zoltan_Map_Create(zz, nprocs, 1, sizeof(int), 1, 0);

  if (destProcMap < 0){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to create second map.");
    ierr = ZOLTAN_FATAL; 
    goto End;
  }

  if (numDestProcs > 0){
    numGIDs = (int *)ZOLTAN_CALLOC(sizeof(int) , numDestProcs);
    if (!numGIDs){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i=0,n=0; i < actual_num_exp; i++){

      ierr = Zoltan_Map_Find(zz, destProcMap, (char *)(actual_exp_procs + i), &dummyPtr);
      if (ierr != ZOLTAN_OK){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to search for proc in second map.");
        goto End;
      }
      
      if (dummyPtr != ZOLTAN_NOT_FOUND){
        dummyPtr++;
      }
      else{
        ierr = Zoltan_Map_Add(zz, destProcMap, (char *)(actual_exp_procs + i), (intptr_t)(numGIDs + n));
        if (ierr != ZOLTAN_OK){
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to add proc to second map.");
          goto End;
        }
        numGIDs[n++] = 1;
      }
    }
  }

  /* Do the application's pre-migrate step */

  if (zz->Migrate.Pre_Migrate_PP != NULL) {
    zz->Migrate.Pre_Migrate_PP(zz->Migrate.Pre_Migrate_PP_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs, import_to_part,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, export_to_part,
                            &ierr);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                      "ZOLTAN_PRE_MIGRATE_PP_FN function.");
      goto End;
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_TRACE_DETAIL){
       printf("%s (%d) pre migrate complete num_import %d num_export %d\n", yo, zz->Proc, 
              num_import, num_export);
    }
  }

  id_size = Zoltan_Align(num_gid_entries * sizeof(ZOLTAN_ID_TYPE));
  aligned_int = Zoltan_Align(sizeof(int));
  tag_size = id_size + aligned_int;

  /*
   * Bucket brigade - at each step, a process receives data from its
   * left and sends data to its right, from and to a different process
   * each time. 
   */

  start_offset = 1;

  if (!(zz->Migrate.Only_Proc_Changes)) {

    ierr = Zoltan_Map_Find(zz, destProcMap, (char *)&rank, &dummyPtr);

    if (ierr != ZOLTAN_OK){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to search second map.");
      goto End;
    }

    if (*dummyPtr != ZOLTAN_NOT_FOUND) {
      /* 
       * send myself objects that are going to a different partition on my process 
       */
      start_offset = 0;
    }
  }

  for (i=start_offset; i < nprocs; i++){ 

    right = (rank + i) % nprocs;
    left = (nprocs + rank - i) % nprocs;

    inCount = outCount = 0;   /* number of objects */
    inBytes = outBytes = 0;   /* number of bytes */
    inBuf = outBuf = NULL;
    mid_migrate_export_gids = mid_migrate_import_gids = NULL; 
    mid_migrate_export_lids = mid_migrate_import_lids = NULL;
    mid_migrate_export_procs = mid_migrate_import_procs = NULL;
    mid_migrate_export_parts = mid_migrate_import_parts = NULL;

    /* Post receive for size of message I'll get from the left */

    MPI_Irecv(&inBytes, 1, MPI_INT, left, sizeTag, comm, &request);

    /* Compute the size of the data I'll send to the right */

    ierr = Zoltan_Map_Find(zz, destProcMap, (char *)&right, &dummyPtr);

    if (ierr != ZOLTAN_OK){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to search for proc in second map.");
      goto End;
    }

    if (dummyPtr != ZOLTAN_NOT_FOUND){   /* "right" was found in the map */
      outCount = dummyPtr;  /* number of objects I send to "right" */

      if (need_export_lists){
        ierr = pack_message_for_proc(zz, right, outCount, 
                  num_export, export_global_ids, export_local_ids, 
                  export_procs, export_to_part,
                  tag_size, id_size, aligned_int,
                  &outBuf, &outBytes,
                  &mid_migrate_export_gids, &mid_migrate_export_lids, 
                  &mid_migrate_export_procs, &mid_migrate_export_parts);
      }
      else{
        ierr = pack_message_for_proc(zz, right, outCount, 
                  num_export, export_global_ids, export_local_ids, 
                  export_procs, export_to_part,
                  tag_size, id_size, aligned_int,
                  &outBuf, &outBytes,
                  NULL, NULL, NULL, NULL);
      }

      if (ierr != ZOLTAN_OK){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from pack_message_for_proc");
        goto End;
      }
    }

    MPI_Send(&outBytes, 1, MPI_INT, right, sizeTag, comm);

    MPI_Wait(&request, &status);

    if (inBytes > 0){
      ierr = get_char_array(zz, inBytes, &_import_buf, &_import_buf_len);
      if (ierr != ZOLTAN_OK){
        goto End;  
      }
      inBuf = _import_buf;
    }
    else{
      inBuf = NULL;
    }

    if (left != right){
      message_receive(comm, left, inBuf, inBytes);   /* post receive */
      message_send(comm, right, outBuf, outBytes);   /* ready send */
    }
    else {
      /* 
       * only happens if nprocs < 3, or if proc sends-to/receives-from self
       */
      if (rank < left){ 
        message_send(comm, right, outBuf, outBytes);
        message_receive(comm, left, inBuf, inBytes);
      }
      else{
        message_receive(comm, left, inBuf, inBytes);
        message_send(comm, right, outBuf, outBytes);
      }
    }

    if (zz->Debug_Level >= ZOLTAN_DEBUG_TRACE_DETAIL){
       printf("%s (%d) send %d to %d, got %d from %d\n", yo, zz->Proc, 
              outBytes, right, inBytes, left);
    }

    /* 
     * We have posted the receive from the left.  We don't need the objects from 
     * the left in order to create the import lists.  So we do this computation 
     * before waiting for the objects.
     */
    if (need_import_lists){
      ierr = create_import_lists(zz, inBytes, left, num_import,
                  import_global_ids, import_local_ids, import_procs, import_to_part,
                  &inCount,
                  &mid_migrate_import_gids, &mid_migrate_import_lids, 
                  (import_procs ? &mid_migrate_import_procs : NULL),
                  (import_to_part ? &mid_migrate_import_parts : NULL));

      if (ierr != ZOLTAN_OK){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from create_import_lists");
        goto End;
      }
    }
    else if (inBytes > 0){   /* still need count of received objects */
      for (k=0; k < num_import; k++){
        if (import_procs[k] == left)
          inCount++;
      }
    }

    message_wait(inBytes);     /* wait for objects from the left */

    /* 
     *  Perform application-specified processing before unpacking the data.
     */
    if (zz->Migrate.Mid_Migrate_PP != NULL) {
      zz->Migrate.Mid_Migrate_PP(zz->Migrate.Mid_Migrate_PP_Data,
                              num_gid_entries, num_lid_entries,
                              inCount,
                              mid_migrate_import_gids, mid_migrate_import_lids, 
                              mid_migrate_import_procs, mid_migrate_import_parts,
                              outCount,
                              mid_migrate_export_gids, mid_migrate_export_lids, 
                              mid_migrate_export_procs, mid_migrate_export_parts,
                              &ierr);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                        "ZOLTAN_MID_MIGRATE_PP_FN function.");
        goto End;
      }
    }

    /*
     *  Unpack the object data.
     */
  
    if (inBytes > 0) {
  
      if (zz->Unpack_Obj_Multi != NULL) {
  
        /* Allocate and fill input arrays for Unpack_Obj_Multi. */

        ierr = get_int_array(zz, inCount, &_sizes, &_sizes_len);
        if (!ierr) goto End;
        sizes = _sizes;

        if (zz->Pack_Obj_Multi != NULL) {
          ierr = get_int_array(zz, inCount, &_idx, &_idx_len);
          if (!ierr) goto End;
          idx = _idx;
        }
        else{
          idx = NULL;
        }

        ierr = get_gid_array(zz, inCount, &_gid_list, &_gid_list_len);
        if (!ierr) goto End;
        tmp_id = _gid_list;

        tmp = inBuf;
        idx_cnt = 0;
        for (k = 0; k < inCount; k++) {
  
          /* Unpack the object's global ID */
          ZOLTAN_SET_GID(zz, &(tmp_id[k*num_gid_entries]), (ZOLTAN_ID_PTR) tmp);
          tmp += id_size;
  
          /* Unpack the object's size */
          sizes[k] = *((int *)tmp);
          tmp += aligned_int;
  
          /* If using ZOLTAN_UNPACK_OBJ_MULTI_FN, build the index array. */
          idx_cnt += tag_size;
          if (idx != NULL) {
            idx[k] = idx_cnt;
          }
  
          tmp += sizes[k];
          idx_cnt += sizes[k];
        }
  
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
          printf("[%1d] DEBUG in %s: Unpacking objects with multi-fn\n",
                 zz->Proc,yo);
        }
        zz->Unpack_Obj_Multi(zz->Unpack_Obj_Multi_Data, num_gid_entries,
                           inCount, tmp_id, sizes, idx, inBuf, &ierr);
        
        if (ierr < 0) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                          "ZOLTAN_UNPACK_OBJ_MULTI_FN.");
          goto End;
        }
      }
      else {
        tmp = inBuf;
        for (k = 0; k < inCount; k++) {
          tmp_size = *((int *)(tmp + id_size));
          if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
            printf("[%1d] DEBUG in %s: Unpacking object with gid ", zz->Proc, yo);
            ZOLTAN_PRINT_GID(zz, (ZOLTAN_ID_PTR)tmp);
            printf("size = %d bytes\n", tmp_size);
          }
  
          /* Unpack the object's data */
         
          zz->Unpack_Obj(zz->Unpack_Obj_Data, num_gid_entries,
                         (ZOLTAN_ID_PTR) tmp, tmp_size,
                         tmp + tag_size, &ierr);
          if (ierr < 0) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                            "ZOLTAN_UNPACK_OBJ_FN.");
            goto End;
          }
          tmp += (tmp_size + tag_size);
        }
      }
    }
  } /* Next left and right process */

  if (zz->Migrate.Post_Migrate_PP != NULL) {
    zz->Migrate.Post_Migrate_PP(zz->Migrate.Post_Migrate_PP_Data,
                             num_gid_entries, num_lid_entries,
                             num_import, import_global_ids,
                             import_local_ids, import_procs, import_to_part,
                             num_export, export_global_ids,
                             export_local_ids, export_procs, export_to_part,
                             &ierr);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                      "ZOLTAN_POST_MIGRATE_PP_FN function.");
      goto End;
    }

    if (zz->Debug_Level >= ZOLTAN_DEBUG_TRACE_DETAIL){
       printf("%s (%d) post migrate complete num_import %d num_export %d\n", yo, zz->Proc, 
              num_import, num_export);
    }
  }

End:

  Zoltan_Map_Destroy(zz, destProcMap);
  ZOLTAN_FREE(&numGIDs);
  reset_arrays();

  if (actual_exp_allocated) {
    Zoltan_Multifree(__FILE__, __LINE__, 4, 
                     &actual_exp_gids, &actual_exp_lids, 
                     &actual_exp_procs, &actual_exp_to_part);
  }
  if (actual_imp_allocated) {
    Zoltan_Multifree(__FILE__, __LINE__, 4, 
                     &actual_imp_gids, &actual_imp_lids, 
                     &actual_imp_procs, &actual_imp_to_part);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Migrate(
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
  int *import_to_part,         /* Array of partition numbers to which imported
                                  objects should be assigned.                */
  int num_export,              /* Number of objs to be exported
                                  to other processors to establish the new
                                  decomposition.                             */
  ZOLTAN_ID_PTR export_global_ids, /* Array of global IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  ZOLTAN_ID_PTR export_local_ids,  /* Array of local IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  int *export_procs,           /* Array of processor IDs
                                  to which objects will be exported 
                                  to establish the new decomposition.        */
  int *export_to_part          /* Array of partition numbers to which exported
                                  objects should be assigned.                */
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

char *yo = "Zoltan_Migrate";
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
int tmp_size;            /* size of a single object's data.                 */
int *idx = NULL;         /* index used for multi-fn packs and unpacks.      */
int idx_cnt = 0;         /* index counter for idx array.                    */
ZOLTAN_ID_PTR tmp_id = NULL; /* pointer to storage for a global ID in comm  
                                buf  */
ZOLTAN_ID_PTR lid;       /* temporary pointer to a local ID; used to pass
                            NULL to query functions when NUM_LID_ENTRIES=0. */
ZOLTAN_COMM_OBJ *imp_plan = NULL; /* Comm obj built from import lists. */
ZOLTAN_COMM_OBJ *exp_plan = NULL; /* Comm obj built from export lists. */
int msgtag, msgtag2;     /* Tags for communication routines                 */
int total_send_size;     /* Total size of outcoming message (in #items)     */
int total_recv_size;     /* Total size of incoming message (in #items)      */
int aligned_int;         /* size of an int padded for alignment             */
int dest;                /* temporary destination partition.                */
int include_parts = 0;   /* flag indicating whether partition info is
                            provided */
int ierr = ZOLTAN_OK;
int actual_num_exp = 0;
int actual_exp_allocated = 0;
ZOLTAN_ID_PTR actual_exp_gids = NULL;    /* Arrays containing only objs to  */
ZOLTAN_ID_PTR actual_exp_lids = NULL;    /* actually be packed.  Objs that  */
int *actual_exp_procs = NULL;            /* are changing partition but not  */
int *actual_exp_to_part = NULL;          /* processor may not be included.  */
int actual_num_imp = 0;
int actual_imp_allocated = 0;
ZOLTAN_ID_PTR actual_imp_gids = NULL;    /* Arrays containing only objs to  */
ZOLTAN_ID_PTR actual_imp_lids = NULL;    /* actually be imported. Objs that  */
int *actual_imp_procs = NULL;            /* are changing partition but not  */
int *actual_imp_to_part = NULL;          /* processor may not be included.  */

  ZOLTAN_TRACE_ENTER(zz, yo);

  /*
   *  Return if this processor is not in the Zoltan structure's
   *  communicator.
   */

  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz)) {
    goto End;
  }


  /*
   *  Check that all procs use the same id types.
   */

  ierr = check_input(zz, 
                    ((num_export >= 0 && export_to_part) || 
                     (num_import >= 0 && import_to_part)),
                     &include_parts);
  if (ierr != ZOLTAN_OK) 
    goto End;

  num_gid_entries = zz->Num_GID;
  num_lid_entries = zz->Num_LID;

  /*
   *  Check that all necessary query functions are available.
   */

  if (zz->Get_Obj_Size == NULL && zz->Get_Obj_Size_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
           "ZOLTAN_OBJ_SIZE_FN or ZOLTAN_OBJ_SIZE_MULTI_FN function "
           "to use the migration-help tools.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->Pack_Obj == NULL && zz->Pack_Obj_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
           "ZOLTAN_PACK_OBJ_FN or ZOLTAN_PACK_OBJ_MULTI_FN function "
           "to use the migration-help tools.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->Unpack_Obj == NULL && zz->Unpack_Obj_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
         "ZOLTAN_UNPACK_OBJ_FN or ZOLTAN_UNPACK_OBJ_MULTI_FN function "
         "to use the migration-help tools.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }


  if (num_export >= 0) {

    /* Build the actual export arrays */
    ierr = actual_arrays(zz, num_gid_entries, num_lid_entries,
                         num_export, export_global_ids, export_local_ids, 
                         export_procs, export_to_part, 
                         &actual_num_exp, &actual_exp_gids, &actual_exp_lids,
                         &actual_exp_procs, &actual_exp_to_part,
                         &actual_exp_allocated);
    if (ierr < 0) 
      goto End;

    /* Compute communication map based on actual exports.  */

    msgtag = 32767;
    ierr = Zoltan_Comm_Create(&exp_plan, actual_num_exp, actual_exp_procs,
                              zz->Communicator, msgtag, &actual_num_imp);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc,yo,"Error returned from Zoltan_Comm_Create.");
      goto End;
    }
  }

  else if (num_import >= 0) {

    /* Build the actual import arrays */
    ierr = actual_arrays(zz, num_gid_entries, num_lid_entries,
                         num_import, import_global_ids, import_local_ids, 
                         import_procs, import_to_part, 
                         &actual_num_imp, &actual_imp_gids, &actual_imp_lids,
                         &actual_imp_procs, &actual_imp_to_part,
                         &actual_imp_allocated);
    if (ierr < 0) 
      goto End;
    
    /* Compute communication map based on imports.  */
    msgtag = 32767;
    ierr = Zoltan_Comm_Create(&imp_plan, actual_num_imp, actual_imp_procs,
                              zz->Communicator, msgtag, &actual_num_exp);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc,yo,"Error returned from Zoltan_Comm_Create.");
      goto End;
    }

    /* Compute actual export lists for packing objects */
    if (actual_num_exp > 0) {
      actual_exp_allocated = 1;
      actual_exp_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, actual_num_exp);
      actual_exp_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, actual_num_exp);
      actual_exp_procs = (int *) ZOLTAN_MALLOC(sizeof(int) * actual_num_exp);
      if (include_parts)
        actual_exp_to_part = (int *) ZOLTAN_MALLOC(sizeof(int)*actual_num_exp);
      if (actual_exp_gids == NULL ||
          (num_lid_entries && actual_exp_lids == NULL) ||
          actual_exp_procs == NULL || 
          (import_to_part != NULL && actual_exp_to_part == NULL)) {
        Zoltan_Multifree(__FILE__, __LINE__, 4, 
                         &actual_exp_gids, &actual_exp_lids, 
                         &actual_exp_procs, &actual_exp_to_part);
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }

    msgtag2 = 32766;
    ierr = Zoltan_Comm_Do(imp_plan, msgtag2, (char *) actual_imp_gids,
                    (int) (sizeof(ZOLTAN_ID_TYPE)*(num_gid_entries)),
                    (char *) actual_exp_gids);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
      goto End;
    }

    if (num_lid_entries) {
      msgtag2--;
      ierr = Zoltan_Comm_Do(imp_plan, msgtag2, (char *) actual_imp_lids,
                      (int) (sizeof(ZOLTAN_ID_TYPE)*num_lid_entries),
                      (char *) actual_exp_lids);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
        goto End;
      }
    }

    Zoltan_Comm_Info(imp_plan, NULL, NULL, NULL, NULL, NULL, NULL,
                     NULL, NULL, NULL, NULL, NULL, actual_exp_procs, NULL);
 
    if (include_parts) {
      msgtag2--;
      ierr = Zoltan_Comm_Do(imp_plan, msgtag2, (char *) actual_imp_to_part,
                      (int) sizeof(int), (char *) actual_exp_to_part);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
        goto End;
      }
    }
   
    /* Create inverse plan (i.e., plan based on exports) so can set 
     * variable sizes. 
     * (Zoltan_Comm_Do_Reverse(imp_plan, ...) allows sending variable 
     * but does not tell how large to allocate receive buffer.
     */
    ierr = Zoltan_Comm_Invert_Plan(&imp_plan);
    exp_plan = imp_plan;
    imp_plan = NULL;
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Import or export lists needed.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->Migrate.Pre_Migrate_PP != NULL) {
    zz->Migrate.Pre_Migrate_PP(zz->Migrate.Pre_Migrate_PP_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs, import_to_part,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, export_to_part,
                            &ierr);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                      "ZOLTAN_PRE_MIGRATE_PP_FN function.");
      goto End;
    }
  }


  if (zz->Migrate.Pre_Migrate != NULL) {
    zz->Migrate.Pre_Migrate(zz->Migrate.Pre_Migrate_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs,
                            &ierr);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                      "ZOLTAN_PRE_MIGRATE_FN function.");
      goto End;
    }
  }


  ZOLTAN_TRACE_DETAIL(zz, yo, "Done pre-migration processing");

  id_size = Zoltan_Align(num_gid_entries * sizeof(ZOLTAN_ID_TYPE));
  /* Note that alignment is not strictly necessary 
     when ZOLTAN_ID_TYPE is int or unsigned int. */
  aligned_int = Zoltan_Align(sizeof(int));
  tag_size = id_size + aligned_int;

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
  if (actual_num_exp > 0) {
    sizes = (int *) ZOLTAN_MALLOC(actual_num_exp * sizeof(int));
    if (!sizes) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    if (zz->Get_Obj_Size_Multi != NULL) {
      zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data, 
                             num_gid_entries, num_lid_entries, actual_num_exp,
                             actual_exp_gids, actual_exp_lids, sizes, &ierr);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                        "ZOLTAN_OBJ_SIZE_MULTI function.");
        goto End;
      }
    }
    else {
      for (i = 0; i < actual_num_exp; i++){
        lid = (num_lid_entries ? &(actual_exp_lids[i*num_lid_entries]) : NULL);
        sizes[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data, 
                       num_gid_entries, num_lid_entries,
                       &(actual_exp_gids[i*num_gid_entries]), 
                       lid, &ierr);
        if (ierr < 0) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                          "ZOLTAN_OBJ_SIZE function.");
          goto End;
        }
      }
    }

    total_send_size = 0;

    for (i = 0; i < actual_num_exp; i++) {
      sizes[i] = Zoltan_Align(sizes[i]);
      total_send_size += sizes[i] + tag_size;
    }
    export_buf = (char *) ZOLTAN_CALLOC(total_send_size, sizeof(char));
    if (!export_buf) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    if (zz->Pack_Obj_Multi != NULL) {
      /* Allocate an index array for ZOLTAN_PACK_OBJ_MULTI_FN. */
      idx = (int *) ZOLTAN_MALLOC(actual_num_exp * sizeof(int));
      if (!idx) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
        ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }

    /*
     *  Pack the objects for export.
     */
  
    idx_cnt = 0;
    tmp = export_buf;
    for (i = 0; i < actual_num_exp; i++) {

      /* Pack the object's global ID */
      tmp_id = (ZOLTAN_ID_PTR) tmp;
      ZOLTAN_SET_GID(zz, tmp_id, &(actual_exp_gids[i*num_gid_entries]));
      tmp += id_size;
    
      /* Pack the object's size */
      *((int *)tmp) = sizes[i];
      tmp += aligned_int;

      /* If using ZOLTAN_PACK_OBJ_MULTI_FN, build the index array. */
      idx_cnt += tag_size;
      if (idx != NULL) {
        idx[i] = idx_cnt;
      }
      tmp += sizes[i];
      idx_cnt += sizes[i];
    }

    if (zz->Pack_Obj_Multi != NULL) {
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] DEBUG in %s: Packing objects with multi-pack\n", 
               zz->Proc, yo);
      }
      zz->Pack_Obj_Multi(zz->Pack_Obj_Multi_Data,
                         num_gid_entries, num_lid_entries, actual_num_exp,
                         actual_exp_gids, actual_exp_lids, 
                         (actual_exp_to_part!=NULL ? actual_exp_to_part 
                                                   : actual_exp_procs),
                         sizes, idx, export_buf, &ierr);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                        "ZOLTAN_PACK_OBJ_MULTI function.");
        goto End;
      }
if (ierr != 0){
 fprintf(stderr,"(%d) error conditions returned from ZOLTAN_PACK_OBJ_MULTI function.",zz->Proc);
}
    }
    else {
      tmp = export_buf + tag_size;
      for (i = 0; i < actual_num_exp; i++) {
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
          printf("[%1d] DEBUG in %s: Packing object with gid ", zz->Proc, yo);
          ZOLTAN_PRINT_GID(zz, &(actual_exp_gids[i*num_gid_entries]));
          printf("size = %d bytes\n", sizes[i]); 
        }

        /* Pack the object's data */
        lid = (num_lid_entries ? &(actual_exp_lids[i*num_lid_entries]) : NULL);
        dest = (actual_exp_to_part != NULL ? actual_exp_to_part[i] 
                                           : actual_exp_procs[i]);

        zz->Pack_Obj(zz->Pack_Obj_Data, 
                           num_gid_entries, num_lid_entries,
                           &(actual_exp_gids[i*num_gid_entries]), lid, dest,
                           sizes[i], tmp, &ierr);
        if (ierr < 0) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                          "ZOLTAN_PACK_OBJ function.");
          goto End;
        }
if (ierr != 0){
 fprintf(stderr,"(%d) error conditions returned from ZOLTAN_PACK_OBJ_MULTI function.",zz->Proc);
}
        tmp += sizes[i] + tag_size;
      }
    }
    ZOLTAN_FREE(&idx);
    tmp_id = NULL;
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done packing objects");


  /* Modify sizes[] to contain message sizes, not object sizes */
  for (i=0; i<actual_num_exp; i++) {
    sizes[i] += tag_size;
  }

  msgtag--;
  ierr = Zoltan_Comm_Resize(exp_plan, sizes, msgtag, &total_recv_size);
  if (ierr < 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Resize.");
    goto End;
  }

  if (actual_num_imp > 0) {
    import_buf = (char *) ZOLTAN_MALLOC(total_recv_size);
    if (!import_buf) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  /*
   *  Send the export data using the communication plan.
   */

  msgtag2 = 32765;
  ierr = Zoltan_Comm_Do(exp_plan, msgtag2, export_buf, 1, import_buf);
  if (ierr < 0) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Comm_Do.");
    goto End;
  }

  /*
   *  Free whatever memory we can.
   */

  Zoltan_Comm_Destroy(&exp_plan);
  ZOLTAN_FREE(&export_buf);
  ZOLTAN_FREE(&sizes);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done communication");

  /* 
   *  Perform application-specified processing before unpacking the data.
   */
  if (zz->Migrate.Mid_Migrate_PP != NULL) {
    zz->Migrate.Mid_Migrate_PP(zz->Migrate.Mid_Migrate_PP_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs, import_to_part,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, export_to_part,
                            &ierr);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                      "ZOLTAN_MID_MIGRATE_PP_FN function.");
      goto End;
    }
  }

  if (zz->Migrate.Mid_Migrate != NULL) {
    zz->Migrate.Mid_Migrate(zz->Migrate.Mid_Migrate_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs,
                            num_export, export_global_ids,
                            export_local_ids, export_procs,
                            &ierr);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                      "ZOLTAN_MID_MIGRATE_FN function.");
      goto End;
    }
  }

  /*
   *  Unpack the object data.
   */

  if (actual_num_imp > 0) {

    if (zz->Unpack_Obj_Multi != NULL) {

      /* Allocate and fill input arrays for Unpack_Obj_Multi. */
      sizes = (int *) ZOLTAN_MALLOC(actual_num_imp * sizeof(int));
      tmp_id = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC_GID_ARRAY(zz, actual_num_imp);
      idx = (int *) ZOLTAN_MALLOC(actual_num_imp * sizeof(int));
      if (!sizes || !tmp_id || !idx) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error.");
        ierr = ZOLTAN_MEMERR;
        goto End;
      }

      tmp = import_buf;
      idx_cnt = 0;
      for (i = 0; i < actual_num_imp; i++) {

        /* Unpack the object's global ID */
        ZOLTAN_SET_GID(zz, &(tmp_id[i*num_gid_entries]), (ZOLTAN_ID_PTR) tmp);
        tmp += id_size;

        /* Unpack the object's size */
        sizes[i] = *((int *)tmp);
        tmp += aligned_int;

        /* If using ZOLTAN_UNPACK_OBJ_MULTI_FN, build the index array. */
        idx_cnt += tag_size;
        if (idx != NULL) {
          idx[i] = idx_cnt;
        }

        tmp += sizes[i];
        idx_cnt += sizes[i];
      }

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] DEBUG in %s: Unpacking objects with multi-fn\n",
               zz->Proc,yo);
      }
      zz->Unpack_Obj_Multi(zz->Unpack_Obj_Multi_Data, num_gid_entries,
                         actual_num_imp, tmp_id, sizes, idx, import_buf, &ierr);
      ZOLTAN_FREE(&import_buf);
      ZOLTAN_FREE(&sizes);
      ZOLTAN_FREE(&tmp_id);
      ZOLTAN_FREE(&idx);
      if (ierr < 0) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                        "ZOLTAN_UNPACK_OBJ_MULTI_FN.");
        goto End;
      }
    }
    else {
      tmp = import_buf;
      for (i = 0; i < actual_num_imp; i++) {
        tmp_size = *((int *)(tmp + id_size));
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
          printf("[%1d] DEBUG in %s: Unpacking object with gid ", zz->Proc, yo);
          ZOLTAN_PRINT_GID(zz, (ZOLTAN_ID_PTR)tmp);
          printf("size = %d bytes\n", tmp_size);
        }

        /* Unpack the object's data */

        zz->Unpack_Obj(zz->Unpack_Obj_Data, num_gid_entries,
                       (ZOLTAN_ID_PTR) tmp, tmp_size,
                       tmp + tag_size, &ierr);
        if (ierr < 0) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                          "ZOLTAN_UNPACK_OBJ_FN.");
          goto End;
        }
        tmp += (tmp_size + tag_size);
      }
      ZOLTAN_FREE(&import_buf);
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done unpacking objects");

  if (zz->Migrate.Post_Migrate_PP != NULL) {
    zz->Migrate.Post_Migrate_PP(zz->Migrate.Post_Migrate_PP_Data,
                             num_gid_entries, num_lid_entries,
                             num_import, import_global_ids,
                             import_local_ids, import_procs, import_to_part,
                             num_export, export_global_ids,
                             export_local_ids, export_procs, export_to_part,
                             &ierr);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                      "ZOLTAN_POST_MIGRATE_PP_FN function.");
      goto End;
    }
  }

  if (zz->Migrate.Post_Migrate != NULL) {
    zz->Migrate.Post_Migrate(zz->Migrate.Post_Migrate_Data,
                             num_gid_entries, num_lid_entries,
                             num_import, import_global_ids,
                             import_local_ids, import_procs,
                             num_export, export_global_ids,
                             export_local_ids, export_procs,
                             &ierr);
    if (ierr < 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                      "ZOLTAN_POST_MIGRATE_FN function.");
      goto End;
    }
  }

End:

  if (actual_exp_allocated) {
    Zoltan_Multifree(__FILE__, __LINE__, 4, 
                     &actual_exp_gids, &actual_exp_lids, 
                     &actual_exp_procs, &actual_exp_to_part);
  }
  if (actual_imp_allocated) {
    Zoltan_Multifree(__FILE__, __LINE__, 4, 
                     &actual_imp_gids, &actual_imp_lids, 
                     &actual_imp_procs, &actual_imp_to_part);
  }

  if (ierr < 0) {
    if (exp_plan) Zoltan_Comm_Destroy(&exp_plan);
    Zoltan_Multifree(__FILE__, __LINE__, 5,
                     &import_buf, &tmp_id, &sizes, &idx, &export_buf);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static int check_input(
  ZZ *zz,
  int parts,
  int *include_parts
)
{
/* 
 * Routine to ensure that all processors have the same values of 
 * zz->Num_GID and zz->Num_LID.
 * Also, check whether partitions are included on any processors; if so,
 * set include_parts to true.
 * All processors return the same error code.
 */
char *yo = "check_input";
char msg[256];
int loc_tmp[6];
int glob[] = {0, 0, 0, 0, 0, 0};
int ierr = ZOLTAN_OK;

  loc_tmp[0] = zz->Num_GID;
  loc_tmp[1] = zz->Num_LID;
  loc_tmp[2] = parts;
  loc_tmp[3] = -(zz->Num_GID);
  loc_tmp[4] = -(zz->Num_LID);
  loc_tmp[5] = -(parts);

  /* 
   * Check both max and min values of IDs so that all processors can 
   * return the same error code. 
   */

  MPI_Allreduce(loc_tmp, glob, 6,
                MPI_INT, MPI_MIN, zz->Communicator);
  *include_parts = -(glob[5]);

  if ((glob[0] != -(glob[3])) ||
      (glob[1] != -(glob[4])))
    ierr = ZOLTAN_FATAL;

  if (zz->Num_GID != -(glob[3])){
    sprintf(msg, "Inconsistent global id sizes: Num_GID = %d "
      "but global max is %d\n", zz->Num_GID, -(glob[3]));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
  }

  if (zz->Num_LID != -(glob[4])){
    sprintf(msg, "Inconsistent local id sizes: Num_LID = %d "
      "but global max is %d\n", zz->Num_LID, -(glob[4]));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
  }

  return ierr;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_Help_Migrate(
  ZZ *zz,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids,
  int *import_procs,
  int num_export,
  ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids,
  int *export_procs
)
{
/*
 *  Wrapper around Zoltan_Migrate with NULL pointers for partition arrays.
 *  Maintained for backward compatibility.
 *  Arguments are same as for Zoltan_Migrate.
 */

char *yo = "Zoltan_Help_Migrate";
int ierr;

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (zz->LB.PartDist != NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Non-uniform distribution of partitions over processors is specified; "
      "use Zoltan_Migrate\n");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  if (zz->Migrate.Pre_Migrate_PP || zz->Migrate.Mid_Migrate_PP || 
      zz->Migrate.Post_Migrate_PP) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "Partition information not available in Zoltan_Help_Migrate for "
      "ZOLTAN_*_MIGRATE_PP_FNs; use ZOLTAN_*_MIGRATE_FNs instead.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /*
   * Wrapper (for backward compatilibity) around Zoltan_Migrate.
   * Passes NULL for partition assignment arrays.
   */
  ierr = Zoltan_Migrate(zz, num_import, import_global_ids, import_local_ids,
                        import_procs, NULL,
                        num_export, export_global_ids, export_local_ids,
                        export_procs, NULL);

End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/****************************************************************************/
static int actual_arrays(
  ZZ *zz,
  int num_gid_entries, 
  int num_lid_entries,
  int num,
  ZOLTAN_ID_PTR gids,
  ZOLTAN_ID_PTR lids,
  int *procs,
  int *to_part,
  int *actual_num,
  ZOLTAN_ID_PTR *actual_gids,
  ZOLTAN_ID_PTR *actual_lids,
  int **actual_procs,
  int **actual_to_part,
  int *actual_allocated 
)
{
char *yo = "actual_arrays";
int i, j;

  /*
   *  Test whether to pack objects that have changed partition
   *  but not changed processor.  
   *  If packing them, the actual objects == objects passed to this function.
   *  If not packing them, build arrays with them stripped out.
   */

  *actual_allocated = 0;
  if (!(zz->Migrate.Only_Proc_Changes)) {
    /* Pack all objects, even if they are not changing processor. */
    *actual_num = num;
    *actual_gids = gids;
    *actual_lids = lids;
    *actual_procs = procs;
    *actual_to_part = to_part;
  }
  else {  /* zz->Migrate.Only_Proc_Changes */
    /* Pack only objects that are actually changing processor. */
    *actual_num = 0;
    for (i = 0; i < num; i++) 
      if (procs[i] != zz->Proc)
        (*actual_num)++;

    if (*actual_num == num) {
      /*  Number of actual objects == number of objects in input arrays. */
      /*  No stripping needed. */
      *actual_gids = gids;
      *actual_lids = lids;
      *actual_procs = procs;
      *actual_to_part = to_part;
    }
    else if (*actual_num != num && *actual_num > 0) {
      /*  Number of actual_num < num.  Build arrays  */
      /*  containing only actual objects. */
      *actual_allocated = 1;
      *actual_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, *actual_num);
      *actual_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, *actual_num);
      *actual_procs = (int *) ZOLTAN_MALLOC(sizeof(int) * (*actual_num));
      if (to_part != NULL)
        *actual_to_part = (int *) ZOLTAN_MALLOC(sizeof(int)*(*actual_num));
      if (*actual_gids == NULL || (num_lid_entries && *actual_lids == NULL) ||
          *actual_procs == NULL || 
          (to_part != NULL && *actual_to_part == NULL)) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory Error.");
        Zoltan_Multifree(__FILE__, __LINE__, 4, 
                         actual_gids, actual_lids, 
                         actual_procs, actual_to_part);
        return (ZOLTAN_MEMERR);
      }
   
      for (j = 0, i = 0; i < num; i++) {
        if (procs[i] != zz->Proc) {
          ZOLTAN_SET_GID(zz, 
                        *actual_gids + j*num_gid_entries,
                        gids + i*num_gid_entries);
          if (num_lid_entries)
            ZOLTAN_SET_LID(zz, 
                          *actual_lids + j*num_lid_entries,
                          lids + i*num_lid_entries);
          (*actual_procs)[j] = procs[i];
          if (to_part) (*actual_to_part)[j] = to_part[i];
          j++;
        }
      }
    }
  }
  return ZOLTAN_OK;
}
/****************************************************************************/
/****************************************************************************/
/* static functions used by Zoltan_Lo_Mem_Migrate                           */
/****************************************************************************/


static int message_receive(MPI_Comm comm, int src, char *buf, int len)
{
  int ack = 1;

  if (len > 0){
    MPI_Irecv(buf, len, MPI_BYTE, src, dataTag, comm, &_request);
    MPI_Send(&ack, 1, MPI_INT, src, ackTag, comm);
  }

  return ZOLTAN_OK;
}

static int message_send(MPI_Comm comm, int dest, char *buf, int len)
{
  int ack;
  MPI_Status s;

  if (len > 0){
    /* we do a ready send to reduce the amount buffering MPI must do */
    MPI_Recv(&ack, 1, MPI_INT, dest, ackTag, comm, &s);
    MPI_Rsend(buf, len, MPI_BYTE, dest, dataTag, comm);
  }

  return 0;
}
static int message_wait(int len)
{
  MPI_Status s;
  if (len > 0){
    MPI_Wait(&_request, &s);
  }
  return 0;
}

static int get_int_array(ZZ *zz, int len, int **buf, int *oldlen)
{
  char *yo = "get_int_array";
  char msg[256];

  if (len > *oldlen){
    *buf = (int *)ZOLTAN_REALLOC(*buf, len * sizeof(int));
    if (*buf) *oldlen = len;
    else{
      sprintf(msg,"can not allocate buffer size %d ints (%ld)\n", len,len*sizeof(int));
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      return ZOLTAN_MEMERR;
    }
  }

  return ZOLTAN_OK;
}

static int get_char_array(ZZ *zz, int len, char **buf, int *oldlen)
{
  char *yo = "get_char_array";
  char msg[256];

  if (len > *oldlen){
    *buf = (char *)ZOLTAN_REALLOC(*buf, len * sizeof(char));
    if (*buf) *oldlen = len;
    else{
      sprintf(msg,"can not allocate buffer size %d\n",len);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      return ZOLTAN_MEMERR;
    }
  }

  return ZOLTAN_OK;
}

static int get_gid_array(ZZ *zz, int len, ZOLTAN_ID_TYPE **buf, int *oldlen)
{
  char *yo = "get_char_array";
  char msg[256];

  if (len > *oldlen){
    *buf = ZOLTAN_MALLOC_GID_ARRAY(zz, len);
    if (*buf) *oldlen = len;
    else{
      sprintf(msg,"can not allocate buffer size %d GIDs (%ld)\n", len,len*sizeof(ZOLTAN_ID_TYPE));
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      return ZOLTAN_MEMERR;
    }
  }

  return ZOLTAN_OK;
}

static void reset_arrays()
{
  if (_sizes_len){ ZOLTAN_FREE(&_sizes); _sizes_len = 0; }
  if (_idx_len){ ZOLTAN_FREE(&_idx); _idx_len = 0; }
  if (_export_buf_len){ ZOLTAN_FREE(&_export_buf); _export_buf_len = 0; }
  if (_import_buf_len){ ZOLTAN_FREE(&_import_buf); _import_buf_len = 0; }
  if (_gid_list_len){ ZOLTAN_FREE(&_gid_list); _gid_list_len = 0;}
  if (_lid_list_len){ ZOLTAN_FREE(&_lid_list); _lid_list_len = 0;}
  if (_proc_list_len){ ZOLTAN_FREE(&_proc_list); _proc_list_len = 0;}
  if (_part_list_len){ ZOLTAN_FREE(&_part_list); _part_list_len = 0;}
  if (_import_gids_len){ ZOLTAN_FREE(&_import_gids); _import_gids_len = 0;}
  if (_import_lids_len){ ZOLTAN_FREE(&_import_lids); _import_lids_len = 0;}
  if (_import_procs_len){ ZOLTAN_FREE(&_import_procs); _import_procs_len = 0;}
  if (_import_parts_len){ ZOLTAN_FREE(&_import_parts); _import_parts_len = 0;}
}

static int pack_message_for_proc(ZZ *zz, int dest_proc, int num,
             int num_export, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, 
             int *to_proc, int *to_part, 
             int tag_size, int id_size, int aligned_int,
             char **send_buf, int *send_size,
             ZOLTAN_ID_PTR *export_gids, ZOLTAN_ID_PTR *export_lids, int **export_procs,
             int **export_parts) 
{
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

char *yo = "pack_message_for_proc";
int ierr = ZOLTAN_OK;
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int total_send_size = 0;
int i, idx_cnt, dest, n;
char *export_buf=NULL, *tmp=NULL;
ZOLTAN_ID_PTR tmp_id=NULL, send_gids=NULL, send_lids=NULL;
int *idx=NULL, *sizes=NULL, *send_parts=NULL, *send_procs=NULL;

  *send_buf = NULL;
  *send_size = 0;    

  if (export_gids){
    *export_gids = *export_lids = NULL;
    *export_procs = *export_parts = NULL;
  }

  if (num == 0){
    goto End;
  }

  ierr = get_gid_array(zz, num, &_gid_list, &_gid_list_len);
  if (ierr != ZOLTAN_OK){
    goto End;
  }

  send_gids = _gid_list;
  if (export_gids) *export_gids = _gid_list;

  ierr = get_int_array(zz, num, &_proc_list, &_proc_list_len);
  if (ierr != ZOLTAN_OK){
    goto End;
  }

  send_procs = _proc_list;
  if (export_procs) *export_procs = _proc_list;

  if (lids){
    ierr = get_gid_array(zz, num, &_lid_list, &_lid_list_len);
    if (ierr != ZOLTAN_OK){
      goto End;
    }

    send_lids = _lid_list;
    if (export_lids) *export_lids = _lid_list;
  }

  if (to_part){
    ierr = get_int_array(zz, num, &_part_list, &_part_list_len);
    if (ierr != ZOLTAN_OK){
      goto End;
    }

    send_parts = _part_list;
    if (export_parts) *export_parts = _part_list;
  }

  for (i=0, n=0; i < num_export; i++){

    if (to_proc[i] == dest_proc){

      ZOLTAN_COPY_GID_ARRAY(send_gids+n, gids+i, zz, 1);
      send_procs[n] = dest_proc;
      if (to_part){
        send_parts[n] = to_part[i];
      }
      if (lids){
        send_lids[n] = lids[i];
      }

      n++;

      if (n == num) break;
    }
  }

  if (num_export > 0) {
    ierr = get_int_array(zz, num, &_sizes, &_sizes_len);
    if (ierr != ZOLTAN_OK){
      goto End;
    }
    sizes = _sizes;

    if (zz->Get_Obj_Size_Multi != NULL) {
      zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data, 
                             num_gid_entries, num_lid_entries, num,
                             send_gids, (send_lids ? send_lids : NULL), 
                             sizes, &ierr);
    }
    else {
      for (i = 0; i < num; i++){
        sizes[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data, 
                       num_gid_entries, num_lid_entries,
                       &(send_gids[i*num_gid_entries]), 
                       (send_lids ? &(send_lids[i*num_lid_entries]) : NULL),
                       &ierr);
        if (ierr != 0) break;
      }
    }

    if (ierr != 0) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from ZOLTAN_OBJ_SIZE function.");
      goto End;
    }

    for (i = 0; i < num; i++) {
      sizes[i] = Zoltan_Align(sizes[i]);
      total_send_size += sizes[i] + tag_size;
    }

    ierr = get_char_array(zz, total_send_size, &_export_buf, &_export_buf_len);

    if (ierr != ZOLTAN_OK){
      goto End;
    }

    export_buf = _export_buf;

    if (zz->Pack_Obj_Multi != NULL) {
      /* Allocate an index array for ZOLTAN_PACK_OBJ_MULTI_FN. */
      ierr = get_int_array(zz, num, &_idx, &_idx_len);

      if (ierr != ZOLTAN_OK) {
        goto End;
      }

      idx = _idx;
    }

    /*
     *  Pack the objects for export.
     */
  
    idx_cnt = 0;
    tmp = export_buf;
    for (i = 0; i < num; i++) {

      /* Pack the object's global ID */
      tmp_id = (ZOLTAN_ID_PTR) tmp;
      ZOLTAN_SET_GID(zz, tmp_id, &(send_gids[i*num_gid_entries]));
      tmp += id_size;
    
      /* Pack the object's size */
      *((int *)tmp) = sizes[i];
      tmp += aligned_int;

      /* If using ZOLTAN_PACK_OBJ_MULTI_FN, build the index array. */
      idx_cnt += tag_size;
      if (idx != NULL) {
        idx[i] = idx_cnt;
      }
      tmp += sizes[i];
      idx_cnt += sizes[i];
    }

    if (zz->Pack_Obj_Multi != NULL) {
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] DEBUG in %s: Packing objects with multi-pack\n", 
               zz->Proc, yo);
      }
      zz->Pack_Obj_Multi(zz->Pack_Obj_Multi_Data,
                         num_gid_entries, num_lid_entries, num,
                         send_gids, (send_lids ? send_lids : NULL), 
                         (send_parts ? send_parts : send_procs),
                         sizes, idx, export_buf, &ierr);
    }
    else {
      tmp = export_buf + tag_size;
      for (i = 0; i < num; i++) {
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
          printf("[%1d] DEBUG in %s: Packing object with gid ", zz->Proc, yo);
          ZOLTAN_PRINT_GID(zz, &(send_gids[i*num_gid_entries]));
          printf("size = %d bytes\n", sizes[i]); 
        }

        /* Pack the object's data */
        dest = (send_parts != NULL ? send_parts[i] : send_procs[i]);

        zz->Pack_Obj(zz->Pack_Obj_Data, 
                           num_gid_entries, num_lid_entries,
                           &(send_gids[i*num_gid_entries]),
                           (send_lids ? &(send_lids[i*num_lid_entries]) : NULL), 
                           dest, sizes[i], tmp, &ierr);

        if (ierr < 0)
          break;

        tmp += sizes[i] + tag_size;
      }
    }

    if (ierr != 0){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from ZOLTAN_PACK_OBJ function.");
      goto End;
    }

    ZOLTAN_FREE(&idx);
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done packing objects");

  *send_buf = export_buf;
  *send_size = total_send_size;    

End:
  return ierr;
}

static int create_import_lists(ZZ *zz, int msglen, int src, int num,
                  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int *procs, int *parts,
                  int *objectCount,
                  ZOLTAN_ID_PTR *inGids, ZOLTAN_ID_PTR *inLids, int **inProcs, int **inParts)
{
  /*
   * List objects imported from process "src".
   */

  int ierr = ZOLTAN_OK;
  int numObj = 0;
  int i, n;
  int gid_size = zz->Num_GID;

  *objectCount = 0;
  *inGids = NULL;
  *inLids = NULL;
  if (inProcs) *inProcs = NULL;
  if (inParts) *inParts = NULL;

  if ((msglen < 1) || (num < 1)) goto End;

  for (i=0; i < num; i++){
    if (procs[i] == src){
      numObj++;
    }
  } 

  ierr = get_gid_array(zz, numObj, &_import_gids, &_import_gids_len);
  if (ierr != ZOLTAN_OK){
    goto End;
  }
  *inGids = _import_gids;

  ierr = get_int_array(zz, numObj, &_import_procs, &_import_procs_len);
  if (ierr != ZOLTAN_OK){
    goto End;
  }
  *inProcs = _import_procs;

  if (lids && inLids){
    ierr = get_gid_array(zz, numObj, &_import_lids, &_import_lids_len);
    if (ierr != ZOLTAN_OK){
      goto End;
    }
    *inLids = _import_lids;
  }

  if (parts && inParts){
    ierr = get_int_array(zz, numObj, &_import_parts, &_import_parts_len);
    if (ierr != ZOLTAN_OK){
      goto End;
    }
    *inParts = _import_parts;
  }

  for (i=0, n=0; i < num; i++){
    if (procs[i] == src){
      ZOLTAN_SET_GID(zz, _import_gids + (n * gid_size), gids + (i * gid_size));
      _import_procs[n] = procs[i];
      if (lids && inLids){
        _import_lids[n] = lids[i];
      }
      if (parts && inParts){
        _import_parts[n] = parts[i];
      }
      n++;
    }
  }

  *objectCount = numObj;

End:

  return ierr;
}

/****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

