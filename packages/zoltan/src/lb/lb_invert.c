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
#include "all_allo_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines implementing the list inverting tools
 *  Zoltan_Invert_Lists and Zoltan_Compute_Destinations.
 *  These functions are all callable by the application.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int check_invert_input(ZZ *, int, int *, int *, int *, int *, int *);


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Invert_Lists(
  ZZ *zz,                       /* Zoltan structure.                  */
  int num_in,                   /* Number of objects in the input lists. */ 
  ZOLTAN_ID_PTR in_global_ids,  /* Array of input global IDs. */
  ZOLTAN_ID_PTR in_local_ids,   /* Array of input local IDs. */
  int *in_procs,                /* Array of processor IDs of processors owning
                                   the input objects. */
  int *in_to_part,              /* Optional:  Array of partition numbers to 
                                   which input objects should be assigned. */
  int *num_out,                 /* Returned value:  Number of output objs. */
  ZOLTAN_ID_PTR *out_global_ids,/* Returned value:  Array of global IDs of
                                   output objects. */
  ZOLTAN_ID_PTR *out_local_ids, /* Returned value:  Array of local IDs of
                                   output objects. */
  int **out_procs,              /* Returned value:  Array of processor IDs
                                   to which output objects are assigned. */
  int **out_to_part             /* Optional:  Returned value:  Array of 
                                   partition numbers to which output
                                   objects should be assigned. */
)
{
/*
 *  Routine to compute the inverse map.  Can be used in two ways:
 *  1.  Given, for each processor, a list of objects to be received by the 
 *  processor, compute the list of objects that the processor needs to send 
 *  to other processors to satisfy their needs.
 *  2.  Given, for each processor, a list of objects to be sent to other 
 *  processors, compute the list of objects that the processor needs to receive
 *  to satisfy its needs.
 */

char *yo = "Zoltan_Invert_Lists";
char msg[256];
ZOLTAN_COMM_OBJ *comm_plan;        /* Object returned communication routines  */
int msgtag, msgtag2;               /* Message tags for communication routines */
int num_gid_entries = zz->Num_GID;
int num_lid_entries = zz->Num_LID;
int include_parts;                 /* Flag indicating whether to compute
                                      inverse list for partitions. */
int ierr, ret_ierr = ZOLTAN_OK;

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

  ierr = check_invert_input(zz, num_in, in_procs, in_to_part, 
                            &num_gid_entries, &num_lid_entries, &include_parts);
  if (ierr != ZOLTAN_OK) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ierr;
  }

  /* Initialize returned arrays. */
  *out_global_ids = NULL;
  *out_local_ids = NULL;
  *out_procs = NULL;
  if (include_parts) *out_to_part = NULL;


  /*
   *  Compute communication map and num_out, the number of objs this
   *  processor has to out to establish the new decomposition.
   */

  msgtag = 32767;
  ierr = Zoltan_Comm_Create(&comm_plan, num_in, in_procs, zz->Communicator, 
                        msgtag, num_out);
  if (ierr != ZOLTAN_OK) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Create.",
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ret_ierr = ierr;
    goto End;
  }
  

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done comm create");

  /*
   *  Allocate space for the object tags that need to be outed.  Communicate
   *  to get the list of objects to be outed.
   */

  if (*num_out > 0) {
    if (!Zoltan_Special_Malloc(zz,(void **)out_global_ids,*num_out,
                           ZOLTAN_SPECIAL_MALLOC_GID)) {
      ret_ierr = ZOLTAN_MEMERR;
      goto End;
    }
    if (num_lid_entries) {
      if (!Zoltan_Special_Malloc(zz,(void **)out_local_ids,*num_out,
                             ZOLTAN_SPECIAL_MALLOC_LID)) {
        ret_ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }
    if (!Zoltan_Special_Malloc(zz,(void **)out_procs,*num_out,
                           ZOLTAN_SPECIAL_MALLOC_INT)) {
      ret_ierr = ZOLTAN_MEMERR;
      goto End;
    }
    if (include_parts) {
      if (!Zoltan_Special_Malloc(zz,(void **)out_to_part,*num_out,
                             ZOLTAN_SPECIAL_MALLOC_INT)) {
        ret_ierr = ZOLTAN_MEMERR;
        goto End;
      }
    }
  }

  /*
   *  Use the communication plan to send global IDs, local IDs, and processor
   *  numbers.  Do in separate communications to avoid a memory copy and to
   *  simplify implementation when a data type is added to the comm. package
   *  (to support heterogeneous computing).
   */

  msgtag2 = 32766;
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) in_global_ids, 
                    (int) (sizeof(ZOLTAN_ID_TYPE)*(num_gid_entries)), 
                    (char *) *out_global_ids);
  if (ierr != ZOLTAN_OK) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Do.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ret_ierr = ierr;
  }

  if (num_lid_entries) {
    msgtag2--;
    ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) in_local_ids, 
                      (int) (sizeof(ZOLTAN_ID_TYPE)*num_lid_entries), 
                      (char *) *out_local_ids);
    if (ierr != ZOLTAN_OK) {
      sprintf(msg, "Error %s returned from Zoltan_Comm_Do.", 
              (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      ret_ierr = ierr;
    }
  }

  Zoltan_Comm_Info(comm_plan, NULL, NULL, NULL, NULL, NULL, NULL,
                   NULL, NULL, NULL, NULL, NULL, *out_procs, NULL);
  
  if (include_parts) {
    msgtag2--;
    ierr = Zoltan_Comm_Do(comm_plan, msgtag2, (char *) in_to_part, 
                      (int) sizeof(int), (char *) *out_to_part);
    if (ierr != ZOLTAN_OK) {
      sprintf(msg, "Error %s returned from Zoltan_Comm_Do.", 
              (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      ret_ierr = ierr;
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done comm_do");

End:

  Zoltan_Comm_Destroy(&comm_plan);

  if (ret_ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    Zoltan_Special_Free(zz,(void**)out_global_ids,ZOLTAN_SPECIAL_MALLOC_GID);
    Zoltan_Special_Free(zz,(void**)out_local_ids,ZOLTAN_SPECIAL_MALLOC_LID);
    Zoltan_Special_Free(zz,(void**)out_procs,ZOLTAN_SPECIAL_MALLOC_INT);
    if (include_parts)
      Zoltan_Special_Free(zz,(void**)out_to_part,ZOLTAN_SPECIAL_MALLOC_INT);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ret_ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int check_invert_input(
  ZZ *zz, 
  int in_num,
  int *in_procs, 
  int *in_to_part, 
  int *num_gid_entries,
  int *num_lid_entries,
  int *include_parts
)
{
/*
 * Routine to ensure that all processors have the same values of
 * zz->Num_GID and zz->Num_LID, and to compute include_parts flag globally.
 * All processors return the same error code.
 */
char *yo = "check_invert_input";
char msg[256];
int loc_tmp[5];
int glob_min[5] = {0,0,0,0,0};
int glob_max[5] = {0,0,0,0,0};
int do_not_include_parts;
int ierr = ZOLTAN_OK;

  loc_tmp[0] = zz->Num_GID;
  loc_tmp[1] = zz->Num_LID;

  /*
   * Check both max and min values of IDs so that all processors can
   * return the same error code.
   */

  MPI_Allreduce(loc_tmp, glob_min, 2,
                MPI_INT, MPI_MIN, zz->Communicator);

  /* 
   * For MPI_MAX operation:
   * include_parts == if any proc has in_num > 0 and specifies in_to_part.
   * do_not_include_parts == if any proc has in_num > 0 and does not specify
   *                         in_to_part.
   * Error:  if include_parts && do_not_include_parts
   * Error:  if any processor has in_num > 0 && in_procs == NULL.
   */
  loc_tmp[2] = (in_num > 0 && in_procs == NULL);
  loc_tmp[3] = (in_num > 0 && in_to_part != NULL);
  loc_tmp[4] = (in_num > 0 && in_to_part == NULL);
  MPI_Allreduce(loc_tmp, glob_max, 5,
                MPI_INT, MPI_MAX, zz->Communicator);

  if (glob_min[0] == glob_max[0])
    *num_gid_entries = glob_max[0];
  else {
    ierr = ZOLTAN_FATAL;
    if (zz->Num_GID != glob_max[0]){
      sprintf(msg, "Inconsistent global id sizes: Num_GID = %d "
        "but global max is %d\n", zz->Num_GID, glob_max[0]);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    }
  }

  if (glob_min[1] == glob_max[1])
    *num_lid_entries = glob_max[1];
  else {
    ierr = ZOLTAN_FATAL;
    if (zz->Num_LID != glob_max[1]){
      sprintf(msg, "Inconsistent local id sizes: Num_LID = %d "
        "but global max is %d\n", zz->Num_LID, glob_max[1]);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    }
  }

  if (glob_max[2]) {
    ierr = ZOLTAN_FATAL;
    if (loc_tmp[2] == glob_max[2]) {
      sprintf(msg, 
              "Inconsistent input: # objects = %d but proc array is NULL\n", 
              in_num);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    }
  }

  *include_parts = glob_max[3];
  do_not_include_parts = glob_max[4];
  if (*include_parts && do_not_include_parts) {
    ierr = ZOLTAN_FATAL;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Inconsistent input; some processors "
      "include partition arrays while others do not.");
  }

  return ierr;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Compute_Destinations(
  ZZ *zz,
  int num_in,
  ZOLTAN_ID_PTR in_global_ids,
  ZOLTAN_ID_PTR in_local_ids,
  int *in_procs,
  int *num_out,
  ZOLTAN_ID_PTR *out_global_ids,
  ZOLTAN_ID_PTR *out_local_ids,
  int **out_procs
)
{
/*  
 *  Wrapper around Zoltan_Invert_Lists, with NULL for excluded partition arrays.
 *  Maintained for backward compatibility.
 *  Arguments are analogous to Zoltan_Invert_Lists.
 */

char *yo = "Zoltan_Compute_Destinations";
int ierr;


  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = Zoltan_Invert_Lists(zz, num_in, 
           in_global_ids, in_local_ids, in_procs, NULL,
           num_out,
           out_global_ids, out_local_ids, out_procs, NULL);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

