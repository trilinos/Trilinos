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
#include "order_const.h"
#include "key_params.h"
#include "params_const.h"
#include "timer_const.h"
#include "ha_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for performing ordering with Zoltan.
 *  These functions are all callable by the application.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/**********  parameters structure for ordering **********/
static PARAM_VARS Order_params[] = {
        { "ORDER_METHOD", NULL, "STRING" },
        { "ORDER_TYPE", NULL, "STRING" },
        { "ORDER_START_INDEX", NULL, "INT" },
        { "REORDER", NULL, "BOOLEAN" },
        { "USE_ORDER_INFO", NULL, "BOOLEAN" },
        { NULL, NULL, NULL } };

int Zoltan_Order(
  ZZ *zz,               /* Zoltan structure */
  int *num_gid_entries, /* # of entries for a global id */
  int *num_lid_entries, /* # of entries for a local id */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
                        /* The application must allocate enough space */
  ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
                        /* The application must allocate enough space */
  int *rank,            /* rank[i] is the rank of gids[i] */
  int *iperm,           /* inverse permutation of rank */
  ZOS *order_info	/* Method-specific ordering info */
)
{
/*
 * Main user-call for ordering.
 * Input:  
 *   zz, a Zoltan structure with appropriate function pointers set.
 *   gids, a list of global ids or enough space to store such a list
 *   lids, a list of local ids or enough space to store such a list
 * Output: 
 *   num_gid_entries
 *   num_lid_entries
 *   gids, a list of global ids (filled in if empty on entry)
 *   lids, a list of local ids (filled in if empty on entry)
 *   rank, rank[i] is the global rank of gids[i]
 *   iperm, inverse permutation of rank
 *   order_info, a Zoltan Ordering Struct with additional info.
 * Return values:
 *   Zoltan error code.
 */

  char *yo = "Zoltan_Order";
  int ierr;
  int *vtxdist;
  double start_time, end_time;
  double order_time[2] = {0.0,0.0};
  char msg[256];
  int comm[2],gcomm[2]; 
  ZOLTAN_ORDER_FN *Order_fn;
  struct Zoltan_Order_Options opt;


  ZOLTAN_TRACE_ENTER(zz, yo);

  if (zz->Proc == zz->Debug_Proc && zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS)
    Zoltan_Print_Key_Params(zz);

  start_time = Zoltan_Time(zz->Timer);

  /* 
   * Compute Max number of array entries per ID over all processors.
   * This is a sanity-maintaining step; we don't want different
   * processors to have different values for these numbers.
   */
  comm[0] = zz->Num_GID;
  comm[1] = zz->Num_LID;
  MPI_Allreduce(comm, gcomm, 2, MPI_INT, MPI_MAX, zz->Communicator);
  zz->Num_GID = *num_gid_entries = gcomm[0];
  zz->Num_LID = *num_lid_entries = gcomm[1];

  /*
   *  Return if this processor is not in the Zoltan structure's
   *  communicator.
   */

  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_OK);
  }

  /*
   *  Get ordering options from parameter list.
   */

  /* Set default parameter values */
  strcpy(opt.method, "PARMETIS");
  strcpy(opt.order_type, "GLOBAL");
  opt.use_order_info = 0;
  opt.start_index = 0;
  opt.reorder = 0;

  Zoltan_Bind_Param(Order_params, "ORDER_METHOD", (void *) opt.method);
  Zoltan_Bind_Param(Order_params, "ORDER_TYPE",   (void *) opt.order_type);
  Zoltan_Bind_Param(Order_params, "ORDER_START_INDEX", (void *) &opt.start_index);
  Zoltan_Bind_Param(Order_params, "REORDER",      (void *) &opt.reorder);
  Zoltan_Bind_Param(Order_params, "USE_ORDER_INFO", (void *) &opt.use_order_info);

  Zoltan_Assign_Param_Vals(zz->Params, Order_params, zz->Debug_Level, 
                           zz->Proc, zz->Debug_Proc);

  if (opt.use_order_info == 0) order_info = NULL;

  /*
   *  Check that the user has allocated space for the return args. 
   */
  if (!(gids && lids && rank && iperm)){
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input argument is NULL. Please allocate all required arrays before calling this routine.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  /*
   *  Find the selected method.
   */

  if (!strcmp(opt.method, "NONE")) {
    if (zz->Proc == zz->Debug_Proc && zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS)
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Ordering method selected == NONE; no ordering performed\n");

    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_WARN);
  }
  else if ((!strcmp(opt.method, "PARMETIS")) 
            || (!strcmp(opt.method, "NODEND"))) {
    Order_fn = Zoltan_ParMetis_Order;
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unknown ordering method");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  /*
   *  Construct the heterogenous machine description.
   */

  ierr = Zoltan_Build_Machine_Desc(zz);

  if (ierr == ZOLTAN_FATAL){
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done machine description");

  /*
   * Call the actual ordering function.
   */

  ierr = (*Order_fn)(zz, gids, lids, rank, iperm, &opt, order_info);

  if (ierr) {
    sprintf(msg, "Ordering routine returned error code %d.", ierr);
    if (ierr == ZOLTAN_WARN){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    } else {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done ordering");

  /* Compute inverse permutation if necessary */
  ierr = Zoltan_Get_Distribution(zz, &vtxdist);
  if (ierr){
    /* Error */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Get_Distribution.\n");
    return (ierr);
  }

  if (!(opt.return_args & RETURN_RANK)){
    /* Compute rank from iperm */
    ZOLTAN_TRACE_DETAIL(zz, yo, "Inverting permutation");
    Zoltan_Inverse_Perm(zz, iperm, rank, vtxdist, opt.order_type, opt.start_index);
  }
  else if (!(opt.return_args & RETURN_IPERM)){
    /* Compute iperm from rank */
    ZOLTAN_TRACE_DETAIL(zz, yo, "Inverting permutation");
    Zoltan_Inverse_Perm(zz, rank, iperm, vtxdist, opt.order_type, opt.start_index);
  }
  ZOLTAN_FREE(&vtxdist);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done ordering");

  end_time = Zoltan_Time(zz->Timer);
  order_time[0] = end_time - start_time;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST) {
    int i, nobjs;
    nobjs = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &i);
    Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
    printf("ZOLTAN: rank for ordering on Proc %d\n", zz->Proc);
    for (i = 0; i < nobjs; i++) {
      printf("GID = ");
      ZOLTAN_PRINT_GID(zz, &(gids[i*(*num_gid_entries)]));
      printf(", rank = %3d\n", rank[i]);
    }
    printf("\n");
    printf("ZOLTAN: inverse permutation on Proc %d\n", zz->Proc);
    for (i = 0; i < nobjs; i++) {
      printf("iperm[%3d] = %3d\n", i, iperm[i]);
    }
    printf("\n");
    Zoltan_Print_Sync_End(zz->Communicator, TRUE);
  }

  
  /* Print timing info */
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ZTIME) {
    if (zz->Proc == zz->Debug_Proc) {
      printf("ZOLTAN Times:  \n");
    }
    Zoltan_Print_Stats (zz->Communicator, zz->Debug_Proc, order_time[0], 
                   "ZOLTAN     Balance:     ");
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  if (ierr)
    return (ierr);
  else
    return (ZOLTAN_OK);
}


/****************** Parameter routine *************************/
int Zoltan_Order_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */

    status = Zoltan_Check_Param(name, val, Order_params, &result, &index);

    return(status);
}

