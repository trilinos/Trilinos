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
        { "USE_ORDER_INFO", NULL, "BOOLEAN" },
        { NULL, NULL, NULL } };

int Zoltan_Order(
  ZZ *zz,               /* Zoltan structure */
  int *num_gid_entries, /* # of entries for a global id */
  int *num_lid_entries, /* # of entries for a local id */
  ZOLTAN_ID_PTR gids,   /* Ordered list of global ids (local to this proc) */
                        /* The application must allocate enough space */
  ZOLTAN_ID_PTR lids,   /* Ordered list of local ids (local to this proc) */
                        /* The application must allocate enough space */
  ZOS *order_info	/* Method-specific ordering info */
)
{
/*
 * Main user-call for ordering.
 * Input:  
 *   zz, a Zoltan structure with appropriate function pointers set.
 * Output: 
 *   num_gid_entries
 *   num_lid_entries
 *   gids, the ordered list of global ids
 *   lids, the ordered list of local ids
 *   order_info, a Zoltan Ordering Struct with additional info.
 * Return values:
 *   Zoltan error code.
 */

  char *yo = "Zoltan_Order";
  int error, use_order_info; 
  double start_time, end_time;
  double order_time[2] = {0.0,0.0};
  char msg[256], method[256];
  int comm[2],gcomm[2]; 
  ZOLTAN_ORDER_FN *Order_fn;


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
   *  Get ordering method from parameter list.
   */

  /* Set default parameter values */
  use_order_info = 0;
  strcpy(method, "PARMETIS");

  Zoltan_Bind_Param(Order_params, "ORDER_METHOD", (void *) &method);
  Zoltan_Bind_Param(Order_params, "USE_ORDER_INFO", (void *) &use_order_info);

  Zoltan_Assign_Param_Vals(zz->Params, Order_params, zz->Debug_Level, 
                           zz->Proc, zz->Debug_Proc);

  if (use_order_info == 0) order_info = NULL;

  /*
   *  Find the selected method.
   */

  if (!strcmp(method, "NONE")) {
    if (zz->Proc == zz->Debug_Proc && zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS)
      printf("%s Ordering method selected == NONE; no ordering performed\n",
              yo);

    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_WARN);
  }
  else if (!strcmp(method, "PARMETIS")) {
    Order_fn = Zoltan_ParMetis_Order;
    /* Set PARMETIS_METHOD to NODEND; needed by Zoltan_ParMetis_Jostle */
    Zoltan_Set_Param(zz, "PARMETIS_METHOD", "NODEND");
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unknown ordering method");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
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
   * Call the actual ordering function.
   */

  error = (*Order_fn)(zz, gids, lids, order_info);

  if (error == ZOLTAN_FATAL){
    sprintf(msg, "Ordering routine returned error code %d.", error);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (error);
  }
  else if (error){
    sprintf(msg, "Ordering routine returned error code %d.", error);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done ordering");

  end_time = Zoltan_Time(zz->Timer);
  order_time[0] = end_time - start_time;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST) {
    int i, nobjs;
    nobjs = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &i);
    Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
    printf("ZOLTAN: GID ordering on Proc %d\n", zz->Proc);
    for (i = 0; i < nobjs; i++) {
      ZOLTAN_PRINT_GID(zz, &gids[i*(zz->Num_GID)]);
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
  if (error)
    return (error);
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

