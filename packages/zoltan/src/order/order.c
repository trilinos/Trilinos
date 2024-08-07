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
#include "key_params.h"
#include "params_const.h"
#include "ha_const.h"
#include "order_const.h"
#include "hsfcOrder.h"
#include "order_params.h"

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

/* New Interface */

int Zoltan_Order (
      struct Zoltan_Struct *zz,
      int num_gid_entries,
      int num_obj,
      ZOLTAN_ID_PTR gids,
      ZOLTAN_ID_PTR permuted_global_ids
    )
{
/*
 * Main user-call for ordering.
 * Input:
 *   zz, a Zoltan structure with appropriate function pointers set.
 *   gids, a list of global ids.
 *   num_gid_entries
 * Output:
 *   permuted_global_ids
 * Return values:
 *   Zoltan error code.
 */

  char *yo = "Zoltan_Order";
  int ierr;
  double start_time, end_time;
  double order_time[2] = {0.0,0.0};
  int comm[2],gcomm[2];
  ZOLTAN_ORDER_FN *Order_fn;
  struct Zoltan_Order_Options opt;
  ZOLTAN_ID_PTR local_gids=NULL, lids=NULL;
  int local_num_obj;
  ZOLTAN_ID_TYPE *local_rank = NULL;  /* global IDs */
  struct Zoltan_DD_Struct *dd = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (zz->Proc == zz->Debug_Proc && zz->Debug_Level >= ZOLTAN_DEBUG_PARAMS) {
    Zoltan_Print_Configuration("  ");
    Zoltan_Print_Key_Params(zz);
  }

  start_time = Zoltan_Time(zz->Timer);

  /*
   * Compute Max number of array entries per ID over all processors.
   * This is a sanity-maintaining step; we don't want different
   * processors to have different values for these numbers.
   */
  comm[0] = zz->Num_GID;
  comm[1] = zz->Num_LID;
  MPI_Allreduce(comm, gcomm, 2, MPI_INT, MPI_MAX, zz->Communicator);
  zz->Num_GID = gcomm[0];

  if (num_gid_entries != zz->Num_GID) {
    char msg[253];
    sprintf(msg, "num_gid_entries=%d is not equal to parameter setting "
                 "NUM_GID_ENTRIES=%d\n", num_gid_entries, zz->Num_GID);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    return (ZOLTAN_FATAL);
  }


  zz->Order.nbr_objects = num_obj;
  zz->Order.start = NULL;
  zz->Order.ancestor = NULL;
  zz->Order.leaves = NULL;
  zz->Order.nbr_leaves = 0;
  zz->Order.nbr_blocks = 0;

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
#ifdef HAVE_MPI
  strncpy(opt.method, "PARMETIS", MAX_PARAM_STRING_LEN);
  strcpy(zz->Order.order_type, "GLOBAL");
#else
  strncpy(opt.method, "METIS", MAX_PARAM_STRING_LEN);
  strcpy(zz->Order.order_type, "LOCAL");
#endif /* HAVE_MPI */

  opt.use_order_info = 0;
  opt.start_index = 0;

  Zoltan_Bind_Param(Order_params, "ORDER_METHOD", (void *) opt.method);
  Zoltan_Bind_Param(Order_params, "USE_ORDER_INFO", (void *) &opt.use_order_info);

  Zoltan_Assign_Param_Vals(zz->Params, Order_params, zz->Debug_Level,
                           zz->Proc, zz->Debug_Proc);

  /*
   *  Check that the user has allocated space for the return args.
   */
  if (num_obj && !(gids && permuted_global_ids)){
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
  else if (!strcmp(opt.method, "LOCAL_HSFC"))
  {
    Order_fn = Zoltan_LocalHSFC_Order;
    strcpy(zz->Order.order_type, "LOCAL"); /*MMW, not sure about this*/
  }
#ifdef ZOLTAN_PARMETIS
  else if (!strcmp(opt.method, "METIS")) {
    Order_fn = Zoltan_ParMetis_Order;
    strcpy(zz->Order.order_type, "LOCAL");
  }
  else if (!strcmp(opt.method, "PARMETIS")) {
    Order_fn = Zoltan_ParMetis_Order;
    strcpy(zz->Order.order_type, "GLOBAL");
  }
#endif /* ZOLTAN_PARMETIS */
#ifdef ZOLTAN_SCOTCH
  else if (!strcmp(opt.method, "SCOTCH")) {
    Order_fn = Zoltan_Scotch_Order;
    strcpy(zz->Order.order_type, "LOCAL");
  }
  else if (!strcmp(opt.method, "PTSCOTCH")) {
    Order_fn = Zoltan_Scotch_Order;
    strcpy(zz->Order.order_type, "GLOBAL");
  }
#endif /* ZOLTAN_SCOTCH */
#ifdef ZOLTAN_HUND
  else if (!strcasecmp(opt.method, "HUND")) {
    ierr = Zoltan_HUND(zz, num_gid_entries, num_obj, gids, permuted_global_ids);
    goto End;
  }
#endif /* ZOLTAN_HUND */
  else {
    fprintf(stderr, "%s\n", opt.method);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unknown ordering method");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  /************************************
   *  Check for required query function
   ************************************/
  if (zz->Get_Num_Obj != NULL) {
    local_num_obj = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_Obj.");
      return (ierr);
    }
  }
  else {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register ZOLTAN_NUM_OBJ_FN.");
    return (ZOLTAN_FATAL);
  }


  /* TODO allocate all this stuff with the graph */
  local_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, local_num_obj);
  local_rank = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC(local_num_obj*sizeof(ZOLTAN_ID_TYPE));

  lids = ZOLTAN_MALLOC_LID_ARRAY(zz, local_num_obj);

  /*
   * Call the actual ordering function.
   * Compute gid according to the local graph.
   */

  ierr = (*Order_fn)(zz, local_num_obj, local_gids, lids, local_rank, NULL, &opt);
  ZOLTAN_FREE(&lids);

  if (ierr) {
    char msg[253];
    sprintf(msg, "Ordering routine returned error code %d.", ierr);
    if (ierr == ZOLTAN_WARN){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    } else {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      Zoltan_Multifree(__FILE__, __LINE__, 2,
                       &local_gids, &local_rank);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done ordering");

  /* TODO: Use directly the "graph" structure to avoid to duplicate things. */

  /* TODO: At this time, I consider rank == permuted_global_ids */

  /* I store : GNO, rank, permuted GID */
  /* MMW: perhaps don't ever use graph here since we need to support geometric orderings, otherwise need if/else */
  ierr = Zoltan_DD_Create (&dd, zz->Communicator, zz->Num_GID, (local_rank==NULL)?0:1, 0, local_num_obj, 0);
  /* Hope a linear assignment will help a little */
  if (local_num_obj)
    Zoltan_DD_Set_Neighbor_Hash_Fn1(dd, local_num_obj);
  /* Associate all the data with our xGNO */

  Zoltan_DD_Update (dd, local_gids, local_rank, NULL, NULL, local_num_obj);

  ZOLTAN_FREE(&local_gids);
  ZOLTAN_FREE(&local_rank);

  Zoltan_DD_Find (dd, gids, permuted_global_ids, NULL, NULL, num_obj, NULL);
  Zoltan_DD_Destroy(&dd);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done Registering results");


  end_time = Zoltan_Time(zz->Timer);
  order_time[0] = end_time - start_time;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST) {
    int i;
    Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
    printf("ZOLTAN: rank for ordering on Proc %d\n", zz->Proc);
    for (i = 0; i < num_obj; i++) {
      printf("GID = ");
      ZOLTAN_PRINT_GID(zz, &(gids[i*(num_gid_entries)]));
      printf(", rank = " ZOLTAN_ID_SPEC "\n", permuted_global_ids[i]);
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

#ifdef ZOLTAN_HUND
 End:
#endif /*ZOLTAN_HUND*/
  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
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


/*****************************************************************************/
/*
 *  Function to return the number of blocks in ordering.
 *  Input:
 *    order_info          --  The ordering computed by Zoltan_Order.
 *  Returned value:       --  The number of blocks in ordering.
 */

int Zoltan_Order_Get_Num_Blocks(
   struct Zoltan_Struct *zz
)
{
  return (zz->TPL_Order.nbr_blocks);
}

/*****************************************************************************/
/*
 *  Function to return the description of an ordering block
 *  Input:
 *    order_info          --  The ordering computed by Zoltan_Order.
 *    block_num           --  The id of the block to take care of.
 *  Output:
 *    first               --  Number of the first element of the block.
 *    last                --  Number of the last element of the block.
 *  For both, number means an indice between 0 and N-1, not in the GID domain.
 *  Returned value:       --  Error code
 */

int Zoltan_Order_Get_Block_Bounds(
   struct Zoltan_Struct *zz,
  int                         block_num,   /* Number of the wanted block */
  int                        *first,       /* First element in block */
  int                        *last         /* Last element in block */
  )
{
  if (block_num >= zz->TPL_Order.nbr_blocks)
    return (ZOLTAN_FATAL);

  *first = zz->TPL_Order.start[block_num];
  *last  = zz->TPL_Order.start[block_num + 1];
  return (ZOLTAN_OK);
}

/*****************************************************************************/
/*
 *  Function to return the number of elements within a block
 *  Input:
 *    order_info          --  The ordering computed by Zoltan_Order.
 *    block_num           --  The id of the block to take care of.
 *  Returned value:       --  Number of elements in the block.
 */

int Zoltan_Order_Get_Block_Size(
  struct Zoltan_Struct *zz,
  int                         block_num   /* Number of the wanted block */
)
{
  if (block_num >= zz->TPL_Order.nbr_blocks)
    return (-1);
  return (zz->TPL_Order.start[block_num+1] - zz->TPL_Order.start[block_num]);
}

/*****************************************************************************/
/*
 *  Function to return the indice of the parent block in the elimination tree.
 *  Input:
 *    order_info          --  The ordering computed by Zoltan_Order.
 *    block_num           --  The id of the block to take care of.
 *  Returned value:       --  Indice of the father, -1 if block is the root.
 */

int Zoltan_Order_Get_Block_Parent(
  struct Zoltan_Struct *zz,
  int                         block_num   /* Number of the wanted block */
)
{
 if (block_num >= zz->TPL_Order.nbr_blocks)
    return (-2);
 return (zz->TPL_Order.ancestor[block_num]);
}

/*****************************************************************************/
/*
 *  Function to return the number of the leaves in the elimination tree
 *  Input:
 *    zz                  --  The ordering computed by Zoltan_Order.
 *  Returned value:       --  Number of leaves in the elimination tree.
 */
int Zoltan_Order_Get_Num_Leaves(
  struct Zoltan_Struct *zz
)
{
  return(zz->TPL_Order.nbr_leaves);
}

/*****************************************************************************/
/*
 *  Function to return the list of the leaves in the elimination tree
 *  Input:
 *    zz                  --  The ordering computed by Zoltan_Order.
 *  Output:
 *    leaves              --  List of block indices that are leaves in the
 *                            elimination tree. -1 marks the end of the list.
 */

void Zoltan_Order_Get_Block_Leaves(
  struct Zoltan_Struct *zz,
  int                        *leaves
)
{
  if (zz->TPL_Order.nbr_leaves > 0)
    memcpy (leaves, zz->TPL_Order.leaves, (zz->TPL_Order.nbr_leaves+1)*sizeof(int));
  else
    *leaves = -1;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
