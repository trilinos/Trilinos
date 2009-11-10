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
#include "params_const.h"
#include "ha_const.h"
#include "order_const.h"
#include "hsfcOrder.h"

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
        { "ORDER_METHOD", NULL, "STRING", 0 },
        { "ORDER_START_INDEX", NULL, "INT", 0 },
        { "USE_ORDER_INFO", NULL, "INT", 0 },
        { NULL, NULL, NULL, 0 } };

int Zoltan_Order(
  ZZ *zz,               /* Zoltan structure */
  int num_gid_entries, /* # of entries for a global id */
  int num_obj,		/* Number of objects to order */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
                        /* The application must allocate enough space */
  int *rank,            /* rank[i] is the rank of gids[i] */
  int *iperm            /* iperm[rank[i]]=i, only for sequential ordering */
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
 * Return values:
 *   Zoltan error code.
 */

  char *yo = "Zoltan_Order";
  int ierr;
  int counti,myrank;
  double start_time, end_time;
  double order_time[2] = {0.0,0.0};
  char msg[256];
  int comm[2],gcomm[2];
  ZOLTAN_ORDER_FN *Order_fn;
  struct Zoltan_Order_Options opt;
  int * vtxdist = NULL;
  ZOLTAN_ID_PTR local_gids=NULL, lids=NULL;
  int local_num_obj;
  int *local_rank = NULL, *local_iperm=NULL;
  struct Zoltan_DD_Struct *dd = NULL;

  /* MMW: needed to temporarily support (ignore) iperm for local hsfc */
  int *tmpIperm;


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
  zz->Num_GID = gcomm[0];

  if (num_gid_entries != zz->Num_GID) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "num_gid_entries doesn't have the good value");
    return (ZOLTAN_FATAL);
  }


  zz->Order.nbr_objects = num_obj;
  zz->Order.rank = rank;
  zz->Order.iperm = iperm;
  zz->Order.gids = gids;
  zz->Order.lids = lids;
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
  Zoltan_Bind_Param(Order_params, "ORDER_START_INDEX", (void *) &opt.start_index);
  Zoltan_Bind_Param(Order_params, "USE_ORDER_INFO", (void *) &opt.use_order_info);

  Zoltan_Assign_Param_Vals(zz->Params, Order_params, zz->Debug_Level,
                           zz->Proc, zz->Debug_Proc);

  zz->Order.start_index = opt.start_index;

  /*
   *  Check that the user has allocated space for the return args.
   */
  if (!(gids && rank)){
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

    /*MMW: tmp code to ignore iperm */
    tmpIperm = iperm;
    iperm = NULL;
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
  else {
    fprintf(stderr, "%s\n", opt.method);
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


  local_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, local_num_obj);
  local_rank = (int*) ZOLTAN_MALLOC(local_num_obj*sizeof(int));
  local_iperm = (int*) ZOLTAN_MALLOC(local_num_obj*sizeof(int));

  lids = ZOLTAN_MALLOC_LID_ARRAY(zz, local_num_obj);

  /*
   * Call the actual ordering function.
   * Compute gid according to the local graph.
   */

  ierr = (*Order_fn)(zz, local_num_obj, local_gids, lids, local_rank, local_iperm, &opt);
  ZOLTAN_FREE(&lids);

  
  if (ierr) {
    sprintf(msg, "Ordering routine returned error code %d.", ierr);
    if (ierr == ZOLTAN_WARN){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    } else {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
      Zoltan_Multifree(__FILE__, __LINE__, 3,
                       &local_gids, &local_rank, &local_iperm);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done ordering");


/*   Compute inverse permutation if necessary */
  if ((!(opt.return_args & RETURN_RANK) && (rank != NULL))
      || (!(opt.return_args & RETURN_IPERM) && (iperm != NULL))) {
    ierr = Zoltan_Get_Distribution(zz, &vtxdist);
    if (ierr){
      /* Error */
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_Get_Distribution.\n");
      return (ierr);
    }

    if (!(opt.return_args & RETURN_RANK) && (rank != NULL)){
      /* Compute rank from iperm */
      ZOLTAN_TRACE_DETAIL(zz, yo, "Inverting permutation");
      Zoltan_Inverse_Perm(zz, local_iperm, local_rank, vtxdist, zz->Order.order_type, opt.start_index); 
    }
    else if (!(opt.return_args & RETURN_IPERM) && (iperm != NULL)){
    /* Compute iperm from rank */
      ZOLTAN_TRACE_DETAIL(zz, yo, "Inverting permutation");
      Zoltan_Inverse_Perm(zz, local_rank, local_iperm, vtxdist, zz->Order.order_type, opt.start_index); 
    }
    ZOLTAN_FREE(&vtxdist);
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done Invert Permutation");


  /* TODO: Use directly the "graph" structure to avoid to duplicate things. */
  /* I store : GNO, rank, iperm */
  /* MMW: perhaps don't ever use graph here since we need to support geometric orderings, otherwise need if/else */
  ierr = Zoltan_DD_Create (&dd, zz->Communicator, zz->Num_GID, (local_rank==NULL)?0:1, (local_iperm==NULL)?0:1, local_num_obj, 0); 
  /* Hope a linear assignment will help a little */
  Zoltan_DD_Set_Neighbor_Hash_Fn1(dd, local_num_obj); 
  /* Associate all the data with our xGNO */

  Zoltan_DD_Update (dd, local_gids, (ZOLTAN_ID_PTR)local_rank, (ZOLTAN_ID_PTR) local_iperm, NULL, local_num_obj);

  ZOLTAN_FREE(&local_gids);
  ZOLTAN_FREE(&local_rank);
  ZOLTAN_FREE(&local_iperm);

  Zoltan_DD_Find (dd, gids, (ZOLTAN_ID_PTR)rank, (ZOLTAN_ID_PTR)iperm, NULL,
		  num_obj, NULL);
  Zoltan_DD_Destroy(&dd);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done Registering results");


  if (!strcmp(opt.method, "LOCAL_HSFC")) 
  {
    /*MMW: tmp code to reset iperm pointer*/
    iperm = tmpIperm;
  }


  end_time = Zoltan_Time(zz->Timer);
  order_time[0] = end_time - start_time;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_LIST) {
    int i;
    Zoltan_Print_Sync_Start(zz->Communicator, TRUE);
    printf("ZOLTAN: rank for ordering on Proc %d\n", zz->Proc);
    for (i = 0; i < num_obj; i++) {
      printf("GID = ");
      ZOLTAN_PRINT_GID(zz, &(gids[i*(num_gid_entries)]));
      printf(", rank = %3d\n", rank[i]);
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
  return (zz->Order.nbr_blocks);
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
  if (block_num >= zz->Order.nbr_blocks)
    return (ZOLTAN_FATAL);

  *first = zz->Order.start[block_num];
  *last  = zz->Order.start[block_num + 1];
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
  if (block_num >= zz->Order.nbr_blocks)
    return (-1);
  return (zz->Order.start[block_num+1] - zz->Order.start[block_num]);
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
 if (block_num >= zz->Order.nbr_blocks)
    return (-2);
 return (zz->Order.ancestor[block_num]);
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
  return(zz->Order.nbr_leaves);
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
  if (zz->Order.nbr_leaves > 0)
    memcpy (leaves, zz->Order.leaves, (zz->Order.nbr_leaves+1)*sizeof(int));
  else
    *leaves = -1;
}

/*****************************************************************************/
/*
 *  Function to return the ordering on the GID
 *  Input:
 *    zz                  --  The Zoltan structure containing
 *                            info for this load-balancing invocation.
 *    order_info          --  The ordering computed by Zoltan_Order.
 *    gids                --  List of global ids.
 *  Ouput:
 *    order_ids           --  New ordering of the gids.
 *  Returned value:       --  Error Code.
 */

int Zoltan_Order_Get_GID_Order(
  struct Zoltan_Struct       *zz,
  ZOLTAN_ID_PTR               global_ids,
  ZOLTAN_ID_PTR               order_ids
  )
{
  int * proctab;
  ZOLTAN_ID_PTR sendgidtab;
  int * sendtab, *recvtab;
  int i;
  int nrecv;
  struct Zoltan_Comm_Obj *plan;
  int ierr = ZOLTAN_OK;
  int *vtxdist;
  int tag = 24061986;
  int offset;

  zz->Order.gidrank = order_ids;

  if (!strcmp(zz->Order.order_type, "LOCAL")) {
    for (i = 0 ; i < zz->Order.nbr_objects ; ++i) {
      order_ids[i] = global_ids[zz->Order.rank[i] - zz->Order.start_index];
    }
    return (ZOLTAN_OK);
  }

  ierr = Zoltan_Get_Distribution(zz, &vtxdist);
  if (ierr != ZOLTAN_OK) return (ierr);

  proctab = (int*) ZOLTAN_MALLOC(3*zz->Order.nbr_objects*sizeof(int));
  sendgidtab = ZOLTAN_MALLOC_GID_ARRAY(zz, 2*zz->Order.nbr_objects);
  if (proctab == NULL)  {
    ZOLTAN_FREE(&vtxdist); return (ZOLTAN_MEMERR);
  }
  if (sendgidtab == NULL) {
    ZOLTAN_FREE(&proctab);
    ZOLTAN_FREE(&vtxdist); return (ZOLTAN_MEMERR);
  }
  sendtab = proctab + zz->Order.nbr_objects;
  recvtab = sendtab + zz->Order.nbr_objects;

  for (i = 0 ; i < zz->Order.nbr_objects ; ++i) {
    proctab[i] = Zoltan_Get_Processor_Graph(vtxdist, zz->Num_Proc, zz->Order.rank[i] - zz->Order.start_index);
    sendtab[i] = zz->Order.rank[i] - zz->Order.start_index;
  }

  ierr = Zoltan_Comm_Create(&plan, zz->Order.nbr_objects, proctab, zz->Communicator, tag++,
                     &nrecv);
  if (ierr != ZOLTAN_OK) {
    ZOLTAN_FREE(&sendgidtab);
    ZOLTAN_FREE(&proctab);
    ZOLTAN_FREE(&vtxdist);
    return (ierr);
  }

  ierr = Zoltan_Comm_Do(plan, tag++, (char *) sendtab, sizeof(int), (char *) recvtab);

  if (ierr != ZOLTAN_OK) {
    ZOLTAN_FREE(&sendgidtab);
    ZOLTAN_FREE(&proctab);
    ZOLTAN_FREE(&vtxdist);
    return (ierr);
  }

  offset = vtxdist[zz->Proc];
  for (i = 0 ; i < zz->Order.nbr_objects ; ++i) {
    int position;
    position = recvtab[i]-offset;
    ZOLTAN_SET_GID(zz, &sendgidtab[i], &global_ids[position]);
  }

  ierr = Zoltan_Comm_Do_Reverse(plan, tag++, (char *) sendgidtab, zz->Num_GID*sizeof(int), NULL,
				(char *) order_ids);

  Zoltan_Comm_Destroy(&plan);
  ZOLTAN_FREE(&sendgidtab);
  ZOLTAN_FREE(&proctab);
  ZOLTAN_FREE(&vtxdist);

  return (ierr);
}




#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
