/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2008 Sandia National Laboratories.                          *
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


#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "all_allo_const.h"
#include "params_const.h"
#include "order_const.h"
#include "third_library.h"
#include "scotch_interface.h"

#ifndef ZOLTAN_PTSCOTCH
#define Scotch_Dgraph Scotch_Graph
#endif /* ZOLTAN_PTSCOTCH */

  /**********  parameters structure for Scotch methods **********/
static PARAM_VARS Scotch_params[] = {
  { "SCOTCH_METHOD", NULL, "STRING", 0 },
  { "SCOTCH_TYPE", NULL, "STRING", 0 },
  { "SCOTCH_STRAT", NULL, "STRING", 0 },
  { "SCOTCH_STRAT_FILE", NULL, "STRING", 0 },
  { NULL, NULL, NULL, 0 } };

static int Zoltan_Scotch_Bind_Param(ZZ * zz, char *alg, char* graph_type, char **strat);

static int
Zoltan_Scotch_Construct_Offset(ZOS *order,
			       indextype *children,
			       int root,
			       indextype* size,
			       indextype* tree,
			       int offset, int *leafnum);

static int compar_int (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

static int
Zoltan_Scotch_Build_Graph(ZOLTAN_Third_Graph * gr,
			  MPI_Comm             comm,
#ifdef ZOLTAN_PTSCOTCH
			  SCOTCH_Dgraph *      dgrafptr,
#endif
			  SCOTCH_Graph *       cgrafptr,
			  SCOTCH_Strat *       stratptr);



/***************************************************************************
 *  The Scotch ordering routine piggy-backs on the Scotch
 *  partitioning routines.
 **************************************************************************/

/*
 * TODO: at this time, only distributed computations are allowed.
 */

int Zoltan_Scotch_Order(
  ZZ *zz,               /* Zoltan structure */
  int num_obj,		/* Number of (local) objects to order. */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
  /* The application must allocate enough space */
  ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
/* The application must allocate enough space */
  int *rank,		/* rank[i] is the rank of gids[i] */
  int *iperm,
  ZOOS *order_opt 	/* Ordering options, parsed by Zoltan_Order */
)
{
  static char *yo = "Zoltan_Scotch_Order";
  int n, ierr;
  ZOLTAN_Output_Order ord;
  ZOLTAN_Third_Graph gr;
  SCOTCH_Strat        stradat;
  SCOTCH_Graph        cgrafdat;
#ifdef ZOLTAN_PTSCOTCH
  SCOTCH_Dgraph       grafdat;
  SCOTCH_Dordering    ordedat;
  /* The following are used to convert elimination tree in Zoltan format */
  indextype          *tree;
  indextype          *size;
  indextype          *children;
  indextype           leafnum;
  indextype           start;
  int root = -1;
#endif /* ZOLTAN_PTSCOTCH */
  int numbloc;

  MPI_Comm comm = zz->Communicator;/* want to risk letting external packages */
  int timer_p = 0;
  int get_times = 0;
  int use_timers = 0;
  double times[5];

  char alg[MAX_PARAM_STRING_LEN+1];
  char graph_type[MAX_PARAM_STRING_LEN+1];
  char *strat = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  memset(&gr, 0, sizeof(ZOLTAN_Third_Graph));
  memset(&ord, 0, sizeof(ZOLTAN_Output_Order));

  strcpy (alg, "NODEND");
  strcpy (graph_type, "GLOBAL");
  ierr = Zoltan_Scotch_Bind_Param(zz, alg, graph_type, &strat);
  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ierr);
  }

  ord.order_info = &(zz->Order);
  ord.order_opt = order_opt;

  if (!order_opt){
    /* If for some reason order_opt is NULL, allocate a new ZOOS here. */
    /* This should really never happen. */
    order_opt = (ZOOS *) ZOLTAN_MALLOC(sizeof(ZOOS));
    strcpy(order_opt->method,"SCOTCH");
  }

  /* Scotch only computes the rank vector */
  order_opt->return_args = RETURN_RANK;

  /* Check that num_obj equals the number of objects on this proc. */
  /* This constraint may be removed in the future. */
  n = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  if ((ierr!= ZOLTAN_OK) && (ierr!= ZOLTAN_WARN)){
    /* Return error code */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Get_Num_Obj returned error.");
    return(ZOLTAN_FATAL);
  }
  if (n != num_obj){
    /* Currently this is a fatal error. */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input num_obj does not equal the number of objects.");
    return(ZOLTAN_FATAL);
  }

  /* Do not use weights for ordering */
/*   gr.obj_wgt_dim = -1; */
/*   gr.edge_wgt_dim = -1; */
  gr.num_obj = num_obj;

  /* Check what ordering type is requested */
  gr.graph_type = - LOCAL_GRAPH;
  if (order_opt){
    if (strcmp(order_opt->order_type, "LOCAL") == 0)
      gr.graph_type = - LOCAL_GRAPH;
#ifdef ZOLTAN_PTSCOTCH
    else if (strcmp(order_opt->order_type, "GLOBAL") == 0)
      gr.graph_type =  - GLOBAL_GRAPH;
#endif /* ZOLTAN_PTSCOTCH */
    else
      gr.graph_type = - NO_GRAPH;
  }
  gr.get_data = 1;

  /* If reorder is true, we already have the id lists. Ignore weights. */
  if ((order_opt && order_opt->reorder))
    gr.id_known = 1;                        /* We already have global_ids and local_ids */

  timer_p = Zoltan_Preprocess_Timer(zz, &use_timers);

    /* Start timer */
  get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
  if (get_times){
    MPI_Barrier(zz->Communicator);
    times[0] = Zoltan_Time(zz->Timer);
  }

  ierr = Zoltan_Preprocess_Graph(zz, &gids, &lids,  &gr, NULL, NULL, NULL);


  if (Zoltan_Scotch_Build_Graph(&gr, comm,
#ifdef ZOLTAN_PTSCOTCH
 &grafdat,
#endif
 &cgrafdat, &stradat) != ZOLTAN_OK) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot construct Scotch Graph");
  }

  if (strat != NULL) {
    if (
#ifdef ZOLTAN_PTSCOTCH
	((gr.graph_type==GLOBAL_GRAPH) && (SCOTCH_stratDgraphOrder (&stradat, strat)) != 0) ||
#endif /* ZOLTAN_PTSCOTCH */
	(SCOTCH_stratGraphOrder (&stradat, strat) != 0)) {
#ifdef ZOLTAN_PTSCOTCH
      if (gr.graph_type == GLOBAL_GRAPH) SCOTCH_dgraphExit(&grafdat);
      else
#endif /* ZOLTAN_PTSCOTCH */
	SCOTCH_graphExit(&cgrafdat);
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      return (ZOLTAN_FATAL);
    }
  }

  /* Allocate space for rank array */
  ord.rank = (indextype *) ZOLTAN_MALLOC(gr.num_obj*sizeof(indextype));
  if (!ord.rank){
    /* Not enough memory */
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
  }

  ord.iperm = NULL;
  if ((gr.graph_type!=GLOBAL_GRAPH) && (iperm != NULL)) {
  /* Allocate space for inverse perm */
    order_opt->return_args |= RETURN_IPERM;
    ord.iperm = (indextype *) ZOLTAN_MALLOC(gr.num_obj*sizeof(indextype));
    if (!ord.iperm){
      /* Not enough memory */
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
  }

  if (gr.graph_type != GLOBAL_GRAPH) { /* Allocate separators tree */
    if (Zoltan_Order_Init_Tree (&zz->Order, gr.num_obj + 1, gr.num_obj) != ZOLTAN_OK) {
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
  }

#ifdef ZOLTAN_PTSCOTCH
  if ((gr.graph_type==GLOBAL_GRAPH) && (SCOTCH_dgraphOrderInit (&grafdat, &ordedat) != 0)) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot construct Scotch graph.");
  }
#endif /* ZOLTAN_PTSCOTCH */

  /* Get a time here */
  if (get_times) times[1] = Zoltan_Time(zz->Timer);

#ifdef ZOLTAN_PTSCOTCH
  if (gr.graph_type==GLOBAL_GRAPH){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the PT-Scotch library");
    ierr = SCOTCH_dgraphOrderCompute (&grafdat, &ordedat, &stradat);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the PT-Scotch library");
  }
  else
#endif /* ZOLTAN_PTSCOTCH */
    {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the Scotch library");
    ierr = SCOTCH_graphOrder (&cgrafdat,  &stradat, ord.rank, ord.iperm,
				     &zz->Order.nbr_blocks, zz->Order.start, zz->Order.ancestor);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the Scotch library");
  }

  if (ierr != 0) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot compute Scotch ordering.");
  }

  /* Get a time here */
  if (get_times) times[2] = Zoltan_Time(zz->Timer);

  if (gr.graph_type != GLOBAL_GRAPH) { /* We already have separator tree, just have to compute the leaves */
    for (numbloc = 0 ; numbloc < zz->Order.nbr_blocks ; ++numbloc) {
      zz->Order.leaves[numbloc] = numbloc;
    }
    for (numbloc = 0 ; numbloc < zz->Order.nbr_blocks ; ++numbloc) {
      if (zz->Order.ancestor[numbloc] < 0)
	continue;
      zz->Order.leaves[zz->Order.ancestor[numbloc]] = zz->Order.nbr_blocks + 1;
    }
    /* TODO : check if there is a normalized sort in Zoltan */
    qsort(zz->Order.leaves, zz->Order.nbr_blocks, sizeof(int), compar_int);
    zz->Order.nbr_leaves = 0;
    for (numbloc = 0 ; numbloc < zz->Order.nbr_blocks ; ++numbloc) {
      if (zz->Order.leaves[numbloc] > zz->Order.nbr_blocks) {
	zz->Order.leaves[numbloc] = -1;
	zz->Order.nbr_leaves = numbloc;
	break;
      }
    }
    if (zz->Order.nbr_leaves == 0) {
      zz->Order.leaves[zz->Order.nbr_blocks] = -1;
      zz->Order.nbr_leaves = zz->Order.nbr_blocks;
    }

  }
#ifdef ZOLTAN_PTSCOTCH
  else{
    /* Compute permutation */
    if (SCOTCH_dgraphOrderPerm (&grafdat, &ordedat, ord.rank) != 0) {
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot compute Scotch rank array");
    }

    /* Construct elimination tree */
    zz->Order.nbr_blocks = SCOTCH_dgraphOrderCblkDist (&grafdat, &ordedat);
    if (zz->Order.nbr_blocks <= 0) {
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot compute Scotch block");
    }

    if (Zoltan_Order_Init_Tree (&zz->Order, 2*zz->Num_Proc, zz->Num_Proc) != ZOLTAN_OK) {
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }

    tree = (indextype *) ZOLTAN_MALLOC((zz->Order.nbr_blocks+1)*sizeof(indextype));
    size = (indextype *) ZOLTAN_MALLOC((zz->Order.nbr_blocks+1)*sizeof(indextype));

    if ((tree == NULL) || (size == NULL)){
      /* Not enough memory */
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }

    if (SCOTCH_dgraphOrderTreeDist (&grafdat, &ordedat, tree, size) != 0) {
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot compute Scotch rank array");
    }

    children = (indextype *) ZOLTAN_MALLOC(3*zz->Order.nbr_blocks*sizeof(indextype));
    for (numbloc = 0 ; numbloc < 3*zz->Order.nbr_blocks ; ++numbloc) {
      children[numbloc] = -2;
    }

    /* Now convert scotch separator tree in Zoltan elimination tree */
    root = -1;
    for (numbloc = 0 ; numbloc < zz->Order.nbr_blocks ; ++numbloc) { /* construct a top-bottom tree */
      indextype tmp;
      int index=0;

      tmp = tree[numbloc];
      if (tmp == -1) {
	root = numbloc;
	continue;
      }
      while ((index<3) && (children[3*tmp+index] > 0))
	index ++;

      if (index >= 3) {
	Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
	ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot compute Scotch tree array");
      }

      children[3*tmp+index] = numbloc;
    }

    leafnum = 0;
    zz->Order.nbr_blocks = Zoltan_Scotch_Construct_Offset(&zz->Order, children, root, size, tree, 0, &leafnum);
    zz->Order.leaves[leafnum] =-1;
    zz->Order.nbr_leaves = leafnum;

    for (numbloc = 0, start=0 ; numbloc < zz->Order.nbr_blocks ; ++numbloc) {
      int tmp;
      tmp = zz->Order.start[numbloc];
      zz->Order.start[numbloc]  = start;
      start += tmp;
      if (zz->Order.ancestor[numbloc] >= 0)
	zz->Order.ancestor[numbloc] = size[zz->Order.ancestor[numbloc]];
    }
    zz->Order.start[zz->Order.nbr_blocks]  = start;

    /* Free temporary tables */
    ZOLTAN_FREE(&tree);
    ZOLTAN_FREE(&size);
    ZOLTAN_FREE(&children);
  }
#endif /* ZOLTAN_PTSCOTCH */

  ierr = Zoltan_Postprocess_Graph (zz, gids, lids, &gr, NULL, NULL, NULL, &ord, NULL);

  /* Get a time here */
  if (get_times) times[3] = Zoltan_Time(zz->Timer);

#ifdef ZOLTAN_PTSCOTCH
  if (gr.graph_type==GLOBAL_GRAPH) {
    SCOTCH_dgraphOrderExit (&grafdat, &ordedat);
    SCOTCH_dgraphExit (&grafdat);
  }
  else
#endif /* ZOLTAN_PTSCOTCH */
  {
    SCOTCH_graphExit (&cgrafdat);
  }
  SCOTCH_stratExit (&stradat);


  if (get_times) Zoltan_Third_DisplayTime(zz, times);

  if (use_timers)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_p, zz->Communicator);

  if ((ord.iperm != NULL) && (iperm != NULL))
    memcpy(iperm, ord.iperm, gr.num_obj*sizeof(indextype));
  if (ord.iperm != NULL)  ZOLTAN_FREE(&ord.iperm);
  if (order_opt->return_args&RETURN_RANK)
    memcpy(rank, ord.rank, gr.num_obj*sizeof(indextype));
  ZOLTAN_FREE(&ord.rank);
  ZOLTAN_FREE(&strat);

  /* Free all other "graph" stuff */
  Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, NULL);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ZOLTAN_OK);
}


static int
Zoltan_Scotch_Construct_Offset(ZOS *order, indextype *children, int root,
			       indextype* size, indextype* tree, int offset, int *leafnum)
{
  int i = 0;
  int childrensize = 0;


  for (i=0 ; i < 2 ; i++) {
    if (children[3*root+i] < 0)
      break;

    childrensize += size[children[3*root+i]];
    offset = Zoltan_Scotch_Construct_Offset(order, children, children[3*root+i], size, tree, offset, leafnum);
  }


  order->start[offset] = size[root] - childrensize;
  order->ancestor[offset] = tree[root];
  size[root] = offset; /* size[root] not used now, can be use to convert indices */
  if (childrensize == 0)  { /* Leaf */
    order->leaves[*leafnum] = offset;
    (*leafnum)++;
  }
  ++offset;
  return (offset);
}



  /**********************************************************/
  /* Interface routine for PT-Scotch: Partitioning          */
  /**********************************************************/

/* TODO: convert multi-weighted vertices to 1 weight */

int Zoltan_Scotch(
		    ZZ *zz,               /* Zoltan structure */
		    float *part_sizes,    /* Input:  Array of size zz->Num_Global_Parts
					     containing the percentage of work to be
					     assigned to each partition.               */
		    int *num_imp,         /* number of objects to be imported */
		    ZOLTAN_ID_PTR *imp_gids,  /* global ids of objects to be imported */
		    ZOLTAN_ID_PTR *imp_lids,  /* local  ids of objects to be imported */
		    int **imp_procs,      /* list of processors to import from */
		    int **imp_to_part,    /* list of partitions to which imported objects are
					     assigned.  */
		    int *num_exp,         /* number of objects to be exported */
		    ZOLTAN_ID_PTR *exp_gids,  /* global ids of objects to be exported */
		    ZOLTAN_ID_PTR *exp_lids,  /* local  ids of objects to be exported */
		    int **exp_procs,      /* list of processors to export to */
		    int **exp_to_part     /* list of partitions to which exported objects are
					     assigned. */
		    )
{
  char *yo = "Zoltan_Scotch";
  int ierr;
  ZOLTAN_Third_Graph gr;
  ZOLTAN_Third_Geom  *geo = NULL;
  ZOLTAN_Third_Vsize vsp;
  ZOLTAN_Third_Part  prt;
  ZOLTAN_Output_Part part;

  ZOLTAN_ID_PTR global_ids = NULL;
  ZOLTAN_ID_PTR local_ids = NULL;

  SCOTCH_Strat        stradat;
#ifdef ZOLTAN_PTSCOTCH
  SCOTCH_Dgraph       grafdat;
#endif
  SCOTCH_Graph        cgrafdat;
  SCOTCH_Arch         archdat;

  int use_timers = 0;
  int timer_p = -1;
  int get_times = 0;
  double times[5];

  char alg[MAX_PARAM_STRING_LEN+1];
  char graph_type[MAX_PARAM_STRING_LEN+1];
  char *strat = NULL;
  weighttype * goal_sizes=NULL;
  weighttype velosum=0;


  int retval=ZOLTAN_OK;
  int i;
  int num_part = zz->LB.Num_Global_Parts;    /* passed to PT-Scotch. Don't             */
  MPI_Comm comm = zz->Communicator;          /* want to risk letting external packages */
                                             /* change our zz struct.                  */


#ifndef ZOLTAN_SCOTCH
  ZOLTAN_PRINT_ERROR(zz->Proc, __func__,
		     "Scotch requested but not compiled into library.");
  return ZOLTAN_FATAL;

#endif /* ZOLTAN_SCOTCH */

  ZOLTAN_TRACE_ENTER(zz, yo);

  Zoltan_Third_Init(&gr, &prt, &vsp, &part,
		    imp_gids, imp_lids, imp_procs, imp_to_part,
		    exp_gids, exp_lids, exp_procs, exp_to_part);

  prt.input_part_sizes = prt.part_sizes = part_sizes;

  strcpy (alg, "RBISECT");
  strcpy (graph_type, "GLOBAL");
  ierr = Zoltan_Scotch_Bind_Param(zz, alg, graph_type, &strat);

  timer_p = Zoltan_Preprocess_Timer(zz, &use_timers);

    /* Start timer */
  get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
  if (get_times){
    MPI_Barrier(zz->Communicator);
    times[0] = Zoltan_Time(zz->Timer);
  }

#ifdef ZOLTAN_PTSCOTCH
  gr.graph_type = - GLOBAL_GRAPH;
  /* Fix type of graph, negative because we impose them */
  if (strcmp (graph_type, "GLOBAL") != 0) {
    gr.graph_type = - LOCAL_GRAPH;
    if (zz->Num_Proc > 1) {
      ZOLTAN_PRINT_ERROR(zz->Proc, __func__, "Distributed graph: cannot call serial Scotch, switching to PT-Scotch");
      gr.graph_type = - GLOBAL_GRAPH;
      retval = ZOLTAN_WARN;
    }
  }
#else
  gr.graph_type = - LOCAL_GRAPH;
#endif /* ZOLTAN_PTSCOTCH */

  /* TODO : take care about multidimensional weights */
  ierr = Zoltan_Preprocess_Graph(zz, &global_ids, &local_ids,  &gr, geo, &prt, &vsp);

  /* Get a time here */
  if (get_times) times[1] = Zoltan_Time(zz->Timer);

  if (gr.obj_wgt_dim > 1) {
    Zoltan_Third_Exit(&gr, NULL, &prt, &vsp, NULL, NULL);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Scotch cannot deal with more than 1 weight.");
  }

  if (Zoltan_Scotch_Build_Graph(&gr, comm,
#ifdef ZOLTAN_PTSCOTCH
 &grafdat,
#endif
 &cgrafdat, &stradat) != ZOLTAN_OK) {
    Zoltan_Third_Exit(&gr, NULL, &prt, &vsp, NULL, NULL);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Cannot construct Scotch graph.");
  }

/*   { */
/*     FILE * graphfile; */
/*     char name[80]; */

/*     sprintf(name, "kdd.%d.src", zz->Proc); */
/*     graphfile = fopen(name, "w+"); */
/*     SCOTCH_dgraphSave(&grafdat, graphfile); */
/*     fclose(graphfile); */
/*   } */

  if (!prt.part_sizes){
    Zoltan_Third_Exit(&gr, NULL, &prt, &vsp, NULL, NULL);
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,"Input parameter part_sizes is NULL.");
  }

  if ((goal_sizes = (weighttype *) ZOLTAN_MALLOC(num_part * sizeof(weighttype))) == NULL) {
    Zoltan_Third_Exit(&gr, NULL, &prt, &vsp, NULL, NULL);
    ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR,"Part_sizes cannot be allocated.");
  }


  if (strat != NULL) {
    if (
#ifdef ZOLTAN_PTSCOTCH
	((gr.graph_type==GLOBAL_GRAPH) && (SCOTCH_stratDgraphMap (&stradat, strat)) != 0) ||
#endif /* ZOLTAN_PTSCOTCH */
	(SCOTCH_stratGraphMap (&stradat, strat) != 0)) {
#ifdef ZOLTAN_PTSCOTCH
      if (gr.graph_type == GLOBAL_GRAPH) SCOTCH_dgraphExit(&grafdat);
      else
#endif /* ZOLTAN_PTSCOTCH */
	SCOTCH_graphExit(&cgrafdat);
      Zoltan_Third_Exit(&gr, NULL, &prt, &vsp, NULL, NULL);
      return (ZOLTAN_FATAL);
    }
  }


  /* Compute size we want for each part */
#ifdef ZOLTAN_PTSCOTCH
  if (gr.graph_type == GLOBAL_GRAPH)
    SCOTCH_dgraphSize (&grafdat, NULL, &velosum, NULL, NULL);
  else
#endif /* ZOLTAN_PTSCOTCH */
    SCOTCH_graphStat (&cgrafdat, NULL, NULL, &velosum, NULL, NULL,
		      NULL, NULL, NULL, NULL,
		      NULL, NULL, NULL, NULL, NULL);

  for (i=0; i < num_part ; ++i) {
    goal_sizes[i] = ceil(prt.part_sizes[i]*velosum);
  }

  /* Construct a complete-graph like */
  SCOTCH_archInit(&archdat);
  SCOTCH_archCmpltw(&archdat, num_part, goal_sizes);

  /* Really Call PT-Scotch or Scotch) */
#ifdef ZOLTAN_PTSCOTCH
  if (gr.graph_type == GLOBAL_GRAPH) {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the PT-Scotch library");
    if (SCOTCH_dgraphMap(&grafdat, &archdat, &stradat, prt.part) != 0) {
      SCOTCH_archExit(&archdat);
      ZOLTAN_FREE(&goal_sizes);
      SCOTCH_dgraphExit(&grafdat);
      Zoltan_Third_Exit(&gr, NULL, &prt, &vsp, NULL, NULL);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,"PT-Scotch partitioning internal error.");
    }
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the PT-Scotch library");
  }
  else
#endif /* ZOLTAN_PTSCOTCH */
 {
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the Scotch library");
    if (SCOTCH_graphMap(&cgrafdat, &archdat, &stradat, prt.part) != 0) {
      SCOTCH_archExit(&archdat);
      ZOLTAN_FREE(&goal_sizes);
      SCOTCH_graphExit(&cgrafdat);
      Zoltan_Third_Exit(&gr, NULL, &prt, &vsp, NULL, NULL);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,"Scotch partitioning internal error.");
    }
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the Scotch library");
  }

  /* Get a time here */
  if (get_times) times[2] = Zoltan_Time(zz->Timer);

#ifdef ZOLTAN_PTSCOTCH
  if (gr.graph_type == GLOBAL_GRAPH) {
    SCOTCH_dgraphExit(&grafdat);
  }
  else
#endif /* ZOLTAN_PTSCOTCH */
  {
    SCOTCH_graphExit(&cgrafdat);
  }
  SCOTCH_archExit(&archdat);
  ZOLTAN_FREE(&goal_sizes);

  ierr = Zoltan_Postprocess_Graph(zz, global_ids, local_ids, &gr, geo, &prt, &vsp, NULL, &part);

  Zoltan_Third_Export_User(&part, num_imp, imp_gids, imp_lids, imp_procs, imp_to_part,
			   num_exp, exp_gids, exp_lids, exp_procs, exp_to_part);

  /* Get a time here */
  if (get_times) times[3] = Zoltan_Time(zz->Timer);

  if (get_times) Zoltan_Third_DisplayTime(zz, times);

  if (use_timers && timer_p >= 0)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_p, zz->Communicator);

  if (gr.final_output) {
    ierr = Zoltan_Postprocess_FinalOutput (zz, &gr, &prt, &vsp,
					   use_timers, 0);
  }

  Zoltan_Third_Exit(&gr, NULL, &prt, &vsp, NULL, NULL);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);

  ZOLTAN_TRACE_EXIT(zz, yo);

  if (ierr == ZOLTAN_OK)
    return (retval);

  return (ierr);
}


/************************
 * Auxiliary function used to parse scotch specific parameters
 ***************************************/

static int Zoltan_Scotch_Bind_Param(ZZ* zz, char *alg, char *graph_type, char **strat)
{
  static char * yo = "Zoltan_Scotch_Bind_Param";
  char stratsmall[MAX_PARAM_STRING_LEN+1];
  char stratfilename[MAX_PARAM_STRING_LEN+1];

  *strat = NULL;
  stratsmall[0] = stratfilename[0] = '\0';
  Zoltan_Bind_Param(Scotch_params, "SCOTCH_METHOD",
		    (void *) alg);
  Zoltan_Bind_Param(Scotch_params, "SCOTCH_TYPE",
		    (void *) graph_type);
  Zoltan_Bind_Param(Scotch_params, "SCOTCH_STRAT",
		    (void *) stratsmall);
  Zoltan_Bind_Param(Scotch_params, "SCOTCH_STRAT_FILE",
		    (void *) stratfilename);
  Zoltan_Assign_Param_Vals(zz->Params, Scotch_params, zz->Debug_Level,
			   zz->Proc, zz->Debug_Proc);

  if ((strlen(stratsmall) > 0) && (strlen(stratfilename) > 0)) {
    ZOLTAN_THIRD_ERROR(ZOLTAN_WARN,
		       "SCOTCH_STRAT and SCOTCH_STRAT_FILE both defined: ignoring\n");
  }

  if (strlen(stratsmall) > 0) {
    *strat = (char *) ZOLTAN_MALLOC((strlen(stratsmall)+1)*sizeof(char));
    strcpy (*strat, stratsmall);
    return (ZOLTAN_OK);
  }

  if (strlen(stratfilename) > 0) {
    long size;
    FILE *stratfile;

    stratfile = fopen(stratfilename, "r");
    if (stratfile == NULL) {
      ZOLTAN_THIRD_ERROR(ZOLTAN_WARN,
			 "Cannot open Scotch strategy file\n");
    }

    fseek(stratfile, (long)0, SEEK_END);
    size = ftell(stratfile);
    *strat = (char *) ZOLTAN_MALLOC((size+2)*sizeof(char));
    if (*strat == NULL) {
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    fseek(stratfile, (long)0, SEEK_SET);
    fgets (*strat, size+1, stratfile);
    fclose(stratfile);


    return (ZOLTAN_OK);
  }

  return (ZOLTAN_OK);
}


/*********************************************************************/
/* Scotch parameter routine                                          */
/*********************************************************************/

int Zoltan_Scotch_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
  int status, i;
  PARAM_UTYPE result;         /* value returned from Check_Param */
  int index;                  /* index returned from Check_Param */
  char *valid_methods[] = {
    "NODEND", /* for nested dissection ordering */
    "RBISECT",
    NULL };

  status = Zoltan_Check_Param(name, val, Scotch_params, &result, &index);
  if (status == 0){
    /* OK so far, do sanity check of parameter values */

    if (strcmp(name, "SCOTCH_METHOD") == 0){
      status = 2;
      for (i=0; valid_methods[i] != NULL; i++){
	if (strcmp(val, valid_methods[i]) == 0){
	  status = 0;
	  break;
	}
      }
    }
  }
  return(status);
}


/** Function to construct scotch graph, in sequential or in parallel
 */
static int
Zoltan_Scotch_Build_Graph(ZOLTAN_Third_Graph * gr,
			  MPI_Comm             comm,
#ifdef ZOLTAN_PTSCOTCH
			  SCOTCH_Dgraph *      dgrafptr,
#endif
			  SCOTCH_Graph *       cgrafptr,
			  SCOTCH_Strat *       stratptr)
{
  int edgelocnbr;

  edgelocnbr =  gr->xadj[gr->num_obj];
#ifdef ZOLTAN_PTSCOTCH
  if (gr->graph_type==GLOBAL_GRAPH){
    if (SCOTCH_dgraphInit (dgrafptr, comm) != 0) {
      return (ZOLTAN_FATAL);
    }

    if (SCOTCH_dgraphBuild (dgrafptr, 0, gr->num_obj, gr->num_obj, gr->xadj, gr->xadj + 1,
		  	    gr->vwgt, NULL,edgelocnbr, edgelocnbr, gr->adjncy, NULL, gr->ewgts) != 0) {
      SCOTCH_dgraphExit(dgrafptr);
      return (ZOLTAN_FATAL);
    }
  }
  else
#endif /* ZOLTAN_PTSCOTCH */
  {/* gr->graph_type==GLOBAL_GRAPH */
    if (SCOTCH_graphInit (cgrafptr) != 0) {
      return (ZOLTAN_FATAL);
    }

    if (SCOTCH_graphBuild (cgrafptr, 0, gr->num_obj, gr->xadj, gr->xadj + 1,
			   gr->vwgt, NULL,edgelocnbr, gr->adjncy, gr->ewgts) != 0) {
      SCOTCH_graphExit(cgrafptr);
      return (ZOLTAN_FATAL);
    }
  }

  SCOTCH_stratInit (stratptr);

  return (ZOLTAN_OK);
}



#ifdef __cplusplus
}
#endif

