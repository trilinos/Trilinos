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
#include "all_allo_const.h"
#include "comm_const.h"
#include "parmetis_jostle_const.h"
#include "params_const.h"
#include "timer_const.h"


/**********  parameters structure for parmetis methods **********/
static PARAM_VARS Parmetis_params[] = {
        { "PARMETIS_METHOD", NULL, "STRING" },
        { "PARMETIS_OUTPUT_LEVEL", NULL, "INT" },
        { "PARMETIS_COARSE_ALG", NULL, "INT" },
        { "PARMETIS_FOLD", NULL, "INT" },
        { NULL, NULL, NULL } };

/**********  parameters structure for jostle methods **********/
static PARAM_VARS Jostle_params[] = {
        { "JOSTLE_OUTPUT_LEVEL", NULL, "INT" },
        { "JOSTLE_THRESHOLD", NULL, "INT" },
        { "JOSTLE_GATHER_THRESHOLD", NULL, "INT" },
        { "JOSTLE_MATCHING", NULL, "STRING" },
        { "JOSTLE_REDUCTION", NULL, "STRING" },
        { "JOSTLE_CONNECT", NULL, "STRING" },
        { "JOSTLE_SCATTER", NULL, "STRING" },
        { NULL, NULL, NULL } };

/**********  parameters structure used by both ParMetis and Jostle **********/
static PARAM_VARS Graph_params[] = {
        { "CHECK_GRAPH", NULL, "INT" },
        { "SCATTER_GRAPH", NULL, "INT" },
        { NULL, NULL, NULL } };

/***************  prototypes for internal functions ********************/

#if (defined(LB_JOSTLE) || defined(LB_PARMETIS))

static int LB_ParMetis_Jostle(LB *lb, int *num_imp, LB_GID **imp_gids,
  LB_LID **imp_lids, int **imp_procs, int *num_exp, LB_GID **exp_gids,
  LB_LID **exp_lids, int **exp_procs, char *alg, int  *options);
static int hash_lookup (struct LB_hash_node **hashtab, LB_GID key, int n);

#endif  /* (defined(LB_JOSTLE) || defined(LB_PARMETIS)) */

/**********************************************************/
/* Interface routine for ParMetis. This is just a simple  */
/* wrapper that sets the options etc and then calls       */
/* LB_ParMetis_Jostle, where the real action is.          */
/**********************************************************/

int LB_ParMetis(
  LB *lb,             /* load balancing structure */
  int *num_imp,       /* number of objects to be imported */
  LB_GID **imp_gids,  /* global ids of objects to be imported */
  LB_LID **imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,    /* list of processors to import from */
  int *num_exp,       /* number of objects to be exported */
  LB_GID **exp_gids,  /* global ids of objects to be exported */
  LB_LID **exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs     /* list of processors to export to */
)
{
#ifndef LB_PARMETIS
  fprintf(stderr, "Zoltan error: ParMetis requested but not compiled into library.\n");
  return LB_FATAL;

#else /* LB_PARMETIS */
  int  i; 
  int  options[MAX_OPTIONS];
  char alg[MAX_PARAM_STRING_LEN+1];

  /* Set parameters */
  strcpy(alg, "REPARTGDIFFUSION");
  for (i=0; i<MAX_OPTIONS; i++)
    options[i] = -1;

  LB_Bind_Param(Parmetis_params, "PARMETIS_METHOD",     
               (void *) alg);
  LB_Bind_Param(Parmetis_params, "PARMETIS_COARSE_ALG", 
                (void *) &(options[OPTION_IPART]));
  LB_Bind_Param(Parmetis_params, "PARMETIS_FOLD",       
                (void *) &(options[OPTION_FOLDF]));
  LB_Bind_Param(Parmetis_params, "PARMETIS_OUTPUT_LEVEL", 
                (void *) &(options[OPTION_DBGLVL]));

  LB_Assign_Param_Vals(lb->Params, Parmetis_params, lb->Debug_Level, lb->Proc,
                       lb->Debug_Proc);

  /* Set options[0] to 1 if any of the low level ParMetis options were set,
     or 0 otherwise. This is required by ParMetis. */
  for (i=1; i<MAX_OPTIONS; i++)
    if (options[i]>0) options[0] = 1;
  if (options[0] == -1) 
    options[0] = 0;
  else {
    /* If one option was set, fill in the others. We need to do this
     * because ParMetis requires all or nothing! 
     * The default values below are consistent with ParMetis 2.0 
     */
    if (options[OPTION_IPART] == -1)
      options[OPTION_IPART] = 2; 
    if (options[OPTION_FOLDF] == -1)
      options[OPTION_FOLDF] = 150; 
    if (options[OPTION_DBGLVL] == -1)
      options[OPTION_DBGLVL] = 0; 
  }

  /* Call the real ParMetis interface */
  return LB_ParMetis_Jostle( lb, num_imp, imp_gids, imp_lids,
            imp_procs, num_exp, exp_gids, exp_lids, exp_procs,
            alg, options);

#endif /* LB_PARMETIS */
}


/**********************************************************/
/* Interface routine for Jostle. This is just a simple    */
/* wrapper that sets the options etc and then calls       */
/* LB_ParMetis_Jostle, where the real action is.          */
/**********************************************************/

int LB_Jostle(
  LB *lb,             /* load balancing structure */
  int *num_imp,       /* number of objects to be imported */
  LB_GID **imp_gids,  /* global ids of objects to be imported */
  LB_LID **imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,    /* list of processors to import from */
  int *num_exp,       /* number of objects to be exported */
  LB_GID **exp_gids,  /* global ids of objects to be exported */
  LB_LID **exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs     /* list of processors to export to */
)
{
#ifndef LB_JOSTLE
  fprintf(stderr, "Zoltan error: Jostle requested but not compiled into library.\n");
  return LB_FATAL;

#else /* LB_JOSTLE */
  static LB *lb_prev = NULL; /* Last lb structure used */
  static char *alg = "JOSTLE";
  char *cptr;
  char str[MAX_PARAM_STRING_LEN+1]; 
  char matching[MAX_PARAM_STRING_LEN+1];
  char reduction[MAX_PARAM_STRING_LEN+1];
  char connect[MAX_PARAM_STRING_LEN+1];
  char scatter[MAX_PARAM_STRING_LEN+1];
  int  i, output_level, threshold, gather_threshold; 
  int num_proc = lb->Num_Proc;     /* Temporary variables whose addresses are */
  int proc = lb->Proc;             /* passed to Jostle. We don't              */
  MPI_Comm comm = lb->Communicator;/* want to risk letting external packages  */
                                   /* change our lb struct.                   */

  /* Initialize Jostle if this is the first call with 
   * this load balancing structure.
   */
  if (lb != lb_prev){
     lb_prev = lb;
     pjostle_init(&num_proc, &proc);
     pjostle_comm(&comm);
  }

  /* Set parameters */
  output_level = 0;
  threshold = 0;
  gather_threshold = 0;
  matching[0] = '\0';
  reduction[0] = '\0';
  connect[0] = '\0';
  scatter[0] = '\0';
  LB_Bind_Param(Jostle_params, "JOSTLE_OUTPUT_LEVEL", 
                (void *) &output_level);
  LB_Bind_Param(Jostle_params, "JOSTLE_THRESHOLD",    
                (void *) &threshold);
  LB_Bind_Param(Jostle_params, "JOSTLE_GATHER_THRESHOLD", 
                (void *) &gather_threshold);
  LB_Bind_Param(Jostle_params, "JOSTLE_MATCHING",     
                (void *) matching);
  LB_Bind_Param(Jostle_params, "JOSTLE_REDUCTION",    
                (void *) reduction);
  LB_Bind_Param(Jostle_params, "JOSTLE_CONNECT",    
                (void *) connect);
  LB_Bind_Param(Jostle_params, "JOSTLE_SCATTER",    
                (void *) scatter);

  LB_Assign_Param_Vals(lb->Params, Jostle_params, lb->Debug_Level, lb->Proc,
                       lb->Debug_Proc); 

  /* Set Jostle parameters using jostle_env() */
  if (threshold){
    sprintf(str, "threshold = %d\0", threshold);
    jostle_env(str);
  }
  if (gather_threshold){
    sprintf(str, "gather threshold = %d\0", gather_threshold);
    jostle_env(str);
  }
  if (matching[0]){
    sprintf(str, "matching = %s\0", matching);
    /* Convert to lower case */
    for (cptr=str; *cptr; cptr++){
      *cptr = tolower(*cptr);
    }
    jostle_env(str);
  }
  if (reduction[0]){
    sprintf(str, "reduction = %s\0", reduction);
    /* Convert to lower case */
    for (cptr=str; *cptr; cptr++){
      *cptr = tolower(*cptr);
    }
    jostle_env(str);
  }
  if (connect[0]){
    sprintf(str, "connect = %s\0", connect);
    /* Convert to lower case */
    for (cptr=str; *cptr; cptr++){
      *cptr = tolower(*cptr);
    }
    jostle_env(str);
  }
  if (scatter[0]){
    sprintf(str, "scatter = %s\0", scatter);
    /* Convert to lower case */
    for (cptr=str; *cptr; cptr++){
      *cptr = tolower(*cptr);
    }
    jostle_env(str);
  }

  /* Set imbalance tolerance */
  sprintf(str, "imbalance = %3d ", (int)(100*(lb->Imbalance_Tol - 1)));
  jostle_env(str);

  /* Multidimensional vertex weights */
  if (lb->Obj_Weight_Dim > 1){
    sprintf(str, "ntypes = %3d ", lb->Obj_Weight_Dim);
    jostle_env(str);
  }

  /* Call the real Jostle/ParMetis interface */
  return LB_ParMetis_Jostle( lb, num_imp, imp_gids, imp_lids,
            imp_procs, num_exp, exp_gids, exp_lids, exp_procs,
            alg, &output_level);

#endif /* LB_JOSTLE */
}

/****************************************************************/
/* This routine constructs the required graph data structures   */
/* and calls ParMetis or Jostle.                                */
/*                                                              */
/* Author: Erik Boman, eboman@cs.sandia.gov (9226)              */
/****************************************************************/

#if (defined(LB_JOSTLE) || defined(LB_PARMETIS))

/* Macro to free all allocated memory */
#define FREE_MY_MEMORY \
  { \
  fprintf(stderr, "[%1d] Error on line %d in %s\n", \
          lb->Proc, __LINE__, __FILE__); \
  LB_FREE(&vtxdist); LB_FREE(&xadj); LB_FREE(&adjncy); \
  LB_FREE(&vwgt); LB_FREE(&adjwgt); LB_FREE(&part); \
  LB_FREE(&float_vwgt); LB_FREE(&xyz); \
  LB_FREE(&sendbuf); LB_FREE(&recvbuf); \
  LB_FREE(&hash_nodes); LB_FREE(&hashtab); \
  LB_FREE(&nbors_proc); LB_FREE(&nbors_global); \
  LB_FREE(&local_ids); LB_FREE(&global_ids); \
  LB_FREE(&proc_list); LB_FREE(&plist); \
  }

static int LB_ParMetis_Jostle(
  LB *lb,             /* load balancing structure */
  int *num_imp,       /* number of objects to be imported */
  LB_GID **imp_gids,  /* global ids of objects to be imported */
  LB_LID **imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,    /* list of processors to import from */
  int *num_exp,       /* number of objects to be exported */
  LB_GID **exp_gids,  /* global ids of objects to be exported */
  LB_LID **exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs,    /* list of processors to export to */
  char *alg,          /* algorithm to use */
  int  *options       /* option array */
)
{
  static char *yo = "LB_ParMetis_Jostle";
  int i, j, jj, ierr, packet_size, offset, tmp, flag, ndims; 
  int obj_wgt_dim, comm_wgt_dim, check_graph, scatter;
  int num_obj, nedges, num_edges, cross_edges, max_edges, edgecut;
  int *nbors_proc, *plist;
  int nsend, nrecv, wgtflag, numflag, num_border, max_proc_list_len;
  int get_graph_data, get_geom_data, get_times; 
  idxtype *vtxdist, *xadj, *adjncy, *vwgt, *adjwgt, *ewgt, *part, *part2;
  int nonint_wgt, tmp_num_obj;
  float *float_vwgt, *xyz, scale, sum_wgt, sum_wgt_local, *imb_tols; 
  double geom_vec[6];
  struct LB_edge_info  *proc_list, *ptr;
  struct LB_hash_node **hashtab, *hash_nodes;
  LB_LID *local_ids;
  LB_GID *global_ids, *nbors_global;
  char *sendbuf, *recvbuf;
  struct Comm_Obj *comm_plan;
  double times[5];
  int num_proc = lb->Num_Proc;     /* Temporary variables whose addresses are */
                                   /* passed to Jostle and ParMETIS.  Don't   */
  MPI_Comm comm = lb->Communicator;/* want to risk letting external packages  */
                                   /* change our lb struct.                   */
  int i99;                         /* Variables used for debugging.           */
#ifdef LB_JOSTLE
  int nnodes;
  int network[4] = {0, 1, 1, 1};
#endif

  LB_TRACE_ENTER(lb, yo);

  /* Set default return values (in case of early exit) */
  *num_exp = 0;
  *num_imp = -1; /* No import data */

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  nbors_proc = NULL;
  vtxdist = xadj = adjncy = vwgt = adjwgt = part = NULL;
  local_ids = global_ids = NULL;
  float_vwgt = xyz = imb_tols = NULL;
  ptr = proc_list = NULL;
  hashtab = NULL;
  hash_nodes = NULL;
  sendbuf = recvbuf = NULL;
  plist = NULL;

  /* Check weight dimensions */
  if (lb->Obj_Weight_Dim<0){
    fprintf(stderr, "ZOLTAN warning: Object weight dimension is %d, "
            "but should be >= 0. Using Obj_Weight_Dim = 0.\n",
            lb->Obj_Weight_Dim);
    obj_wgt_dim = 0;
  }
  else {
    obj_wgt_dim = lb->Obj_Weight_Dim;
  }
  if (lb->Comm_Weight_Dim<0){
    fprintf(stderr, "ZOLTAN warning: Communication weight dimension is %d, "
            "but should be >= 0. Using Comm_Weight_Dim = 0.\n",
            lb->Comm_Weight_Dim);
    comm_wgt_dim = 0;
  }
  else if (lb->Comm_Weight_Dim>1){
    fprintf(stderr, "ZOLTAN warning: This method does not support "
        "multidimensional communication weights. Using Comm_Weight_Dim = 1.\n");
    comm_wgt_dim = 1;
  }
  else {
    comm_wgt_dim = lb->Comm_Weight_Dim;
  }

  if (lb->Debug_Level >= LB_DEBUG_ALL) {
    printf("[%1d] Debug: alg=%s, Obj_Weight_Dim=%d, Comm_Weight_Dim=%d\n", 
      lb->Proc, alg, obj_wgt_dim, comm_wgt_dim);
    printf("[%1d] Debug: ParMetis options = %d, %d, %d, %d\n", lb->Proc,
      options[0], options[1], options[2], options[3]);
  }

  /* Start timer */
  get_times = (lb->Debug_Level >= LB_DEBUG_ATIME);
  if (get_times){
    MPI_Barrier(lb->Communicator);
    times[0] = LB_Time();
  }

  /* Get parameter options shared by ParMetis and Jostle */
  check_graph = 1;          /* default */
  scatter = 1;              /* default */
  LB_Bind_Param(Graph_params, "CHECK_GRAPH", (void *) &check_graph);
  LB_Bind_Param(Graph_params, "SCATTER_GRAPH", (void *) &scatter);
  LB_Assign_Param_Vals(lb->Params, Graph_params, lb->Debug_Level, lb->Proc,
                       lb->Debug_Proc);

  /* Most ParMetis methods use only graph data */
  get_graph_data = 1;
  get_geom_data = 0;
  ndims = 0;

  /* Some algorithms use geometry data */
  if (strcmp(alg, "PARTGEOMKWAY") == 0){
    get_graph_data = 1;
    get_geom_data = 1;
  }
  else if (strcmp(alg, "PARTGEOM") == 0){
    get_graph_data = 0;
    get_geom_data = 1;
  }

  /* Set up ParMETIS data structures */
  num_obj = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr){
    /* Return error code */
    FREE_MY_MEMORY;
    LB_TRACE_EXIT(lb, yo);
    return (ierr);
  }
  
  if (lb->Debug_Level >= LB_DEBUG_ALL)
    printf("[%1d] Debug: num_obj =%d\n", lb->Proc, num_obj);

  
  vtxdist = (idxtype *)LB_MALLOC((lb->Num_Proc+1)* sizeof(idxtype));
  if (num_obj>0){
    global_ids = (LB_GID *) LB_MALLOC(num_obj * sizeof(LB_GID) );
    local_ids = (LB_LID *) LB_MALLOC(num_obj * sizeof(LB_LID) );
    if (obj_wgt_dim)
      float_vwgt = (float *)LB_MALLOC(obj_wgt_dim*num_obj * sizeof(float));
    else {
      float_vwgt = NULL;
      vwgt = NULL;
    }
    if (!vtxdist || !global_ids || !local_ids || (obj_wgt_dim && !float_vwgt)){
      /* Not enough memory */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    LB_Get_Obj_List(lb, global_ids, local_ids, obj_wgt_dim, float_vwgt, &ierr);
    if (ierr){
      /* Return error */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_FATAL;
    }
  
    if (lb->Debug_Level >= LB_DEBUG_ALL) {
      printf("[%1d] Debug: Global ids = ", lb->Proc);
      for (i99=0; i99<num_obj; i99++) printf("%d ", global_ids[i99]);
      printf("\n");
    }
  }
  
  /* Construct vtxdist[i] = the number of objects on all procs < i. */
  /* Scan to compute partial sums of the number of objs */
  MPI_Scan (&num_obj, vtxdist, 1, IDX_DATATYPE, MPI_SUM, lb->Communicator);
  /* Gather data from all procs */
  MPI_Allgather (&vtxdist[0], 1, IDX_DATATYPE, 
                 &vtxdist[1], 1, IDX_DATATYPE, lb->Communicator);
  vtxdist[0] = 0;
  
  if (lb->Debug_Level >= LB_DEBUG_ALL) {
    printf("[%1d] Debug: vtxdist = ", lb->Proc);
    for (i99=0; i99<=lb->Num_Proc; i99++)
      printf("%d ", vtxdist[i99]);
    printf("\n");
  }
  
  if (check_graph >= 1){
     if ((lb->Proc ==0) && (vtxdist[lb->Num_Proc]==0))
        fprintf(stderr, "Zoltan warning: No objects to balance in %s\n", yo);
  }

  if (get_graph_data){

    num_edges = 0;
    max_edges = 0;
    for (i=0; i< num_obj; i++){
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], &ierr);
      if (ierr){
        /* Return error */
        FREE_MY_MEMORY;
        LB_TRACE_EXIT(lb, yo);
        return (ierr);
      }
      num_edges += nedges;
      if (nedges>max_edges) max_edges = nedges;
    }
    if (lb->Debug_Level >= LB_DEBUG_ALL)
      printf("[%1d] Debug: num_edges = %d\n", lb->Proc, num_edges);

    if (check_graph >= 1){
       ierr = MPI_Reduce(&nedges, &tmp, 1, MPI_INT, MPI_SUM, 0, lb->Communicator);
       if ((lb->Proc ==0) && (tmp==0))
          fprintf(stderr, "Zoltan warning: No edges in graph.\n");
    }

  
    /* Allocate space for ParMETIS data structs */
    xadj   = (idxtype *)LB_MALLOC((num_obj+1) * sizeof(idxtype));
    adjncy = (idxtype *)LB_MALLOC(num_edges * sizeof(idxtype));
    if (comm_wgt_dim) 
      adjwgt = (idxtype *)LB_MALLOC(comm_wgt_dim 
                * num_edges * sizeof(idxtype));
  
    if (!xadj || (num_edges && !adjncy) || (num_edges && comm_wgt_dim && !adjwgt)){
      /* Not enough memory */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    if (lb->Debug_Level >= LB_DEBUG_ALL)
      printf("[%1d] Debug: Successfully allocated ParMetis space\n", lb->Proc);
  
    /* Construct ParMETIS graph */
    /* First compute a global dense numbering of the objects/vertices */
  
    /* Construct local hash table */
    hash_nodes = (struct LB_hash_node *)LB_MALLOC(num_obj *
      sizeof(struct LB_hash_node));
    hashtab = (struct LB_hash_node **) LB_MALLOC(num_obj *
      sizeof(struct LB_hash_node *) );
    if (num_obj && ((!hash_nodes) || (!hashtab))){
      /* Not enough memory */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    
    for (i=0; i< num_obj; i++){
      hashtab[i] = NULL;
      hash_nodes[i].gid = global_ids[i];
      hash_nodes[i].gno = vtxdist[lb->Proc]+i;
    }
  
    for (i=0; i< num_obj; i++){
      /* insert hashed elements into hash table */
      j = LB_Hash(global_ids[i], num_obj);
      hash_nodes[i].next = hashtab[j];
      hashtab[j] = &hash_nodes[i];
    }
  
    /* Estimate the number of inter-proc edges. 
     * First estimate the number of border objects 
     * based on a 2d regular grid. (We could alternatively
     * use Get_Num_Border_Obj but this could be expensive
     * and the function may not be reqistered. )
     */
    num_border = 4*sqrt((double) num_obj);
    if (num_border > num_obj) num_border = num_obj;
     
    /* Assume that the edges are approx. evenly distributed among the objs. */
    if (num_obj>0){
       max_proc_list_len = num_edges * num_border / num_obj;
       if (max_proc_list_len < CHUNKSIZE)
          max_proc_list_len = CHUNKSIZE;
    }
    else
       max_proc_list_len = 0;
    
    /* Allocate edge list data */
    nbors_global = (LB_GID *)LB_MALLOC(max_edges * sizeof(LB_GID));
    nbors_proc = (int *)LB_MALLOC(max_edges * sizeof(int));
    plist = (int *)LB_MALLOC(lb->Num_Proc * sizeof(int));

    if ((max_edges && ((!nbors_global) || (!nbors_proc))) || 
        (!plist)){
      /* Not enough memory */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    for (i=0; i<lb->Num_Proc; i++)
      plist[i] = -1;
  
    if (max_proc_list_len){
      /* Allocate space for processor list */
      while ((proc_list==NULL) && (max_proc_list_len>=CHUNKSIZE)){
        proc_list = (struct LB_edge_info *) LB_MALLOC(max_proc_list_len *
          sizeof(struct LB_edge_info) );
        if (!proc_list){
          /* Not enough memory, try shorter list */
          if (lb->Debug_Level >= LB_DEBUG_ALL) {
            printf("[%1d] Debug: Could not allocate %d list nodes, "
                   "trying %d instead.\n", lb->Proc,
                   max_proc_list_len, max_proc_list_len/2);
          }
          max_proc_list_len /= 2;
        }
      }
      if (!proc_list){
        /* Not enough memory */
        FREE_MY_MEMORY;
        LB_TRACE_EXIT(lb, yo);
        return LB_MEMERR;
      }
    }

    /* proc_list[i] will contain a struct with data to send
     * to another processor. We don't know yet the total number
     * of inter-proc edges so we may have to adjust the size of 
     * proc_list through REALLOC.  We add a chunk of space at a time.
     * The motivation for this design is to reduce the number of calls
     * to REALLOC.  We could have reduced it into a single MALLOC 
     * if we were willing to call the query function get_edge_list 
     * twice for every object, but that is probably
     * even more expensive.
     *
     * Note that proc_list may contain duplicate elements.
     * These are eliminated before the communication phase.
     */
    
    nsend = 0;      /* Number of objects we need to send info about */
    offset = 0;     /* Index of next available node in proc_list */
    ptr = proc_list;
    ewgt = NULL;
    xadj[0] = 0;
    jj = 0;
  
    for (i=0; i< num_obj; i++){
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], &ierr);
      xadj[i+1] = xadj[i] + nedges;
      if (comm_wgt_dim)
          ewgt = &adjwgt[jj*comm_wgt_dim];
      lb->Get_Edge_List(lb->Get_Edge_List_Data, global_ids[i], local_ids[i],
          nbors_global, nbors_proc, comm_wgt_dim, ewgt, &ierr);
      if (ierr){
        /* Return error */
        FREE_MY_MEMORY;
        LB_TRACE_EXIT(lb, yo);
        return (ierr);
      }
  
      if (lb->Debug_Level >= LB_DEBUG_ALL) {
        printf("[%1d] Debug: i=%d, gid=%d, lid=%d, nedges=%d\n", lb->Proc, i, 
          global_ids[i], local_ids[i], nedges);
      }

      /* Separate inter-processor edges from the local ones */
      for (j=0; j<nedges; j++, jj++){
        if (nbors_proc[j] == lb->Proc){
          /* local edge */
          adjncy[jj] = hash_lookup(hashtab, nbors_global[j], num_obj);
        } else {
          /* Inter-processor edge. */
          /* Check if we already have gid[i] in proc_list with */
          /* the same destination.                             */
          flag = 0;
          if (plist[nbors_proc[j]] < i){
            /* We need to send info about this edge */
            flag = 1;
            nsend++; 
            plist[nbors_proc[j]] = i;
          }

          /* Check if we need to allocate more space for proc_list.*/
          if (offset == max_proc_list_len){
            if (lb->Debug_Level >= LB_DEBUG_ALL)
              printf("[%1d] Debug: Allocating more list space, "
                     "max_proc_list_len = %d, increasing by %d\n", 
                     lb->Proc, max_proc_list_len, CHUNKSIZE);

            max_proc_list_len += CHUNKSIZE;
            proc_list = (struct LB_edge_info *) LB_REALLOC(proc_list,
                         max_proc_list_len*sizeof(struct LB_edge_info));
            if (!proc_list){
              /* Not enough memory */
              FREE_MY_MEMORY;
              LB_TRACE_EXIT(lb, yo);
              return LB_MEMERR;
            }
          }
          ptr = &proc_list[offset];
          ptr->my_gid = global_ids[i];
          ptr->my_gno = hash_lookup(hashtab, global_ids[i], num_obj);
          ptr->nbor_gid = nbors_global[j];
          if (flag)
            ptr->nbor_proc = nbors_proc[j];
          else
            ptr->nbor_proc = -1;
          ptr->adj = &adjncy[jj];

          if (lb->Debug_Level >= LB_DEBUG_ALL)
            printf("[%1d] Debug: proc_list[%1d] my_gid=%d, my_gno=%d, "
                   "nbor_proc=%d\n",
                   lb->Proc, offset, ptr->my_gid, ptr->my_gno, ptr->nbor_proc);

          adjncy[jj] = -1; /* We need to come back here later */
          offset++;
        }
      }
    }
    cross_edges = offset;

    LB_FREE(&plist);
    LB_FREE(&nbors_global);
    LB_FREE(&nbors_proc);

    /* Sanity check */
    if ((check_graph >= 1) && (jj != xadj[num_obj])){
      fprintf(stderr, "[%1d] ZOLTAN ERROR: Internal error in %s. "
              "Something may be wrong with the edges in the graph\n", 
               lb->Proc, yo);
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_FATAL;
    }
  
    /* Exchange info between processors to resolve global number 
     * for objects that are off-proc.
     */

    /* Allocate send buffer */
    packet_size = sizeof(LB_GID) + sizeof(int);
    sendbuf = (char *) LB_MALLOC(nsend * packet_size);
    plist = (int *) LB_MALLOC(nsend * sizeof(int));

    if (nsend && (!sendbuf || !plist) ){
      /* Not enough space */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }

    /* Pack the data to send */
    if (lb->Debug_Level >= LB_DEBUG_ALL)
      printf("[%1d] Debug: %d messages to send.\n", lb->Proc, nsend);

    offset = 0;
    j = 0;
    for (i=0, ptr=proc_list; i<cross_edges; i++, ptr++){
      if (ptr->nbor_proc >= 0){
        if (lb->Debug_Level >= LB_DEBUG_ALL)
          printf("[%1d] Debug: Sending (%d,%d) to proc %d\n", lb->Proc, 
            ptr->my_gid, ptr->my_gno, ptr->nbor_proc);
        memcpy(&sendbuf[offset], (char *) &(ptr->my_gid), sizeof(LB_GID)); 
        offset += sizeof(LB_GID);
        memcpy(&sendbuf[offset], (char *) &(ptr->my_gno), sizeof(int)); 
        offset += sizeof(int);
        plist[j++] = ptr->nbor_proc;
      }
    }

    /* Create the communication plan */
    if (lb->Debug_Level >= LB_DEBUG_ALL)
      printf("[%1d] Debug: Calling LB_Comm_Create with %d packets to send.\n",
             lb->Proc, nsend);

    ierr = LB_Comm_Create( &comm_plan, nsend, plist, comm, TAG1, 
                           lb->Deterministic, &nrecv);
    if (ierr){
      /* Return error code */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return (ierr);
    }

    /* Allocate recv buffer */
    recvbuf = (char *) LB_MALLOC(nrecv * packet_size);
    if (nrecv && (!sendbuf || !plist) ){
      /* Not enough space */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    if (lb->Debug_Level >= LB_DEBUG_ALL)
      printf("[%1d] Debug: Ready to receive %d packets.\n", 
        lb->Proc, nrecv);

    /* Do the communication */
    ierr = LB_Comm_Do( comm_plan, TAG2, sendbuf, packet_size, recvbuf);
    if (ierr){
      /* Return error code */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return (ierr);
    }

    /* Destroy the comm. plan */
    LB_Comm_Destroy( &comm_plan);
  
    /* Unpack data into the ParMETIS/Jostle struct.
     * Resolve off-proc global numbers by:
     * 1) Insert global ids from recvbuf into hash table;
     * 2) Look up missing references in proc_list 
     *    using the hash table.
     */

    /* Change hash table to contain only border objects */
    /* that we received from other procs.               */
    hash_nodes = (struct LB_hash_node *)LB_REALLOC(hash_nodes,
      nrecv * sizeof(struct LB_hash_node));
    hashtab = (struct LB_hash_node **) LB_REALLOC(hashtab,
      nrecv * sizeof(struct LB_hash_node *) );
    if (nrecv && ((!hash_nodes) || (!hashtab))){
      /* Not enough memory */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    
    /* Copy data from recvbuf into hash table nodes */
    for (i=0; i< nrecv; i++){
      hashtab[i] = NULL;
      hash_nodes[i].gid = *((LB_GID *)&recvbuf[i*packet_size]);
      hash_nodes[i].gno = *((int *)&recvbuf[i*packet_size+sizeof(LB_GID)]);
      /* Do we need to pad for byte alignment? */
    }
  
    /* Insert nodes into hash table */
    for (i=0; i< nrecv; i++){
      j = LB_Hash(hash_nodes[i].gid, nrecv);
      hash_nodes[i].next = hashtab[j];
      hashtab[j] = &hash_nodes[i];
      if (lb->Debug_Level >= LB_DEBUG_ALL)
        printf("[%1d] Debug: Hashed GID %d to %d, gno = %d\n",
               lb->Proc, hash_nodes[i].gid, j, hash_nodes[i].gno);
    }

    for (i=0; i<cross_edges; i++){
      /* Look up unresolved global_ids */
      if ((tmp=hash_lookup(hashtab, proc_list[i].nbor_gid, nrecv)) <0){
        /* Error. This should never happen! */
        fprintf(stderr, "[%1d] ZOLTAN ERROR: Internal error in %s. "
               "Off-proc global ID %d is not in hash table.\n", 
               lb->Proc, yo, proc_list[i].nbor_gid);
        FREE_MY_MEMORY;
        LB_TRACE_EXIT(lb, yo);
        return LB_FATAL;
      }
      else{
        /* Insert the global number into adjncy vector */
        *(proc_list[i].adj) = tmp;
        if (lb->Debug_Level >= LB_DEBUG_ALL)
          printf("[%1d] Debug: GID %d has global number %d\n",
            lb->Proc, proc_list[i].nbor_gid, tmp);
      }
      
    }

    /* Free space */
    LB_FREE(&sendbuf);
    LB_FREE(&recvbuf);
    LB_FREE(&proc_list);
    LB_FREE(&plist);
    LB_FREE(&hash_nodes);
    LB_FREE(&hashtab);
  
    /* Get vertex weights if needed */
    if (obj_wgt_dim){
      vwgt = (idxtype *)LB_MALLOC(obj_wgt_dim*num_obj
                          * sizeof(idxtype));

      /* Compute local sum of the obj weights (over all dimensions) */
      /* Check if all weights are integers */
      flag = 0;
      sum_wgt_local = 0;
      for (i=0; i<num_obj*obj_wgt_dim; i++){
        if (!flag){ 
          tmp = float_vwgt[i]; /* Converts to nearest int */
          if (fabs((double)tmp-float_vwgt[i]) > .001){
            flag = 1;
          }
        }
        sum_wgt_local += float_vwgt[i];
      }
      /* Compute global sum of the obj weights */
      MPI_Allreduce(&flag, &nonint_wgt, 1, 
          MPI_INT, MPI_LOR, lb->Communicator);
      MPI_Allreduce(&sum_wgt_local, &sum_wgt, 1, 
          MPI_FLOAT, MPI_SUM, lb->Communicator);

      if (sum_wgt == 0){
        fprintf(stderr, "ZOLTAN ERROR: All object weights are zero.\n");
        FREE_MY_MEMORY;
        LB_TRACE_EXIT(lb, yo);
        return LB_FATAL;
      }

      /* Convert weights to integers between 1 and MAX_WGT_SUM/sum_wgt */
      scale = 1;
      /* Scale unless all weights are integers */
      if (nonint_wgt || (sum_wgt > MAX_WGT_SUM)){
        scale = MAX_WGT_SUM/sum_wgt;
      }
      if (lb->Debug_Level >= LB_DEBUG_ALL)
        printf("[%1d] Debug: Converting vertex (object) weights, nonint_wgt = %d "
               "sum_wgt = %g, scale = %g\n", lb->Proc, nonint_wgt, sum_wgt, scale);

      for (i=0; i<num_obj*obj_wgt_dim; i++){
        if (scale==1)
           vwgt[i] = (int) float_vwgt[i];
        else
           vwgt[i] = (int) ceil(float_vwgt[i]*scale);
      }
      LB_FREE(&float_vwgt);
    }

  } /* end get_graph_data */

  if (get_geom_data){
    /* Determine how many dimensions the data have */
    ndims = lb->Get_Num_Geom(lb->Get_Num_Geom_Data, &ierr);
    if (ierr){
      /* Return error */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return (ierr);
    }
    /* Allocate space for the geometry data */
    xyz = (float *) LB_MALLOC(ndims*num_obj * sizeof(float));
    if (ndims && num_obj && !xyz){
      /* Not enough space */
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    /* Get the geometry data */
    for (i=0; i<num_obj; i++){
      lb->Get_Geom(lb->Get_Geom_Data, global_ids[i], local_ids[i], 
        geom_vec, &ierr);
      if (ierr) {
        /* Return error code */
        FREE_MY_MEMORY;
        LB_TRACE_EXIT(lb, yo);
        return (ierr);
      }
      for (j=0; j<ndims; j++)
        xyz[i*ndims+j] = geom_vec[j];
    }
  }

  tmp_num_obj = num_obj;
  if ((scatter>0) && (scatter<3)){
    /* Decide if the data imbalance is so bad that we should scatter the graph. */
    /* scatter==1: Scatter if all the objects are on a single processor.        */
    /* scatter==2: Scatter if any processor has no objects.                     */
    if (num_obj==0)
      j = 1;
    else 
      j = 0;
    MPI_Allreduce(&j, &tmp, 1, MPI_INT, MPI_SUM, lb->Communicator);
    if (scatter == 1){
      if (tmp < lb->Num_Proc-1)
        scatter = 0;
    }
    else if (scatter==2){
      if (tmp == 0)
        scatter = 0;
    }
  }

  if (scatter){
    ierr = LB_scatter_graph(get_graph_data, &vtxdist, &xadj, &adjncy, &vwgt, &adjwgt, &xyz, ndims, 
              lb, &comm_plan);
    if ((ierr == LB_FATAL) || (ierr == LB_MEMERR)){
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return ierr;
    }
    tmp_num_obj = vtxdist[lb->Proc+1]-vtxdist[lb->Proc];
  }

  /* Get ready to call ParMETIS or Jostle */
  edgecut = -1; 
  wgtflag = 2*(obj_wgt_dim>0) + (comm_wgt_dim>0); /* Multidim wgts not supported yet */
  numflag = 0;
  part = (idxtype *)LB_MALLOC((tmp_num_obj+1) * sizeof(idxtype));
  if (!part){
    /* Not enough memory */
    FREE_MY_MEMORY;
    LB_TRACE_EXIT(lb, yo);
    return LB_MEMERR;
  }
  if (obj_wgt_dim>0){
    /* Set Imbalance Tolerance for each component. For now, they are all the same. */
    imb_tols = (float *) LB_MALLOC(obj_wgt_dim * sizeof(float));
    for (i=0; i<obj_wgt_dim; i++)
      imb_tols[i] = lb->Imbalance_Tol;
  }

  /* Verify that graph is correct */
  if (get_graph_data){
     ierr = LB_verify_graph(lb->Communicator, vtxdist, xadj, adjncy, vwgt, 
               adjwgt, obj_wgt_dim, comm_wgt_dim, check_graph);
  }
  
  /* Get a time here */
  if (get_times) times[1] = LB_Time();

  /* Select the desired ParMetis or Jostle function */

  if (strcmp(alg, "JOSTLE") == 0){
#ifdef LB_JOSTLE
    offset = 0;            /* Index of the first object/node. */
    j = 0;                 /* Dummy variable for Jostle */
    nnodes = vtxdist[lb->Num_Proc]; /* Global number of objects */ 
    ndims = 1;             /* Topological dimension of the computer */
    network[1] = lb->Num_Proc;
    /* Convert xadj array to Jostle format (degree) */
    for (i=0; i<tmp_num_obj; i++){
      xadj[i] = xadj[i+1] - xadj[i];
    }
    /* Convert vtxdist array to Jostle format (degree) */
    for (i=0; i<lb->Num_Proc; i++){
      vtxdist[i] = vtxdist[i+1] - vtxdist[i];
    }
    LB_TRACE_DETAIL(lb, yo, "Calling the Jostle library");
    jostle_env("format = contiguous");
    if (check_graph >= 2){
      jostle_env("check = on");
    }
    pjostle(&nnodes, &offset, &(vtxdist[lb->Proc]), &j, vtxdist, 
       xadj, vwgt, part, &num_edges, adjncy, adjwgt, network,
       NULL, options, &ndims, NULL); 
    LB_TRACE_DETAIL(lb, yo, "Returned from the Jostle library");
#else
    /* We don't have Jostle */
    fprintf(stderr, "Sorry, Jostle is not available on this system.\n");
    FREE_MY_MEMORY;
    LB_TRACE_EXIT(lb, yo);
    return LB_FATAL;
#endif
  }
  else if (strcmp(alg, "PARTKWAY") == 0){
    LB_TRACE_DETAIL(lb, yo, "Calling the ParMETIS library");
    if (obj_wgt_dim <= 1)
      ParMETIS_PartKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
        &numflag, &num_proc, options, &edgecut, part, &comm);
    else /* obj_wgt_dim >= 2 */
      /* Beta version of multiconstraint ParMetis. Interface may change! */
      Moc_ParMETIS_PartKway (&obj_wgt_dim, vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
        &numflag, &num_proc, imb_tols, options, &edgecut, part, &comm, (lb->Proc +1));
    LB_TRACE_DETAIL(lb, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "PARTGEOMKWAY") == 0){
    LB_TRACE_DETAIL(lb, yo, "Calling the ParMETIS library");
    ParMETIS_PartGeomKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag,
      &numflag, &ndims, xyz, &num_proc, options, &edgecut, 
      part, &comm);
    LB_TRACE_DETAIL(lb, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "PARTGEOM") == 0){
    LB_TRACE_DETAIL(lb, yo, "Calling the ParMETIS library");
    ParMETIS_PartGeom (vtxdist, &ndims, xyz, part, &comm);
    LB_TRACE_DETAIL(lb, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REPARTLDIFFUSION") == 0){
    LB_TRACE_DETAIL(lb, yo, "Calling the ParMETIS library");
    ParMETIS_RepartLDiffusion (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    LB_TRACE_DETAIL(lb, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REPARTGDIFFUSION") == 0){
    LB_TRACE_DETAIL(lb, yo, "Calling the ParMETIS library");
    ParMETIS_RepartGDiffusion (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    LB_TRACE_DETAIL(lb, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REPARTREMAP") == 0){
    LB_TRACE_DETAIL(lb, yo, "Calling the ParMETIS library");
    if (obj_wgt_dim <= 1)
      ParMETIS_RepartRemap (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
        &numflag, options, &edgecut, part, &comm);
    else { /* obj_wgt_dim >= 2 */
      /* Beta version of multiconstraint ParMetis. Interface may change! */
      Moc_ParMETIS_SR (&obj_wgt_dim, vtxdist, xadj, adjncy, vwgt, adjwgt, NULL,
        &wgtflag, &numflag, &num_proc, imb_tols, options, &edgecut, part, &comm, (lb->Proc+1));
    }
    LB_TRACE_DETAIL(lb, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REPARTMLREMAP") == 0){
    LB_TRACE_DETAIL(lb, yo, "Calling the ParMETIS library");
    ParMETIS_RepartMLRemap (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    LB_TRACE_DETAIL(lb, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REFINEKWAY") == 0){
    LB_TRACE_DETAIL(lb, yo, "Calling the ParMETIS library");
    ParMETIS_RefineKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    LB_TRACE_DETAIL(lb, yo, "Returned from the ParMETIS library");
  }
  else {
    /* This should never happen! */
    fprintf(stderr, "ZOLTAN Error: Unknown ParMetis or Jostle algorithm %s\n", alg);
    FREE_MY_MEMORY;
    LB_TRACE_EXIT(lb, yo);
    return LB_FATAL;
  }

  /* Get a time here */
  if (get_times) times[2] = LB_Time();

  if (lb->Debug_Level >= LB_DEBUG_ALL)
    printf("[%1d] Debug: Returned from partitioner with edgecut= %d\n", 
      lb->Proc, edgecut);

  /* Free weights; they are no longer needed */
  if (obj_wgt_dim) {
    LB_FREE(&vwgt);
    LB_FREE(&imb_tols);
  }
  if (comm_wgt_dim) LB_FREE(&adjwgt);

  /* If we have been using a scattered graph, convert partition result back to original distribution */
  if (scatter){
    /* Allocate space for partition array under original distribution */
    part2 = (idxtype *) LB_MALLOC(num_obj*sizeof(idxtype)); 
    if (num_obj && !part2){
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    /* Use reverse communication to compute the partition array under the original distribution */
    ierr = LB_Comm_Do_Reverse(comm_plan, TAG3, (char *) part, sizeof(idxtype), (char *) part2);
    if ((ierr == LB_FATAL) || (ierr == LB_MEMERR)){
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return ierr;
    }
    LB_Comm_Destroy(&comm_plan); /* Destroy the comm. plan */
    LB_FREE(&part); /* We don't need the partition array with the scattered distribution any more */
    part = part2;   /* part is now the partition array under the original distribution */
  }
 
  /* Determine number of objects to export */
  nsend = 0;
  for (i=0; i<num_obj; i++)
    if (part[i] != lb->Proc) nsend++;
  (*num_exp) = nsend;

  /* Create export lists */
  if (nsend>0){
    if (!LB_Special_Malloc(lb,(void **)exp_gids,nsend,LB_SPECIAL_MALLOC_GID)) {
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    if (!LB_Special_Malloc(lb,(void **)exp_lids,nsend,LB_SPECIAL_MALLOC_LID)) {
      LB_Special_Free(lb,(void **)exp_gids,LB_SPECIAL_MALLOC_GID);
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    if (!LB_Special_Malloc(lb,(void **)exp_procs,nsend,LB_SPECIAL_MALLOC_INT)) {
      LB_Special_Free(lb,(void **)exp_lids,LB_SPECIAL_MALLOC_LID);
      LB_Special_Free(lb,(void **)exp_gids,LB_SPECIAL_MALLOC_GID);
      FREE_MY_MEMORY;
      LB_TRACE_EXIT(lb, yo);
      return LB_MEMERR;
    }
    j = 0;
    for (i=0; i<num_obj; i++){
      if (part[i] != lb->Proc){
        (*exp_gids)[j] = global_ids[i];
        (*exp_lids)[j] = local_ids[i];
        (*exp_procs)[j] = part[i];
        j++;
      }
    }
  }

  /* Free space */
  LB_FREE(&part);
  LB_FREE(&local_ids);
  LB_FREE(&global_ids);
  LB_FREE(&vtxdist);
  LB_FREE(&xadj);
  LB_FREE(&adjncy);
  LB_FREE(&xyz);

  /* Get a time here */
  if (get_times) times[3] = LB_Time();

  /* Output timing results if desired */
  if (get_times){
    if (lb->Proc==0) printf("\nZOLTAN timing statistics:\n");
    LB_Print_Stats(lb, times[1]-times[0], " Partitioner Pre-processing time  ");
    LB_Print_Stats(lb, times[2]-times[1], " Partitioner Library time         ");
    LB_Print_Stats(lb, times[3]-times[2], " Partitioner Post-processing time ");
    LB_Print_Stats(lb, times[3]-times[0], " Partitioner Total time           ");
    if (lb->Proc==0) printf("\n");
  }

  LB_TRACE_EXIT(lb, yo);
  return LB_OK;
}

#endif /* defined (LB_JOSTLE) || defined (LB_PARMETIS) */


/*********************************************************************/
/* ParMetis parameter routine                                        */
/*********************************************************************/

int LB_Set_ParMetis_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status, i;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */
    char *valid_methods[] = {
        "PARTKWAY", "PARTGEOMKWAY", "PARTGEOM", 
        "REPARTLDIFFUSION", "REPARTGDIFFUSION",
        "REPARTREMAP", "REPARTMLREMAP",
        "REFINEKWAY",
         NULL };

    status = LB_Check_Param(name, val, Parmetis_params, &result, &index);
    if (status == 1)
       status = LB_Check_Param(name, val, Graph_params, &result, &index);

    if (status == 0){
      /* OK so far, do sanity check of parameter values */

      if (strcmp(name, "PARMETIS_METHOD") == 0){
        status = 2;
        for (i=0; valid_methods[i] != NULL; i++){
          if (strcmp(val, valid_methods[i]) == 0){
            status = 0; 
            break;
          }
        }
      }
      else {
        /* All integer parameters should be non-negative */
        if (result.ival < 0)
          status = 2; 
      }
    }

    return(status);
}

/*********************************************************************/
/* Jostle parameter routine                                          */
/*********************************************************************/

int LB_Set_Jostle_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */

    status = LB_Check_Param(name, val, Jostle_params, &result, &index);
    if (status == 1)
       status = LB_Check_Param(name, val, Graph_params, &result, &index);

    if (status == 0){
      /* OK so far, do sanity check of parameter values */

      if (strcmp(name, "JOSTLE_OUTPUT_LEVEL") == 0){
        if (result.ival < 0)
          status = 2; /* output level must be non-negative */
      }
      /* Pass the other parameters to Jostle without checking the values */
    }

    return(status);
}

#if (defined(LB_JOSTLE) || defined(LB_PARMETIS))

/*******************************************************************
 * hash_lookup uses LB_Hash to lookup a key 
 *
 * Input:
 *   hashtab, pointer to the hash table
 *   key, a key to look up of type LB_GID (any data type)
 *   n,   dimension of the hash table
 *
 * Return value:
 *   the global number of the element with key key,
 *   or -1 if the key is not in the hash table
 *
 *******************************************************************/

static int hash_lookup (struct LB_hash_node **hashtab, LB_GID key, int n)
{
  int i;
  struct LB_hash_node *ptr;

  i = LB_Hash(key, n);
  for (ptr=hashtab[i]; ptr != NULL; ptr = ptr->next){
    if (LB_EQ_GID(ptr->gid, key))
      return (ptr->gno);
  }
  /* Key not in hash table */
  return -1;
}

#endif /* (defined(LB_JOSTLE) || defined(LB_PARMETIS)) */

