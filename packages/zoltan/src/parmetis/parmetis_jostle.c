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

#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "all_allo_const.h"
#include "parmetis_jostle.h"
#include "params_const.h"
#include "timer_const.h"

/* Data structures used in ParMetis interface routines */
/* An array of this data structure works with a parallel array of
 * ZOLTAN_ID_PTR called proc_list_nbor containing the global IDs of the
 * neighboring object.
 * This separate array is needed to prevent individual mallocs of
 * neighboring global IDs.
 */
struct Edge_Info {
  ZOLTAN_ID_PTR my_gid;  /* Pointer to the Global id of local vtx */
  int my_gno;        /* Global number of local vtx */
  int nbor_proc;     /* Proc id for the neighboring proc */
  int *adj;          /* Pointer to adjcny array */
};

struct Hash_Node {
  ZOLTAN_ID_PTR gid;     /* Pointer to a Global id */
  int gno;           /* Global number */
  struct Hash_Node * next;
};


/**********  parameters structure for parmetis methods **********/
static PARAM_VARS Parmetis_params[] = {
        { "PARMETIS_METHOD", NULL, "STRING" },
        { "PARMETIS_OUTPUT_LEVEL", NULL, "INT" },
        { "PARMETIS_SEED", NULL, "INT" },
        { "PARMETIS_IPC2REDIST", NULL, "FLOAT" },
        { "PARMETIS_USE_OBJ_SIZE", NULL, "INT" },
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

#if (defined(ZOLTAN_JOSTLE) || defined(ZOLTAN_PARMETIS))

static int Zoltan_ParMetis_Jostle(ZZ *zz, int *num_imp, ZOLTAN_ID_PTR *imp_gids,
  ZOLTAN_ID_PTR *imp_lids, int **imp_procs, int *num_exp, ZOLTAN_ID_PTR *exp_gids,
  ZOLTAN_ID_PTR *exp_lids, int **exp_procs, char *alg, int  *options, float *ipc2r);
static int hash_lookup (ZZ *, struct Hash_Node **, ZOLTAN_ID_PTR, int);
static int scale_round_weights(float *fwgts, idxtype *iwgts, int n, int dim, 
                 int mode, int max_wgt_sum, int debug_level, MPI_Comm comm);
#if PARMETIS_VERSION >= 3
static int pmv3method(char *alg);
#endif

#endif  /* (defined(ZOLTAN_JOSTLE) || defined(ZOLTAN_PARMETIS)) */

/**********************************************************/
/* Interface routine for ParMetis. This is just a simple  */
/* wrapper that sets the options etc and then calls       */
/* Zoltan_ParMetis_Jostle, where the real action is.          */
/**********************************************************/

int Zoltan_ParMetis(
  ZZ *zz,               /* Zoltan structure */
  int *num_imp,         /* number of objects to be imported */
  ZOLTAN_ID_PTR *imp_gids,  /* global ids of objects to be imported */
  ZOLTAN_ID_PTR *imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,      /* list of processors to import from */
  int *num_exp,         /* number of objects to be exported */
  ZOLTAN_ID_PTR *exp_gids,  /* global ids of objects to be exported */
  ZOLTAN_ID_PTR *exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs       /* list of processors to export to */
)
{
#ifndef ZOLTAN_PARMETIS
  char *yo="Zoltan_ParMetis";
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
     "ParMetis requested but not compiled into library.");
  return ZOLTAN_FATAL;

#else /* ZOLTAN_PARMETIS */
  int  i; 
  int  options[MAX_OPTIONS];
  char alg[MAX_PARAM_STRING_LEN+1];
  int output_level, seed, coarse_alg, fold, use_obj_size;
  float ipc2redist;

  /* Set parameters */
#if PARMETIS_VERSION >= 3
  strcpy(alg, "ADAPTIVEREPART");
  ipc2redist = 100.0;  /* 100 gives about the same partition quality as GDiffusion */
#else
  strcpy(alg, "REPARTGDIFFUSION");
#endif
 
  /* Always use ParMetis option array because Zoltan by default 
     produces no output (silent mode). ParMetis requires options[0]=1
     when options array is to be used. */
  options[0] = 1;
  for (i = 1; i < MAX_OPTIONS; i++)
    options[i] = 0;

  /* Set the default option values. */
  output_level = 0;
  coarse_alg = 2;
  use_obj_size = 0;
  fold = 0;
  seed = GLOBAL_SEED;

  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_METHOD",     
                    (void *) alg);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_OUTPUT_LEVEL", 
                    (void *) &output_level);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_SEED", 
                    (void *) &seed);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_IPC2REDIST",       
                    (void *) &ipc2redist);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_USE_OBJ_SIZE",       
                    (void *) &use_obj_size);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_COARSE_ALG", 
                    (void *) &coarse_alg);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_FOLD",       
                    (void *) &fold);

  Zoltan_Assign_Param_Vals(zz->Params, Parmetis_params, zz->Debug_Level, 
                           zz->Proc, zz->Debug_Proc);

  /* Copy option values to ParMetis options array */
#if PARMETIS_VERSION >= 3
  if (pmv3method(alg)){
    /* ParMetis 3.0 options */
    options[PMV3_OPTION_DBGLVL] = output_level; 
    options[PMV3_OPTION_SEED] = seed; 
    options[PMV3_OPT_USE_OBJ_SIZE] = use_obj_size; 
  }
  else
#endif
  {
    /* ParMetis 2.0 options */
    options[OPTION_IPART] = coarse_alg; 
    options[OPTION_FOLDF] = fold; 
    options[OPTION_DBGLVL] = output_level; 
  }

  /* Call the real ParMetis interface */
  return Zoltan_ParMetis_Jostle( zz, num_imp, imp_gids, imp_lids,
            imp_procs, num_exp, exp_gids, exp_lids, exp_procs,
            alg, options, &ipc2redist);

#endif /* ZOLTAN_PARMETIS */
}

#if PARMETIS_VERSION >= 3
static int pmv3method( char *alg)
{
   /* Check if alg is a ParMetis 3.0 method */
   return ((!strcmp(alg, "PARTKWAY")) || (!strcmp(alg, "PARTGEOMKWAY"))
          || (!strcmp(alg, "ADAPTIVEREPART")));
}
#endif

/**********************************************************/
/* Interface routine for Jostle. This is just a simple    */
/* wrapper that sets the options etc and then calls       */
/* Zoltan_ParMetis_Jostle, where the real action is.          */
/**********************************************************/

int Zoltan_Jostle(
  ZZ *zz,               /* Zoltan structure */
  int *num_imp,         /* number of objects to be imported */
  ZOLTAN_ID_PTR *imp_gids,  /* global ids of objects to be imported */
  ZOLTAN_ID_PTR *imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,      /* list of processors to import from */
  int *num_exp,         /* number of objects to be exported */
  ZOLTAN_ID_PTR *exp_gids,  /* global ids of objects to be exported */
  ZOLTAN_ID_PTR *exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs       /* list of processors to export to */
)
{
#ifndef ZOLTAN_JOSTLE
  char *yo = "Zoltan_Jostle";
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
     "Jostle requested but not compiled into library.");
  return ZOLTAN_FATAL;

#else /* ZOLTAN_JOSTLE */
  static ZZ *zz_prev = NULL; /* Last zz structure used */
  static char *alg = "JOSTLE";
  char *cptr;
  char str[MAX_PARAM_STRING_LEN+1]; 
  char matching[MAX_PARAM_STRING_LEN+1];
  char reduction[MAX_PARAM_STRING_LEN+1];
  char connect[MAX_PARAM_STRING_LEN+1];
  char scatter[MAX_PARAM_STRING_LEN+1];
  int  output_level, threshold, gather_threshold; 
  int num_proc = zz->Num_Proc;     /* Temporary variables whose addresses are */
  int proc = zz->Proc;             /* passed to Jostle. We don't              */
  MPI_Comm comm = zz->Communicator;/* want to risk letting external packages  */
                                   /* change our zz struct.                   */

  /* Initialize Jostle if this is the first call with 
   * this Zoltan structure.
   */
  if (zz != zz_prev){
     zz_prev = zz;
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
  Zoltan_Bind_Param(Jostle_params, "JOSTLE_OUTPUT_LEVEL", 
                (void *) &output_level);
  Zoltan_Bind_Param(Jostle_params, "JOSTLE_THRESHOLD",    
                (void *) &threshold);
  Zoltan_Bind_Param(Jostle_params, "JOSTLE_GATHER_THRESHOLD", 
                (void *) &gather_threshold);
  Zoltan_Bind_Param(Jostle_params, "JOSTLE_MATCHING",     
                (void *) matching);
  Zoltan_Bind_Param(Jostle_params, "JOSTLE_REDUCTION",    
                (void *) reduction);
  Zoltan_Bind_Param(Jostle_params, "JOSTLE_CONNECT",    
                (void *) connect);
  Zoltan_Bind_Param(Jostle_params, "JOSTLE_SCATTER",    
                (void *) scatter);

  Zoltan_Assign_Param_Vals(zz->Params, Jostle_params, zz->Debug_Level, zz->Proc,
                       zz->Debug_Proc); 

  /* Set Jostle parameters using jostle_env() */
  if (threshold){
    sprintf(str, "threshold = %d", threshold);
    jostle_env(str);
  }
  if (gather_threshold){
    sprintf(str, "gather threshold = %d", gather_threshold);
    jostle_env(str);
  }
  if (matching[0]){
    sprintf(str, "matching = %s", matching);
    /* Convert to lower case */
    for (cptr=str; *cptr; cptr++){
      *cptr = tolower(*cptr);
    }
    jostle_env(str);
  }
  if (reduction[0]){
    sprintf(str, "reduction = %s", reduction);
    /* Convert to lower case */
    for (cptr=str; *cptr; cptr++){
      *cptr = tolower(*cptr);
    }
    jostle_env(str);
  }
  if (connect[0]){
    sprintf(str, "connect = %s", connect);
    /* Convert to lower case */
    for (cptr=str; *cptr; cptr++){
      *cptr = tolower(*cptr);
    }
    jostle_env(str);
  }
  if (scatter[0]){
    sprintf(str, "scatter = %s", scatter);
    /* Convert to lower case */
    for (cptr=str; *cptr; cptr++){
      *cptr = tolower(*cptr);
    }
    jostle_env(str);
  }

  /* Set imbalance tolerance */
  sprintf(str, "imbalance = %3d ", (int)(100*(zz->LB.Imbalance_Tol - 1)));
  jostle_env(str);

  /* Multidimensional vertex weights */
  if (zz->Obj_Weight_Dim > 1){
    sprintf(str, "ntypes = %3d ", zz->Obj_Weight_Dim);
    jostle_env(str);
  }

  /* Call the real Jostle/ParMetis interface */
  return Zoltan_ParMetis_Jostle( zz, num_imp, imp_gids, imp_lids,
            imp_procs, num_exp, exp_gids, exp_lids, exp_procs,
            alg, &output_level, NULL);

#endif /* ZOLTAN_JOSTLE */
}

/****************************************************************/
/* This routine constructs the required graph data structures   */
/* and calls ParMetis or Jostle.                                */
/*                                                              */
/* Author: Erik Boman, eboman@cs.sandia.gov (9226)              */
/****************************************************************/

#if (defined(ZOLTAN_JOSTLE) || defined(ZOLTAN_PARMETIS))

/* Macro to free all allocated memory */
#define FREE_MY_MEMORY \
  { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "ParMETIS/Jostle error."); \
  ZOLTAN_FREE(&vtxdist); ZOLTAN_FREE(&xadj); ZOLTAN_FREE(&adjncy); \
  ZOLTAN_FREE(&vwgt); ZOLTAN_FREE(&adjwgt); ZOLTAN_FREE(&part); \
  ZOLTAN_FREE(&float_vwgt); ZOLTAN_FREE(&ewgts); ZOLTAN_FREE(&xyz); \
  ZOLTAN_FREE(&sendbuf); ZOLTAN_FREE(&recvbuf); \
  ZOLTAN_FREE(&hash_nodes); ZOLTAN_FREE(&hashtab); \
  ZOLTAN_FREE(&nbors_proc); ZOLTAN_FREE(&nbors_global); \
  ZOLTAN_FREE(&local_ids); ZOLTAN_FREE(&global_ids); \
  ZOLTAN_FREE(&proc_list); ZOLTAN_FREE(&proc_list_nbor); ZOLTAN_FREE(&plist); \
  ZOLTAN_FREE(&tmp_ewgts); ZOLTAN_FREE(&imb_tols);\
  ZOLTAN_FREE(&vsize); ZOLTAN_FREE(&tpwgt);\
  }

static int Zoltan_ParMetis_Jostle(
  ZZ *zz,               /* Zoltan structure */
  int *num_imp,         /* number of objects to be imported */
  ZOLTAN_ID_PTR *imp_gids,  /* global ids of objects to be imported */
  ZOLTAN_ID_PTR *imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,      /* list of processors to import from */
  int *num_exp,         /* number of objects to be exported */
  ZOLTAN_ID_PTR *exp_gids,  /* global ids of objects to be exported */
  ZOLTAN_ID_PTR *exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs,      /* list of processors to export to */
  char *alg,            /* algorithm to use */
  int  *options,        /* option array */
  float *ipc2redist     /* ParMetis 3.0 parameter */
)
{
  static char *yo = "Zoltan_ParMetis_Jostle";
  int i, j, jj, k, ierr, packet_size, offset, tmp, flag, ndims; 
  int obj_wgt_dim, comm_wgt_dim, check_graph, scatter;
  int num_obj, nedges, num_edges, cross_edges, max_edges, edgecut;
  int *nbors_proc, *plist;
  int nsend, nrecv, wgtflag, numflag, num_border, max_proc_list_len;
  int get_graph_data, get_geom_data, get_times; 
  idxtype *vtxdist, *xadj, *adjncy, *vwgt, *adjwgt, *part, *part2, *vsize;
  int tmp_num_obj, nself, ncon;
  float *float_vwgt, *ewgts, *xyz, *imb_tols, *tmp_ewgts, *tpwgt; 
  double geom_vec[6];
  struct Edge_Info *ptr;
  struct Edge_Info *proc_list;   /* Edge information; contains global IDs
                                       of objects with off-processor nbors. */
  ZOLTAN_ID_PTR proc_list_nbor;         /* Global IDs of neighbors of proc_list 
                                       entries.  This array is separate from
                                       proc_list to prevent individual mallocs
                                       for nbor global IDs.   */
  struct Hash_Node **hashtab, *hash_nodes;
  ZOLTAN_ID_PTR local_ids;
  ZOLTAN_ID_PTR global_ids;     /* Do not deallocate while still using the hash
                               table with num_obj (hash_node) or 
                               proc_list (Edge_Info); these data structures
                               point to global IDs in this array. */
  ZOLTAN_ID_PTR nbors_global;
  ZOLTAN_ID_PTR lid;            /* Temporary pointer to a local id; used to pass
                               NULL to query fns when NUM_LID_ENTRIES == 0. */
  char *sendbuf; 
  char *recvbuf;            /* Do not deallocate while still using the hash
                               table with nrecv; the hash_nodes point to 
                               global IDs in this array. */
  ZOLTAN_COMM_OBJ *comm_plan;
  double times[5];
  char msg[256];
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int gid_size = num_gid_entries * sizeof(ZOLTAN_ID_TYPE);
  int gid_off, lid_off;
  int num_proc = zz->Num_Proc;     /* Temporary variables whose addresses are */
                                   /* passed to Jostle and ParMETIS.  Don't   */
  MPI_Comm comm = zz->Communicator;/* want to risk letting external packages  */
                                   /* change our zz struct.                   */
  int i99;                         /* Variables used for debugging.           */
#ifdef ZOLTAN_JOSTLE
  int nnodes;
  int network[4] = {0, 1, 1, 1};
#endif

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Set default return values (in case of early exit) */
  *num_exp = -1;
  *num_imp = -1; /* No import data */

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  vtxdist = xadj = adjncy = vwgt = adjwgt = part = vsize = NULL;
  nbors_proc = NULL;
  nbors_global = NULL;
  local_ids = NULL;
  global_ids = NULL;
  float_vwgt = ewgts = tmp_ewgts = xyz = imb_tols = tpwgt = NULL;
  ptr = proc_list = NULL;
  proc_list_nbor = NULL;
  hashtab = NULL;
  hash_nodes = NULL;
  sendbuf = recvbuf = NULL;
  plist = NULL;

  /* Check weight dimensions */
  if (zz->Obj_Weight_Dim<0){
    sprintf(msg, "Object weight dimension is %d, "
            "but should be >= 0. Using Obj_Weight_Dim = 0.",
            zz->Obj_Weight_Dim);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    obj_wgt_dim = 0;
  }
  else {
    obj_wgt_dim = zz->Obj_Weight_Dim;
  }
  if (zz->Edge_Weight_Dim<0){
    sprintf(msg, "Communication weight dimension is %d, "
            "but should be >= 0. Using Edge_Weight_Dim = 0.",
            zz->Edge_Weight_Dim);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    comm_wgt_dim = 0;
  }
  else if (zz->Edge_Weight_Dim>1){
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "This method does not support "
        "multidimensional communication weights. Using Edge_Weight_Dim = 1.");
    comm_wgt_dim = 1;
  }
  else {
    comm_wgt_dim = zz->Edge_Weight_Dim;
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    printf("[%1d] Debug: alg=%s, Obj_Weight_Dim=%d, Edge_Weight_Dim=%d\n", 
      zz->Proc, alg, obj_wgt_dim, comm_wgt_dim);
    printf("[%1d] Debug: ParMetis options = %d, %d, %d, %d\n", zz->Proc,
      options[0], options[1], options[2], options[3]);
  }

  /* Start timer */
  get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
  if (get_times){
    MPI_Barrier(zz->Communicator);
    times[0] = Zoltan_Time(zz->Timer);
  }

  /* Get parameter options shared by ParMetis and Jostle */
  check_graph = 1;          /* default */
  scatter = 1;              /* default */
  Zoltan_Bind_Param(Graph_params, "CHECK_GRAPH", (void *) &check_graph);
  Zoltan_Bind_Param(Graph_params, "SCATTER_GRAPH", (void *) &scatter);
  Zoltan_Assign_Param_Vals(zz->Params, Graph_params, zz->Debug_Level, zz->Proc,
                       zz->Debug_Proc);

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
  num_obj = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);
  if (ierr){
    /* Return error code */
    FREE_MY_MEMORY;
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }
  
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("[%1d] Debug: num_obj =%d\n", zz->Proc, num_obj);

  
  vtxdist = (idxtype *)ZOLTAN_MALLOC((zz->Num_Proc+1)* sizeof(idxtype));
  if (num_obj>0){
    global_ids = ZOLTAN_MALLOC_GID_ARRAY(zz, num_obj);
    local_ids =  ZOLTAN_MALLOC_LID_ARRAY(zz, num_obj);
    if (obj_wgt_dim)
      float_vwgt = (float *)ZOLTAN_MALLOC(obj_wgt_dim*num_obj * sizeof(float));
    if (!vtxdist || !global_ids || (num_lid_entries && !local_ids) || 
        (obj_wgt_dim && !float_vwgt)){
      /* Not enough memory */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    Zoltan_Get_Obj_List(zz, global_ids, local_ids, obj_wgt_dim, float_vwgt, &ierr);
    if (ierr){
      /* Return error */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_FATAL;
    }
  
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
      printf("[%1d] Debug: Global ids = ", zz->Proc);
      for (i99=0; i99<num_obj; i99++) {
        printf("    ");
        ZOLTAN_PRINT_GID(zz, &(global_ids[i99*num_gid_entries]));
        printf("\n");
      }
    }
  }
  
  /* Construct vtxdist[i] = the number of objects on all procs < i. */
  /* Scan to compute partial sums of the number of objs */
  MPI_Scan (&num_obj, vtxdist, 1, IDX_DATATYPE, MPI_SUM, zz->Communicator);
  /* Gather data from all procs */
  MPI_Allgather (&vtxdist[0], 1, IDX_DATATYPE, 
                 &vtxdist[1], 1, IDX_DATATYPE, zz->Communicator);
  vtxdist[0] = 0;
  
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    printf("[%1d] Debug: vtxdist = ", zz->Proc);
    for (i99=0; i99<=zz->Num_Proc; i99++)
      printf("%d ", vtxdist[i99]);
    printf("\n");
  }
  
  if (zz->Debug_Level){
     if ((zz->Proc ==0) && (vtxdist[zz->Num_Proc]==0))
        ZOLTAN_PRINT_WARN(zz->Proc, yo, "No objects to balance.");
  }

  if (get_graph_data){

    num_edges = 0;
    max_edges = 0;
    for (i=0; i< num_obj; i++){
      lid = (num_lid_entries ? &(local_ids[i*num_lid_entries]) : NULL);
      nedges = zz->Get_Num_Edges(zz->Get_Edge_List_Data, 
                                 num_gid_entries, num_lid_entries,
                                 &(global_ids[i*num_gid_entries]), 
                                 lid, &ierr);
      if (ierr){
        /* Return error */
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ierr);
      }
      num_edges += nedges;
      if (nedges>max_edges) max_edges = nedges;
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: num_edges = %d\n", zz->Proc, num_edges);

    /* Allocate space for ParMETIS data structs */
    xadj   = (idxtype *)ZOLTAN_MALLOC((num_obj+1) * sizeof(idxtype));
    adjncy = (idxtype *)ZOLTAN_MALLOC(num_edges * sizeof(idxtype));
    if (comm_wgt_dim) 
      adjwgt = (idxtype *)ZOLTAN_MALLOC(comm_wgt_dim * num_edges 
                 * sizeof(idxtype));
  
    if (!xadj || (num_edges && !adjncy) || 
        (num_edges && comm_wgt_dim && !adjwgt)){
      /* Not enough memory */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: Successfully allocated ParMetis space\n", zz->Proc);
  
    /* Construct ParMETIS graph */
    /* First compute a global dense numbering of the objects/vertices */
  
    /* Construct local hash table */
    hash_nodes = (struct Hash_Node *)ZOLTAN_MALLOC(num_obj *
      sizeof(struct Hash_Node));
    hashtab = (struct Hash_Node **) ZOLTAN_MALLOC(num_obj *
      sizeof(struct Hash_Node *) );
    if (num_obj && ((!hash_nodes) || (!hashtab))){
      /* Not enough memory */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    
    for (i=0; i< num_obj; i++){
      hashtab[i] = NULL;
      hash_nodes[i].gid = &(global_ids[i*num_gid_entries]);
      hash_nodes[i].gno = vtxdist[zz->Proc]+i;
    }

    for (i=0; i< num_obj; i++){
      /* insert hashed elements into hash table */
      j = Zoltan_Hash(&(global_ids[i*num_gid_entries]), num_gid_entries, (unsigned int)num_obj);
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
    nbors_global = ZOLTAN_MALLOC_GID_ARRAY(zz, max_edges);
    nbors_proc = (int *)ZOLTAN_MALLOC(max_edges * sizeof(int));
    plist = (int *)ZOLTAN_MALLOC(zz->Num_Proc * sizeof(int));
    if (comm_wgt_dim && max_edges){
      tmp_ewgts = (float *)ZOLTAN_MALLOC(comm_wgt_dim * max_edges * sizeof(float));
      ewgts = (float *)ZOLTAN_MALLOC(comm_wgt_dim * num_edges * sizeof(float));
    }

    if ((max_edges && ((!nbors_global) || (!nbors_proc) ||
                       (comm_wgt_dim && !ewgts) || 
                       (comm_wgt_dim && !tmp_ewgts))) || (!plist)){
      /* Not enough memory */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    for (i=0; i<zz->Num_Proc; i++)
      plist[i] = -1;
  
    if (max_proc_list_len){
      /* Allocate space for processor list */
      while ((proc_list==NULL || proc_list_nbor == NULL)
          && (max_proc_list_len>=CHUNKSIZE)){
        proc_list = (struct Edge_Info *) ZOLTAN_MALLOC(max_proc_list_len *
          sizeof(struct Edge_Info) );
        proc_list_nbor = ZOLTAN_MALLOC_GID_ARRAY(zz, max_proc_list_len);
        if (!proc_list || !proc_list_nbor){
          /* Not enough memory, try shorter list */
          ZOLTAN_FREE(&proc_list);
          ZOLTAN_FREE(&proc_list_nbor);
          if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
            printf("[%1d] Debug: Could not allocate %d list nodes, "
                   "trying %d instead.\n", zz->Proc,
                   max_proc_list_len, max_proc_list_len/2);
          }
          max_proc_list_len /= 2;
        }
      }
      if (!proc_list || !proc_list_nbor){
        /* Not enough memory */
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
    }

    /* proc_list[i] will contain a struct with data to send
     * to another processor.  proc_list_nbor[i*num_gid_entries] contains
     * the global ID of the neighboring object.
     * We don't know yet the total number
     * of inter-proc edges so we may have to adjust the size of proc_list and
     * proc_list_nbor through REALLOC.  We add a chunk of space at a time.
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
    xadj[0] = 0;
    jj = 0;
    nself = 0;      /* Number of self-edges in the graph */
  
    for (i=0; i< num_obj; i++){
      gid_off = i * num_gid_entries;
      lid_off = i * num_lid_entries;
      lid = (num_lid_entries ? &(local_ids[lid_off]) : NULL);

      nedges = zz->Get_Num_Edges(zz->Get_Edge_List_Data, 
                                 num_gid_entries, num_lid_entries,
                                 &(global_ids[gid_off]), lid, 
                                 &ierr);
      zz->Get_Edge_List(zz->Get_Edge_List_Data,
                        num_gid_entries, num_lid_entries,
                        &(global_ids[gid_off]), lid, 
                        nbors_global, nbors_proc, comm_wgt_dim, 
                        tmp_ewgts, &ierr);
      if (ierr){
        /* Return error */
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ierr);
      }
  
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
        printf("[%1d] Debug: i=%d, gid=", zz->Proc, i);
        ZOLTAN_PRINT_GID(zz, &(global_ids[gid_off]));
        printf("lid=");
        ZOLTAN_PRINT_LID(zz, lid);
        printf("nedges=%d\n", nedges);
      }

      /* Separate inter-processor edges from the local ones */
      /* Also remove any self-edges */
      for (j=0; j<nedges; j++){
        if (nbors_proc[j] == zz->Proc){
          /* Local edge */
          tmp = hash_lookup(zz, hashtab, 
                            &(nbors_global[j*num_gid_entries]), num_obj);
          if (tmp == i+vtxdist[zz->Proc]){
            /* Self-edge! */
            nself++;
          }
          else{
            /* Copy over edge weights. */
            for (k=0; k<comm_wgt_dim; k++)
              ewgts[jj*comm_wgt_dim+k] = tmp_ewgts[j*comm_wgt_dim+k];
            /* Put the global number into the adjacency array */
            adjncy[jj++] = tmp;
          }
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
            if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
              printf("[%1d] Debug: Allocating more list space, "
                     "max_proc_list_len = %d, increasing by %d\n", 
                     zz->Proc, max_proc_list_len, CHUNKSIZE);

            max_proc_list_len += CHUNKSIZE;
            proc_list = (struct Edge_Info *) ZOLTAN_REALLOC(proc_list,
                         max_proc_list_len*sizeof(struct Edge_Info));
            proc_list_nbor = ZOLTAN_REALLOC_GID_ARRAY(zz, proc_list_nbor,
                              max_proc_list_len);
            if (!proc_list){
              /* Not enough memory */
              FREE_MY_MEMORY;
              ZOLTAN_TRACE_EXIT(zz, yo);
              return ZOLTAN_MEMERR;
            }
          }
          ptr = &proc_list[offset];
          ptr->my_gid = &(global_ids[gid_off]);
          ptr->my_gno = hash_lookup(zz, hashtab, 
                                    &(global_ids[gid_off]), num_obj);
          ZOLTAN_SET_GID(zz, &(proc_list_nbor[offset*num_gid_entries]),
                         &(nbors_global[j*num_gid_entries]));
          if (flag)
            ptr->nbor_proc = nbors_proc[j];
          else
            ptr->nbor_proc = -1;
          ptr->adj = &adjncy[jj];

          if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
            printf("[%1d] Debug: proc_list[%1d] my_gid=", zz->Proc, offset);
            ZOLTAN_PRINT_GID(zz, ptr->my_gid);
            printf(", my_gno=%d, nbor_proc=%d\n", ptr->my_gno, ptr->nbor_proc);
          }

          /* Copy over edge weights. */
          for (k=0; k<comm_wgt_dim; k++)
            ewgts[jj*comm_wgt_dim+k] = tmp_ewgts[j*comm_wgt_dim+k];

          /* Still don't know the global number, need to come back here later */
          adjncy[jj++] = -1; 

          offset++;
        }
      }
      xadj[i+1] = jj; /* NB: We do not count self-edges. */
    }
    cross_edges = offset;

    ZOLTAN_FREE(&plist);
    ZOLTAN_FREE(&nbors_global);
    ZOLTAN_FREE(&nbors_proc);
    ZOLTAN_FREE(&tmp_ewgts);

    /* Warn if we removed any self-edges */
    if (check_graph >= 1){
      if (nself>0) ierr = ZOLTAN_WARN;
      MPI_Reduce(&nself, &tmp, 1, MPI_INT, MPI_SUM, 0, zz->Communicator);
      if ((zz->Proc==0) && (tmp>0) && (zz->Debug_Level>0)){
          sprintf(msg, "Found and removed %d self edges in the graph.\n", 
                  tmp);
          ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
      }
    }

    /* Sanity check */
    if ((check_graph >= 1) && (xadj[num_obj] + nself != num_edges)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid graph. "
              "Something may be wrong with the edges in the graph, "
              "or perhaps you found a bug in Zoltan.\n"); 
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_FATAL;
    }
  
    /* Exchange info between processors to resolve global number 
     * for objects that are off-proc.
     */

    /* Allocate send buffer */
    packet_size = gid_size + sizeof(int);
    sendbuf = (char *) ZOLTAN_MALLOC(nsend * packet_size);
    plist = (int *) ZOLTAN_MALLOC(nsend * sizeof(int));

    if (nsend && (!sendbuf || !plist) ){
      /* Not enough space */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }

    /* Pack the data to send */
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: %d messages to send.\n", zz->Proc, nsend);

    offset = 0;
    j = 0;
    for (i=0, ptr=proc_list; i<cross_edges; i++, ptr++){
      if (ptr->nbor_proc >= 0){
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
          printf("[%1d] Debug: Sending gid = ", zz->Proc);
          ZOLTAN_PRINT_GID(zz, ptr->my_gid);
          printf(", gno = %d to proc %d\n", ptr->my_gno, ptr->nbor_proc);
        }
        memcpy(&sendbuf[offset], (char *) (ptr->my_gid), gid_size); 
        offset += gid_size;
        memcpy(&sendbuf[offset], (char *) &(ptr->my_gno), sizeof(int)); 
        offset += sizeof(int);
        plist[j++] = ptr->nbor_proc;
      }
    }

    /* Create the communication plan */
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: Calling Zoltan_Comm_Create with %d packets to send.\n",
             zz->Proc, nsend);

    ierr = Zoltan_Comm_Create(&comm_plan, nsend, plist, comm, TAG1, &nrecv);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      /* Return error code */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }

    /* Allocate recv buffer */
    recvbuf = (char *) ZOLTAN_MALLOC(nrecv * packet_size);
    if (nrecv && (!sendbuf || !plist) ){
      /* Not enough space */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("[%1d] Debug: Ready to receive %d packets.\n", 
        zz->Proc, nrecv);

    /* Do the communication */
    ierr = Zoltan_Comm_Do( comm_plan, TAG2, sendbuf, packet_size, recvbuf);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      /* Return error code */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }

    /* Destroy the comm. plan */
    Zoltan_Comm_Destroy( &comm_plan);
  
    /* Unpack data into the ParMETIS/Jostle struct.
     * Resolve off-proc global numbers by:
     * 1) Insert global ids from recvbuf into hash table;
     * 2) Look up missing references in proc_list 
     *    using the hash table.
     */

    /* Change hash table to contain only border objects */
    /* that we received from other procs.               */
    hash_nodes = (struct Hash_Node *)ZOLTAN_REALLOC(hash_nodes,
      nrecv * sizeof(struct Hash_Node));
    hashtab = (struct Hash_Node **) ZOLTAN_REALLOC(hashtab,
      nrecv * sizeof(struct Hash_Node *) );
    if (nrecv && ((!hash_nodes) || (!hashtab))){
      /* Not enough memory */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    
    /* Copy data from recvbuf into hash table nodes */
    for (i=0; i< nrecv; i++){
      hashtab[i] = NULL;
      hash_nodes[i].gid = (ZOLTAN_ID_PTR) &(recvbuf[i*packet_size]);
      hash_nodes[i].gno = *((int *)&recvbuf[i*packet_size+gid_size]);
      /* Do we need to pad for byte alignment? */
    }
  
    /* Insert nodes into hash table */
    for (i=0; i< nrecv; i++){
      j = Zoltan_Hash(hash_nodes[i].gid, num_gid_entries, (unsigned int)nrecv);
      hash_nodes[i].next = hashtab[j];
      hashtab[j] = &hash_nodes[i];
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
        printf("[%1d] Debug: Hashed GID ", zz->Proc);
        ZOLTAN_PRINT_GID(zz, hash_nodes[i].gid);
        printf(" to %d, gno = %d\n", j, hash_nodes[i].gno);
      }
    }

    for (i=0; i<cross_edges; i++){
      /* Look up unresolved global_ids */
      if ((tmp=hash_lookup(zz, hashtab, &(proc_list_nbor[i*num_gid_entries]), 
                           nrecv)) <0){
        /* Error: Global ID is not in hash table. 
           This only happens if the graph is invalid. */
        ZOLTAN_PRINT_ERROR(zz->Proc, yo,"Invalid graph. Please check that "
           "your graph query functions are correct.\n");
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_FATAL;
      }
      else{
        /* Insert the global number into adjncy vector */
        *(proc_list[i].adj) = tmp;
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
          printf("[%1d] Debug: GID ", zz->Proc);
          ZOLTAN_PRINT_GID(zz, &(proc_list_nbor[i*num_gid_entries]));
          printf(" has global number %d\n", tmp);
        }
      }
      
    }

    /* Free space */
    ZOLTAN_FREE(&sendbuf);
    ZOLTAN_FREE(&recvbuf);
    ZOLTAN_FREE(&proc_list);
    ZOLTAN_FREE(&proc_list_nbor);
    ZOLTAN_FREE(&plist);
    ZOLTAN_FREE(&hash_nodes);
    ZOLTAN_FREE(&hashtab);
  
    /* Get vertex weights if needed */
    if (obj_wgt_dim){
      vwgt = (idxtype *)ZOLTAN_MALLOC(obj_wgt_dim*num_obj
                          * sizeof(idxtype));
      ierr = scale_round_weights(float_vwgt, vwgt, num_obj, obj_wgt_dim, 1,
                                 (int) MAX_WGT_SUM, zz->Debug_Level, zz->Communicator);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        /* Return error code */
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ierr;
      }

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
        for (i99=0; i99 < (num_obj<3 ? num_obj : 3); i99++){
          for (k=0; k<obj_wgt_dim; k++)
            sprintf(msg+10*k, "%9d ", vwgt[i99*obj_wgt_dim+k]);
          printf("[%1d] Debug: scaled weights for vertex %d = %s\n", 
                 zz->Proc, vtxdist[zz->Proc]+i99, msg);
        }
      ZOLTAN_FREE(&float_vwgt);
    }

    /* Get edge weights if needed */
    if (comm_wgt_dim){
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] Debug: Edge weights are (before scaling) = \n", zz->Proc);
        for (j=0; j<num_edges*comm_wgt_dim; j++)
          printf("%f ", ewgts[j]);
        printf("\n");
      }
      /* Scale and round edge weights to integers */
      ierr = scale_round_weights(ewgts, adjwgt, num_edges, comm_wgt_dim, 1,
                                 (int) MAX_WGT_SUM, zz->Debug_Level, zz->Communicator);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        /* Return error code */
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ierr;
      }

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] Debug: Edge weights are (after scaling) = \n", zz->Proc);
        for (j=0; j<num_edges; j++)
          printf("%d ", adjwgt[j]);
        printf("\n");
      }
      ZOLTAN_FREE(&ewgts);
    }

  } /* end get_graph_data */

  if (get_geom_data){
    /* Determine how many dimensions the data have */
    ndims = zz->Get_Num_Geom(zz->Get_Num_Geom_Data, &ierr);
    if (ierr){
      /* Return error */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }
    /* Allocate space for the geometry data */
    xyz = (float *) ZOLTAN_MALLOC(ndims*num_obj * sizeof(float));
    if (ndims && num_obj && !xyz){
      /* Not enough space */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    /* Get the geometry data */
    for (i=0; i<num_obj; i++){
      lid = (num_lid_entries ? &(local_ids[i*num_lid_entries]) : NULL);
      zz->Get_Geom(zz->Get_Geom_Data, 
                   num_gid_entries, num_lid_entries,
                   &(global_ids[i*num_gid_entries]), 
                   lid, geom_vec, &ierr);
      if (ierr) {
        /* Return error code */
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ierr);
      }
      for (j=0; j<ndims; j++)
        xyz[i*ndims+j] = geom_vec[j];
    }
  }

#if (PARMETIS_VERSION >= 3)
  /* Get object sizes if requested */ 
  if (options[PMV3_OPT_USE_OBJ_SIZE] && zz->Migrate.Get_Obj_Size){
    vsize = (idxtype *) ZOLTAN_MALLOC(num_obj*sizeof(idxtype));
    if (num_obj && !vsize){
      /* Not enough space */
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    for (i=0; i<num_obj; i++){
      lid = (num_lid_entries ? &(local_ids[i*num_lid_entries]) : NULL);
      vsize[i] = zz->Migrate.Get_Obj_Size(zz->Migrate.Get_Obj_Size_Data, 
                     num_gid_entries, num_lid_entries,
                     &(global_ids[i*num_gid_entries]), 
                     lid, &ierr);
    }
  }
#endif

  /* Scatter graph? */
  tmp_num_obj = num_obj;
  if ((scatter>0) && (scatter<3)){
    /* Decide if the data imbalance is so bad that we should scatter the graph. */
    /* scatter==1: Scatter if all the objects are on a single processor.        */
    /* scatter==2: Scatter if any processor has no objects.                     */
    if (num_obj==0)
      j = 1;
    else 
      j = 0;
    MPI_Allreduce(&j, &tmp, 1, MPI_INT, MPI_SUM, zz->Communicator);
    if (scatter == 1){
      if (tmp < zz->Num_Proc-1)
        scatter = 0;
    }
    else if (scatter==2){
      if (tmp == 0)
        scatter = 0;
    }
  }

  if (scatter){
    ierr = Zoltan_Scatter_Graph(&vtxdist, &xadj, &adjncy, &vwgt, &vsize,
              &adjwgt, &xyz, ndims, zz, &comm_plan);
    if ((ierr == ZOLTAN_FATAL) || (ierr == ZOLTAN_MEMERR)){
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ierr;
    }
    tmp_num_obj = vtxdist[zz->Proc+1]-vtxdist[zz->Proc];
  }

  /* Get ready to call ParMETIS or Jostle */
  edgecut = -1; 
  wgtflag = 2*(obj_wgt_dim>0) + (comm_wgt_dim>0); 
  numflag = 0;
  part = (idxtype *)ZOLTAN_MALLOC((tmp_num_obj+1) * sizeof(idxtype));
  if (!part){
    /* Not enough memory */
    FREE_MY_MEMORY;
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }
  ncon = (obj_wgt_dim > 0 ? obj_wgt_dim : 1);

  /* Set Imbalance Tolerance for each weight component. For now, they are all the same. */
  imb_tols = (float *) ZOLTAN_MALLOC(ncon * sizeof(float));
  if (!imb_tols){
    /* Not enough memory */
    FREE_MY_MEMORY;
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<ncon; i++)
    imb_tols[i] = zz->LB.Imbalance_Tol;

  /* Verify that graph is correct */
  if (get_graph_data){
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) 
      flag = 2; /* Maximum output level */
    else
      flag = 1; /* Medium output level */
    ierr = Zoltan_Verify_Graph(zz->Communicator, vtxdist, xadj, adjncy, vwgt, 
              adjwgt, obj_wgt_dim, comm_wgt_dim, check_graph, flag);
  
    /* Special error checks to avoid certain death in ParMETIS 2.0 */
    if (xadj[tmp_num_obj] == 0){
      /* No edges on a proc is a fatal error in ParMETIS 2.0
       * and in Jostle 1.2. This error test should be removed
       * when the bugs in ParMETIS and Jostle have been fixed.
       */
      if (strcmp(alg, "JOSTLE") == 0){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No edges on this processor; "
                      "Jostle will likely crash. "
                      "Please use a different load balancing method.");
      } else { /* ParMETIS */
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "No edges on this processor; "
                      "ParMETIS will likely crash. "
                      "Please use a different load balancing method.");
      }
      ierr = ZOLTAN_FATAL;
      FREE_MY_MEMORY; 
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }
  }
  
  /* Get a time here */
  if (get_times) times[1] = Zoltan_Time(zz->Timer);

  /* Select the desired ParMetis or Jostle function */

#if (PARMETIS_VERSION >= 3)
  tpwgt = (float *) ZOLTAN_MALLOC(ncon*num_proc * sizeof(float));
  if (!tpwgt){
    /* Not enough memory */
    FREE_MY_MEMORY;
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_MEMERR;
  }
  /* For now, all desired partition sizes are equal */
  for (i=0; i<ncon*num_proc; i++)
    tpwgt[i] = 1.0/num_proc;
#else /* PARMETIS_VERSION < 3 */
  if ((ncon >= 2) && strcmp(alg, "JOSTLE")) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "You need ParMETIS version 3.0 or higher to use multi-weights.\n"
      "If you have this installed, please make sure you include the appropriate\n"
      "version of parmetis.h and recompile Zoltan.");
    FREE_MY_MEMORY;
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }
#endif 
  if (strcmp(alg, "JOSTLE") == 0){
#ifdef ZOLTAN_JOSTLE
    offset = 0;            /* Index of the first object/node. */
    j = 0;                 /* Dummy variable for Jostle */
    nnodes = vtxdist[zz->Num_Proc]; /* Global number of objects */ 
    ndims = 1;             /* Topological dimension of the computer */
    network[1] = zz->Num_Proc;
    /* Convert xadj array to Jostle format (degree) */
    for (i=0; i<tmp_num_obj; i++){
      xadj[i] = xadj[i+1] - xadj[i];
    }
    /* Convert vtxdist array to Jostle format (degree) */
    for (i=0; i<zz->Num_Proc; i++){
      vtxdist[i] = vtxdist[i+1] - vtxdist[i];
    }
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the Jostle library");
    jostle_env("format = contiguous");
    if (check_graph >= 2){
      jostle_env("check = on");
    }
    pjostle(&nnodes, &offset, &(vtxdist[zz->Proc]), &j, vtxdist, 
       xadj, vwgt, part, &num_edges, adjncy, adjwgt, network,
       NULL, options, &ndims, NULL); 
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the Jostle library");
#else
    /* We don't know about Jostle */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "Sorry, Jostle is not available on this system.\n"
      "If you have Jostle, please set the JOSTLE_XXXPATHs appropriately "
      "in the Zoltan configuration files and recompile Zoltan. Otherwise, "
      "use a different method, for example ParMETIS.");
    FREE_MY_MEMORY;
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
#endif /* ZOLTAN_JOSTLE */
  }
#if PARMETIS_VERSION >= 3
  else if (strcmp(alg, "PARTKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
    ParMETIS_V3_PartKway (vtxdist, xadj, adjncy, vwgt, adjwgt, 
      &wgtflag, &numflag, &ncon, &num_proc, tpwgt,
      imb_tols, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "PARTGEOMKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
    ParMETIS_V3_PartGeomKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag,
      &numflag, &ndims, xyz, &ncon, &num_proc, tpwgt,
      imb_tols, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "ADAPTIVEREPART") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
    ParMETIS_V3_AdaptiveRepart (vtxdist, xadj, adjncy, vwgt, vsize, adjwgt, 
      &wgtflag, &numflag, &ncon, &num_proc, tpwgt, imb_tols, 
      ipc2redist, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
#endif  /* PARMETIS_VERSION >= 3 */
  /* Check for ParMetis 2.0 routines */
  else if (strcmp(alg, "PARTKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_PartKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, &num_proc, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "PARTGEOMKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_PartGeomKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag,
      &numflag, &ndims, xyz, &num_proc, options, &edgecut, 
      part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "PARTGEOM") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_PartGeom (vtxdist, &ndims, xyz, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REPARTLDIFFUSION") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_RepartLDiffusion (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REPARTGDIFFUSION") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_RepartGDiffusion (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REPARTREMAP") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_RepartRemap (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REPARTMLREMAP") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_RepartMLRemap (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REFINEKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_RefineKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else {
    /* This should never happen! */
    sprintf(msg, "Unknown ParMetis or Jostle algorithm %s.", alg);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    FREE_MY_MEMORY;
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ZOLTAN_FATAL;
  }

  /* Get a time here */
  if (get_times) times[2] = Zoltan_Time(zz->Timer);

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("[%1d] Debug: Returned from partitioner with edgecut= %d\n", 
      zz->Proc, edgecut);

  /* Free weights; they are no longer needed */
  if (obj_wgt_dim) {
    ZOLTAN_FREE(&vwgt);
    if (float_vwgt) ZOLTAN_FREE(&float_vwgt);
  }
  if (comm_wgt_dim){
    ZOLTAN_FREE(&adjwgt);
  }
  /* Also free temp arrays for ParMetis 3.0 */
  ZOLTAN_FREE(&tpwgt);
  ZOLTAN_FREE(&imb_tols);
  if (vsize) ZOLTAN_FREE(&vsize);

  /* If we have been using a scattered graph, convert partition result back to 
   * original distribution 
   */
  if (scatter){
    /* Allocate space for partition array under original distribution */
    part2 = (idxtype *) ZOLTAN_MALLOC(num_obj*sizeof(idxtype)); 
    if (num_obj && !part2){
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
    /* Use reverse communication to compute the partition array under the 
     * original distribution 
     */
    ierr = Zoltan_Comm_Do_Reverse(comm_plan, TAG3, (char *) part, sizeof(idxtype), 
                              NULL, (char *) part2);
    if ((ierr == ZOLTAN_FATAL) || (ierr == ZOLTAN_MEMERR)){
      FREE_MY_MEMORY;
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ierr);
    }
    Zoltan_Comm_Destroy(&comm_plan); /* Destroy the comm. plan */
    /* We don't need the partition array with the scattered distribution 
     * any more */
    ZOLTAN_FREE(&part); 
    /* part is now the partition array under the original distribution */
    part = part2;   
  }
 
  /* Determine number of objects to export */
  nsend = 0;
  for (i=0; i<num_obj; i++)
    if (part[i] != zz->Proc) nsend++;

  /* Create export lists */
  if (zz->LB.Return_Lists){
    (*num_exp) = nsend;
    if (nsend > 0) {
      if (!Zoltan_Special_Malloc(zz,(void **)exp_gids,nsend,ZOLTAN_SPECIAL_MALLOC_GID)) {
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      if (!Zoltan_Special_Malloc(zz,(void **)exp_lids,nsend,ZOLTAN_SPECIAL_MALLOC_LID)) {
        Zoltan_Special_Free(zz,(void **)exp_gids,ZOLTAN_SPECIAL_MALLOC_GID);
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      if (!Zoltan_Special_Malloc(zz,(void **)exp_procs,nsend,ZOLTAN_SPECIAL_MALLOC_INT)) {
        Zoltan_Special_Free(zz,(void **)exp_lids,ZOLTAN_SPECIAL_MALLOC_LID);
        Zoltan_Special_Free(zz,(void **)exp_gids,ZOLTAN_SPECIAL_MALLOC_GID);
        FREE_MY_MEMORY;
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ZOLTAN_MEMERR;
      }
      j = 0;
      for (i=0; i<num_obj; i++){
        if (part[i] != zz->Proc){
          ZOLTAN_SET_GID(zz, &((*exp_gids)[j*num_gid_entries]),
                         &(global_ids[i*num_gid_entries]));
          ZOLTAN_SET_LID(zz, &((*exp_lids)[j*num_lid_entries]),
                         &(local_ids[i*num_lid_entries]));
          (*exp_procs)[j] = part[i];
          j++;
        }
      }
    }
  }

  /* Free space */
  ZOLTAN_FREE(&part);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&vtxdist);
  ZOLTAN_FREE(&xadj);
  ZOLTAN_FREE(&adjncy);
  ZOLTAN_FREE(&xyz);

  /* Get a time here */
  if (get_times) times[3] = Zoltan_Time(zz->Timer);

  /* Output timing results if desired */
  if (get_times){
    if (zz->Proc==0) printf("\nZOLTAN timing statistics:\n");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[1]-times[0], 
                   " Partitioner Pre-processing time  ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[2]-times[1], 
                   " Partitioner Library time         ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[3]-times[2], 
                   " Partitioner Post-processing time ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[3]-times[0], 
                   " Partitioner Total time           ");
    if (zz->Proc==0) printf("\n");
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ZOLTAN_OK;
}

/* 
 * Scale and round float weights to integer (non-negative) weights 
 * subject to sum(weights) <= max_wgt_sum.
 * Only scale if deemed necessary. 
 *
 *   mode == 0 : no scaling, just round to int
 *   mode == 1 : scale each weight dimension separately
 *   mode == 2 : use same scale factor for all weights
 * 
 */

static int scale_round_weights(float *fwgts, idxtype *iwgts, int n, int dim, 
                 int mode, int max_wgt_sum, int debug_level, MPI_Comm comm)
{
  int i, j, tmp, ierr, proc; 
  int *nonint, *nonint_local; 
  float *scale, *sum_wgt_local, *sum_wgt, *max_wgt_local, *max_wgt;
  char msg[256];
  static char *yo = "scale_round_weights";

  ierr = ZOLTAN_OK;
  MPI_Comm_rank(comm, &proc);

  if (mode == 0) {
    /* No scaling; just convert to int */
    for (i=0; i<n*dim; i++){
      iwgts[i] = (int) fwgts[i];
    }
  }
  else{
      /* Allocate local arrays */
      nonint = (int *)ZOLTAN_MALLOC(dim*sizeof(int));
      nonint_local = (int *)ZOLTAN_MALLOC(dim*sizeof(int));
      scale = (float *)ZOLTAN_MALLOC(dim*sizeof(float));
      sum_wgt = (float *)ZOLTAN_MALLOC(dim*sizeof(float));
      sum_wgt_local = (float *)ZOLTAN_MALLOC(dim*sizeof(float));
      max_wgt = (float *)ZOLTAN_MALLOC(dim*sizeof(float));
      max_wgt_local = (float *)ZOLTAN_MALLOC(dim*sizeof(float));
      if (!(nonint && nonint_local && scale && sum_wgt && sum_wgt_local
           && max_wgt && max_wgt_local)){
        ZOLTAN_PRINT_ERROR(proc, yo, "Out of memory.");
        ZOLTAN_FREE(&nonint);
        ZOLTAN_FREE(&nonint_local);
        ZOLTAN_FREE(&scale);
        ZOLTAN_FREE(&sum_wgt);
        ZOLTAN_FREE(&sum_wgt_local);
        ZOLTAN_FREE(&max_wgt);
        ZOLTAN_FREE(&max_wgt_local);
        return ZOLTAN_MEMERR;
      }
      /* Initialize */
      for (j=0; j<dim; j++){
        nonint_local[j] = 0;
        sum_wgt_local[j] = 0;
        max_wgt_local[j] = 0;
      }

      /* Compute local sums of the weights */
      /* Check if all weights are integers */
#define EPSILON (1e-5)
      for (i=0; i<n; i++){
        for (j=0; j<dim; j++){
          if (!nonint_local[j]){ 
            tmp = fwgts[i*dim+j]; /* Converts to nearest int */
            if (fabs((double)tmp-fwgts[i*dim+j]) > EPSILON){
              nonint_local[j] = 1;
            }
          }
          sum_wgt_local[j] += fwgts[i*dim+j];
          if (fwgts[i*dim+j] > max_wgt_local[j])
            max_wgt_local[j] = fwgts[i*dim+j]; 
        }
      }
      /* Compute global sum of the weights */
      MPI_Allreduce(nonint_local, nonint, dim, 
          MPI_INT, MPI_LOR, comm);
      MPI_Allreduce(sum_wgt_local, sum_wgt, dim, 
          MPI_FLOAT, MPI_SUM, comm);
      MPI_Allreduce(max_wgt_local, max_wgt, dim, 
          MPI_FLOAT, MPI_MAX, comm);

      /* Calculate scale factor */
      for (j=0; j<dim; j++){
        scale[j] = 1.;
        /* Scale unless all weights are integers (not all zero) */
        if (nonint[j] || (max_wgt[j] <= EPSILON) || (sum_wgt[j] > max_wgt_sum)){
          if (sum_wgt[j] == 0){
            ierr = ZOLTAN_WARN;
            if (proc == 0){
              sprintf(msg, "All weights are zero in component %1d", j);
              ZOLTAN_PRINT_WARN(proc, yo, msg);
            }
          }
          else /* sum_wgt[j] != 0) */
            scale[j] = max_wgt_sum/sum_wgt[j];
        }
      }

      /* If mode==2, let the scale factor be the same for all weights */
      if (mode==2){
        for (j=1; j<dim; j++){
          if (scale[j]<scale[0])
            scale[0] = scale[j];
        }
        for (j=1; j<dim; j++){
          scale[j] = scale[0];
        }
      }

      if ((debug_level >= ZOLTAN_DEBUG_ALL) && (proc==0)){
        printf("ZOLTAN DEBUG in %s: scaling weights with scale factors = ", yo);
        for (j=0; j<dim; j++)
          printf("%f ", scale[j]);
        printf("\n");
      }

      /* Convert weights to positive integers using the computed scale factor */
      for (i=0; i<n; i++){
        for (j=0; j<dim; j++){
          if (scale[j] == 1)
             iwgts[i*dim+j] = (int) fwgts[i*dim+j];
          else
             iwgts[i*dim+j] = (int) ceil(fwgts[i*dim+j]*scale[j]);
        }
      }

    ZOLTAN_FREE(&nonint);
    ZOLTAN_FREE(&nonint_local);
    ZOLTAN_FREE(&scale);
    ZOLTAN_FREE(&sum_wgt);
    ZOLTAN_FREE(&sum_wgt_local);
    ZOLTAN_FREE(&max_wgt);
    ZOLTAN_FREE(&max_wgt_local);
  }
  return ierr;
}
#undef EPSILON

#endif /* defined (ZOLTAN_JOSTLE) || defined (ZOLTAN_PARMETIS) */


/*********************************************************************/
/* ParMetis parameter routine                                        */
/*********************************************************************/

int Zoltan_ParMetis_Set_Param(
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

    status = Zoltan_Check_Param(name, val, Parmetis_params, &result, &index);
    if (status == 1)
       status = Zoltan_Check_Param(name, val, Graph_params, &result, &index);

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

int Zoltan_Jostle_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */

    status = Zoltan_Check_Param(name, val, Jostle_params, &result, &index);
    if (status == 1)
       status = Zoltan_Check_Param(name, val, Graph_params, &result, &index);

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

#if (defined(ZOLTAN_JOSTLE) || defined(ZOLTAN_PARMETIS))

/*******************************************************************
 * hash_lookup uses Zoltan_Hash to lookup a key 
 *
 * Input:
 *   hashtab, pointer to the hash table
 *   key, a key to look up of type ZOLTAN_ID_PTR (any data type)
 *   n,   dimension of the hash table
 *
 * Return value:
 *   the global number of the element with key key,
 *   or -1 if the key is not in the hash table
 *
 *******************************************************************/

static int hash_lookup (ZZ *zz, struct Hash_Node **hashtab, ZOLTAN_ID_PTR key,
                        int n)
{
  int i;
  struct Hash_Node *ptr;

  i = Zoltan_Hash(key, zz->Num_GID, (unsigned int)n);
  for (ptr=hashtab[i]; ptr != NULL; ptr = ptr->next){
    if (ZOLTAN_EQ_GID(zz, ptr->gid, key))
      return (ptr->gno);
  }
  /* Key not in hash table */
  return -1;
}

#endif /* (defined(ZOLTAN_JOSTLE) || defined(ZOLTAN_PARMETIS)) */

