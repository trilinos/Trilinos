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


#include <ctype.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "all_allo_const.h"
#include "parmetis_jostle.h"
#include "params_const.h"
#include "timer_const.h"
#include "order_const.h"

/**********  parameters structure for parmetis methods **********/
static PARAM_VARS Parmetis_params[] = {
        { "PARMETIS_METHOD", NULL, "STRING", 0 },
        { "PARMETIS_OUTPUT_LEVEL", NULL, "INT", 0 },
        { "PARMETIS_SEED", NULL, "INT", 0 },
        { "PARMETIS_ITR", NULL, "DOUBLE", 0 },
        { "PARMETIS_USE_OBJ_SIZE", NULL, "INT", 0 },
        { "PARMETIS_COARSE_ALG", NULL, "INT", 0 },
        { "PARMETIS_FOLD", NULL, "INT", 0 },
        { NULL, NULL, NULL, 0 } };

/**********  parameters structure for jostle methods **********/
static PARAM_VARS Jostle_params[] = {
        { "JOSTLE_OUTPUT_LEVEL", NULL, "INT", 0 },
        { "JOSTLE_THRESHOLD", NULL, "INT", 0 },
        { "JOSTLE_GATHER_THRESHOLD", NULL, "INT", 0 },
        { "JOSTLE_MATCHING", NULL, "STRING", 0 },
        { "JOSTLE_REDUCTION", NULL, "STRING", 0 },
        { "JOSTLE_CONNECT", NULL, "STRING", 0 },
        { "JOSTLE_SCATTER", NULL, "STRING", 0 },
        { NULL, NULL, NULL, 0 } };

/**********  parameters structure used by both ParMetis and Jostle **********/
static PARAM_VARS Graph_params[] = {
        { "CHECK_GRAPH", NULL, "INT", 0 },
        { "SCATTER_GRAPH", NULL, "INT", 0 },
        { NULL, NULL, NULL, 0 } };

/***************  prototypes for internal functions ********************/

/* Zoltan_ParMetis_Shared should be compiled even when ZOLTAN_PARMETIS 
   is not defined so it can return an error.
 */
static int Zoltan_ParMetis_Shared(
  ZZ *zz,
  float *part_sizes,
  int *num_imp,
  ZOLTAN_ID_PTR *imp_gids,
  ZOLTAN_ID_PTR *imp_lids,
  int **imp_procs,
  int **imp_to_part,
  int *num_exp,
  ZOLTAN_ID_PTR *exp_gids,
  ZOLTAN_ID_PTR *exp_lids,
  int **exp_procs,
  int **exp_to_part,
  int *rank,
  int *iperm,
  ZOOS *order_opt,
  ZOS *order_info
);

#if (defined(ZOLTAN_JOSTLE) || defined(ZOLTAN_PARMETIS))

static int Zoltan_ParMetis_Jostle(ZZ *zz, float *part_sizes,
  int *num_imp, ZOLTAN_ID_PTR *imp_gids,
  ZOLTAN_ID_PTR *imp_lids, int **imp_procs, int **imp_to_part,
  int *num_exp, ZOLTAN_ID_PTR *exp_gids, ZOLTAN_ID_PTR *exp_lids, 
  int **exp_procs, int **exp_to_part, char *alg, int  *options, 
  float *itr, int *rank, int *iperm, ZOOS *order_opt, ZOS *order_info);
static int scale_round_weights(float *fwgts, idxtype *iwgts, int n, int dim, 
                 int mode, int max_wgt_sum, int debug_level, MPI_Comm comm);

#if PARMETIS_MAJOR_VERSION >= 3
static int pmv3method(char *alg);
#endif

#endif  /* (defined(ZOLTAN_JOSTLE) || defined(ZOLTAN_PARMETIS)) */

/**********************************************************/
/* Interface routine for ParMetis. This is just a simple  */
/* wrapper that sets the options etc and then calls       */
/* Zoltan_ParMetis_Jostle, where the real action is.      */
/**********************************************************/

int Zoltan_ParMetis(
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
  return Zoltan_ParMetis_Shared(zz, part_sizes, num_imp, imp_gids, imp_lids,
         imp_procs, imp_to_part, num_exp, exp_gids, exp_lids, exp_procs, 
         exp_to_part, NULL, NULL, NULL, NULL);
}

/* Zoltan_ParMetis_Shared is shared by Zoltan_ParMetis
   and Zoltan_ParMetis_Order */

static int Zoltan_ParMetis_Shared(
  ZZ *zz,               /* Zoltan structure */
  float *part_sizes,    /* Input:  Array of size zz->Num_Global_Parts
                           containing the percentage of work to be
                           assigned to each partition.               */
  int *num_imp,         /* number of objects to be imported */
  ZOLTAN_ID_PTR *imp_gids,  /* global ids of objects to be imported */
  ZOLTAN_ID_PTR *imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,      /* list of processors to import from */
  int **imp_to_part,    /* list of partitions to import to */
  int *num_exp,         /* number of objects to be exported */
  ZOLTAN_ID_PTR *exp_gids,  /* global ids of objects to be exported */
  ZOLTAN_ID_PTR *exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs,      /* list of processors to export to */
  int **exp_to_part,    /* list of partitions to export to */
  int *rank,            /* rank[i] is the rank of gids[i] */
  int *iperm,           /* inverse permutation of rank */
  ZOOS *order_opt,	/* ordering options */
  ZOS *order_info	/* ordering info */
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
  float itr = 0.0;
  double pmv3_itr = 0.0;

  /* Set parameters */
#if PARMETIS_MAJOR_VERSION >= 3
  strcpy(alg, "ADAPTIVEREPART");
  pmv3_itr = 100.; /* Ratio of inter-proc comm. time to data redist. time;
                      100 gives similar partition quality to GDiffusion */
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
  use_obj_size = 1;
  fold = 0;
  seed = GLOBAL_SEED;

  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_METHOD",     
                    (void *) alg);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_OUTPUT_LEVEL", 
                    (void *) &output_level);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_SEED", 
                    (void *) &seed);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_ITR",       
                    (void *) &pmv3_itr);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_USE_OBJ_SIZE",       
                    (void *) &use_obj_size);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_COARSE_ALG", 
                    (void *) &coarse_alg);
  Zoltan_Bind_Param(Parmetis_params, "PARMETIS_FOLD",       
                    (void *) &fold);

  Zoltan_Assign_Param_Vals(zz->Params, Parmetis_params, zz->Debug_Level, 
                           zz->Proc, zz->Debug_Proc);

  /* Copy option values to ParMetis options array */

#if PARMETIS_MAJOR_VERSION >= 3
  /* In this version of Zoltan, processors and partitions are coupled. */
  /* This will likely change in future releases, and then the options  */
  /* value should change to DISCOUPLED.                                */
  options[PMV3_OPTION_PSR] = COUPLED; 

  if (pmv3method(alg)){
    /* ParMetis 3.0 options */
    options[PMV3_OPTION_DBGLVL] = output_level; 
    options[PMV3_OPTION_SEED] = seed; 
    options[PMV3_OPT_USE_OBJ_SIZE] = use_obj_size; 
    itr = pmv3_itr;
  }
  else
#endif
  {
    /* ParMetis 2.0 options */
    options[OPTION_IPART] = coarse_alg; 
    options[OPTION_FOLDF] = fold; 
    options[OPTION_DBGLVL] = output_level; 
  }

  /* If ordering, use ordering method instead of load-balancing method */
  if (order_opt && order_opt->method){
    strcpy(alg, order_opt->method);
  }

  /* Call the real ParMetis interface */
  return Zoltan_ParMetis_Jostle(zz, part_sizes,
            num_imp, imp_gids, imp_lids, imp_procs, imp_to_part, 
            num_exp, exp_gids, exp_lids, exp_procs, exp_to_part,
            alg, options, &itr, rank, iperm, order_opt, order_info);

#endif /* ZOLTAN_PARMETIS */
}

#if PARMETIS_MAJOR_VERSION >= 3
static int pmv3method( char *alg)
{
   /* Check if alg is a supported ParMetis 3.0 method */
   return ((!strcmp(alg, "PARTKWAY")) 
           || (!strcmp(alg, "PARTGEOMKWAY"))
           || (!strcmp(alg, "ADAPTIVEREPART")) 
           || (!strcmp(alg, "REFINEKWAY"))
           || (!strcmp(alg, "NODEND"))
          );
}
#endif

/***************************************************************************
 *  The ParMetis ordering routine piggy-backs on the ParMetis 
 *  partitioning routines. 
 **************************************************************************/

int Zoltan_ParMetis_Order(
  ZZ *zz,               /* Zoltan structure */
  int num_obj,		/* Number of (local) objects to order. */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
                        /* The application must allocate enough space */
  ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
                        /* The application must allocate enough space */
  int *rank,		/* rank[i] is the rank of gids[i] */
  int *iperm,		/* inverse permutation of rank */
  ZOOS *order_opt, 	/* Ordering options, parsed by Zoltan_Order */
  ZOS *order_info       /* Ordering info for this particular ordering */
)
{
  static char *yo = "Zoltan_ParMetis_Order";
  int n, ierr;

  if (!order_opt){
    /* If for some reason order_opt is NULL, allocate a new ZOOS here. */
    /* This should really never happen. */
    order_opt = (ZOOS *) ZOLTAN_MALLOC(sizeof(ZOOS));
    strcpy(order_opt->method,"PARMETIS");
  }

  /* ParMetis only computes the rank vector */
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
   
  /* Call ParMetis_Shared */
  return Zoltan_ParMetis_Shared(zz, NULL, NULL, &gids, &lids, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, rank, iperm, order_opt, order_info);
}

/**********************************************************/
/* Interface routine for Jostle. This is just a simple    */
/* wrapper that sets the options etc and then calls       */
/* Zoltan_ParMetis_Jostle, where the real action is.      */
/**********************************************************/

int Zoltan_Jostle(
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
  sprintf(str, "imbalance = %3d ", (int)(100*(zz->LB.Imbalance_Tol[0] - 1)));
  jostle_env(str);

  /* Multidimensional vertex weights */
  if (zz->Obj_Weight_Dim > 1){
    sprintf(str, "ntypes = %3d ", zz->Obj_Weight_Dim);
    jostle_env(str);
  }

  /* Call the real Jostle/ParMetis interface */
  return Zoltan_ParMetis_Jostle(zz, part_sizes,
            num_imp, imp_gids, imp_lids, imp_procs, imp_to_part,
            num_exp, exp_gids, exp_lids, exp_procs, exp_to_part,
            alg, &output_level, NULL, NULL, NULL, NULL, NULL);

#endif /* ZOLTAN_JOSTLE */
}

/****************************************************************/
/* This routine constructs the required graph data structures   */
/* and calls ParMetis or Jostle.                                */
/*                                                              */
/* Author: Erik Boman, eboman@cs.sandia.gov (9226)              */
/****************************************************************/

#if (defined(ZOLTAN_JOSTLE) || defined(ZOLTAN_PARMETIS))

static int Zoltan_ParMetis_Jostle(
  ZZ *zz,               /* Zoltan structure */
  float *part_sizes,    /* partition sizes */
  int *num_imp,         /* number of objects to be imported */
  ZOLTAN_ID_PTR *imp_gids,  /* global ids of objects to be imported */
  ZOLTAN_ID_PTR *imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,      /* list of processors to import from */
  int **imp_to_part,    /* list of partitions to import to */
  int *num_exp,         /* number of objects to be exported */
  ZOLTAN_ID_PTR *exp_gids,  /* global ids of objects to be exported */
  ZOLTAN_ID_PTR *exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs,      /* list of processors to export to */
  int **exp_to_part,    /* list of partitions to export to */
  char *alg,            /* algorithm to use */
  int  *options,        /* ParMetis option array */
  float *itr,    	/* ParMetis 3.0 parameter for adaptive repart. */
  int *rank,            /* Ordering only: rank[i] is the rank of gids[i] */
  int *iperm,           /* Ordering only: inverse permutation of rank */
  ZOOS *order_opt, 	/* Ordering only: options */
  ZOS *order_info	/* Ordering only: Zoltan ordering struct */
)
{
  static char *yo = "Zoltan_ParMetis_Jostle";
  int i, j, k, ierr, tmp, flag, ndims;
  int obj_wgt_dim, edge_wgt_dim, check_graph, scatter;
  int num_obj=0, num_edges, edgecut;
  int nsend, wgtflag, numflag, graph_type; 
  int get_graph_data, get_geom_data, get_times; 
  int compute_only_part_changes=0; /* EBEB: Read parameter when implemented. */
  idxtype *vtxdist, *xadj, *adjncy, *vwgt, *adjwgt, *part, *vsize;
  idxtype *sep_sizes, *part_orig;
  int num_obj_orig, ncon, start_index, compute_order=0;
  float *float_vwgt, *ewgts, *xyz, *imb_tols; 
  double *geom_vec;
  ZOLTAN_ID_PTR local_ids;
  ZOLTAN_ID_PTR global_ids;    
#ifdef PARMETIS_V3_1_MEMORY_ERROR_FIXED
#if (PARMETIS_MAJOR_VERSION >= 3) 
  ZOLTAN_ID_PTR lid;        /* Temporary pointer to a local id; used to pass
                               NULL to query fns when NUM_LID_ENTRIES == 0. */
#endif
#endif /* PARMETIS_V3_1_MEMORY_ERROR_FIXED */
  int *input_parts;         /* Initial partitions for objects. */
  int *newproc;             /* New processor for each object. */
  ZOLTAN_COMM_OBJ *comm_plan;
  double times[5];
  char msg[256];
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int num_proc = zz->Num_Proc;     /* Temporary variables whose addresses are*/
  int num_part = zz->LB.Num_Global_Parts;/* passed to Jostle/ParMETIS. Don't */
  int new_map;          /* flag indicating whether parts were remapped */
#ifdef ZOLTAN_PARMETIS
  MPI_Comm comm = zz->Communicator;/* want to risk letting external packages */
                                   /* change our zz struct.                  */
#endif /* ZOLTAN_PARMETIS */
  int i99;                         /* Variables used for debugging.          */
#ifdef ZOLTAN_JOSTLE
  int nnodes;
  int network[4] = {0, 1, 1, 1};
#endif

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Set default return values (in case of early exit) */
  if (num_exp) *num_exp = -1;
  if (num_imp) *num_imp = -1; /* No import data */

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  vtxdist = xadj = adjncy = vwgt = adjwgt = part = NULL;
  vsize = sep_sizes = NULL;
  float_vwgt = ewgts = xyz = imb_tols = NULL;
  geom_vec = NULL;
  local_ids = NULL;
  global_ids = NULL;
  input_parts = part_orig = NULL;
  newproc = NULL;

  /* Start timer */
  get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
  if (get_times){
    MPI_Barrier(zz->Communicator);
    times[0] = Zoltan_Time(zz->Timer);
  }

  /* Check for outdated/unsupported ParMetis versions. */
#if (PARMETIS_MAJOR_VERSION < 3) 
  if (zz->Proc == 0)
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "ParMetis 2.0 is obsolete. Zoltan currently works with this version, but please upgrade to ParMetis 3.1 (or later) soon.");
  ierr = ZOLTAN_WARN;
#elif (PARMETIS_MAJOR_VERSION == 3) && (PARMETIS_MINOR_VERSION == 0)
  if (zz->Proc == 0)
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "ParMetis 3.0 is no longer supported by Zoltan. Please upgrade to ParMetis 3.1 (or later).");
  ierr = ZOLTAN_WARN;
#endif

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
    sprintf(msg, "Edge weight dimension is %d, "
            "but should be >= 0. Using Edge_Weight_Dim = 0.",
            zz->Edge_Weight_Dim);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    edge_wgt_dim = 0;
  }
  else if (zz->Edge_Weight_Dim>1){
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "This method does not support "
        "multidimensional edge weights. Using Edge_Weight_Dim = 1.");
    edge_wgt_dim = 1;
  }
  else {
    edge_wgt_dim = zz->Edge_Weight_Dim;
  }

  /* Default graph type is GLOBAL. */
  graph_type = GLOBAL_GRAPH;

  /* Get parameter options shared by ParMetis and Jostle */
  check_graph = 1;          /* default */
  scatter = 1;              /* default */
  Zoltan_Bind_Param(Graph_params, "CHECK_GRAPH", (void *) &check_graph);
  Zoltan_Bind_Param(Graph_params, "SCATTER_GRAPH", (void *) &scatter);
  Zoltan_Assign_Param_Vals(zz->Params, Graph_params, zz->Debug_Level, zz->Proc,
                       zz->Debug_Proc);

  /* Are we doing ordering or partitioning? */
  /* Note: If ORDER_METHOD=PARMETIS, then PARMETIS_METHOD=NODEND */
  if (strcmp(alg, "NODEND") == 0)
    compute_order = 1;
  else 
    compute_order = 0;

  if (compute_order){
    /* Do not use weights for ordering */
    obj_wgt_dim = 0;
    edge_wgt_dim = 0;

    /* Check what ordering type is requested */
    if (order_opt){
       if (strcmp(order_opt->order_type, "LOCAL") == 0)
          graph_type = LOCAL_GRAPH;
       else if (strcmp(order_opt->order_type, "GLOBAL") == 0)
          graph_type = GLOBAL_GRAPH;
       else
          graph_type = NO_GRAPH;
    }

    /* Allocate space for separator sizes */
    sep_sizes = (idxtype *) ZOLTAN_MALLOC(2*num_part*sizeof(idxtype));
    if (!sep_sizes){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    printf("[%1d] Debug: alg=%s, Obj_Weight_Dim=%d, Edge_Weight_Dim=%d\n", 
      zz->Proc, alg, obj_wgt_dim, edge_wgt_dim);
    printf("[%1d] Debug: ParMetis options = %d, %d, %d, %d\n", zz->Proc,
      options[0], options[1], options[2], options[3]);
  }
  
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
    graph_type = NO_GRAPH;
  }

  /* For ordering, the ids are passed in through the *imp_ids
     parameters. These id lists are only populated if
     reorder=True, otherwise we need to initialize them. 
  */
  if (compute_order){
    global_ids = *imp_gids;
    local_ids =  *imp_lids;
  }

  /* If reorder is true, we already have the id lists. Ignore weights. */
  if (!(order_opt && order_opt->reorder)){
    ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids,
                               obj_wgt_dim, &float_vwgt, &input_parts);
    if (ierr){
      /* Return error */
      ZOLTAN_PARMETIS_ERROR(ierr, "Get_Obj_List returned error.");
    }
#if (PARMETIS_MAJOR_VERSION >= 3) && (PARMETIS_MINOR_VERSION == 0)
    /* Special error checks to avoid incorrect results from ParMetis 3.0.
     * ParMETIS 3.0 Partkway ignores partition sizes for problems with 
     * less than 10000 objects.
     */
    if (!strcmp(alg, "PARTKWAY") && !(zz->LB.Uniform_Parts) 
                                 && (zz->Obj_Weight_Dim <= 1)) {
      int gsum;
      MPI_Allreduce(&num_obj, &gsum, 1, MPI_INT, MPI_SUM, comm);
      if (gsum < 10000) {
        char str[256];
        sprintf(str, "Total objects %d < 10000 causes ParMETIS 3.0 PARTKWAY "
                "to ignore partition sizes; uniform partition sizes will be "
                "produced. Please try a different load-balancing method.\n", 
                gsum);
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, str);
      }
    }
#endif

  }
  
  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    printf("[%1d] Debug: num_obj =%d\n", zz->Proc, num_obj);
    printf("[%1d] Debug: Global ids = ", zz->Proc);
    for (i99=0; i99<num_obj; i99++) {
      printf("    ");
      ZOLTAN_PRINT_GID(zz, &(global_ids[i99*num_gid_entries]));
      printf("\n");
    }
  }
  
  /* Build ParMetis data structures, or just get vtxdist. */
  ierr = Zoltan_Build_Graph(zz, graph_type, check_graph, num_obj,
         global_ids, local_ids, obj_wgt_dim, edge_wgt_dim,
         &vtxdist, &xadj, &adjncy, &ewgts);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      ZOLTAN_PARMETIS_ERROR(ierr, "Zoltan_Build_Graph returned error.");
  }
  part = (idxtype *)ZOLTAN_MALLOC((num_obj+1) * sizeof(idxtype));
  if (!part){
    /* Not enough memory */
    ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
  }
  /* Copy input_parts array to part, in case ParMetis needs it. */
  for (i=0; i<num_obj; i++) 
    part[i] = input_parts[i];

  /* Special error checks to avoid certain death in ParMETIS 3.0.
   * AdaptiveRepart uses input partition number to index into an array 
   * of size num_part.  Thus, it segfaults if an input 
   * partition number >= num_part. 
   * Problem fixed in ParMETIS 3.1.
   */
#if (PARMETIS_MAJOR_VERSION == 3) && (PARMETIS_MINOR_VERSION == 0)
  if (strcmp(alg, "ADAPTIVEREPART") == 0) {
    int gmax, maxpart = -1;
    for (i = 0; i < num_obj; i++)
      if (part[i] > maxpart) maxpart = part[i];
    MPI_Allreduce(&maxpart, &gmax, 1, MPI_INT, MPI_MAX, zz->Communicator);
    if (gmax >= num_part) {
      sprintf(msg, "Partition number %1d >= number of partitions %1d.\n"
        "ParMETIS 3.0 with %s will fail, please upgrade to 3.1 or later.",
        gmax, num_part, alg);
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, msg);
    }
  }
#endif

  /* ParMetis 2 and 3 require integer weights. Convert from float. */

    /* Get vertex weights if needed */
    if (obj_wgt_dim){
      vwgt = (idxtype *)ZOLTAN_MALLOC(obj_wgt_dim*num_obj
                          * sizeof(idxtype));
      ierr = scale_round_weights(float_vwgt, vwgt, num_obj, obj_wgt_dim, 1,
                                 (int) MAX_WGT_SUM, zz->Debug_Level, zz->Communicator);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        /* Return error code */
        ZOLTAN_PARMETIS_ERROR(ierr, "Error in scaling of weights.");
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
    if (get_graph_data)
       num_edges = xadj[num_obj];
    else
       num_edges = 0;

    if (edge_wgt_dim){
      adjwgt = (idxtype *)ZOLTAN_MALLOC(edge_wgt_dim * num_edges
                 * sizeof(idxtype));
      if (num_edges && (!adjwgt)){
        /* Not enough memory */
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
      }
  

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] Debug: Edge weights are (before scaling) = \n", zz->Proc);        for (j=0; j<num_edges*edge_wgt_dim; j++)
          printf("%f ", ewgts[j]);
        printf("\n");
      }
      /* Scale and round edge weights to integers */
      ierr = scale_round_weights(ewgts, adjwgt, num_edges, edge_wgt_dim, 1,
                                 (int) MAX_WGT_SUM, zz->Debug_Level, zz->Communicator);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        /* Return error code */
        ZOLTAN_PARMETIS_ERROR(ierr, "Error in scaling of weights.");
      }

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] Debug: Edge weights are (after scaling) = \n", zz->Proc);
        for (j=0; j<num_edges; j++)
          printf("%d ", adjwgt[j]);

        printf("\n");
      }
      ZOLTAN_FREE(&ewgts);
    }

  if (get_geom_data){

    /* Get coordinate information */
    ierr = Zoltan_Get_Coordinates(zz, num_obj, global_ids, local_ids, 
                                  &ndims, &geom_vec);
    if (ierr) {
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, 
        "Error returned from Zoltan_Get_Coordinates");
    }
    /* Convert geometry info from double to float for ParMETIS */
    if (num_obj && ndims) {
      xyz = (float *) ZOLTAN_MALLOC(num_obj * ndims * sizeof(float));
      if (xyz == NULL)  {
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Memory error.");
      }
      for (i = 0; i < num_obj * ndims; i++) xyz[i] = (float) geom_vec[i];
      ZOLTAN_FREE(&geom_vec);
    }
  }

#ifdef PARMETIS_V3_1_MEMORY_ERROR_FIXED
#if (PARMETIS_MAJOR_VERSION >= 3) 
  /* Get object sizes if requested */ 
  if (options[PMV3_OPT_USE_OBJ_SIZE] && 
      (zz->Get_Obj_Size || zz->Get_Obj_Size_Multi)) {
    vsize = (idxtype *) ZOLTAN_MALLOC(num_obj*sizeof(idxtype));
    if (num_obj && !vsize){
      /* Not enough space */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    if (zz->Get_Obj_Size_Multi) {
      zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data, 
                       num_gid_entries, num_lid_entries, num_obj,
                       global_ids, local_ids, vsize, &ierr);
    }
    else {
      for (i=0; i<num_obj; i++){
        lid = (num_lid_entries ? &(local_ids[i*num_lid_entries]) : NULL);
        vsize[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data, 
                       num_gid_entries, num_lid_entries,
                       &(global_ids[i*num_gid_entries]), 
                       lid, &ierr);
      }
    }
  }
#endif
#endif /* PARMETIS_V3_1_MEMORY_ERROR_FIXED */

  if (get_geom_data){
    /* ParMETIS will crash if geometric method and some procs have no nodes. */
    /* Avoid fatal crash by setting scatter to level 2 or higher. */
    if (scatter<2) scatter = 2;
  }

  /* Scatter graph? 
   * If data distribution is highly imbalanced, it is better to
   * redistribute the graph data structure before calling ParMetis.
   * After partitioning, the results must be mapped back.
   */
  if ((scatter>0) && (scatter<3)){
    /* Decide if the data imbalance is so bad that we should scatter the graph. */
    /* scatter==1: Scatter if all the objects are on a single processor.        */
    /* scatter==2: Scatter if any processor has no objects.                     */
    num_obj_orig = num_obj; /* Save value for later. */
    if (num_obj==0)
      j = 1;
    else 
      j = 0;
    MPI_Allreduce(&j, &tmp, 1, MPI_INT, MPI_SUM, zz->Communicator);
    if (scatter == 1){
      if (tmp < num_proc-1)
        scatter = 0;
    }
    else if (scatter==2){
      if (tmp == 0)
        scatter = 0;
    }
  }

  /* We need to make sure we don't scatter the graph 
   * if graph_type = LOCAL_GRAPH, i.e. METIS is used. 
   */
  if (scatter && (graph_type == LOCAL_GRAPH)){
    scatter = 0;
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Setting scatter_graph=0 since the graph"
      " is local on each proc");
    ierr = ZOLTAN_WARN;
  }

  if (scatter){
    ierr = Zoltan_Scatter_Graph(&vtxdist, &xadj, &adjncy, &vwgt, &vsize,
              &adjwgt, &xyz, ndims, zz, &comm_plan);
    if ((ierr == ZOLTAN_FATAL) || (ierr == ZOLTAN_MEMERR)){
      ZOLTAN_PARMETIS_ERROR(ierr, "Zoltan_Scatter_Graph returned error.");
    }
    num_obj = vtxdist[zz->Proc+1]-vtxdist[zz->Proc];
    part_orig = part;
    part = (idxtype *) ZOLTAN_MALLOC((num_obj+1) * sizeof(int));
    Zoltan_Comm_Do(comm_plan, TAG1, (char *) part_orig, sizeof(idxtype), 
                   (char *) part);
  }

  /* Get ready to call ParMETIS or Jostle */
  edgecut = -1; 
  wgtflag = 2*(obj_wgt_dim>0) + (edge_wgt_dim>0); 
  numflag = 0;
  ncon = (obj_wgt_dim > 0 ? obj_wgt_dim : 1);

  /* Set Imbalance Tolerance for each weight component. */
  imb_tols = (float *) ZOLTAN_MALLOC(ncon * sizeof(float));
  if (!imb_tols){
    /* Not enough memory */
    ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
  }
  for (i=0; i<ncon; i++)
    imb_tols[i] = zz->LB.Imbalance_Tol[i];

  /* Verify that graph is correct */
  if (get_graph_data){
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) 
      flag = 2; /* Maximum output level */
    else
      flag = 1; /* Medium output level */
    ierr = Zoltan_Verify_Graph(zz->Communicator, vtxdist, xadj, adjncy, vwgt, 
              adjwgt, obj_wgt_dim, edge_wgt_dim, graph_type, check_graph, flag);
  
    /* Special error checks to avoid certain death in ParMETIS */
    if (xadj[num_obj] == 0){
      /* No edges on a proc is a fatal error in ParMETIS 2.0/3.0
       * and in Jostle 1.2. This error test should be removed
       * when the bugs in ParMETIS and Jostle have been fixed.
       */
      if (strcmp(alg, "JOSTLE") == 0){
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, "No edges on this proc. Jostle "
                              "will likely fail. Please try a different "
                              "load-balancing method.");
      } else { /* ParMETIS */
#if (PARMETIS_MAJOR_VERSION == 2) 
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, "No edges on this proc. "
                              "ParMETIS 2.0 will likely fail. Please "
                              "upgrade to version 3.1 or later.");
#elif (PARMETIS_MAJOR_VERSION == 3) && (PARMETIS_MINOR_VERSION == 0)
        if (strcmp(alg, "ADAPTIVEREPART") == 0)
          ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, "No edges on this proc. "
                                "ParMETIS 3.0 will likely fail with method "
                                "AdaptiveRepart. Please upgrade to 3.1 or later.");
#endif
      }
    }
  }
  
  /* Get a time here */
  if (get_times) times[1] = Zoltan_Time(zz->Timer);

  /* Select the desired ParMetis or Jostle function */

#if (PARMETIS_MAJOR_VERSION >= 3)
  if (!compute_order){ /* ordering does not need partsizes */
    /* Assume partition sizes are given as input. */
    if (!part_sizes){
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL,"Input parameter part_sizes is NULL.");
    }
    if ((zz->Proc == 0) && (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)) {
      for (i=0; i<num_part; i++){
        printf("Debug: Size(s) for partition %1d = ", i);
        for (j=0; j<ncon; j++)
          printf("%f ", part_sizes[i*ncon+j]);
        printf("\n");
      }
    }

    /* Special error checks to avoid certain death in ParMETIS */
    /* ParMETIS 3.x divides by partition sizes so it           */
    /* will produce a FPE error when part_sizes[i] == 0        */

    /* if (strcmp(alg, "ADAPTIVEREPART") == 0) */
      for (i = 0; i < num_part*ncon; i++)
        if (part_sizes[i] == 0) 
          ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, "Zero-sized partition(s) requested! "
            "ParMETIS 3.x will likely fail. Please use a different method, or remove the zero-sized partitions from the problem.");

  }
#else /* PARMETIS_MAJOR_VERSION < 3 */
  if ((ncon >= 2) && strcmp(alg, "JOSTLE")) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "You need ParMETIS version 3.0 or higher to use multi-weights.\n"
      "If you have this installed, please make sure you include the appropriate\n"
      "version of parmetis.h and recompile Zoltan.");
    ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, "Multi-constraint error.");
  }
  /* Error if non-uniform partition sizes are input. */
  if (!(zz->LB.Uniform_Parts)) {
    ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL,
      "You need ParMETIS v3.0 or higher to non-uniform partition sizes.\n"
      "If you have this version installed, please include the appropriate\n"
      "version of parmetis.h and recompile Zoltan.");
  }
#endif 

  if (strcmp(alg, "JOSTLE") == 0){
#ifdef ZOLTAN_JOSTLE
    k = 0;                 /* Index of the first object/node. */
    j = 0;                 /* Dummy variable for Jostle */
    nnodes = vtxdist[zz->Num_Proc]; /* Global number of objects */ 
    ndims = 1;             /* Topological dimension of the computer */
    network[1] = zz->Num_Proc;
    /* Convert xadj array to Jostle format (degree) */
    for (i=0; i<num_obj; i++){
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
    pjostle(&nnodes, &k, &(vtxdist[zz->Proc]), &j, vtxdist, 
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
    ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, "Jostle not available.");
#endif /* ZOLTAN_JOSTLE */
  }
#ifdef ZOLTAN_PARMETIS
#if PARMETIS_MAJOR_VERSION >= 3
  /* First check for ParMetis 3 routines */
  else if (strcmp(alg, "PARTKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
    ParMETIS_V3_PartKway (vtxdist, xadj, adjncy, vwgt, adjwgt, 
      &wgtflag, &numflag, &ncon, &num_part, part_sizes,
      imb_tols, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "PARTGEOMKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
    ParMETIS_V3_PartGeomKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag,
      &numflag, &ndims, xyz, &ncon, &num_part, part_sizes,
      imb_tols, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "PARTGEOM") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
    ParMETIS_V3_PartGeom (vtxdist, &ndims, xyz, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "ADAPTIVEREPART") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
    ParMETIS_V3_AdaptiveRepart (vtxdist, xadj, adjncy, vwgt, vsize, adjwgt, 
      &wgtflag, &numflag, &ncon, &num_part, part_sizes, imb_tols, 
      itr, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "REFINEKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
    ParMETIS_V3_RefineKway (vtxdist, xadj, adjncy, vwgt, adjwgt, 
      &wgtflag, &numflag, &ncon, &num_part, part_sizes, imb_tols,
      options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "NODEND") == 0){
    if (graph_type==GLOBAL_GRAPH){
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");
      ParMETIS_V3_NodeND (vtxdist, xadj, adjncy, 
        &numflag, options, part, sep_sizes, &comm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
    }
    else {
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the METIS library");
      options[0] = 0;  /* Use default options for METIS. */
      METIS_NodeND (&num_obj, xadj, adjncy, 
        &numflag, options, part, iperm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the METIS library");
    }
  }
#endif  /* PARMETIS_MAJOR_VERSION >= 3 */
  /* Check for ParMetis 2.0 routines */
  else if (strcmp(alg, "PARTKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_PartKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag, 
      &numflag, &num_part, options, &edgecut, part, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
  }
  else if (strcmp(alg, "PARTGEOMKWAY") == 0){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
    ParMETIS_PartGeomKway (vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag,
      &numflag, &ndims, xyz, &num_part, options, &edgecut, 
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
  else if (strcmp(alg, "NODEND") == 0){
    if (graph_type==GLOBAL_GRAPH){
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 2 library");
      ParMETIS_NodeND (vtxdist, xadj, adjncy, 
        &numflag, options, part, sep_sizes, &comm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
    }
    else {
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the METIS library");
      options[0] = 0;  /* Use default options for METIS. */
      METIS_NodeND (&num_obj, xadj, adjncy, 
        &numflag, options, part, iperm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the METIS library");
    }
  }
#endif /* ZOLTAN_PARMETIS */
  else {
    /* Sanity check: This should never happen! */
    sprintf(msg, "Unknown ParMetis or Jostle algorithm %s.", alg);
    ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, msg);
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
  if (edge_wgt_dim){
    ZOLTAN_FREE(&adjwgt);
  }
  /* Also free temp arrays for ParMetis 3.0 */
  ZOLTAN_FREE(&imb_tols);
  if (vsize) ZOLTAN_FREE(&vsize);

  /* If we have been using a scattered graph, convert partition result back to 
   * original distribution 
   */
  if (scatter){
    num_obj = num_obj_orig; /* Restore value from original distribution. */

    /* Use reverse communication to compute the partition array under the 
     * original distribution.
     */
    ierr = Zoltan_Comm_Do_Reverse(comm_plan, TAG2, (char *) part, 
           sizeof(idxtype), NULL, (char *) part_orig);
    if ((ierr == ZOLTAN_FATAL) || (ierr == ZOLTAN_MEMERR)){
      ZOLTAN_PARMETIS_ERROR(ierr, "Zoltan_Comm_Do_Reverse returned error.");
    }
    Zoltan_Comm_Destroy(&comm_plan); /* Destroy the comm. plan */
    /* We don't need the partition array with the scattered distribution 
     * any more */
    ZOLTAN_FREE(&part); 
    /* part is now the new partition array under the original distribution */
    part = part_orig;   
    part_orig = NULL;
  }
 
  if (compute_order){
    /* Ordering */
    /* ParMetis produces the rank vector in Zoltan lingo */

    /* Check if start_index != 0 */
    if (order_opt && order_opt->start_index)
       start_index = order_opt->start_index;
    else
       start_index = 0;
    /* Copy from 'part' array to rank */
    if (rank != NULL){
      for (i=0; i<num_obj; i++){
        rank[i] = part[i] + start_index;
      }
    }
    else {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "rank is NULL, no data returned");
      ierr = ZOLTAN_WARN;
    }
    /* If we did local ordering via METIS, then we also have the inv. perm. */
    if (graph_type == LOCAL_GRAPH){
      for (i=0; i<num_obj; i++){
        iperm[i] = iperm[i] + start_index;
      }
      /* EBEB: Return parameter that says we have computed both return args? */
    }

    /* Fill in the Zoltan Order Struct */
    /* EBEB: For now, discard separator info */ 
    ZOLTAN_FREE(&sep_sizes); 
  }
  else{
    /* Partitioning */
    /* Determine new processor and number of objects to export */
    nsend = 0;
    newproc = (int *) ZOLTAN_MALLOC(num_obj * sizeof(int));
    if (num_obj && !newproc){
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory. ");
    }
    for (i=0; i<num_obj; i++){
      newproc[i] = Zoltan_LB_Part_To_Proc(zz, part[i],
                                          &(global_ids[i*num_gid_entries]));
      if (newproc[i]<0){
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL, 
         "Zoltan_LB_Part_To_Proc returned invalid processor number.");
      }
    }
  
    if (zz->LB.Remap_Flag) {
      ierr = Zoltan_LB_Remap(zz, &new_map, num_obj, newproc, input_parts, 
                             part, 1);
      if (ierr < 0) {
        ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL,
                              "Error returned from Zoltan_LB_Remap");
      }
    }

    for (i=0; i<num_obj; i++){
      if ((part[i] != input_parts[i]) || ((!compute_only_part_changes) && 
                                          (newproc[i] != zz->Proc))) 
        nsend++;
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
        printf("[%1d] DEBUG: local object %1d: old part = %1d, new part = %1d\n", 
        zz->Proc, i, input_parts[i], part[i]);
    }
  
    /* Create export lists */
    if (zz->LB.Return_Lists){
      (*num_exp) = nsend;
      if (nsend > 0) {
        if (!Zoltan_Special_Malloc(zz,(void **)exp_gids,nsend,ZOLTAN_SPECIAL_MALLOC_GID)) {
          ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
        }
        if (!Zoltan_Special_Malloc(zz,(void **)exp_lids,nsend,ZOLTAN_SPECIAL_MALLOC_LID)) {
          Zoltan_Special_Free(zz,(void **)exp_gids,ZOLTAN_SPECIAL_MALLOC_GID);
          ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
        }
        if (!Zoltan_Special_Malloc(zz,(void **)exp_procs,nsend,ZOLTAN_SPECIAL_MALLOC_INT)) {
          Zoltan_Special_Free(zz,(void **)exp_lids,ZOLTAN_SPECIAL_MALLOC_LID);
          Zoltan_Special_Free(zz,(void **)exp_gids,ZOLTAN_SPECIAL_MALLOC_GID);
          ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
        }
        if (!Zoltan_Special_Malloc(zz,(void **)exp_to_part,nsend,ZOLTAN_SPECIAL_MALLOC_INT)) {
          Zoltan_Special_Free(zz,(void **)exp_lids,ZOLTAN_SPECIAL_MALLOC_LID);
          Zoltan_Special_Free(zz,(void **)exp_gids,ZOLTAN_SPECIAL_MALLOC_GID);
          Zoltan_Special_Free(zz,(void **)exp_procs,ZOLTAN_SPECIAL_MALLOC_INT);
          ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Not enough memory.");
        }
        j = 0;
        for (i=0; i<num_obj; i++){
          if ((part[i] != input_parts[i]) || ((!compute_only_part_changes) 
               && (newproc[i] != zz->Proc))){ 
            /* Object should move to new partition or processor */
            ZOLTAN_SET_GID(zz, &((*exp_gids)[j*num_gid_entries]),
                           &(global_ids[i*num_gid_entries]));
            if (num_lid_entries)
              ZOLTAN_SET_LID(zz, &((*exp_lids)[j*num_lid_entries]),
                             &(local_ids[i*num_lid_entries]));
            (*exp_to_part)[j] = part[i];  
            (*exp_procs)[j] = newproc[i];
            /* printf("[%1d] Debug: Move object %1d to part %1d, proc %1d\n",
                zz->Proc, i, part[i], newproc[i]); */
            j++;
          }
        }
      }
    }
  }

  /* Get a time here */
  if (get_times) times[3] = Zoltan_Time(zz->Timer);

  /* Output timing results if desired */
  if (get_times){
    if (zz->Proc == zz->Debug_Proc) printf("\nZOLTAN timing statistics:\n");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[1]-times[0], 
                   " Partitioner Pre-processing time  ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[2]-times[1], 
                   " Partitioner Library time         ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[3]-times[2], 
                   " Partitioner Post-processing time ");
    Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[3]-times[0], 
                   " Partitioner Total time           ");
    if (zz->Proc==zz->Debug_Proc) printf("\n");
  }

  /* Successful finish */
  ierr = ZOLTAN_OK;

End:
  /* If an error was encountered, the following data may need to be freed */
  if (vwgt)      ZOLTAN_FREE(&vwgt); 
  if (adjwgt)    ZOLTAN_FREE(&adjwgt); 
  if (float_vwgt)ZOLTAN_FREE(&float_vwgt); 
  if (ewgts)     ZOLTAN_FREE(&ewgts); 
  if (imb_tols)  ZOLTAN_FREE(&imb_tols);
  if (vsize)     ZOLTAN_FREE(&vsize); 
  if (sep_sizes) ZOLTAN_FREE(&sep_sizes);
  if (newproc)   ZOLTAN_FREE(&newproc);
  if (part_orig) ZOLTAN_FREE(&part_orig);

  /* Free local_ids and global_ids if they were allocated here */
  if (!compute_order){
    ZOLTAN_FREE(&local_ids);
    ZOLTAN_FREE(&global_ids);
  }

  /* Always free these arrays */
  ZOLTAN_FREE(&vtxdist);
  ZOLTAN_FREE(&xadj);
  ZOLTAN_FREE(&adjncy);
  ZOLTAN_FREE(&xyz);
  ZOLTAN_FREE(&part);
  ZOLTAN_FREE(&input_parts);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
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
 * Note that we use ceil() instead of round() to avoid
 * rounding to zero weights.
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
      iwgts[i] = (int) ceil((double) fwgts[i]);
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
      for (i=0; i<n; i++){
        for (j=0; j<dim; j++){
          if (!nonint_local[j]){ 
            /* tmp = (int) roundf(fwgts[i]);  EB: Valid C99, but not C89 */
            tmp = (int) floor((double) fwgts[i] + .5); /* Nearest int */
            if (fabs((double)tmp-fwgts[i*dim+j]) > INT_EPSILON){
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
        if (nonint[j] || (max_wgt[j] <= INT_EPSILON) || (sum_wgt[j] > max_wgt_sum)){
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
          iwgts[i*dim+j] = (int) ceil((double) fwgts[i*dim+j]*scale[j]);
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
        "REFINEKWAY", "ADAPTIVEREPART",
        "NODEND", /* for nested dissection ordering */
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


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
