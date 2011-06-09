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
#include "zz_rand.h"
#include "params_const.h"
#include "all_allo_const.h"
#include "hier.h"
#include "zoltan_comm.h"
#include "key_params.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>

static PARAM_VARS graph_type_params[] = {
  { "GRAPH_PACKAGE", NULL, "STRING", 0 },
  { "ORDER_TYPE", NULL, "STRING", 0 },
  { "GRAPH_SYMMETRIZE", NULL, "STRING", 0 },
  { "GRAPH_SYM_WEIGHT", NULL, "STRING", 0 },
  { "GRAPH_BIPARTITE_TYPE", NULL, "STRING", 0},
  { "GRAPH_BUILD_TYPE", NULL, "STRING"},       
  { "GRAPH_FAST_BUILD_BASE", NULL, "INTEGER"},       
  { "ORDER_METHOD", NULL, "STRING", 0 },
  { "USE_ORDER_INFO", NULL, "INT", 0 },
  {"HYPERGRAPH_PACKAGE",              NULL,  "STRING", 0},
  {"PHG_MULTILEVEL",                  NULL,  "INT", 0},
  {"PHG_CUT_OBJECTIVE",               NULL,  "STRING", 0},
  {"PHG_OUTPUT_LEVEL",                NULL,  "INT",    0},
  {"CHECK_GRAPH",                     NULL,  "INT",    0},
  {"CHECK_HYPERGRAPH",                NULL,  "INT",    0},
  {"PHG_NPROC_VERTEX",                NULL,  "INT",    0},
  {"PHG_NPROC_EDGE",                  NULL,  "INT",    0},
  {"PHG_COARSENING_LIMIT",            NULL,  "INT",    0},
  {"PHG_COARSENING_NCANDIDATE",       NULL,  "INT",    0},
  {"PHG_COARSENING_METHOD",           NULL,  "STRING", 0},
  {"PHG_COARSENING_METHOD_FAST",      NULL,  "STRING", 0},
  {"PHG_VERTEX_VISIT_ORDER",          NULL,  "INT",    0},
  {"PHG_EDGE_SCALING",                NULL,  "INT",    0},
  {"PHG_VERTEX_SCALING",              NULL,  "INT",    0},
  {"PHG_COARSEPARTITION_METHOD",      NULL,  "STRING", 0},
  {"PHG_REFINEMENT_METHOD",           NULL,  "STRING", 0},
  {"PHG_DIRECT_KWAY",                 NULL,  "INT",    0},
  {"PHG_REFINEMENT_LOOP_LIMIT",       NULL,  "INT",    0},
  {"PHG_REFINEMENT_MAX_NEG_MOVE",     NULL,  "INT",    0},
  {"PHG_REFINEMENT_QUALITY",          NULL,  "FLOAT",  0},
  {"PHG_USE_TIMERS",                  NULL,  "INT",    0},
  {"USE_TIMERS",                      NULL,  "INT",    0},
  {"PHG_EDGE_SIZE_THRESHOLD",         NULL,  "FLOAT",  0},
  {"PHG_MATCH_EDGE_SIZE_THRESHOLD",   NULL,  "INT",    0},
  {"PHG_BAL_TOL_ADJUSTMENT",          NULL,  "FLOAT",  0},
  {"PHG_EDGE_WEIGHT_OPERATION",       NULL,  "STRING",  0},
  {"PARKWAY_SERPART",                 NULL,  "STRING", 0},
  {"PHG_RANDOMIZE_INPUT",             NULL,  "INT",    0},
  {"PHG_PROCESSOR_REDUCTION_LIMIT",   NULL,  "FLOAT",  0},
  {"PHG_REPART_MULTIPLIER",           NULL,  "FLOAT",  0},
  {"PATOH_ALLOC_POOL0",               NULL,  "INT",    0},
  {"PATOH_ALLOC_POOL1",               NULL,  "INT",    0},
  {"PHG_KEEP_TREE",                   NULL,  "INT",    0},
  { "PARMETIS_METHOD", NULL, "STRING", 0 },
  { "PARMETIS_OUTPUT_LEVEL", NULL, "INT", 0 },
  { "PARMETIS_SEED", NULL, "INT", 0 },
  { "PARMETIS_ITR", NULL, "DOUBLE", 0 },
  { "PARMETIS_COARSE_ALG", NULL, "INT", 0 },
  { "PARMETIS_FOLD", NULL, "INT", 0 },
  { NULL, NULL, NULL, 0 } };

static PARAM_VARS geometric_type_params[] = {
        { "RCB_OVERALLOC", NULL, "DOUBLE", 0 },
        { "CHECK_GEOM", NULL, "INT", 0 },
        { "RCB_OUTPUT_LEVEL", NULL, "INT", 0 },
        { "RCB_LOCK_DIRECTIONS", NULL, "INT", 0 },
        { "RCB_SET_DIRECTIONS", NULL, "INT", 0 },
        { "RCB_RECTILINEAR_BLOCKS", NULL, "INT", 0 },
        { "OBJ_WEIGHTS_COMPARABLE", NULL, "INT", 0 },
        { "RCB_MULTICRITERIA_NORM", NULL, "INT", 0 },
        { "RCB_MAX_ASPECT_RATIO", NULL, "DOUBLE", 0 },
        { "AVERAGE_CUTS", NULL, "INT", 0 },
        { "RANDOM_PIVOTS", NULL, "INT", 0 },
        { "RCB_RECOMPUTE_BOX", NULL, "INT", 0 },
        { "REDUCE_DIMENSIONS", NULL, "INT", 0 },
        { "DEGENERATE_RATIO", NULL, "DOUBLE", 0 },
        { "RIB_OVERALLOC", NULL, "DOUBLE", 0 },
        { "RIB_OUTPUT_LEVEL", NULL, "INT", 0 },
        { "AVERAGE_CUTS", NULL, "INT", 0 },
        { "REDUCE_DIMENSIONS", NULL, "INT", 0 },
        { "DEGENERATE_RATIO", NULL, "DOUBLE", 0 },
        {"FINAL_OUTPUT", NULL,  "INT",    0},
        { NULL, NULL, NULL, 0 } };

static char *make_platform_name_string();
static void view_hierarchy_specification(zoltan_platform_specification *spec, int rank,
    int verbose);

#define ZOLTAN_MAX_SIBLINGS  32768      /* A sanity check, valid in 2010 AD */

/* 
 * The topology can be provided by giving a platform name in the PLATFORM_NAME parameter.  
 *
 * Alternatively the topology can be indicated in a string with the TOPOLOGY parameter.
 *
 * Quad-socket six-core:  TOPOLOGY="4,6".
 * Quad-socket six-core , every core pair shares a cache: TOPOLOGY="4,3,2"
 * We assume that every node has the same structure.
 *
 * If you add a node topology to the zoltan_hier_platform_specs, include a "1" before the
 * node topology numbers, and include it in the level count.  Process are partitioned
 * across nodes before they are partitioned within nodes.
 *
 * The machine name should be lower case.
 */

zoltan_platform_specification zoltan_hier_platform_specs[ZOLTAN_HIER_LAST_PLATFORM]={

{"glory",         /* machine named Glory */
  3,              /* 3-level hierarchy */
  {1, 4, 4}},     /* 1 node, 4 sockets, 4 cpus */

{"redsky",       /* machine named RedSky */
  3,             /* 3-level hierarchy */
  {1, 2, 4}},     /* 1 node, 2 sockets, 4 cpus */

{"ctx",          /* machine named CTX */
  3,             /* 3-level hierarchy */
  {1, 2, 6}},    /* 1 node, 2 sockets, 6 cpus */

{"odin",         /* machine named odin */
  3,             /* 3-level hierarchy */
  {1, 2, 4}},    /* 1 node, 2 sockets, 4 cpus */

{"octopi",        /* eight-core machine named octopi */
  2,             /* 2-level hierarchy */
  {2, 4}},       /* 2 sockets, 4 cpus */

{"s861036",      /* dual-core machine named s861036 */
  1,             /* 1-level hierarchy */
  {2}}           /* 2 cpus */
};

static int Zoltan_Hier_Assist_Num_Levels(void *data, int *ierr)
{
  zoltan_platform_specification *spec = (zoltan_platform_specification *)data;
  *ierr = ZOLTAN_OK;

  if (spec == NULL){
    *ierr = ZOLTAN_FATAL;
    return 0;
  }
  else{
    return spec->numLevels;
  }
}

static int Zoltan_Hier_Assist_Part_Number(void *data, int level, int *ierr)
{
  zoltan_platform_specification *spec = (zoltan_platform_specification *)data;
  *ierr = ZOLTAN_OK;

  return spec->my_part[level];
}

static void Zoltan_Hier_Assist_Method(void *data, int level, struct Zoltan_Struct *to, int *ierr)
{
  struct Zoltan_Struct *from;
  *ierr = ZOLTAN_OK;

  from = (struct Zoltan_Struct *)data;

  if (from->LB.Imb_Tol_Len > 0){
    memcpy(to->LB.Imbalance_Tol, from->LB.Imbalance_Tol, sizeof(float) * from->LB.Imb_Tol_Len);
  }

  to->Debug_Proc = 0;

  strcpy(to->LB.Approach, from->LB.Approach);

  if (from->Seed > 0){
    to->Seed = from->Seed;
    Zoltan_Srand(from->Seed, NULL);
  }

  if ((from->Get_Num_Edges != NULL || from->Get_Num_Edges_Multi != NULL) &&
           (from->Get_Edge_List != NULL || from->Get_Edge_List_Multi != NULL)) {

    Zoltan_Filter_Params(to, from, graph_type_params, from->Debug_Level, to->Proc, 0);
    Zoltan_Set_Param(to, "LB_METHOD", "GRAPH");
  }
  else if (from->Get_Num_Geom != NULL &&
      (from->Get_Geom != NULL || from->Get_Geom_Multi != NULL)) {

    Zoltan_Filter_Params(to, from, geometric_type_params, from->Debug_Level, to->Proc, 0);
    Zoltan_Set_Param(to, "LB_METHOD", "RIB");   /* TODO figure out RIB, RCB or HSFC? */
  }
  else{
    *ierr = ZOLTAN_FATAL;
  }

  return;
}
/*****************************************************************************/
/* parameters for the hierarchical balancing.  Used in  */
/* Zoltan_Hier_Set_Param and Zoltan_Hier          */
static PARAM_VARS Hier_params[] = {
  {  "HIER_DEBUG_LEVEL", NULL, "INT", 0},
  {  "HIER_GENERATE_FILES", NULL, "INT", 0},
  {  "HIER_CHECKS", NULL, "INT" , 0},

  {  "HIER_ASSIST", NULL, "INT", 0},  /* If "1", Zoltan determines the hierarchy */

     /* If HIER_ASSIST is "1", define either PLATFORM_NAME or TOPOLOGY */

  {  "PLATFORM_NAME", NULL, "STRING", 0},  /* a name from zoltan_hier_platform_specs above */

     /* topology: for example
          double socket, quad core: "2, 4"
          dual processor work station: "2"
          quad socket, each with 2 L3 caches, 3 cores per cache: "4,2,3" 
      */
  {  "TOPOLOGY", NULL, "STRING", 0},

  {  NULL,              NULL,  NULL, 0 }};

/* prototypes for static functions: */
static int Zoltan_Hier_Initialize_Params(ZZ*, HierPartParams*);
static int Zoltan_Hier_Num_Obj_Fn(void *data, int *ierr);
static void Zoltan_Hier_Obj_List_Fn(void *data, int num_gid_entries,
				    int num_lid_entries, 
				    ZOLTAN_ID_TYPE * global_ids, 
				    ZOLTAN_ID_TYPE * local_ids, 
				    int wgt_dim, float *obj_wgts, int *ierr);
static int Zoltan_Hier_Num_Geom_Fn(void *data, int *ierr);
static void Zoltan_Hier_Geom_Fn(void *data, int num_gid_entries, 
				int num_lid_entries, 
				ZOLTAN_ID_TYPE * global_id, 
				ZOLTAN_ID_TYPE * local_id,
				double *coor, int *ierr);
static int Zoltan_Hier_Num_Edges_Fn(void *data, int num_gid_entries, 
				    int num_lid_entries, 
				    ZOLTAN_ID_TYPE * global_id, 
				    ZOLTAN_ID_TYPE * local_id,
				    int *ierr);
static void Zoltan_Hier_Edge_List_Fn(void *data, int num_gid_entries, 
				     int num_lid_entries, 
				     ZOLTAN_ID_TYPE * global_id,
				     ZOLTAN_ID_TYPE * local_id, 
				     ZOLTAN_ID_TYPE * nbor_global_id, 
				     int *nbor_procs,
				     int wgt_dim, float *ewgts, int *ierr);
static int find_needed_gno_procs(HierPartParams *);
static int Zoltan_Hier_Obj_Size_Fn(void *data, 
				   int num_gid_entries, int num_lid_entries, 
				   ZOLTAN_ID_TYPE * global_id, 
				   ZOLTAN_ID_TYPE * local_id, int *ierr);
static void Zoltan_Hier_Pack_Obj_Fn(void *data, 
				    int num_gid_entries, int num_lid_entries,
				    ZOLTAN_ID_TYPE * global_id, 
				    ZOLTAN_ID_TYPE * local_id, int dest_proc,
				    int size, char *buf, int *ierr);
static void Zoltan_Hier_Unpack_Obj_Fn(void *data, int num_gid_entries,
				      ZOLTAN_ID_TYPE * global_id, 
				      int size, char *buf, int *ierr);
static void Zoltan_Hier_Mid_Migrate_Fn(void *data, int num_gid_entries,
				       int num_lid_entries, int num_import,
				       ZOLTAN_ID_TYPE * import_gids,
				       ZOLTAN_ID_TYPE * import_lids,
				       int *import_procs, int *import_parts,
				       int num_export,
				       ZOLTAN_ID_TYPE * export_gids,
				       ZOLTAN_ID_TYPE * export_lids,
				       int *export_procs, int *export_parts,
				       int *ierr);
static void free_hier_mig_data(void *);
static void Zoltan_Hier_Check_Data(HierPartParams *, int *);

/* for debugging of internal GIDs for small meshes.  define only if 
   the number of objects is small and debugging sanity checks are desired */
#if 0
#define HIER_CHECK_GNO_RANGE(x) {if ((x)<0 || (x)>hpp->vtxdist[hpp->origzz->Num_Proc]) printf("[%d] suspicious GNO %ld on line %d\n", hpp->origzz->Proc, (x), __LINE__);}
#else
#define HIER_CHECK_GNO_RANGE(x)
#endif

/*****************************************************************************/

int Zoltan_Hier_Set_Param(
  char *name,                 /* name of variable */
  char *val                   /* value of variable */
)
{
  int status;
  PARAM_UTYPE result;        /* value returned from Zoltan_Check_Param */
  int index;                 /* index returned from Zoltan_Check_Param */

  status = Zoltan_Check_Param(name, val, Hier_params, &result, &index);
  return(status);
}

/*****************************************************************************/
/* static helper function for Zoltan_Hier */
/* compute new groups, rank maps, communicators for next level */
/* returns error condition */
static int split_comm(HierPartParams *hpp) {
  int ierr = ZOLTAN_OK;

  MPI_Group orig_group, hier_group;
  MPI_Comm hier_comm_copy;
  int *origranks;
  int i;
  
  if (hpp->num_parts > 1){
    /* split the communicator based on the partition it just participated in */
    /* not sure if this much swapping around communicators is needed */
    MPI_Comm_dup(hpp->hier_comm, &hier_comm_copy);
    MPI_Comm_free(&hpp->hier_comm);  
    MPI_Comm_split(hier_comm_copy, hpp->part_to_compute, 0, &hpp->hier_comm);
    MPI_Comm_free(&hier_comm_copy);
  
    /* compute correspondence for hier_ranks_of_orig array */
    MPI_Comm_group(hpp->origzz->Communicator, &orig_group);
    MPI_Comm_group(hpp->hier_comm, &hier_group);
  
    origranks = (int *)ZOLTAN_MALLOC(hpp->origzz->Num_Proc*sizeof(int));
    if (!origranks) {
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, "split_comm", "Out of memory");
      return ZOLTAN_MEMERR;
    }
    for (i=0; i<hpp->origzz->Num_Proc; i++) {
      origranks[i] = i;
    }
  
    MPI_Group_translate_ranks(orig_group, hpp->origzz->Num_Proc, origranks, 
  			    hier_group, hpp->hier_ranks_of_orig);
    MPI_Group_free(&orig_group);
    MPI_Group_free(&hier_group);
  
    ZOLTAN_FREE(&origranks);
  }
  
  return ierr;
}

/*****************************************************************************/
/* static helper function for Zoltan_Hier */
/* compute partition sizes for the current level and set them in hierzz */
/* part_sizes: input array of size
   hpp->origzz->Num_Global_Parts * hpp->origzz->Obj_Weight_Dim
   containing the percentage of work to be
   assigned to each final global partition.               */
/* returns error condition */
static int set_hier_part_sizes(HierPartParams *hpp, float *part_sizes) {
  int ierr = ZOLTAN_OK;
  float *my_level_part_sizes=NULL, *level_part_sizes=NULL;
  int *part_ids=NULL, *wgt_idx=NULL;
  int i;
  char msg[256];
  int part_weight_dim = hpp->origzz->Obj_Weight_Dim;

  /* when this is called, hpp->num_parts contains the number of
     partitions to be computed at this level, and hpp->part_to_compute
     contains the partition id to be computed by this process at this
     level, hpp->hier_comm is a communicator for all procs participating
     at this level. */

  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] set_hier_part_sizes at level %d, computing %d partitions\n",
	   hpp->origzz->Proc, hpp->level, hpp->num_parts);
  }

  /* careful of part_weight_dim of 0 for variable partition sizes */
  if (part_weight_dim == 0) part_weight_dim = 1;

  /* allocate an array for input to reduction to compute
     partition sizes for this level */
  my_level_part_sizes = (float *)ZOLTAN_MALLOC(hpp->num_parts *
					       part_weight_dim * 
					       sizeof(float));
  if (!my_level_part_sizes) {
    sprintf(msg, "Out of memory, tried to alloc %u bytes", 
	(unsigned int)(hpp->num_parts * part_weight_dim * sizeof(float)));
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, "set_hier_part_sizes", msg);
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  for (i=0; i<hpp->num_parts * part_weight_dim; i++) {
    my_level_part_sizes[i] = 0;
  }

  /* put in my part_sizes for the partition I'll be computing */
  for (i=0; i<part_weight_dim; i++) {
    my_level_part_sizes[hpp->part_to_compute*part_weight_dim+i] =
      part_sizes[hpp->origzz->Proc*part_weight_dim+i];
  }

  /* allocate an array for result of reduction of
     partition sizes for this level */
  level_part_sizes = (float *)ZOLTAN_MALLOC(hpp->num_parts *
					    part_weight_dim *
					    sizeof(float));
  if (!level_part_sizes) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, "set_hier_part_sizes",
		       "Out of memory");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /* do the reduction to get global array on each proc */
  MPI_Allreduce(my_level_part_sizes, level_part_sizes, 
		hpp->num_parts * part_weight_dim,
		MPI_FLOAT, MPI_SUM, hpp->hier_comm);

  /* allocate and populate extra args to set_part_sizes) */
  part_ids = (int *)ZOLTAN_MALLOC(hpp->num_parts *
				  part_weight_dim *
				  sizeof(int));
  wgt_idx = (int *)ZOLTAN_MALLOC(hpp->num_parts *
				 part_weight_dim *
				 sizeof(int));
  if (!part_ids || !wgt_idx) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, "set_hier_part_sizes",
		       "Out of memory");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  for (i=0; i<hpp->num_parts * part_weight_dim; i++) {
    part_ids[i] = i/part_weight_dim;
    wgt_idx[i] = i%part_weight_dim;
  }

  if (hpp->output_level >= HIER_DEBUG_ALL) {
    for (i=0; i<hpp->num_parts * part_weight_dim; i++) {
      printf("[%d] setting part_size[%d] to %.3f\n", hpp->origzz->Proc,
	     i, level_part_sizes[i]);
    }
  }

  /* set the partition sizes in hpp->hierzz */
  Zoltan_LB_Set_Part_Sizes(hpp->hierzz, 1, 
			   hpp->num_parts * part_weight_dim,
			   part_ids, wgt_idx, level_part_sizes);

End:
  if (my_level_part_sizes) ZOLTAN_FREE(&my_level_part_sizes);
  if (level_part_sizes) ZOLTAN_FREE(&level_part_sizes);
  if (part_ids) ZOLTAN_FREE(&part_ids);
  if (wgt_idx) ZOLTAN_FREE(&wgt_idx);

  return ierr;
}

/*****************************************************************************/
/*****************************************************************************/
/** Zoltan_Hier: main routine for hierarchical balancing *********************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Hier(
  ZZ *zz,                 /* Zoltan structure */
  float *part_sizes,    /* Input:  Array of size
			   zz->LB.Num_Global_Parts * zz->Obj_Weight_Dim
                           containing the percentage of work to be
                           assigned to each partition.               */
  int *num_imp,         /* number of objects to be imported */
  ZOLTAN_ID_TYPE **imp_gids,  /* global ids of objects to be imported */
  ZOLTAN_ID_TYPE **imp_lids,  /* local  ids of objects to be imported */
  int **imp_procs,      /* list of processors to import from */
  int **imp_to_part,    /* list of partitions to which imported objects are 
                           assigned.  */
  int *num_exp,         /* number of objects to be exported */
  ZOLTAN_ID_TYPE **exp_gids,  /* global ids of objects to be exported */
  ZOLTAN_ID_TYPE **exp_lids,  /* local  ids of objects to be exported */
  int **exp_procs,      /* list of processors to export to */
  int **exp_to_part     /* list of partitions to which exported objects are
                           assigned. */
) 
{
  int ierr = ZOLTAN_OK;   /* error flag for initialization checks */
  HierPartParams hpp;     /* hierarchical partitioning parameters,
			     mainly things that will be needed in the 
			     callbacks (a pointer to this is passed as 
			     the user data) */
  int i99, i, exp_index, last_level;
  char msg[256];
  char *yo = "Zoltan_Hier";
  /* internal Zoltan_LB_Partition call parameters */
  int hier_changes=0, hier_num_gid_entries=0, 
    hier_num_lid_entries=0, hier_num_import_objs=0;
  ZOLTAN_ID_TYPE * hier_import_gids=NULL, *hier_import_lids=NULL;
  int *hier_import_procs=NULL, *hier_import_to_part=NULL;
  int hier_num_export_objs=0;
  ZOLTAN_ID_TYPE * hier_export_gids=NULL, *hier_export_lids=NULL;
  int *hier_export_procs=NULL, *hier_export_to_part=NULL;
  int gno_size_for_dd, graph_type = 0;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialize return arguments */
  *num_imp   = *num_exp   = -1;
  *imp_gids  = *exp_gids  = NULL;
  *imp_lids  = *exp_lids  = NULL;
  *imp_procs = *exp_procs = NULL;
  *imp_to_part = *exp_to_part = NULL;

  /* Initialize hpp structure */
  hpp.local_ids = NULL;
  hpp.global_ids = NULL;
  hpp.gnos_of_interest = NULL;
  hpp.gnos_of_interest_procs = NULL;
  hpp.migrated = NULL;
  hpp.vwgt = NULL;
  hpp.input_parts = NULL;
  hpp.vtxdist = NULL;
  hpp.xadj = NULL;
  hpp.adjncy = NULL;
  hpp.ewgts = NULL;
  hpp.adjproc = NULL;
  hpp.geom_vec = NULL;
  hpp.hier_ranks_of_orig = NULL;
  hpp.migrated_in_gnos = NULL;
  hpp.migrated_in_data = NULL;
  hpp.dd = NULL;
  hpp.hierzz = NULL;
  hpp.hier_comm = MPI_COMM_NULL;
  hpp.spec = NULL;
  hpp.gen_files = 0;

  gno_size_for_dd = sizeof(ZOLTAN_GNO_TYPE) / sizeof(ZOLTAN_ID_TYPE);

  /* Cannot currently do hierarchical balancing for num_parts != num_procs */
  if ((zz->Num_Proc != zz->LB.Num_Global_Parts) ||
      (!zz->LB.Single_Proc_Per_Part)) {
    ZOLTAN_HIER_ERROR(ZOLTAN_FATAL, "number_partitions != number_processes not yet supported by LB_METHOD HIER");
  }

  /* Initialize hierarchical partitioning parameters. */
  ierr = Zoltan_Hier_Initialize_Params(zz, &hpp);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_HIER_ERROR(ierr, "Zoltan_Hier_Initialize_Params returned error");
  }

  if (!hpp.spec){

    /* Make sure we have the callbacks we need */

    if (zz->Get_Hier_Num_Levels == NULL) {
      ZOLTAN_HIER_ERROR(ZOLTAN_FATAL, "Must register ZOLTAN_HIER_NUM_LEVELS_FN");
    }
    if (zz->Get_Hier_Part == NULL) {
      ZOLTAN_HIER_ERROR(ZOLTAN_FATAL, "Must register ZOLTAN_HIER_PART_FN");
    }
    if (zz->Get_Hier_Method == NULL) {
      ZOLTAN_HIER_ERROR(ZOLTAN_FATAL, "Must register ZOLTAN_HIER_METHOD_FN");
    }
  }
  else{

    /* Zoltan defines the callbacks based on the network topology */

    Zoltan_Set_Hier_Num_Levels_Fn(zz, Zoltan_Hier_Assist_Num_Levels, (void *)hpp.spec);
    Zoltan_Set_Hier_Part_Fn(zz, Zoltan_Hier_Assist_Part_Number, (void *)hpp.spec);
    Zoltan_Set_Hier_Method_Fn(zz, Zoltan_Hier_Assist_Method, (void *)zz);
  }

  /* do we have callbacks to get geometric and/or graph information? */
  hpp.use_geom = ((zz->Get_Num_Geom != NULL) ||
		  (zz->Get_Geom_Multi != NULL));
  hpp.use_graph = ((zz->Get_Num_Edges != NULL) || 
		   (zz->Get_Num_Edges_Multi != NULL));

  if (!hpp.use_geom){
    i = 0;
#ifdef ZOLTAN_PARMETIS
    i = 1;
#endif
#ifdef ZOLTAN_PTSCOTCH
    i = 1;
#endif

    if (i == 0){
      ZOLTAN_HIER_ERROR(ZOLTAN_FATAL, 
         "LB_METHOD HIER requires a third party graph library "
         "if only graph callbacks are defined.");
    }
  }

  /* Check weight dimensions */
  if (zz->Obj_Weight_Dim<0){
    sprintf(msg, "Object weight dimension is %d, "
            "but should be >= 0. Using Obj_Weight_Dim = 0.",
            zz->Obj_Weight_Dim);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    hpp.obj_wgt_dim = 0;
  }
  else {
    hpp.obj_wgt_dim = zz->Obj_Weight_Dim;
  }
  if (zz->Edge_Weight_Dim<0){
    sprintf(msg, "Edge weight dimension is %d, "
            "but should be >= 0. Using Edge_Weight_Dim = 0.",
            zz->Edge_Weight_Dim);
    ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    hpp.edge_wgt_dim = 0;
  }
  else {
    hpp.edge_wgt_dim = zz->Edge_Weight_Dim;
  }

  /* build our initial intermediate representation.  at least for now,
     this is the same graph that gets fed to parmetis, so we just call
     that function */

  ierr = Zoltan_Get_Obj_List(zz, &hpp.init_num_obj,
			     &hpp.global_ids, &hpp.local_ids,
			     hpp.obj_wgt_dim, &hpp.vwgt, &hpp.input_parts);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      ZOLTAN_HIER_ERROR(ierr, "Zoltan_Get_Obj_List returned error.");
  }

  /* current number of objects starts as initial number of objects */
  hpp.hier_num_obj = hpp.init_num_obj;

  if (hpp.init_num_obj) {
    hpp.migrated = (char *)ZOLTAN_MALLOC(sizeof(char)*hpp.init_num_obj);
    if (!hpp.migrated) {
      ZOLTAN_HIER_ERROR(ZOLTAN_MEMERR, "Out of memory");
    }
    for (i=0; i<hpp.init_num_obj; i++) {
      hpp.migrated[i] = 0;
    }
  }

  /* use a distributed data manager to track gids as they are
     migrated.  Objects without an entry are in their original
     location as determined by the gid value (index into global
     ordering) and the vtxdist array */

  ierr = Zoltan_DD_Create(&hpp.dd, zz->Communicator, gno_size_for_dd, 0, 0, 0, 0);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      ZOLTAN_HIER_ERROR(ierr, "Zoltan_DD_Create returned error.");
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
    printf("[%1d] Debug: init_num_obj =%d\n", zz->Proc, hpp.init_num_obj);
    printf("[%1d] Debug: Global ids = ", zz->Proc);
    for (i99=0; i99<hpp.init_num_obj; i99++) {
      printf("    ");
      ZOLTAN_PRINT_GID(zz, &(hpp.global_ids[i99*zz->Num_GID]));
      printf("\n");
    }
  }

  /* build a graph (only vtxdist if we don't need a graph) */
  if (hpp.use_graph)
    SET_GLOBAL_GRAPH(&graph_type);
  else
    graph_type |= (1 << NO_GRAPH);
  ierr = Zoltan_Build_Graph(zz, &graph_type,
			    hpp.checks, hpp.init_num_obj,
			    hpp.global_ids, hpp.local_ids,       /* input ZOLTAN_ID_TYPEs */
			    hpp.obj_wgt_dim, &hpp.edge_wgt_dim,
			    &hpp.vtxdist, &hpp.xadj, &hpp.adjncy, /* internal ZOLTAN_GNO_TYPEs */
			    &hpp.ewgts, &hpp.adjproc);

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
    ZOLTAN_HIER_ERROR(ierr, "Zoltan_Build_Graph returned error.");
  }
  
  hpp.num_edges = (hpp.use_graph ? hpp.xadj[hpp.init_num_obj] : 0);

  /* if we're going to need coordinates */
  if (hpp.use_geom){

    /* Get coordinate information */
    ierr = Zoltan_Get_Coordinates(zz, hpp.init_num_obj, hpp.global_ids, 
				  hpp.local_ids, &hpp.ndims, &hpp.geom_vec);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_HIER_ERROR(ierr, "Error returned from Zoltan_Get_Coordinates");
    }
  }

  /* find out how many levels of hierarchy this proc will participate in */
  hpp.num_levels = zz->Get_Hier_Num_Levels(zz->Get_Hier_Num_Levels_Data,
					   &ierr);

  if (hpp.output_level >= HIER_DEBUG_ALL) {
    printf("HIER: Proc %d to compute %d levels\n", zz->Proc, hpp.num_levels);
  }

  MPI_Allreduce(&hpp.num_levels, &hpp.global_num_levels, 1, MPI_INT, MPI_MAX,
		zz->Communicator);

  last_level = hpp.global_num_levels - 1;

  /* initialize our communicator to the "world" as seen by Zoltan */
  MPI_Comm_dup(zz->Communicator, &hpp.hier_comm);

  /* initialize hier_ranks array of correspondence between original
     comm and hier_comm - starts out the same */
  hpp.hier_ranks_of_orig = (int *)ZOLTAN_MALLOC(sizeof(int)*zz->Num_Proc);
  
  if (!hpp.hier_ranks_of_orig ) {
    ZOLTAN_HIER_ERROR(ZOLTAN_MEMERR, "Out of memory");
  }
  for (i=0; i<zz->Num_Proc; i++) {
    hpp.hier_ranks_of_orig[i]=i;
  }

  /* loop over levels of hierarchical balancing to be done */
  for (hpp.level = 0; hpp.level < hpp.num_levels; hpp.level++) {

    /* determine partitions to compute at this level */
    hpp.part_to_compute = 
      zz->Get_Hier_Part(zz->Get_Hier_Part_Data, hpp.level, &ierr);
    /* number of partitions is one more than the highest partition id
       specified on procs in the current hier_comm */
    MPI_Allreduce(&hpp.part_to_compute, &hpp.num_parts, 1, MPI_INT,
		  MPI_MAX, hpp.hier_comm);
    hpp.num_parts++;

    if (hpp.num_parts == 1){
      /* 
       * If there is only one part, and we've not done any partitioning
       * yet, we can skip this step.
       */
      MPI_Comm_size(hpp.hier_comm, &i);
      if (i == zz->Num_Proc)
        continue;
    }

    if (hpp.output_level >= HIER_DEBUG_ALL || 
	zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
      printf("HIER: Proc %d computing part %d of %d at level %d\n",
	     zz->Proc, hpp.part_to_compute, hpp.num_parts, hpp.level);
    }

    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL) {
      printf("[%1d] Debug: hier_num_obj =%d\n", zz->Proc, hpp.hier_num_obj);
      printf("[%1d] Debug: Internal gids (I:imported,E:exported)=", zz->Proc);
      for (i99=0; i99<hpp.init_num_obj; i99++) {
	printf(" %d",i99);
	if (hpp.migrated[i99]) printf("E");	
      }
      for (i99=0; i99<hpp.num_migrated_in_gnos; i99++) {
	printf(ZOLTAN_GNO_SPEC "I", hpp.migrated_in_gnos[i99]);
      }
      printf("\n");
    }

    /* should make sure we have reasonable partitions to compute */

    hpp.hierzz = NULL;

    if (hpp.num_parts > 1){

      /* construct appropriate ZZ and input arrays */
      /* create a brand new one */
      hpp.hierzz = Zoltan_Create(hpp.hier_comm);
  
      /* and copy in some specified params from zz where appropriate */
  
      /* just copy debug level to child Zoltan_Struct, use can override
         by setting params of the hierzz in the Get_Hier_Method callback */
  
      hpp.hierzz->Debug_Level = zz->Debug_Level;
      hpp.hierzz->Timer = zz->Timer;
      hpp.hierzz->Deterministic = zz->Deterministic;
      hpp.hierzz->Obj_Weight_Dim = zz->Obj_Weight_Dim;
      hpp.hierzz->Edge_Weight_Dim = zz->Edge_Weight_Dim;
  
      /* remapping does not make sense for internal steps, only at the end */
      hpp.hierzz->LB.Remap_Flag = 0;
  
      /* let the application specify any balancing params for this level */
      zz->Get_Hier_Method(zz->Get_Hier_Method_Data, hpp.level,
  			hpp.hierzz, &ierr);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Get_Hier_Method callback returned error.");
      }
  
      /* set the numbers of partitions */
      sprintf(msg, "%d", hpp.num_parts);
      Zoltan_Set_Param(hpp.hierzz, "NUM_GLOBAL_PARTS", msg);

      /* specify the callbacks */

      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_NUM_OBJ_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Num_Obj_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
   
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_OBJ_LIST_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Obj_List_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_NUM_GEOM_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Num_Geom_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_GEOM_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Geom_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_NUM_EDGES_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Num_Edges_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_EDGE_LIST_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Edge_List_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_OBJ_SIZE_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Obj_Size_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_PACK_OBJ_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Pack_Obj_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_UNPACK_OBJ_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Unpack_Obj_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      ierr = Zoltan_Set_Fn(hpp.hierzz, ZOLTAN_MID_MIGRATE_PP_FN_TYPE,
  			 (void (*)()) Zoltan_Hier_Mid_Migrate_Fn,
  			 (void *) &hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_Set_Fn returned error.");
      }
  
      /* specify the GIDs (just the global numbering) */  
      Zoltan_Set_Param(hpp.hierzz, "NUM_GID_ENTRIES", "1");
      Zoltan_Set_Param(hpp.hierzz, "NUM_LID_ENTRIES", "0");
  
      /* specify that we need import and export lists (at least for now) */
      Zoltan_Set_Param(hpp.hierzz, "RETURN_LISTS", "ALL");
  
      /* we want to have Zoltan's migration routines called automatically
         by Zoltan_LB_Partition, except in the last round, where we can
         save some work by processing result arrays manually */
      if (hpp.level == last_level) {
        Zoltan_Set_Param(hpp.hierzz, "AUTO_MIGRATE", "0");
      }
      else {
        Zoltan_Set_Param(hpp.hierzz, "AUTO_MIGRATE", "1");
      }
  
      /* deal with partition sizes, etc */
      /* we have the assumption here that the final result is one
         partition per process */
      ierr = set_hier_part_sizes(&hpp, part_sizes);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_HIER_ERROR(ierr, "set_hier_part_sizes returned error");
      }
    }

    /* fill in array of proc assignments for migrated objects 
     *
     * Zoltan_DD_* is called here, so all processes must participate,
     *  even if they computing partition 0 out of a total of 1.
     */

    if (hpp.use_graph && hpp.level) {     /* if nlocalprocs < nglobalprocs */
      ierr = find_needed_gno_procs(&hpp);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
	ZOLTAN_HIER_ERROR(ierr, "find_needed_gno_procs returned error.");
      }
    }

    if (hpp.num_parts > 1){

      if (hpp.gen_files){
        sprintf(msg,"level_%d",hpp.level);
        ierr = Zoltan_Generate_Files(hpp.hierzz, msg, zz->Proc, 0, 1, 0);
        if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
          ZOLTAN_HIER_ERROR(ierr, "Zoltan_Generate_Files returned error.");
        }
      }

      /* call partitioning method to compute the partitions at this level */
      ierr = Zoltan_LB_Partition(hpp.hierzz, &hier_changes, 
  			       &hier_num_gid_entries, &hier_num_lid_entries, 
  			       &hier_num_import_objs,
  			       &hier_import_gids, &hier_import_lids,
  			       &hier_import_procs, &hier_import_to_part,
  			       &hier_num_export_objs,
  			       &hier_export_gids, &hier_export_lids,
  			       &hier_export_procs, &hier_export_to_part);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
        ZOLTAN_HIER_ERROR(ierr, "Zoltan_LB_Partition returned error.");
      }
    }
    else{
      hier_changes=0;
      hier_num_gid_entries = 1; hier_num_lid_entries = 0;
      hier_num_import_objs = hier_num_export_objs = 0;
      hier_import_gids = hier_import_lids=NULL;
      hier_import_procs = hier_import_to_part=NULL;
      hier_export_gids = hier_export_lids=NULL;
      hier_export_procs = hier_export_to_part=NULL;
    }

    /* clean up array of proc assignments for migrated objects */
    if (hpp.gnos_of_interest) ZOLTAN_FREE(&hpp.gnos_of_interest);
    if (hpp.gnos_of_interest_procs) ZOLTAN_FREE(&hpp.gnos_of_interest_procs);
    hpp.num_gnos_of_interest = 0;
    hpp.allocsize_gnos_of_interest = 0;

    /* processing of output arrays is done by migration callbacks 
       except in the last round when we need to call just the
       mid-migrate function explicitly */

    /* If hpp.num_parts == 1, then we didn't call Zoltan_LB_Partition, which is
     *   is going to call Zoltan_Migrate which will call Mid_Migrate.  Mid_Migrate
     *   does global communication, so all procs need to call it.
     */
       
    if ((hpp.level == last_level) || (hpp.num_parts == 1)) {
      Zoltan_Hier_Mid_Migrate_Fn((void *)&hpp, hier_num_gid_entries,
				 hier_num_lid_entries,
				 hier_num_import_objs,
				 hier_import_gids, hier_import_lids,
				 hier_import_procs, hier_import_to_part,
				 hier_num_export_objs,
				 hier_export_gids, hier_export_lids,
				 hier_export_procs, hier_export_to_part,
				 &ierr);
      if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
	ZOLTAN_HIER_ERROR(ierr, "Zoltan_Hier_Mid_Migrate_Fn returned error");
      }
    }

    /* free memory from this Zoltan_LB_Partition call */
    ierr = Zoltan_LB_Free_Part(&hier_import_gids, &hier_import_lids,
			       &hier_import_procs, &hier_import_to_part);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      ZOLTAN_HIER_ERROR(ierr, "Zoltan_LB_Free_Part returned error.");
    }
    
    ierr = Zoltan_LB_Free_Part(&hier_export_gids, &hier_export_lids,
			       &hier_export_procs, &hier_export_to_part);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      ZOLTAN_HIER_ERROR(ierr, "Zoltan_LB_Free_Part returned error.");
    }

    /* clean up hierzz */
    Zoltan_Destroy(&hpp.hierzz);

    /* create communicators to do next level */
    /* group procs that worked on the same partition at this level to
       work in a group at the next level
     */
    ierr = split_comm(&hpp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_HIER_ERROR(ierr, "split_comm returned error");
    }
      
    if (hpp.level != last_level) {
      if (hpp.checks) {
	/* check intermediate structure correctness when HIER_CHECKS is 1 
	   since the structure is not fully updated on the last level,
	   we skip the check if it's the last level */
	Zoltan_Hier_Check_Data(&hpp, &ierr);
	if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
	  ZOLTAN_HIER_ERROR(ierr, "Zoltan_Hier_Check_Data returned error");
	}
      }
    }
  }

  /* figure out what all that did and construct return arrays for
     original hierarchical balancing call -- for now at least we can
     get this information from the DD -- all migrated out objects
     should be in there.  All migrated out objects are marked as such
     in the local migrated array */
  /* first compute the number of exports */
  *num_exp = 0;
  for (i=0; i<hpp.init_num_obj; i++) {
    if (hpp.migrated[i]) (*num_exp)++;
  }

  /* allocate export arrays and gnos_of_interest array to find
     destination procs */
  if (*num_exp > 0) {
    hpp.gnos_of_interest = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(*num_exp * sizeof(ZOLTAN_GNO_TYPE));
    if (!hpp.gnos_of_interest) {
      ZOLTAN_HIER_ERROR(ZOLTAN_MEMERR, "Out of memory");
    }
    
    if (!Zoltan_Special_Malloc(zz, (void **)exp_gids, *num_exp,
                               ZOLTAN_SPECIAL_MALLOC_GID) ||
	!Zoltan_Special_Malloc(zz, (void **)exp_lids, *num_exp,
			       ZOLTAN_SPECIAL_MALLOC_LID) ||
	!Zoltan_Special_Malloc(zz, (void **)exp_procs, *num_exp,
			       ZOLTAN_SPECIAL_MALLOC_INT) ||
	!Zoltan_Special_Malloc(zz, (void **)exp_to_part, *num_exp,
			       ZOLTAN_SPECIAL_MALLOC_INT)) {
      ZOLTAN_HIER_ERROR(ZOLTAN_MEMERR, 
			"Memory error allocating export arrays.");
    }

    /* fill in GID and LID return arrays and gnos_of_interest */
    exp_index = 0;
    for (i=0; i<hpp.init_num_obj; i++) {
      if (hpp.migrated[i]) {
	/* internal gid is global numbering */
	hpp.gnos_of_interest[exp_index] = i + hpp.vtxdist[zz->Proc];
	ZOLTAN_SET_GID(zz, &((*exp_gids)[exp_index*zz->Num_GID]),
		       &(hpp.global_ids[i*zz->Num_GID]));
	if (zz->Num_LID) {
	  ZOLTAN_SET_LID(zz, &((*exp_lids)[exp_index*zz->Num_LID]),
			 &(hpp.local_ids[i*zz->Num_LID]));
	}
	exp_index++;
      }
    }

  }
  /* query the DD for procs - should put them right where we need them */
  ierr = Zoltan_DD_Find(hpp.dd, (ZOLTAN_ID_PTR)hpp.gnos_of_interest, NULL, NULL, NULL,   
			*num_exp, *exp_procs);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
    ZOLTAN_HIER_ERROR(ierr, "Zoltan_DD_Find returned error");
  }

  /* copy processor assignments to partition assignments.  works here since
     we restrict to num_parts == num_procs for now */
  for (i=0; i<*num_exp; i++) {
    (*exp_to_part)[i] = (*exp_procs)[i];
  }

  /* all errors as well as normal termination go here, so we can do
     all cleanup in one place - we check for NULL pointers since not
     all of these are ever allocated in all situations, and we could
     be exiting early because of an error condition */
End:
  if (hpp.local_ids) ZOLTAN_FREE(&hpp.local_ids);
  if (hpp.global_ids) ZOLTAN_FREE(&hpp.global_ids);
  if (hpp.gnos_of_interest) ZOLTAN_FREE(&hpp.gnos_of_interest);
  if (hpp.gnos_of_interest_procs) ZOLTAN_FREE(&hpp.gnos_of_interest_procs);
  /*if (hpp.migrated_to) ZOLTAN_FREE(&hpp.migrated_to);*/
  if (hpp.migrated) ZOLTAN_FREE(&hpp.migrated);
  if (hpp.vwgt) ZOLTAN_FREE(&hpp.vwgt);
  if (hpp.input_parts) ZOLTAN_FREE(&hpp.input_parts);
  if (hpp.vtxdist) ZOLTAN_FREE(&hpp.vtxdist);
  if (hpp.xadj) ZOLTAN_FREE(&hpp.xadj);
  if (hpp.adjncy) ZOLTAN_FREE(&hpp.adjncy);
  if (hpp.ewgts) ZOLTAN_FREE(&hpp.ewgts);
  if (hpp.adjproc) ZOLTAN_FREE(&hpp.adjproc);
  if (hpp.geom_vec) ZOLTAN_FREE(&hpp.geom_vec);
  if (hpp.hier_ranks_of_orig) ZOLTAN_FREE(&hpp.hier_ranks_of_orig);
  /*if (hpp.orig_ranks_of_hier) ZOLTAN_FREE(&hpp.orig_ranks_of_hier);*/
  if (hpp.migrated_in_gnos) ZOLTAN_FREE(&hpp.migrated_in_gnos);
  if (hpp.migrated_in_data) {
    for (i=0; i<hpp.num_migrated_in_gnos; i++) {
      if (hpp.migrated_in_data[i]) {
	free_hier_mig_data(hpp.migrated_in_data[i]);
      }
    }
    ZOLTAN_FREE(&hpp.migrated_in_data);
  }

  if (hpp.spec && (hpp.spec->platform_name == NULL)){
    ZOLTAN_FREE(&hpp.spec);
  }

  if (hpp.dd) Zoltan_DD_Destroy(&hpp.dd);

  if (hpp.hierzz) Zoltan_Destroy(&hpp.hierzz);
  if (hpp.hier_comm != MPI_COMM_NULL) MPI_Comm_free(&hpp.hier_comm);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/
/* Initialize the parameter structure for hierarchical */
static int Zoltan_Hier_Initialize_Params(ZZ *zz, HierPartParams *hpp) {

char *yo = "Zoltan_Hier_Initialize_Params";
int assist, i=0, j, len;
int num_cpus, num_siblings;
char platform[MAX_PARAM_STRING_LEN+1];
char topology[MAX_PARAM_STRING_LEN+1];
char *msg=NULL, *pnames=NULL, *c=NULL;
div_t result;

  Zoltan_Bind_Param(Hier_params, "HIER_DEBUG_LEVEL", (void *) &hpp->output_level);
  Zoltan_Bind_Param(Hier_params, "HIER_GENERATE_FILES", (void *) &hpp->gen_files);
  Zoltan_Bind_Param(Hier_params, "HIER_CHECKS", (void *) &hpp->checks);
  Zoltan_Bind_Param(Hier_params, "HIER_ASSIST", (void *) &assist);
  Zoltan_Bind_Param(Hier_params, "PLATFORM_NAME", (void *) platform);
  Zoltan_Bind_Param(Hier_params, "TOPOLOGY", (void *) topology);

  /* set default values */
  hpp->output_level = HIER_DEBUG_NONE;
  hpp->checks = 0;
  assist = 0;
  platform[0] = topology[0] = 0;

  /* Get application values of parameters. */
  Zoltan_Assign_Param_Vals(zz->Params, Hier_params, zz->Debug_Level, zz->Proc,
                           zz->Debug_Proc);

  /* initialize other members of hpp */
  hpp->init_num_obj = 0;
  hpp->hier_num_obj = 0;
  hpp->local_ids = NULL;
  hpp->global_ids = NULL;
  /*hpp->migrated_to = NULL;*/
  hpp->migrated = NULL;
  hpp->allocsize_gnos_of_interest = 0;
  hpp->num_gnos_of_interest = 0;
  hpp->gnos_of_interest = NULL;
  hpp->gnos_of_interest_procs = NULL;
  hpp->vwgt = NULL;
  hpp->input_parts = NULL;
  hpp->vtxdist = NULL;
  hpp->xadj = NULL;
  hpp->adjncy = NULL;
  hpp->ewgts = NULL;
  hpp->adjproc = NULL;
  hpp->ndims = 0;
  hpp->num_edges = 0;
  hpp->geom_vec = NULL;
  hpp->origzz = zz;
  hpp->hier_ranks_of_orig = NULL;
  /*hpp->orig_ranks_of_hier = NULL;*/
  hpp->num_migrated_in_gnos = 0;
  hpp->alloc_migrated_in_gnos = 0;
  hpp->migrated_in_gnos = NULL;
  hpp->migrated_in_data = NULL;

  hpp->spec = NULL;

  if (!assist)
    return ZOLTAN_OK;

  if (platform[0]){ 
    len = strlen(platform);
    for (i=0; i < len; i++){
      if (isupper((int)platform[i]))
        platform[i] = (char)tolower((int)platform[i]);
    }

    for (i=0; i < ZOLTAN_HIER_LAST_PLATFORM; i++){
      if (strcmp(platform, zoltan_hier_platform_specs[i].platform_name)) continue;
      hpp->spec = zoltan_hier_platform_specs + i;
      break;
    }
  }

  if (!hpp->spec && topology[0]){
    hpp->spec = 
      (zoltan_platform_specification *)ZOLTAN_CALLOC(sizeof(zoltan_platform_specification), 1);

    if (!hpp->spec){
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
      return ZOLTAN_MEMERR;
    }

    hpp->spec->platform_name = NULL; 

    if (topology[0]){
      hpp->spec->num_siblings[0] = 1;            /* the node or machine itself is the first level */
      i = 1;
      c = topology;
    }

    j = 0;

    while (*c){
      while (*c && !isdigit(*c)) c++;

      if (*c){
        if (i == PLATFORM_MAX_LEVELS){
          ZOLTAN_FREE(&(hpp->spec)); 
          break;
        }

        sscanf(c, "%d", hpp->spec->num_siblings +  i);

        if ((hpp->spec->num_siblings[i] < 1) || (hpp->spec->num_siblings[i] > ZOLTAN_MAX_SIBLINGS)){
          ZOLTAN_FREE(&(hpp->spec)); 
          break;
        }
        i++;
        j++;
      } 

      while (*c && isdigit(*c)) c++;

    }

    hpp->spec->numLevels = i;

    if (j == 0){
      ZOLTAN_FREE(&(hpp->spec));
    }
  }

  if (!hpp->spec){
    if (zz->Proc == 0){
      pnames = make_platform_name_string();
      if (pnames == NULL){
        ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
        return ZOLTAN_MEMERR;
      }
      i = strlen(pnames) + 2048;
      msg = (char *)ZOLTAN_MALLOC(i);
      if (!msg){
        ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
        return ZOLTAN_MEMERR;
      }
      strcpy(msg,"Error:\n");
      strcat(msg, "HIER_ASSIST requested but insufficient topology information provided.\n\n" 
        "Specify PLATFORM_NAME or TOPOLOGY.\n\n");
   
      strcat(msg,"TOPOLOGY is the number of hardware siblings at each level in a topology.\n"
        "  Ex. TOPOLOGY=\"2, 4\" describes a dual-socket quad-core computing cluster.\n"
        "  Ex. TOPOLOGY=\"4\" describes a quad-core desktop computer.\n\n");
  
      strcat(msg,"Zoltan assumes the run-time system has pinned each process to a CPU.\n");
      strcat(msg,"It assumes MPI process ranks map to the topology.  (In the 2,4 example,\n");
      strcat(msg,"this means ranks 0-7 are on the same node, and 0-3 on the same socket.)\n\n");
      
      strcat(msg, "PLATFORM_NAME can be one of the following:\n");
      strcat(msg, pnames);
  
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
  
      ZOLTAN_FREE(&pnames);
      ZOLTAN_FREE(&msg);
    }
  
    return ZOLTAN_FATAL;
  }

  /* 
   * Compute which part my process has at each level.  We are assuming that
   * MPI laid out the process ranks with respect to the topology.  This
   * may not be true and eventually we want a way to determine the
   * topological rank of each process.
   */
  
  num_cpus = 1;
  for (i = 0; i < hpp->spec->numLevels; i++)
    num_cpus *= hpp->spec->num_siblings[i];
  
  result = div(zz->Num_Proc, num_cpus);

  hpp->spec->num_siblings[0] = result.quot; 

  if (result.rem > 0)
    hpp->spec->num_siblings[0]++;  /* number of nodes */

  for (i=0; i < hpp->spec->numLevels; i++){

    /* total number of objects at this level */
    num_siblings = hpp->spec->num_siblings[i];

    /* total number of cpus within an object at this level */
    num_cpus = 1;

    for (j = hpp->spec->numLevels-1; j > i; j--)
      num_cpus *= hpp->spec->num_siblings[j];

    result = div(zz->Proc, num_cpus);

    result = div(result.quot, num_siblings);

    hpp->spec->my_part[i] = result.rem;
  }

  if (hpp->output_level >= HIER_DEBUG_LIST){
    MPI_Barrier(MPI_COMM_WORLD);
    for (i=0; i < zz->Num_Proc; i++){
      if (i == zz->Proc){
        view_hierarchy_specification(hpp->spec, i, (i==0));
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return ZOLTAN_OK;
}

/*****************************************************************************/
/****************** Static helper functions **********************************/
/*****************************************************************************/

/* do a binary search on the migrated in array and return the index
   into that array of the given gno */
static int migrated_in_index_of_gno(HierPartParams *hpp, ZOLTAN_GNO_TYPE gno) {
  int low, mid, high;

  HIER_CHECK_GNO_RANGE(gno);

  if (hpp->num_migrated_in_gnos == 0) return -1;

  low = 0;
  high = hpp->num_migrated_in_gnos;
  mid = high/2;

  while (low!=mid && hpp->migrated_in_gnos[mid] != gno) {
    if (hpp->migrated_in_gnos[mid] > gno) {
      high = mid;
    }
    else {
      low = mid;
    }
    mid = (high+low)/2;
  }

  /* found? */
  if (hpp->migrated_in_gnos[mid] == gno) return mid;

  /* nope */
  return -1;  
}

/* data structure for migrated in GID info */
/* maybe this can be reworked as a union or a few separate structures
   to reduce overhead when some info is not needed */
/* now, the wgts array has obj_wgt_dim entries for vertex
   weights followed by edge_wgt_dim edge weights for each adjacent edge
   and the coords array has ndims entries for geometric coordinates */
struct HierGNOInfo {
  float *wgts;
  double *coords;
  int num_adj;
  ZOLTAN_GNO_TYPE *adj;
};

/* free a HierGNOInfo structure */
static void free_hier_mig_data(void *voiddata) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;
  if (data->wgts) ZOLTAN_FREE(&data->wgts);
  if (data->coords) ZOLTAN_FREE(&data->coords);
  if (data->adj) ZOLTAN_FREE(&data->adj);
  ZOLTAN_FREE(&data);		       
}

/* get the edge list from an element of the migrated in list */
/* we return a pointer into the structure since this is static
   function and we no one unauthorized will get their hands on a
   pointer that could go stale */
static void get_hier_mig_adj_info(void *voiddata, int *count, 
				  ZOLTAN_GNO_TYPE **edgelist) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;
  
  *count = data->num_adj;
  *edgelist = data->adj;
}

/* get vertex weights from an element of the migrated list */
static float *get_hier_mig_vwgts(void *voiddata) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;

  return data->wgts;
}

/* get geometric coordinates an element of the migrated list */
static double *get_hier_mig_coords(void *voiddata, HierPartParams *hpp) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;

  return data->coords;
}

/* get number of adjacencies of an element of the migrated list */
static int get_hier_mig_num_adj(void *voiddata) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;

  return data->num_adj;
}

/* get adjacent gno list from an element of the migrated list */
static ZOLTAN_GNO_TYPE *get_hier_mig_adj(void *voiddata) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;

  return data->adj;
}

/* get ith adjacent gid from an element of the migrated list */
static ZOLTAN_GNO_TYPE get_hier_mig_ith_adj(void *voiddata, int adj_index) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;

  return data->adj[adj_index];
}

/* get edge weights of ith adjacent gid from an element of the migrated list */
static void get_hier_mig_ith_adj_wgts(void *voiddata, HierPartParams *hpp, 
				      int adj_index, float *weights) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;
  int i;

  for (i=0; i<hpp->edge_wgt_dim; i++) {
    weights[i] = 
      data->wgts[hpp->obj_wgt_dim + adj_index*hpp->edge_wgt_dim + i];
  }
}

/* get list of edge weights of adjacent gids from an element of the
   migrated list */
static float *get_hier_mig_adj_wgts(void *voiddata, HierPartParams *hpp) {
  struct HierGNOInfo *data = (struct HierGNOInfo *)voiddata;

  return &(data->wgts[hpp->obj_wgt_dim]);
}

/* is a gno located here? */
static int is_gno_local(HierPartParams *hpp, ZOLTAN_GNO_TYPE gno) {

  /* to be local it must be either
     1) in our range of original gnos and not migrated
     or
     2) migrated to here 
  */
  int index;

  HIER_CHECK_GNO_RANGE(gno);

  /* possibility 1 */
  if ((gno >= hpp->vtxdist[hpp->origzz->Proc]) &&
      (gno < hpp->vtxdist[hpp->origzz->Proc+1]) &&
      !hpp->migrated[gno-hpp->vtxdist[hpp->origzz->Proc]])
    return 1;

  /* possibility 2 */
  index = migrated_in_index_of_gno(hpp, gno);
  if (index != -1) return 1;

  return 0;
}

/* insert a gid uniquely into the (ordered) gnos_of_interest array */
/* it might be more efficient to make a big unsorted array with duplicates 
   then remove duplicates and sort */
static int insert_unique_gnos_of_interest(HierPartParams *hpp,
					   ZOLTAN_GNO_TYPE gno) {
  int low, mid, high, index;
  char *yo = "insert_unique_gnos_of_interest";

  HIER_CHECK_GNO_RANGE(gno);

  /* we will not always insert, but we realloc first just in case */
  if (hpp->allocsize_gnos_of_interest == hpp->num_gnos_of_interest) {
    /* reasonable increment? */
    hpp->allocsize_gnos_of_interest += hpp->hier_num_obj + 1;
    hpp->gnos_of_interest = 
      (ZOLTAN_GNO_TYPE *)ZOLTAN_REALLOC(hpp->gnos_of_interest, 
          sizeof(ZOLTAN_GNO_TYPE)* hpp->allocsize_gnos_of_interest);

    if (!hpp->gnos_of_interest){
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
      return ZOLTAN_MEMERR;
    }
  }

  /* find it (or not) */
  if (hpp->num_gnos_of_interest) {
    low=0;
    mid=hpp->num_gnos_of_interest/2;
    high=hpp->num_gnos_of_interest;
    while (low!=mid && hpp->gnos_of_interest[mid] != gno) {
      if (hpp->gnos_of_interest[mid] > gno) {
	high = mid;
      }
      else {
	low = mid;
      }
      mid = (high+low)/2;
    }
    /* did we find it? */
    if (hpp->gnos_of_interest[mid] == gno) return ZOLTAN_OK;
  }

  /* no, insert it */
  index = hpp->num_gnos_of_interest;
  while (index && hpp->gnos_of_interest[index-1] > gno) {
    hpp->gnos_of_interest[index] = hpp->gnos_of_interest[index-1];
    index--;
  }
  hpp->gnos_of_interest[index] = gno;
  hpp->num_gnos_of_interest++;

  return ZOLTAN_OK;
}

/* request locations of migrated gids of interest to graph callbacks */
static int find_needed_gno_procs(HierPartParams *hpp) {
  int i, adjindex;
  ZOLTAN_GNO_TYPE adjgid;
  int ierr=ZOLTAN_OK;
  int numadj;
  ZOLTAN_GNO_TYPE *adjlist;
  char *yo = "find_needed_gno_procs";

  if (hpp->hier_num_obj) {
    /* build a list of remote gids whose locations we will need later */
    hpp->allocsize_gnos_of_interest = hpp->hier_num_obj; /* initial estimate */
    hpp->gnos_of_interest = 
      (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE)*
				   hpp->allocsize_gnos_of_interest);
    if (!hpp->gnos_of_interest) {
      ierr = ZOLTAN_MEMERR; 
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
      goto End;
    }
    hpp->num_gnos_of_interest = 0;
    
    /* traverse adjacency array of locally-owned gids */
    for (i=0; i<hpp->init_num_obj; i++) {
      if (!hpp->migrated[i]) {
	/* insert each non-local adjacent gid to be located */
	for (adjindex=hpp->xadj[i]; adjindex<hpp->xadj[i+1]; adjindex++) {
	  /* is the adjacent object remote? */
	  adjgid = hpp->adjncy[adjindex];
	  HIER_CHECK_GNO_RANGE(adjgid);
	  if (!is_gno_local(hpp, adjgid)) {
	    ierr = insert_unique_gnos_of_interest(hpp, adjgid);
            if (ierr != ZOLTAN_OK)
              goto End;
	  }
	}
      }
    }
    
    /* also traverse migrated-in adjacency lists */
    for (i=0; i<hpp->num_migrated_in_gnos; i++) {
      if (!hpp->migrated_in_data[i]) {
	char msg[256];
	ierr = ZOLTAN_FATAL;
	sprintf(msg, "num_migrated_in_gnos=%d but hpp->migrated_in_data[%d] is NULL\n", hpp->num_migrated_in_gnos, i);
	ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, msg);
	goto End;
      }
      
      get_hier_mig_adj_info(hpp->migrated_in_data[i], &numadj, &adjlist);
      for (adjindex=0; adjindex<numadj; adjindex++) {
	adjgid = adjlist[adjindex];
	HIER_CHECK_GNO_RANGE(adjgid);
	if (!is_gno_local(hpp, adjgid)) {
	  ierr = insert_unique_gnos_of_interest(hpp, adjgid);
          if (ierr != ZOLTAN_OK)
            goto End;
	}
      }
    }
    
    if (hpp->num_gnos_of_interest) {
      hpp->gnos_of_interest_procs =
	(int *)ZOLTAN_MALLOC(sizeof(int)*hpp->num_gnos_of_interest);
      if (!hpp->gnos_of_interest_procs) {
	ierr = ZOLTAN_MEMERR; 
	ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
	goto End;
      }
    } 
  }

  /* do the communication - need to make sure all procs participate! */
  /* note that hpp->gnos_of_interest_procs will be NULL if this proc
     has no gids of interest */
  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] calling DD_Find on %d GNOs\n", hpp->origzz->Proc,
	   hpp->num_gnos_of_interest);
    for (i=0; i<hpp->num_gnos_of_interest; i++) {
      printf("[%d] DD_Find slot %d GNO " ZOLTAN_GNO_SPEC "\n",
	     hpp->origzz->Proc, i, hpp->gnos_of_interest[i]);
    }
  }
  ierr = Zoltan_DD_Find(hpp->dd, (ZOLTAN_ID_PTR)hpp->gnos_of_interest, NULL, NULL, NULL, 
			hpp->num_gnos_of_interest, 
			hpp->gnos_of_interest_procs);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
		       "Zoltan_DD_Find returned error.");
    goto End;
  }
  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] DD_Find on %d GNOs returned\n", hpp->origzz->Proc,
	   hpp->num_gnos_of_interest);
    for (i=0; i<hpp->num_gnos_of_interest; i++) {
      printf("[%d] DD_Find slot %d GNO " ZOLTAN_GNO_SPEC " proc %d\n",
	     hpp->origzz->Proc, i, hpp->gnos_of_interest[i],
	     hpp->gnos_of_interest_procs[i]);
    }
  }

  End:
  return ierr;

}

/* look up the proc of a gid in the gnos_of_interest.  Return proc if
   found, return -1 otherwise */
static int gno_of_interest_proc(HierPartParams *hpp, ZOLTAN_GNO_TYPE gno) {
  int high, low, mid;

  HIER_CHECK_GNO_RANGE(gno);

  if (hpp->num_gnos_of_interest == 0) return -1;

  low = 0;
  high = hpp->num_gnos_of_interest;
  mid = high/2;

  while (low!=mid && hpp->gnos_of_interest[mid] != gno) {
    if (hpp->gnos_of_interest[mid] > gno) {
      high = mid;
    }
    else {
      low = mid;
    }
    mid = (high+low)/2;
  }
  /* did we find it? */
  if (hpp->gnos_of_interest[mid] == gno) 
    return hpp->gnos_of_interest_procs[mid];

  return -1;
}

/*****************************************************************************/
/******** Internal hierarchical balancing structure callback functions *******/
/*****************************************************************************/

/* helper function to look up local_index of a gid */
/* return index into local arrays if gid is in the local range and has
   not been migrated, return -1 otherwise */
static int get_local_index(HierPartParams *hpp, ZOLTAN_GNO_TYPE gno) {

  HIER_CHECK_GNO_RANGE(gno);
  if ((gno >= hpp->vtxdist[hpp->origzz->Proc]) &&
      (gno < hpp->vtxdist[hpp->origzz->Proc+1])) {
    /* it's in our range, just check that it hasn't been migrated */
    if (!hpp->migrated[gno - hpp->vtxdist[hpp->origzz->Proc]]) {
      return (int)(gno - hpp->vtxdist[hpp->origzz->Proc]);
    }
  }

  /* it was migrated or we just don't know anything about it */
  return -1;
}

/* return index into local arrays if gno is in the local range even if
   it has been migrated, return -1 otherwise */
static int get_starting_local_index(HierPartParams *hpp, ZOLTAN_GNO_TYPE gno) {

  HIER_CHECK_GNO_RANGE(gno);
  if ((gno >= hpp->vtxdist[hpp->origzz->Proc]) &&
      (gno < hpp->vtxdist[hpp->origzz->Proc+1])) {
    return (int)(gno - hpp->vtxdist[hpp->origzz->Proc]);
  }
  
  /* it didn't start here */
  return -1;
}

/* helper function to determine current proc assignment of a gno - gno
   must either be local or be in the gnos_of_interest as computed by
   find_needed_gno_procs */
static int current_proc_of_gno(HierPartParams *hpp, ZOLTAN_GNO_TYPE gno) {
  int proc;
  int dd_proc;
  char msg[64];
  
  HIER_CHECK_GNO_RANGE(gno);

  for (proc=0; proc<hpp->origzz->Num_Proc; proc++) {
    if ((gno >= hpp->vtxdist[proc]) && (gno < hpp->vtxdist[proc+1])) {
      /* if the gno started here, and has not migrated, return proc */
      if ((proc == hpp->origzz->Proc) && 
	  !hpp->migrated[gno - hpp->vtxdist[proc]]) {
	return proc;
      }
      /* is it migrated in? */
      else if (migrated_in_index_of_gno(hpp, gno) != -1) {
	return hpp->origzz->Proc;
      }
      else {
	/* check for it in the list computed from the dd */
	dd_proc = gno_of_interest_proc(hpp, gno);
	/* if we got -1, it was not in the list and must be in its
	   original location */
	if (dd_proc == -1) return proc;
	return dd_proc;
      }
    }
  }
  /* should never get here */
  sprintf(msg, "GNO " ZOLTAN_GNO_SPEC " not found", gno);
  ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, "current_proc_of_gno", msg);
  return -1;	       
}

/* helper function to convert a rank in origzz's communicator to the
   corresponding rank in hierzz's communicator

   returns -1 if proc is not in hierzz's communicator 
*/
static int orig_rank_to_hier_rank(HierPartParams *hpp, int origrank) {

  if (hpp->hier_ranks_of_orig[origrank] != MPI_UNDEFINED) 
    return hpp->hier_ranks_of_orig[origrank];
  return -1;
}


/* the actual callbacks */
static int Zoltan_Hier_Num_Obj_Fn(void *data, int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;

  *ierr = ZOLTAN_OK;
  return hpp->hier_num_obj;
}

static void Zoltan_Hier_Obj_List_Fn(void *data, int num_gid_entries,
				   int num_lid_entries, 
				   ZOLTAN_ID_TYPE * global_ids, 
				   ZOLTAN_ID_TYPE  * local_ids, 
				   int wgt_dim, float *obj_wgts, int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;
  int objindex, j, k;
  float *vwgts;
  char *yo = "Zoltan_Hier_Obj_List_Fn";

  if ((num_gid_entries != 1) || (num_lid_entries != 0)) {
    ZOLTAN_PRINT_ERROR(hpp->hierzz->Proc, yo, "Expected 1 GID, 0 LID");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  /* fill in the GIDs and vertex weights */
  objindex=0;
  for (j=0; j<hpp->init_num_obj; j++) {
    if (!hpp->migrated[j]) {
      global_ids[objindex] = (ZOLTAN_ID_TYPE)hpp->vtxdist[hpp->origzz->Proc]+j;
      for (k=0; k<wgt_dim; k++) 
	obj_wgts[wgt_dim*objindex+k] = 
	  hpp->vwgt[wgt_dim*objindex+k];
      objindex++;
    }
  }

  /* add in objects that have been migrated to here */
  for (j=0; j<hpp->num_migrated_in_gnos; j++) {
    global_ids[objindex] = (ZOLTAN_ID_TYPE)hpp->migrated_in_gnos[j];
    vwgts = get_hier_mig_vwgts(hpp->migrated_in_data[j]);
    for (k=0; k<wgt_dim; k++) {
      obj_wgts[wgt_dim*objindex+k] = vwgts[k];
    }
    objindex++;
  }
  
  /* make sure things make sense - we should have added exactly
     hpp.hier_num_obj objects to the GID array */
  if (objindex != hpp->hier_num_obj) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
		       "Logic error: inconsistency in number of objects");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  *ierr = ZOLTAN_OK;
}

static int Zoltan_Hier_Num_Geom_Fn(void *data, int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;

  *ierr = ZOLTAN_OK;
  return hpp->ndims;
}

static void Zoltan_Hier_Geom_Fn(void *data, int num_gid_entries, 
				int num_lid_entries, 
				ZOLTAN_ID_TYPE * global_id, 
				ZOLTAN_ID_TYPE * local_id,
				double *coord, int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;
  char *yo = "Zoltan_Hier_Geom_Fn";
  ZOLTAN_ID_TYPE gid = global_id[0];
  int local_index;
  int i;
  double *coord_ptr;

  HIER_CHECK_GNO_RANGE((ZOLTAN_GNO_TYPE)gid);

  /* everything is OK unless we don't find it */
  *ierr = ZOLTAN_OK;

  /* make sure we have geometric info */
  if (!hpp->use_geom) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
		       "GEOM callbacks not specified, geom info not available");
    *ierr = ZOLTAN_FATAL;
    return;
  }

  /* check for an object that started locally */
  local_index = get_local_index(hpp, (ZOLTAN_GNO_TYPE)gid);
  if (local_index != -1) {
    for (i=0; i<hpp->ndims; i++) {
      coord[i] = hpp->geom_vec[hpp->ndims*local_index+i];
    }
    return;
  }

  /* check for an object that has been migrated in in a previous step */
  local_index = migrated_in_index_of_gno(hpp, (ZOLTAN_GNO_TYPE)gid);
  if (local_index != -1) {
    coord_ptr = get_hier_mig_coords(hpp->migrated_in_data[local_index], hpp);
    for (i=0; i<hpp->ndims; i++) {
      coord[i] = coord_ptr[i];
    }
    return;
  }

  ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "GID not found");
  *ierr = ZOLTAN_FATAL;
}

static int Zoltan_Hier_Num_Edges_Fn(void *data, int num_gid_entries, 
				    int num_lid_entries, 
				    ZOLTAN_ID_TYPE * global_id, 
				    ZOLTAN_ID_TYPE  * local_id,
				    int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;
  ZOLTAN_ID_TYPE gid = global_id[0];
  int local_index;
  char *yo = "Zoltan_Hier_Num_Edges_Fn";
  int edge, count_edges, edge_proc_orig_rank, num_adj; 
  ZOLTAN_GNO_TYPE next_adj;
  ZOLTAN_GNO_TYPE gno;

  gno = (ZOLTAN_GNO_TYPE)gid;

  HIER_CHECK_GNO_RANGE(gno);

  /* everything is OK unless we don't find it */
  *ierr = ZOLTAN_OK;

  /* make sure we have graph info */
  if (!hpp->use_graph) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
		       "Graph callbacks not specified, edges not available");
    *ierr = ZOLTAN_FATAL;
    return 0;
  }

  /* check for an object that started locally */
  local_index = get_local_index(hpp, gno);
  if (local_index != -1) {
    /* need to be careful that the edges are to objects located on a
       proc in the current hierzz->Comm communicator */
    count_edges = 0;
    for (edge = 0; edge < hpp->xadj[local_index+1] - hpp->xadj[local_index]; 
	 edge++) {
      edge_proc_orig_rank = 
	current_proc_of_gno(hpp, hpp->adjncy[hpp->xadj[local_index]+edge]);
      if (orig_rank_to_hier_rank(hpp, edge_proc_orig_rank) != -1) {
	count_edges++;
      }
    }
    return count_edges;
  }

  /* check for an object that has been migrated in in a previous step */
  local_index = migrated_in_index_of_gno(hpp, gno);
  if (local_index != -1) {
    count_edges = 0;
    num_adj = get_hier_mig_num_adj(hpp->migrated_in_data[local_index]);
    for (edge = 0; edge < num_adj; edge++) {
      next_adj = get_hier_mig_ith_adj(hpp->migrated_in_data[local_index], edge);
      edge_proc_orig_rank = current_proc_of_gno(hpp, next_adj);
      if (orig_rank_to_hier_rank(hpp, edge_proc_orig_rank) != -1) {
	count_edges++;
      }
    }
    return count_edges;
  }

  ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "GID not found");
  *ierr = ZOLTAN_FATAL;
  return 0;
}

static void Zoltan_Hier_Edge_List_Fn(void *data, int num_gid_entries, 
				     int num_lid_entries, 
				     ZOLTAN_ID_TYPE * global_id,
				     ZOLTAN_ID_TYPE * local_id, 
				     ZOLTAN_ID_TYPE * nbor_global_id, 
				     int *nbor_procs,
				     int wgt_dim, float *ewgts, int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;
  ZOLTAN_ID_TYPE gid = global_id[0];
  int local_index, offset, nadj;
  int i, j, edge, edge_proc_orig_rank, edge_proc_hier_rank, num_adj;
  ZOLTAN_GNO_TYPE next_adj, gno;
  char *yo = "Zoltan_Hier_Edge_List_Fn";

  gno = (ZOLTAN_GNO_TYPE)gid;

  HIER_CHECK_GNO_RANGE(gno);

  /* everything is OK unless we don't find it */
  *ierr = ZOLTAN_OK;
  
  /* make sure we have graph info */
  if (!hpp->use_graph) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
		       "Graph callbacks not specified, edges not available");
    *ierr = ZOLTAN_FATAL;
    return;
  }
  
  /* check for an object that started locally */
  local_index = get_local_index(hpp, gno);
  if (local_index != -1) {
    edge = 0;
    offset = hpp->xadj[local_index];
    nadj = hpp->xadj[local_index+1] - offset;
    for (i=0; i<nadj; i++,offset++) {
      edge_proc_orig_rank = current_proc_of_gno(hpp, hpp->adjncy[offset]);
      edge_proc_hier_rank = orig_rank_to_hier_rank(hpp, edge_proc_orig_rank);
      if (edge_proc_hier_rank != -1) {
	nbor_global_id[edge] = (ZOLTAN_ID_TYPE)hpp->adjncy[offset];
	nbor_procs[edge] = edge_proc_hier_rank;
	for (j=0; j<wgt_dim; j++) {
	  ewgts[edge*wgt_dim+j] = hpp->ewgts[offset*wgt_dim+j];
	}
	edge++;
      }
      /* would be good to verify that we found the same number of edges 
	 as we expected */
    }

    return;
  }

  /* check for an object that has been migrated in in a previous step */
  /* it would be nice to factor out some code common to this and the 
     loop above */
  local_index = migrated_in_index_of_gno(hpp, gno);
  if (local_index != -1) {
    edge = 0;
    num_adj = get_hier_mig_num_adj(hpp->migrated_in_data[local_index]);
    for (i=0; i<num_adj; i++) {
      next_adj = get_hier_mig_ith_adj(hpp->migrated_in_data[local_index], i);
      edge_proc_orig_rank = current_proc_of_gno(hpp, next_adj);
      edge_proc_hier_rank = orig_rank_to_hier_rank(hpp, edge_proc_orig_rank);
      if (edge_proc_hier_rank != -1) {
	nbor_global_id[edge] = (ZOLTAN_ID_TYPE)next_adj;
	nbor_procs[edge] = edge_proc_hier_rank;
	get_hier_mig_ith_adj_wgts(hpp->migrated_in_data[local_index],
				  hpp, i, &(ewgts[edge*wgt_dim]));
	edge++;
      }
    }

    return;
  }

  ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "GID not found");
  *ierr = ZOLTAN_FATAL;
}

/* migration callbacks */

/* migration helpers */

/* find the number of edges of a gid that we know is local */
static int num_graph_edges_of_gno(HierPartParams *hpp, ZOLTAN_GNO_TYPE gno,
				  int *count) {
  int local_index;
  char *yo = "num_graph_edges_of_gno";
  char msg[256];

  HIER_CHECK_GNO_RANGE(gno);

  *count = 0;

  local_index = get_starting_local_index(hpp, gno);
  if (local_index != -1) {
    /* it started here so its edges are here */
    *count = hpp->xadj[local_index+1] - hpp->xadj[local_index];
    return ZOLTAN_OK;
  }

  /* look it up in our migrated in lists */
  local_index = migrated_in_index_of_gno(hpp, gno);
  if (local_index != -1) {
    *count = get_hier_mig_num_adj(hpp->migrated_in_data[local_index]);
    return ZOLTAN_OK;
  }

  sprintf(msg, "could not find gno " ZOLTAN_GNO_SPEC "\n", gno);
  ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, msg);
  return ZOLTAN_FATAL;
}

static void Zoltan_Hier_Check_Data(HierPartParams *hpp, int *ierr) {

  char *yo="Zoltan_Hier_Check_Data";
  int i;
  struct Zoltan_Comm_Obj *plan=NULL;
  int nreturn=0;
  int *proclist=NULL;
  int proc;
  ZOLTAN_GNO_TYPE *sendbuf=NULL, *recvbuf=NULL;
  int local_index;
  int *owners=NULL;
  ZOLTAN_GNO_TYPE *ddlookup=NULL;
  char msg[256];

  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] %s called\n", hpp->origzz->Proc, yo);
  }

  /* first make sure the migrated_in_gnos array is sorted and
     that each entry has a corresponding data entry */
  /* check entry 0 */
  if (hpp->num_migrated_in_gnos && !hpp->migrated_in_data[0]) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo,
		       "migrated_in_data[0] is NULL");
    *ierr=ZOLTAN_FATAL;
    goto End;
  }
  /* check the rest */
  for (i=1; i<hpp->num_migrated_in_gnos; i++) {
    if (hpp->migrated_in_gnos[i] <= hpp->migrated_in_gnos[i-1]) {
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo,
			 "migrated_in_gnos not sorted properly");
      *ierr=ZOLTAN_FATAL;
      goto End;
    }
    if (!hpp->migrated_in_data[i]) {
      sprintf(msg, "migrated_in_data[%d] is NULL\n", i);
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, msg);
      *ierr=ZOLTAN_FATAL;
      goto End;
    }
  }

  /* send data about each migrated_in_gid to its original owner to make
     sure everything is consistent */

  /* first populate a proclist to create a communication plan */
  if (hpp->num_migrated_in_gnos) {
    proclist = (int *)ZOLTAN_MALLOC(hpp->num_migrated_in_gnos*sizeof(int));
    if (!proclist){
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "out of memory");
      *ierr=ZOLTAN_MEMERR;
    }
    for (i=0; i<hpp->num_migrated_in_gnos; i++) {
      for (proc=0; proc<hpp->origzz->Num_Proc; proc++) {
	if ((hpp->migrated_in_gnos[i] >= hpp->vtxdist[proc]) && 
	    (hpp->migrated_in_gnos[i] < hpp->vtxdist[proc+1])) {
	  proclist[i] = proc;
	  break;
	}
      }
    }
  }

  /* create the plan */
  *ierr = Zoltan_Comm_Create(&plan, hpp->num_migrated_in_gnos, proclist,
			     hpp->origzz->Communicator, 0, &nreturn);

  /* allocate space for the information we'll send and receive */
  if (hpp->num_migrated_in_gnos) {
    sendbuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(2*hpp->num_migrated_in_gnos*sizeof(ZOLTAN_GNO_TYPE));
    if (!sendbuf){
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "out of memory");
      *ierr=ZOLTAN_MEMERR;
    }
    for (i=0; i<hpp->num_migrated_in_gnos; i++) {
      sendbuf[2*i] = hpp->migrated_in_gnos[i];
      sendbuf[2*i+1] = (ZOLTAN_GNO_TYPE)hpp->origzz->Proc;  /* assumes proc ID can fit in an ZOLTAN_GNO_TYPE */
    }
  }

  if (nreturn) {
    recvbuf = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(2*nreturn*sizeof(ZOLTAN_GNO_TYPE));
    if (!recvbuf){
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "out of memory");
      *ierr=ZOLTAN_MEMERR;
    }
  }

  /* do the communication */
  *ierr = Zoltan_Comm_Do(plan, 0, (char *)sendbuf, 2*sizeof(ZOLTAN_GNO_TYPE), 
			 (char *)recvbuf);

  /* check to make sure each received gid lives here, and that we know
     it's been migrated.  We make sure each migrated gid was sent to
     us somewhere by temporarily making the migrated flag a 2 instead
     of a 1 */
  for (i=0; i<nreturn; i++) {
    local_index = get_starting_local_index(hpp, recvbuf[2*i]);
    if (hpp->migrated[local_index] == 0) {
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo,
			 "received gid that was not migrated!");
      *ierr = ZOLTAN_FATAL;
      goto End;
    }
    else if (hpp->migrated[local_index] == 1) {
      /* mark as migrated and reported back in */
      hpp->migrated[local_index] = 2;
    } else {
      /* we already marked it! */
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo,
			 "received gid twice!");
      *ierr = ZOLTAN_FATAL;
      goto End;
    }
  }

  /* check that no 1's remain in migrated array, put 2's back to 1's while
     we're in the loop */
  for (i=0; i<hpp->init_num_obj; i++) {
    if (hpp->migrated[i] == 1) {
      sprintf(msg, "gid " ZOLTAN_GNO_SPEC " at pos %d marked as migrated, not reported back!",
	      hpp->vtxdist[hpp->origzz->Proc]+i, i);
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, msg);
      *ierr = ZOLTAN_FATAL;
      goto End;
    }
    if (hpp->migrated[i] == 2) hpp->migrated[i] = 1;
  }

  /* check that the information we received about current locations is
     consistent with what is in the DD */
  if (nreturn) {
    ddlookup = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(nreturn*sizeof(ZOLTAN_GNO_TYPE));   
    owners = (int *)ZOLTAN_MALLOC(nreturn*sizeof(int));
    if (!ddlookup || !owners){
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "out of memory");
      *ierr=ZOLTAN_MEMERR;
    }
    for (i=0; i<nreturn; i++) {
      ddlookup[i]=(ZOLTAN_ID_TYPE)recvbuf[2*i];
    }
  }

  *ierr = Zoltan_DD_Find(hpp->dd, (ZOLTAN_ID_TYPE *)ddlookup, NULL, NULL, NULL, nreturn, owners);

  for (i=0; i<nreturn; i++) {
    if (owners[i] != (int)recvbuf[2*i+1]) {
      sprintf(msg, "Owner mismatch for GNO " ZOLTAN_GNO_SPEC ": DD has %d, msg came from " ZOLTAN_GNO_SPEC "\n",
	      recvbuf[2*i], owners[i], recvbuf[2*i+1]);
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, msg);
      *ierr = ZOLTAN_FATAL;
      goto End;
    }
  }

End:
  /* clean up memory */
  if (sendbuf) ZOLTAN_FREE(&sendbuf);
  if (recvbuf) ZOLTAN_FREE(&recvbuf);
  if (proclist) ZOLTAN_FREE(&proclist);
  if (ddlookup) ZOLTAN_FREE(&ddlookup);
  if (owners) ZOLTAN_FREE(&owners);

  Zoltan_Comm_Destroy(&plan);
}

/* callback to determine migration size */
static int Zoltan_Hier_Obj_Size_Fn(void *data, 
				   int num_gid_entries, int num_lid_entries, 
				   ZOLTAN_ID_TYPE * global_id, 
				   ZOLTAN_ID_TYPE  * local_id, int *ierr) {
  HierPartParams *hpp = (HierPartParams *)data;
  ZOLTAN_ID_TYPE gid = global_id[0];
  ZOLTAN_GNO_TYPE gno;
  char *yo = "Zoltan_Hier_Obj_Size_Fn";
  int num_bytes = 0;
  int num_edges;

  gno = (ZOLTAN_GNO_TYPE)gid;

  HIER_CHECK_GNO_RANGE(gno);

  /* everything is OK unless we don't find it */
  *ierr = ZOLTAN_OK;

  /* figure out what we need to include in the packed buffer */
  /* object weights */
  num_bytes += hpp->obj_wgt_dim*sizeof(float);
  /* geometry */
  num_bytes += hpp->ndims*sizeof(double);
  /* adjacent edges */
  if (hpp->use_graph) {
    *ierr = num_graph_edges_of_gno(hpp, gno, &num_edges);
    if (*ierr == ZOLTAN_FATAL) {
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
			 "num_graph_edges_of_gno returned error");
      return 0;
    }
    /* 1 int for edge count, 1 int and 1 float per edge */
    num_bytes += sizeof(int) + num_edges*(sizeof(ZOLTAN_GNO_TYPE)+ hpp->edge_wgt_dim*sizeof(float));
  }

  /*
  printf("[%d] pack size of GID " ZOLTAN_GNO_SPEC " is %d with %d edges, %d obj_wgt_dim, %d ndims, %d edge_wgt_dim\n", 
	 hpp->origzz->Proc, gno, num_bytes, num_edges, hpp->obj_wgt_dim, 
	 hpp->ndims,  hpp->edge_wgt_dim);
  */

  return num_bytes;
}

/* callback to pack the buffer with needed information */
/* note that we could avoid all object migration when the object is
   going back to its starting processor, since the original
   intermediate structure remains there.  However,
   Zoltan_Hier_Obj_Size_Fn has already allocated space and that
   callback does not take a desintation processor as an argument. */
static void Zoltan_Hier_Pack_Obj_Fn(void *data, 
				    int num_gid_entries, int num_lid_entries,
				    ZOLTAN_ID_TYPE * global_id, 
				    ZOLTAN_ID_TYPE  * local_id, int dest_proc,
				    int size, char *buf, int *ierr) {

  HierPartParams *hpp = (HierPartParams *)data;
  ZOLTAN_ID_TYPE gid = global_id[0];
  char *yo = "Zoltan_Hier_Pack_Obj_Fn";
  int local_index;
  /* pointers to data in local structures */
  float *vwgts = NULL;
  double *coords = NULL;
  int num_adj = 0;
  ZOLTAN_GNO_TYPE *adj = NULL;
  ZOLTAN_GNO_TYPE gno;
  float *ewgts = NULL;
  float *buf_float_ptr;
  ZOLTAN_GNO_TYPE * buf_gno_ptr;
  double *buf_double_ptr;
  int *buf_int_ptr;
  int buf_index, i;
 
  gno = (ZOLTAN_GNO_TYPE)gid;

  HIER_CHECK_GNO_RANGE(gno);

  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] packing GID " ZOLTAN_GNO_SPEC " to send to %d, size=%d\n", hpp->origzz->Proc, 
	   gno, dest_proc, size);
  }

  /* assume failure -- set to OK when we find the gid */
  *ierr = ZOLTAN_FATAL;

  /* check for an object that started locally */
  local_index = get_local_index(hpp, gno);
  if (local_index != -1) {
    *ierr = ZOLTAN_OK;
    if (hpp->obj_wgt_dim) {
      vwgts = &(hpp->vwgt[local_index*hpp->obj_wgt_dim]);
    }
    if (hpp->use_geom) {
      coords = &(hpp->geom_vec[local_index*hpp->ndims]);
    }
    if (hpp->use_graph) {
      num_adj = hpp->xadj[local_index+1] - hpp->xadj[local_index];
      adj = &(hpp->adjncy[hpp->xadj[local_index]]);
      if (hpp->edge_wgt_dim) {
	ewgts = &(hpp->ewgts[hpp->xadj[local_index]*hpp->edge_wgt_dim]);
      }
    }
  }
  else {
    /* check for an object that has been migrated in in a previous step */
    local_index = migrated_in_index_of_gno(hpp, gno);
    if (local_index != -1) {
      *ierr = ZOLTAN_OK;
      if (hpp->obj_wgt_dim) {
	vwgts = get_hier_mig_vwgts(hpp->migrated_in_data[local_index]);
      }
      if (hpp->use_geom) {
	coords = get_hier_mig_coords(hpp->migrated_in_data[local_index], hpp);
      }
      if (hpp->use_graph) {
	num_adj = get_hier_mig_num_adj(hpp->migrated_in_data[local_index]);
	adj = get_hier_mig_adj(hpp->migrated_in_data[local_index]);
	if (hpp->edge_wgt_dim) {
	  ewgts = get_hier_mig_adj_wgts(hpp->migrated_in_data[local_index],
					hpp);
	}
      }
    }
  }

  if (*ierr == ZOLTAN_OK) {
    /* we found the GID and have pointers to its info - pack into the buffer */
    /* this is not pretty but hopefully it is safe */
    /* we pack in num_adj, all of the floats we have, all doubles, then all of
       the ints */
    if (hpp->use_graph) {
      buf_int_ptr = (int *)buf;
      buf_int_ptr[0] = num_adj;
      
      buf_float_ptr = (float *)&buf_int_ptr[1];
    }
    else {
      buf_float_ptr = (float *)buf;
    }
    buf_index = 0;

    /* pack in object weights */
    for (i=0; i<hpp->obj_wgt_dim; i++) {
      buf_float_ptr[buf_index++] = vwgts[i];
    }

    /* pack in edge weights */
    for (i=0; i<hpp->edge_wgt_dim*num_adj; i++) {
      buf_float_ptr[buf_index++] = ewgts[i];
    }

    /* pack in coordinates */
    buf_double_ptr = (double *)&buf_float_ptr[buf_index];
    buf_index = 0;
    for (i=0; i<hpp->ndims; i++) {
      buf_double_ptr[buf_index++] = coords[i];
    }

    if (hpp->use_graph) {
      buf_gno_ptr = (ZOLTAN_GNO_TYPE *)&buf_double_ptr[buf_index];
      buf_index = 0;

      /* pack in adjacent GIDs */
      for (i=0; i<num_adj; i++) {
	HIER_CHECK_GNO_RANGE(adj[i]);
	buf_gno_ptr[buf_index++] = adj[i];
      }
      
    }

    /* some sort of sanity check on the size seems appropriate here */
    /* it seems that these pointers should be the same if everything is OK */
    /* once this code is trusted, this should probably be able to go away */
    /*if (Zoltan_Align((int)&buf_gno_ptr[buf_index]) != 
	Zoltan_Align((int)&buf[size])) {
      *ierr = ZOLTAN_FATAL;
      sprintf(msg, "data inserted (%d) not same size as buffer allocated (%d)",
	      (void *)&buf_gno_ptr[buf_index], (void *)&buf[size]);
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, msg);
      return;
    }*/
    
  } else {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "GID not found");
  }
}


/* callback to unpack an arriving buffer */
static void Zoltan_Hier_Unpack_Obj_Fn(void *data, 
				      int num_gid_entries,
				      ZOLTAN_ID_TYPE * global_id, 
				      int size, char *buf, int *ierr) {

  HierPartParams *hpp = (HierPartParams *)data;
  int local_index;
  float *buf_float_ptr;
  ZOLTAN_GNO_TYPE * buf_gno_ptr;
  ZOLTAN_GNO_TYPE gno;
  double *buf_double_ptr;
  int *buf_int_ptr;
  int buf_index, i;
  struct HierGNOInfo *info = NULL;
  char *yo = "Zoltan_Hier_Unpack_Obj_Fn";

  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] unpacking GID " ZOLTAN_ID_SPEC "\n", hpp->origzz->Proc, global_id[0]);
  }

  gno = (ZOLTAN_GNO_TYPE)global_id[0];

  /* everything is OK unless we don't find it */
  *ierr = ZOLTAN_OK;
  
  /* is this an object returning home to where it started?  If so, we
     can ignore the buffer.  The mid-migration callback should have
     already unmarked its migrated flag */
  if (get_starting_local_index(hpp, gno) != -1) {
    if (hpp->output_level >= HIER_DEBUG_ALL) {
      printf("[%d] GID " ZOLTAN_GNO_SPEC " returns home, no unpack needed\n",
	     hpp->origzz->Proc, gno);
    }
    return;
  }

  /* find the slot for our gid in the migrated_in array */
  local_index = migrated_in_index_of_gno(hpp, gno);
  if (local_index != -1) {
    if (hpp->migrated_in_data[local_index]) {
      *ierr = ZOLTAN_FATAL;
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo,
			 "Data slot for unpack GID occupied!");
      return;
    }

    /* allocate a struct HierGNOInfo to hold GID info */
    info = (struct HierGNOInfo *)ZOLTAN_MALLOC(sizeof(struct HierGNOInfo));
    if (!info) {
      *ierr = ZOLTAN_MEMERR;
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory.");
      goto End;
    }
    if (hpp->output_level >= HIER_DEBUG_ALL) {
      printf("[%d] unpacking GID " ZOLTAN_GNO_SPEC " into slot %d\n", hpp->origzz->Proc, gno,
	     local_index);
    }
    hpp->migrated_in_data[local_index] = (void *)info;
    

    /* allocate and unpack the things we need */

    /*    info->num_adj = 
      (size - hpp->obj_wgt_dim*sizeof(float) - hpp->ndims*sizeof(double)) /
      (hpp->edge_wgt_dim*sizeof(float) + sizeof(ZOLTAN_ID_TYPE)); */

    /* num_adj */

    if (hpp->use_graph) {
      buf_int_ptr = (int *)buf;
      info->num_adj = buf_int_ptr[0];

      buf_float_ptr = (float *)&buf_int_ptr[1];
    }
    else {
      info->num_adj = 0;
      buf_float_ptr = (float *)buf;
    }

    /* float array */
    buf_index = 0;

    if (hpp->obj_wgt_dim+hpp->edge_wgt_dim*info->num_adj) {
      info->wgts = 
	(float *)ZOLTAN_MALLOC((hpp->obj_wgt_dim+
				hpp->edge_wgt_dim*info->num_adj)*
			       sizeof(float));
      if (!info->wgts) {
	*ierr = ZOLTAN_MEMERR;
	ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory.");
	goto End;
      }

      for (i=0; i<hpp->obj_wgt_dim; i++) {
	info->wgts[i] = buf_float_ptr[buf_index++];
      }
      for (i=0; i<hpp->edge_wgt_dim*info->num_adj; i++) {
	info->wgts[i+hpp->obj_wgt_dim] = buf_float_ptr[buf_index++];
      }
    }
    else info->wgts = NULL;
    
    /* coordinates? */
    buf_double_ptr = (double *)&buf_float_ptr[buf_index];
    buf_index = 0;
    if (hpp->ndims) {
      info->coords = (double *)ZOLTAN_MALLOC(hpp->ndims*sizeof(double));
      if (!info->coords) {
	*ierr = ZOLTAN_MEMERR;
	ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory.");
	goto End;
      }

      for (i=0; i<hpp->ndims; i++) {
	info->coords[i] = buf_double_ptr[buf_index++];
      }
    }
    else info->coords = NULL;

    /* adjacent gids? */
    buf_gno_ptr = (ZOLTAN_GNO_TYPE *)&buf_double_ptr[buf_index];
    buf_index = 0;
    if (info->num_adj) {
      info->adj = (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(info->num_adj*sizeof(ZOLTAN_GNO_TYPE));
      if (!info->adj) {
	*ierr = ZOLTAN_MEMERR;
	ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory.");
	goto End;
      }

      for (i=0; i<info->num_adj; i++) {
	info->adj[i] = buf_gno_ptr[buf_index++];
	HIER_CHECK_GNO_RANGE(info->adj[i]);
      }
    }
    else info->adj = NULL;

    /* we unpacked everything */
    return;
  }

  /* if we're here, we didn't find a slot for this GID.  This should
     not happen */
  *ierr = ZOLTAN_FATAL;
  ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Unexpected GID to unpack");
End:
  if (info) free_hier_mig_data(info);
}

/* static helper function used by Zoltan_Hier_Mid_Migrate_Fn */
static int next_insert_index(HierPartParams *hpp, int start) {

  /* always bump up one slot, then skip over occupied spots */
  start++;
  while (hpp->migrated_in_data[start] != NULL) start++;
  return start;
}

/* static helper for quicksort of ints with parallel array of ptrs */
/* based on Zoltan_Comm_Sort_Ints, but that can only do two arrays of
   int.  This could be moved somewhere else if it is more generally
   useful */
static int sort_gno_ptrs(
ZOLTAN_GNO_TYPE * vals_sort,	/* values to be sorted */
void    **vals_other,		/* other array to be reordered w/ sort */
int       nvals)		/* length of these two arrays */
{
  ZOLTAN_GNO_TYPE temp_int;	/* swapping value */
  void     *temp_ptr;		/* swapping value */
  int       lo, hi;		/* counters from bottom and top of array */
  ZOLTAN_GNO_TYPE pivot;		/* value to partition with */
  
  if (nvals <= 1) return(ZOLTAN_OK);
  
  /* Partition */
  lo = nvals/2;
  if (lo == nvals - 1) --lo;
  pivot = vals_sort[lo];
  
  lo = -1;
  hi = nvals;
  
  while (lo < hi) {
    do {
      hi--;
    } while (vals_sort[hi] > pivot);
    
    do {
      lo++;
    } while (vals_sort[lo] < pivot);
    
    if (lo < hi) {	/* Swap low and high items */
      temp_int = vals_sort[lo];
      vals_sort[lo] = vals_sort[hi];
      vals_sort[hi] = temp_int;

      temp_ptr = vals_other[lo];
      vals_other[lo] = vals_other[hi];
      vals_other[hi] = temp_ptr;
    }
  } 
  
  /* Recurse */
  if (hi + 1 > 1) sort_gno_ptrs(vals_sort, vals_other, hi + 1);
  if (nvals - hi - 1 > 1)
    sort_gno_ptrs(&vals_sort[hi + 1], &vals_other[hi + 1], nvals - hi - 1);
  
  return(ZOLTAN_OK);
}

/* callback to be done between packing for exports and unpacking of imports */
/* this is not called for the last global level - at that time, a subset of
   the functionality here is done in Zoltan_Hier just after the call to
   Zoltan_LB_Partition for the last level */
static void Zoltan_Hier_Mid_Migrate_Fn(void *data, int num_gid_entries,
				       int num_lid_entries, int num_import,
				       ZOLTAN_ID_TYPE * import_gids,
				       ZOLTAN_ID_TYPE * import_lids,
				       int *import_procs, int *import_parts,
				       int num_export,
				       ZOLTAN_ID_TYPE * export_gids,
				       ZOLTAN_ID_TYPE * export_lids,
				       int *export_procs, int *export_parts,
				       int *ierr) {
  int i, local_index;
  HierPartParams *hpp = (HierPartParams *)data;
  char *yo = "Zoltan_Hier_Mid_Migrate_Fn";
  int removed_from_migrated_in;
  int add_to_migrated_in;
  /*int export_partition_change_only;*/
  int import_partition_change_only;
  int insert_index = 0, trail_index;
  ZOLTAN_GNO_TYPE * dd_updates=NULL;
  int remove_count = 0, update_count = 0;
  int last_level, idx;

  /* everything is OK until we find otherwise */
  *ierr = ZOLTAN_OK;

  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] mid-migrate with %d imports, %d exports\n", 
	   hpp->origzz->Proc, num_import, num_export);
  }

  /* check if this is the last level, if so, we only update the
     migrated arrays for GIDs that started locally and make any
     needed updates in the DD to find final locations of migrated
     GIDs */
  last_level = (hpp->level == hpp->global_num_levels - 1);

  /* do exports first, since we can make some space in the migrated in
     arrays and avoid some memory usage */

  /* remove stuff from the migrated in arrays, mark as migrated in 
     orig arrays */
  removed_from_migrated_in = 0;
  /*export_partition_change_only = 0;*/
  for (i=0; i<num_export; i++) {
    HIER_CHECK_GNO_RANGE((ZOLTAN_GNO_TYPE)export_gids[i]);
    /* is this just a partition reassignment and not a real migration? */
    if (export_procs[i] == hpp->hierzz->Proc) {
      /*export_partition_change_only++;*/
      if (hpp->output_level >= HIER_DEBUG_ALL) {
	printf("[%d] mid-migrate exporting GID " ZOLTAN_ID_SPEC " - partition change only\n", 
	       hpp->origzz->Proc, export_gids[i]);
      }
      /* skip to next export */
      continue;
    }

    /* see if the exported GID is originally local */
    local_index = get_local_index(hpp,(ZOLTAN_GNO_TYPE)export_gids[i]);
    if (local_index != -1) {
      /* it started here, so we just need to mark it as migrated */
      if (hpp->output_level >= HIER_DEBUG_ALL) {
	printf("[%d] mid-migrate exporting GID " ZOLTAN_ID_SPEC " to %d, started here\n", 
	       hpp->origzz->Proc, export_gids[i], export_procs[i]);
      }
      hpp->migrated[local_index] = 1;
    }
    else {
      if (!last_level) {
	local_index = migrated_in_index_of_gno(hpp, (ZOLTAN_GNO_TYPE)export_gids[i]);
	if (local_index != -1) {
	  /* it's in our migrated in arrays, so we can get rid of it.
	     To do so, we tag the entry by freeing the migrated_in_data
	     entry and setting it to NULL.  The actual gid entries in
	     migrated_in_gnos need to remain until we are done with this
	     loop to make sure the binary search for
	     migrated_in_index_of_gno doesn't break */
	  if (hpp->output_level >= HIER_DEBUG_ALL) {
	    printf("[%d] mid-migrate exporting GID " ZOLTAN_ID_SPEC " to %d, was imported and was in slot %d\n",
		   hpp->origzz->Proc, export_gids[i], export_procs[i],
		   local_index);
	  }
	  free_hier_mig_data(hpp->migrated_in_data[local_index]);
	  hpp->migrated_in_data[local_index] = NULL;
	  removed_from_migrated_in++;
	}
	else {
	  *ierr = ZOLTAN_FATAL;
	  ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
			     "cannot find gid for export");
	  return;
	}
      }
    }
  }

  /* find out how many entries to migrated in arrays we will have */
  /* count imports that we need to store, ignore those that are coming
     back here to their original location */
  add_to_migrated_in = 0;
  import_partition_change_only = 0;
  for (i=0; i<num_import; i++) {
    HIER_CHECK_GNO_RANGE((ZOLTAN_GNO_TYPE)import_gids[i]);
    /* check for a partition change only */
    if (import_procs[i] == hpp->hierzz->Proc) {
      import_partition_change_only++;
      if (hpp->output_level >= HIER_DEBUG_ALL) {
	printf("[%d] mid-migrate importing GID " ZOLTAN_ID_SPEC " - partition change only\n", 
	       hpp->origzz->Proc, import_gids[i]);
      }
      /* skip to next import */
      continue;
    }
    if (!last_level) {
      local_index = get_starting_local_index(hpp, (ZOLTAN_GNO_TYPE)import_gids[i]);
      if (local_index == -1) add_to_migrated_in++;
    }
  }

  if (!last_level) {
    /* is there enough space using the holes we created from items being
       migrated out and other allocated space? */
    if (hpp->num_migrated_in_gnos - removed_from_migrated_in + 
	add_to_migrated_in > hpp->alloc_migrated_in_gnos) {
      /* reallocate arrays */
      hpp->alloc_migrated_in_gnos = hpp->num_migrated_in_gnos - 
	removed_from_migrated_in + add_to_migrated_in;
      if (hpp->migrated_in_gnos) {
	hpp->migrated_in_gnos = 
	  (ZOLTAN_GNO_TYPE *)ZOLTAN_REALLOC(hpp->migrated_in_gnos,
					sizeof(ZOLTAN_GNO_TYPE)* hpp->alloc_migrated_in_gnos);
	hpp->migrated_in_data = 
	  (void **)ZOLTAN_REALLOC(hpp->migrated_in_data,
				  sizeof(void *)*hpp->alloc_migrated_in_gnos);
      }
      else {
	hpp->migrated_in_gnos =
	  (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE)*
				       hpp->alloc_migrated_in_gnos);
	hpp->migrated_in_data = 
	  (void **)ZOLTAN_MALLOC(sizeof(void *)*hpp->alloc_migrated_in_gnos);
      }
      if (!hpp->migrated_in_gnos || !hpp->migrated_in_data) {
	*ierr = ZOLTAN_MEMERR;
	ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
	return;
      }
    }

    /* make sure we have NULLs in any expansion of hpp->migrated_in_data */
    for (i=hpp->num_migrated_in_gnos; i<hpp->alloc_migrated_in_gnos; i++) {
      hpp->migrated_in_data[i] = NULL;
    }
    
    insert_index = -1;
  }

  if (num_import - import_partition_change_only) {
    dd_updates = 
      (ZOLTAN_GNO_TYPE *)ZOLTAN_MALLOC(sizeof(ZOLTAN_GNO_TYPE)*
				   (num_import-import_partition_change_only));
    if (!dd_updates) {
      *ierr = ZOLTAN_MEMERR;
      ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, "Out of memory");
      return;
    }
  }

  idx = num_import - import_partition_change_only;

  remove_count = 0; update_count = 0;

  for (i=0; i<num_import; i++) {
    /* check for a partition change only */
    if (import_procs[i] == hpp->hierzz->Proc) {
      continue;
    }
    local_index = get_starting_local_index(hpp, (ZOLTAN_GNO_TYPE)import_gids[i]);
    /* did it start its life here? */
    if (local_index != -1) {
      if (hpp->output_level >= HIER_DEBUG_ALL) {
	printf("[%d] mid-migrate importing GID " ZOLTAN_ID_SPEC " from %d, started here\n",
	       hpp->origzz->Proc, import_gids[i], import_procs[i]);
      }
      hpp->migrated[local_index] = 0;
      /* remove from DD - it's back home */

      dd_updates[--idx] = (ZOLTAN_GNO_TYPE)import_gids[i];
   
      remove_count++;
    }
    else {
      /* update location in the DD */
      dd_updates[update_count] = (ZOLTAN_GNO_TYPE)import_gids[i];
      update_count++;
      if (!last_level) {
	/* put gid in the next open slot in the migrated_in array */
	insert_index = next_insert_index(hpp, insert_index);
	hpp->migrated_in_gnos[insert_index] = (ZOLTAN_GNO_TYPE)import_gids[i];
	if (hpp->output_level >= HIER_DEBUG_ALL) {
	  printf("[%d] mid-migrate importing GID " ZOLTAN_ID_SPEC " from %d, new arrival into slot %d\n",
		 hpp->origzz->Proc, import_gids[i], import_procs[i], insert_index);
	}
      }
    }
  }

  /* do DD updates */
  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] calling DD_Remove on %d GIDs\n", hpp->origzz->Proc, 
	   remove_count);
    for (i=0; i<remove_count; i++) {
      printf("[%d] DD_Remove slot %d GID " ZOLTAN_GNO_SPEC "\n",
	     hpp->origzz->Proc, i, dd_updates[update_count+i]);
    }
  }
  *ierr = Zoltan_DD_Remove(hpp->dd, 
			   (ZOLTAN_ID_PTR)(remove_count ? &(dd_updates[idx]) : NULL),
			   remove_count);
  if (*ierr != ZOLTAN_OK  && *ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
		       "Zoltan_DD_Remove returned error");
    goto End;
  }
  
  if (hpp->output_level >= HIER_DEBUG_ALL) {
    printf("[%d] calling DD_Update on %d GIDs\n", hpp->origzz->Proc, 
	   update_count);
    for (i=0; i<update_count; i++) {
      printf("[%d] DD_Update slot %d GID " ZOLTAN_GNO_SPEC "\n",
	     hpp->origzz->Proc, i, dd_updates[i]);
    }
  }
  *ierr = Zoltan_DD_Update(hpp->dd, (ZOLTAN_ID_PTR)(update_count ? dd_updates : NULL),
			   NULL, NULL, NULL, update_count);
  if (*ierr != ZOLTAN_OK  && *ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
		       "Zoltan_DD_Update returned error");
    goto End;
  }

  if (!last_level) {
    /* if we removed more from than we added to the migrated in arrays,
       we need to backfill the holes before we sort */
    if (removed_from_migrated_in > add_to_migrated_in) {
      if (hpp->output_level >= HIER_DEBUG_ALL ) {
	printf("[%d] backfilling migrated_in, removed=%d, add=%d\n",
	       hpp->origzz->Proc, removed_from_migrated_in, add_to_migrated_in);
      }
      
      /* first candidate to be moved up to backfill is the original last item */
      trail_index = hpp->num_migrated_in_gnos-1;
      
      for (i=0; i<removed_from_migrated_in-add_to_migrated_in; i++) {
	
	/* if trail_index entry is not already a null, move it to
	   backfill a null, otherwise, we can skip it and move on */
	if (hpp->migrated_in_data[trail_index]) {
	  insert_index = next_insert_index(hpp, insert_index);
	  if (hpp->output_level >= HIER_DEBUG_ALL ) {
	    printf("[%d] backfill from trail=%d to insert=%d, GID " ZOLTAN_GNO_SPEC ", data=%p\n",
		   hpp->origzz->Proc, trail_index, insert_index, 
		   hpp->migrated_in_gnos[trail_index], 
		   hpp->migrated_in_data[trail_index]);
	  }
	  hpp->migrated_in_gnos[insert_index] = 
	    hpp->migrated_in_gnos[trail_index];
	  hpp->migrated_in_data[insert_index] = 
	    hpp->migrated_in_data[trail_index];
	}
	trail_index--;
      }
    }
    
    /* we now should have a packed but unsorted array of gids, update
       the size of the migrated in arrays */
    hpp->num_migrated_in_gnos = 
      hpp->num_migrated_in_gnos - removed_from_migrated_in + add_to_migrated_in;
    
    
    if (hpp->output_level >= HIER_DEBUG_ALL) {
      for (i=0; i<hpp->num_migrated_in_gnos; i++) { 
	printf("[%d] GID " ZOLTAN_GNO_SPEC " in slot %d unsorted\n", hpp->origzz->Proc, 
	       hpp->migrated_in_gnos[i], i);
      }
    }
    
    /* re-sort migrated in arrays */
    if (hpp->num_migrated_in_gnos) {
      *ierr = sort_gno_ptrs(hpp->migrated_in_gnos, hpp->migrated_in_data,
			     hpp->num_migrated_in_gnos);
      if (*ierr != ZOLTAN_OK  && *ierr != ZOLTAN_WARN) {
	ZOLTAN_PRINT_ERROR(hpp->origzz->Proc, yo, 
			   "sort_ints_ptrs returned error");
	goto End;
      }
    }
    if (hpp->output_level >= HIER_DEBUG_ALL) {
      for (i=0; i<hpp->num_migrated_in_gnos; i++) { 
	printf("[%d] GID " ZOLTAN_GNO_SPEC " in slot %d sorted\n", hpp->origzz->Proc,
	       hpp->migrated_in_gnos[i], i);
      }
    }
  }
  hpp->hier_num_obj = hpp->hier_num_obj + num_import - num_export;

End:
  if (dd_updates) ZOLTAN_FREE(&dd_updates);
}

static void view_hierarchy_specification(zoltan_platform_specification *spec, int rank, int verbose)
{
int i;

  if (verbose){
    if (spec->platform_name){
      printf("%s\n",spec->platform_name);
    }
    printf("Number of siblings at each level: ");
    for (i=0; i < spec->numLevels; i++){
      printf("%d ",spec->num_siblings[i]);
    }
    printf("\n");
  }

  printf("Part for MPI rank %d at each level: ", rank);
  for (i=0; i < spec->numLevels; i++){
    printf("%d ",spec->my_part[i]);
  }
  printf("\n");

  fflush(stdout);
}

static char *make_platform_name_string()
{
int i;
int len;
char *msg;
char *yo = "make_platform_name_string";


  for (i=0, len=0; i < ZOLTAN_HIER_LAST_PLATFORM; i++){
    len += strlen(zoltan_hier_platform_specs[i].platform_name);
  }

  len += ((ZOLTAN_HIER_LAST_PLATFORM * 3) + 64);

  msg = (char *)ZOLTAN_MALLOC(len);
  if (!msg){
    ZOLTAN_PRINT_ERROR(-1, yo, "Out of memory");
    return NULL;
  }
  msg[0] = 0;

  for (i=1; i <= ZOLTAN_HIER_LAST_PLATFORM; i++){
    strcat(msg, zoltan_hier_platform_specs[i].platform_name);
    strcat(msg, " ");
    if (i % 5  == 0) strcat(msg, "\n");
  }

  if (ZOLTAN_HIER_LAST_PLATFORM % 5)
    strcat(msg, "\n");

  return msg;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
