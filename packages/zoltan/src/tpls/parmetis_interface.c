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
#include "params_const.h"
#include "order_const.h"
#include "third_library.h"
#include "parmetis_interface.h"


/*********** COMPATIBILITY CHECKING AT COMPILE TIME ************/
#if (PARMETIS_MAJOR_VERSION < 3)
#error "Specified version of ParMETIS is not compatible with Zoltan; upgrade to ParMETIS v3.1 or later, or build Zoltan without ParMETIS."
#endif

#if (PARMETIS_MAJOR_VERSION == 3) && (PARMETIS_MINOR_VERSION < 1)
#error "Specified version of ParMETIS is not compatible with Zoltan; upgrade to ParMETIS v3.1 or later, or build Zoltan without ParMETIS."
#endif

/********** Workaround for memory bug in ParMETIS 3.1.0 **********/
#ifndef PARMETIS_SUBMINOR_VERSION
#define PARMETIS_SUBMINOR_VERSION 0
#endif

#if (PARMETIS_MAJOR_VERSION == 3) && (PARMETIS_MINOR_VERSION == 1) && (PARMETIS_SUBMINOR_VERSION == 0)
#define  PARMETIS31_ALWAYS_FREES_VSIZE
#endif

/**********  parameters structure for parmetis methods **********/
static PARAM_VARS Parmetis_params[] = {
  { "PARMETIS_METHOD", NULL, "STRING", 0 },
  { "PARMETIS_OUTPUT_LEVEL", NULL, "INT", 0 },
  { "PARMETIS_SEED", NULL, "INT", 0 },
  { "PARMETIS_ITR", NULL, "DOUBLE", 0 },
  { "PARMETIS_COARSE_ALG", NULL, "INT", 0 },
  { "PARMETIS_FOLD", NULL, "INT", 0 },
  { NULL, NULL, NULL, 0 } };

static int pmv3method(char *alg);

static int mylog2(int x)
{
  int i = 0;

  for (i=0 ; (1<<i) <= x ; ++i);
  return (i-1);
}

static int Zoltan_Parmetis_Parse(ZZ*, indextype *, char*, realtype*, double *, 
                                 ZOLTAN_Output_Order*);


  /**********************************************************/
  /* Interface routine for ParMetis: Partitioning           */
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
  char *yo = "Zoltan_ParMetis";
  int ierr;
  ZOLTAN_Third_Graph gr;
  ZOLTAN_Third_Geom  *geo = NULL;
  ZOLTAN_Third_Vsize vsp;
  ZOLTAN_Third_Part  prt;
  ZOLTAN_Output_Part part;

  ZOLTAN_ID_PTR global_ids = NULL;
  ZOLTAN_ID_PTR local_ids = NULL;

  int use_timers = 0;
  int timer_p = -1;
  int get_times = 0;
  double times[5];

  double pmv3_itr = 0.0;
  realtype itr = 0.0;
  indextype options[MAX_OPTIONS];
  char alg[MAX_PARAM_STRING_LEN+1];

#ifdef ZOLTAN_PARMETIS
  MPI_Comm comm = zz->Communicator;/* don't risk letting external packages */
                                   /* change our zz struct.                  */
#endif

  indextype i;
  realtype *imb_tols;
  indextype ncon;
  indextype edgecut;
  indextype wgtflag;
  indextype numflag = 0;
  indextype num_part = zz->LB.Num_Global_Parts; /* passed to ParMETIS. */

#ifndef ZOLTAN_PARMETIS
  ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                     "ParMetis requested but not compiled into library.");
  return ZOLTAN_FATAL;

#endif /* ZOLTAN_PARMETIS */

  ZOLTAN_TRACE_ENTER(zz, yo);

  Zoltan_Third_Init(&gr, &prt, &vsp, &part,
                    imp_gids, imp_lids, imp_procs, imp_to_part,
                    exp_gids, exp_lids, exp_procs, exp_to_part);

  if (sizeof(realtype) != sizeof(float)) {
    prt.input_part_sizes = (realtype *)
                   ZOLTAN_MALLOC(zz->LB.Num_Global_Parts * sizeof(realtype));
    for (i = 0; i < zz->LB.Num_Global_Parts; i++) 
      prt.input_part_sizes[i] = (realtype) part_sizes[i];
    prt.part_sizes = prt.input_part_sizes;
  }
  else
    prt.input_part_sizes = prt.part_sizes = (realtype *) part_sizes;


  ierr = Zoltan_Parmetis_Parse(zz, options, alg, &itr, &pmv3_itr, NULL);
  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)) {
    Zoltan_Third_Exit(&gr, geo, &prt, &vsp, &part, NULL);
    return (ierr);
  }

  gr.graph_type = 0;

#ifdef ZOLTAN_PARMETIS
  SET_GLOBAL_GRAPH(&gr.graph_type);
  /* Select type of graph, negative because we impose them */
  /* TODO: add a parameter to select the type, shared with Scotch */
/*   if (strcmp (graph_type, "GLOBAL") != 0) { */
/*     gr.graph_type = - LOCAL_GRAPH; */
/*     if (zz->Num_Proc > 1) { */
/*       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Distributed graph: cannot call METIS, switching to ParMetis"); */
/*       gr.graph_type = - GLOBAL_GRAPH; */
/*       retval = ZOLTAN_WARN; */
/*     } */
/*   } */
#else /* graph is local */
  SET_LOCAL_GRAPH(&gr.graph_type);
#endif /* ZOLTAN_PARMETIS */


  /* Some algorithms use geometry data */
  if (strncmp(alg, "PARTGEOM", 8) == 0){          /* PARTGEOM & PARTGEOMKWAY */
    geo = (ZOLTAN_Third_Geom*) ZOLTAN_MALLOC(sizeof(ZOLTAN_Third_Geom));
    memset (geo, 0, sizeof(ZOLTAN_Third_Geom));
    /* ParMETIS will crash if geometric method and some procs have no nodes. */
    /* Avoid fatal crash by setting scatter to level 2 or higher. */
    gr.scatter_min = 2;
    if (geo == NULL) {
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Out of memory.");
      return (ZOLTAN_MEMERR);
    }
    if (strcmp(alg, "PARTGEOM") == 0) {
      gr.get_data = 0;
    }
  }

  timer_p = Zoltan_Preprocess_Timer(zz, &use_timers);

  /* Start timer */
  get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
  if (get_times){
    MPI_Barrier(zz->Communicator);
    times[0] = Zoltan_Time(zz->Timer);
  }

  vsp.vsize_malloc = 0;
#ifdef PARMETIS31_ALWAYS_FREES_VSIZE
  if (!strcmp(alg, "ADAPTIVEREPART") && (zz->Num_Proc > 1)) {
    /* ParMETIS will free this memory; use malloc to allocate so
       ZOLTAN_MALLOC counters don't show an error. */
    vsp.vsize_malloc = 1 ;
  }
#endif /* PARMETIS31_ALWAYS_FREES_VSIZE */


  ierr = Zoltan_Preprocess_Graph(zz, &global_ids, &local_ids,  &gr, 
                                 geo, &prt, &vsp);
  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)) {
    Zoltan_Third_Exit(&gr, geo, &prt, &vsp, &part, NULL);
    return (ierr);
  }

  /* Get object sizes if requested */
  if (options[PMV3_OPT_USE_OBJ_SIZE] &&
      (zz->Get_Obj_Size || zz->Get_Obj_Size_Multi) &&
      (!strcmp(alg, "ADAPTIVEREPART") || gr.final_output))
    gr.showMoveVol = 1;


  /* Get a time here */
  if (get_times) times[1] = Zoltan_Time(zz->Timer);

  /* Get ready to call ParMETIS */
  edgecut = -1;
  wgtflag = 2*(gr.obj_wgt_dim>0) + (gr.edge_wgt_dim>0);
  numflag = 0;
  ncon = (gr.obj_wgt_dim > 0 ? gr.obj_wgt_dim : 1);

  if (!prt.part_sizes){
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,"Input parameter part_sizes is NULL.");
  }
  if ((zz->Proc == 0) && (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)) {
    for (i=0; i<num_part; i++){
      indextype j;

      printf("Debug: Size(s) for part " TPL_IDX_SPEC " = ", i);
      for (j=0; j<ncon; j++)
        printf("%f ", prt.part_sizes[i*ncon+j]);
      printf("\n");
    }
  }

  /* if (strcmp(alg, "ADAPTIVEREPART") == 0) */
  for (i = 0; i < num_part*ncon; i++)
    if (prt.part_sizes[i] == 0) 
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, "Zero-sized part(s) requested! "
                            "ParMETIS 3.x will likely fail. Please use a "
                            "different method, or remove the zero-sized "
                            "parts from the problem.");


  /* Set Imbalance Tolerance for each weight component. */
  imb_tols = (realtype *) ZOLTAN_MALLOC(ncon * sizeof(realtype));
  if (!imb_tols){
    /* Not enough memory */
    ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
  }
  for (i=0; i<ncon; i++)
    imb_tols[i] = (realtype) (zz->LB.Imbalance_Tol[i]);

  /* Now we can call ParMetis */

  /* Zoltan_Third_Graph_Print(zz, &gr, "Before calling parmetis"); */


#ifdef ZOLTAN_PARMETIS
  if (!IS_LOCAL_GRAPH(gr.graph_type)) { /* May be GLOBAL or NO GRAPH */

    /* First check for ParMetis 3 routines */
    if (strcmp(alg, "PARTKWAY") == 0){
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library "
                                  "ParMETIS_V3_PartKway");
      ParMETIS_V3_PartKway(gr.vtxdist, gr.xadj, gr.adjncy, gr.vwgt, gr.ewgts,
                           &wgtflag, &numflag, &ncon, &num_part, prt.part_sizes,
                           imb_tols, options, &edgecut, prt.part, &comm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
    }
    else if (strcmp(alg, "PARTGEOMKWAY") == 0){
      indextype ndims = geo->ndims;
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library "
                                  "ParMETIS_V3_PartGeomKway");
      ParMETIS_V3_PartGeomKway(gr.vtxdist, gr.xadj, gr.adjncy, gr.vwgt,gr.ewgts,
                               &wgtflag, &numflag, &ndims, geo->xyz, &ncon,
                               &num_part, prt.part_sizes,
                               imb_tols, options, &edgecut, prt.part, &comm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
    }
    else if (strcmp(alg, "PARTGEOM") == 0){
      indextype ndims = geo->ndims;
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library "
                                  "ParMETIS_V3_PartGeom");
      ParMETIS_V3_PartGeom(gr.vtxdist, &ndims, geo->xyz, prt.part, &comm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
    }
    else if (strcmp(alg, "ADAPTIVEREPART") == 0){
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library "
                                  "ParMETIS_V3_AdaptiveRepart");
      ParMETIS_V3_AdaptiveRepart(gr.vtxdist, gr.xadj, gr.adjncy, gr.vwgt,
                                 vsp.vsize, gr.ewgts, &wgtflag, &numflag, &ncon,
                                 &num_part, prt.part_sizes, imb_tols,
                                 &itr, options, &edgecut, prt.part, &comm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
    }
    else if (strcmp(alg, "REFINEKWAY") == 0){
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library "
                                  "ParMETIS_V3_RefineKway");
      ParMETIS_V3_RefineKway(gr.vtxdist, gr.xadj, gr.adjncy, gr.vwgt, gr.ewgts,
                             &wgtflag, &numflag, &ncon, &num_part,
                             prt.part_sizes, imb_tols,
                             options, &edgecut, prt.part, &comm);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");
    }
    else {
      /* Sanity check: This should never happen! */
      char msg[256];
      sprintf(msg, "Unknown ParMetis algorithm %s.", alg);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, msg);
    }
  }
#endif /* ZOLTAN_PARMETIS */
#ifdef ZOLTAN_METIS
  /* TODO: I don't know how to set balance ! */
  if (IS_LOCAL_GRAPH(gr.graph_type)) {
    /* Check for Metis routines */
    if (strcmp(alg, "PARTKWAY") == 0){
      ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the METIS 4 library "
                                  "METIS_WPartGraphKway");
      METIS_WPartGraphKway (gr.vtxdist+1, gr.xadj, gr.adjncy, 
                            gr.vwgt, gr.ewgts, &wgtflag,
                            &numflag, &num_part, prt.part_sizes, 
                            options, &edgecut, prt.part);
      ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the METIS library");
    }
    else {
      /* Sanity check: This should never happen! */
      char msg[256];
      sprintf(msg, "Unknown Metis algorithm %s.", alg);
      ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL, msg);
    }
  }
#endif /* ZOLTAN_METIS */


  /* Get a time here */
  if (get_times) times[2] = Zoltan_Time(zz->Timer);


  if (gr.final_output) { 
    /* Do final output now because after the data will not be coherent:
       unscatter only unscatter part data, not graph */
    ierr = Zoltan_Postprocess_FinalOutput (zz, &gr, &prt, &vsp, use_timers, itr);
  }
  /* Ignore the timings of Final Ouput */
  if (get_times) times[3] = Zoltan_Time(zz->Timer);

  ierr = Zoltan_Postprocess_Graph(zz, global_ids, local_ids, &gr, 
                                  geo, &prt, &vsp, NULL, &part);

  Zoltan_Third_Export_User(&part, 
                           num_imp, imp_gids, imp_lids, imp_procs, imp_to_part,
                           num_exp, exp_gids, exp_lids, exp_procs, exp_to_part);

  /* Get a time here */
  if (get_times) times[4] = Zoltan_Time(zz->Timer);

  if (get_times) Zoltan_Third_DisplayTime(zz, times);

  if (use_timers && timer_p >= 0)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_p, zz->Communicator);

  Zoltan_Third_Exit(&gr, geo, &prt, &vsp, NULL, NULL);
  if (imb_tols != NULL) ZOLTAN_FREE(&imb_tols);
  if (geo != NULL) ZOLTAN_FREE(&geo);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}


static int Zoltan_Parmetis_Parse(
  ZZ* zz, 
  indextype *options, 
  char* alg,
  realtype* itr, 
  double *pmv3_itr,
  ZOLTAN_Output_Order *ord
)
{
  static char * yo = "Zoltan_Parmetis_Parse";

  int  i;
  int output_level, seed, coarse_alg, fold, use_obj_size;

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

    if(ord == NULL) {
    /* Map LB_APPROACH to suitable PARMETIS_METHOD */
      if (!strcasecmp(zz->LB.Approach, "partition")){
        strcpy(alg, "PARTKWAY");
      }
      else if (!strcasecmp(zz->LB.Approach, "repartition")){
        strcpy(alg, "ADAPTIVEREPART");
        *pmv3_itr = 100.; /* Ratio of inter-proc comm. time to data redist. time;
                          100 gives similar partition quality to GDiffusion */
      }
      else if (!strcasecmp(zz->LB.Approach, "refine")){
        strcpy(alg, "REFINEKWAY");
      }
      else { /* If no LB_APPROACH is set, use repartition */
        strcpy(alg, "ADAPTIVEREPART");
        *pmv3_itr = 100.; /* Ratio of inter-proc comm. time to data redist. time;
                          100 gives similar partition quality to GDiffusion */
      }
    }
    else {
      strcpy(alg, "NODEND");
    }

    Zoltan_Bind_Param(Parmetis_params, "PARMETIS_METHOD",
                      (void *) alg);
    Zoltan_Bind_Param(Parmetis_params, "PARMETIS_OUTPUT_LEVEL",
                      (void *) &output_level);
    Zoltan_Bind_Param(Parmetis_params, "PARMETIS_SEED",
                      (void *) &seed);
    Zoltan_Bind_Param(Parmetis_params, "PARMETIS_ITR",
                      (void *) pmv3_itr);
    Zoltan_Bind_Param(Parmetis_params, "PARMETIS_COARSE_ALG",
                      (void *) &coarse_alg);
    Zoltan_Bind_Param(Parmetis_params, "PARMETIS_FOLD",
                      (void *) &fold);

    Zoltan_Assign_Param_Vals(zz->Params, Parmetis_params, zz->Debug_Level,
                             zz->Proc, zz->Debug_Proc);

    /* Copy option values to ParMetis options array */

    /* In this version of Zoltan, processors and partitions are coupled. */
    /* This will likely change in future releases, and then the options  */
    /* value should change to DISCOUPLED.                                */
    options[PMV3_OPTION_PSR] = COUPLED;

    if (pmv3method(alg)){
      /* ParMetis 3.0 options */
      options[PMV3_OPTION_DBGLVL] = output_level;
      options[PMV3_OPTION_SEED] = seed;
      options[PMV3_OPT_USE_OBJ_SIZE] = use_obj_size;
      if (ord == NULL)
        *itr = (realtype)*pmv3_itr;
    }

    /* If ordering, use ordering method instead of load-balancing method */
    if (ord && ord->order_opt && ord->order_opt->method){
      strcpy(alg, ord->order_opt->method);
    }

    if ((zz->Num_Proc == 1) &&
        (!strcmp(alg, "ADAPTIVEREPART") ||
         !strcmp(alg, "REPARTLDIFFUSION") ||
         !strcmp(alg, "REPARTGDIFFUSION") ||
         !strcmp(alg, "REPARTREMAP") ||
         !strcmp(alg, "REPARTMLREMAP"))) {
      /* These ParMETIS methods fail on one processor; an MPI command assumes
         at least two processors. */
      char str[256];
      sprintf(str, "ParMETIS method %s fails on one processor due to a bug"
                   " in ParMETIS v3.x; please select another method.", alg);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, str);
      return (ZOLTAN_FATAL);
    }

    return(ZOLTAN_OK);
}

/***************************************************************************/

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


/***************************************************************************
 *  The ParMetis ordering routine piggy-backs on the ParMetis
 *  partitioning routines.
 **************************************************************************/

int Zoltan_ParMetis_Order(
  ZZ *zz,               /* Zoltan structure */
  int num_obj,          /* Number of (local) objects to order. */
  ZOLTAN_ID_PTR gids,   /* List of global ids (local to this proc) */
                        /* The application must allocate enough space */
  ZOLTAN_ID_PTR lids,   /* List of local ids (local to this proc) */
                        /* The application must allocate enough space */
  ZOLTAN_ID_PTR rank,   /* rank[i] is the rank of gids[i] */
  int *iperm,
  ZOOS *order_opt       /* Ordering options, parsed by Zoltan_Order */
)
{
  static char *yo = "Zoltan_ParMetis_Order";
  int i, n, ierr;
  ZOLTAN_Output_Order ord;
  ZOLTAN_Third_Graph gr;

#ifdef ZOLTAN_PARMETIS
  MPI_Comm comm = zz->Communicator;/* don't want to risk letting external 
                                      packages changing our communicator */
  indextype numflag = 0;
#endif

  int timer_p = 0;
  int get_times = 0;
  int use_timers = 0;
  double times[5];

  ZOLTAN_ID_PTR       l_gids = NULL;
  ZOLTAN_ID_PTR       l_lids = NULL;

  indextype options[MAX_OPTIONS];
  char alg[MAX_PARAM_STRING_LEN+1];

  ZOLTAN_TRACE_ENTER(zz, yo);

#ifdef ZOLTAN_PARMETIS
#if TPL_USE_DATATYPE != TPL_METIS_DATATYPES

#ifdef TPL_FLOAT_WEIGHT
  i = 1;
#else
  i = 0;
#endif

  if ((sizeof(indextype) != sizeof(idxtype)) ||
      (sizeof(weighttype) != sizeof(idxtype)) || i){

    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,
          "Not supported: Multiple 3rd party libraries with incompatible "
          "data types.");
    return ZOLTAN_FATAL;
  }
#endif
#endif

  memset(&gr, 0, sizeof(ZOLTAN_Third_Graph));
  memset(&ord, 0, sizeof(ZOLTAN_Output_Order));
  memset(times, 0, sizeof(times));

  ord.order_opt = order_opt;

  if (!order_opt){
    /* If for some reason order_opt is NULL, allocate a new ZOOS here. */
    /* This should really never happen. */
    order_opt = (ZOOS *) ZOLTAN_MALLOC(sizeof(ZOOS));
    strcpy(order_opt->method,"PARMETIS");
  }

  ierr = Zoltan_Parmetis_Parse(zz, options, alg, NULL, NULL, &ord);
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
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input num_obj does not equal the "
                                     "number of objects.");
    return(ZOLTAN_FATAL);
  }

  /* Do not use weights for ordering */
  gr.obj_wgt_dim = -1;
  gr.edge_wgt_dim = -1;
  gr.num_obj = num_obj;

  /* Check what ordering type is requested */
  if (order_opt){
      SET_GLOBAL_GRAPH(&gr.graph_type); /* GLOBAL by default */

#ifdef ZOLTAN_PARMETIS
      if ((strcmp(order_opt->method, "METIS") == 0))
#endif /* ZOLTAN_PARMETIS */
      SET_LOCAL_GRAPH(&gr.graph_type);
  }
  gr.get_data = 1;

  if (IS_LOCAL_GRAPH(gr.graph_type) && zz->Num_Proc > 1) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Serial ordering on more than 1 process: "
                                     "set ParMetis instead.");
    return(ZOLTAN_FATAL);
  }

  timer_p = Zoltan_Preprocess_Timer(zz, &use_timers);

    /* Start timer */
  get_times = (zz->Debug_Level >= ZOLTAN_DEBUG_ATIME);
  if (get_times){
    MPI_Barrier(zz->Communicator);
    times[0] = Zoltan_Time(zz->Timer);
  }

  ierr = Zoltan_Preprocess_Graph(zz, &l_gids, &l_lids,  &gr, NULL, NULL, NULL);
  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)) {
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, NULL);
    return (ierr);
  }

  /* Allocate space for separator sizes */

  if (IS_GLOBAL_GRAPH(gr.graph_type)) {
    if (Zoltan_TPL_Order_Init_Tree(&zz->TPL_Order, 2*zz->Num_Proc, zz->Num_Proc) != ZOLTAN_OK) {
      /* Not enough memory */
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    ord.sep_sizes = (indextype*)ZOLTAN_MALLOC((2*zz->Num_Proc+1)*sizeof(indextype));
    if (ord.sep_sizes == NULL) {
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    memset(ord.sep_sizes, 0, (2*zz->Num_Proc+1)*sizeof(int)); /* It seems parmetis don't initialize correctly */
  }

  /* Allocate space for direct perm */
  ord.rank = (indextype *) ZOLTAN_MALLOC(gr.num_obj*sizeof(indextype));
  if (!ord.rank){
    /* Not enough memory */
    Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
    ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
  }
  if (IS_LOCAL_GRAPH(gr.graph_type)){
  /* Allocate space for inverse perm */
    ord.iperm = (indextype *) ZOLTAN_MALLOC(gr.num_obj*sizeof(indextype));
    if (!ord.iperm){
      /* Not enough memory */
      Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, &ord);
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
  }
  else
    ord.iperm = NULL;

  /* Get a time here */
  if (get_times) times[1] = Zoltan_Time(zz->Timer);

#ifdef ZOLTAN_PARMETIS
  if (IS_GLOBAL_GRAPH(gr.graph_type)){
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the ParMETIS 3 library");

    ParMETIS_V3_NodeND (gr.vtxdist, gr.xadj, gr.adjncy, 
                        &numflag, options, ord.rank, ord.sep_sizes, &comm);
    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the ParMETIS library");

  }
  else
#endif /* ZOLTAN_PARMETIS */
#if defined(ZOLTAN_METIS) || defined(ZOLTAN_PARMETIS)
 if (IS_LOCAL_GRAPH(gr.graph_type)) { /* Be careful : permutation parameters are in the opposite order */
    indextype numobj = gr.num_obj;
    ZOLTAN_TRACE_DETAIL(zz, yo, "Calling the METIS library");
    options[0] = 0;  /* Use default options for METIS. */
    order_opt->return_args = RETURN_RANK|RETURN_IPERM; /* We provide directly all the permutations */

    METIS_NodeND(&numobj, gr.xadj, gr.adjncy, &numflag, options, 
                 ord.iperm, ord.rank);

    ZOLTAN_TRACE_DETAIL(zz, yo, "Returned from the METIS library");
  }
#endif /* ZOLTAN_METIS */

  /* Get a time here */
  if (get_times) times[2] = Zoltan_Time(zz->Timer);

  if (IS_GLOBAL_GRAPH(gr.graph_type)){ /* Update Elimination tree */
    int numbloc;
    int start;
    int leaf;
    int *converttab;
    int levelmax;

    levelmax = mylog2(zz->Num_Proc) + 1;
    converttab = (int*)ZOLTAN_MALLOC(zz->Num_Proc*2*sizeof(int));

    memset(converttab, 0, zz->Num_Proc*2*sizeof(int));
     /* Determine the first node in each separator, store it in zz->TPL_Order.start */
    for (numbloc = 0, start=0, leaf=0; numbloc < zz->Num_Proc /2; numbloc++) {
      int father;

      father = zz->Num_Proc + numbloc;
      converttab[start] = 2*numbloc;
      zz->TPL_Order.leaves[leaf++]=start;
      zz->TPL_Order.ancestor[start] = start + 2;
      converttab[start+1] = 2*numbloc+1;
      zz->TPL_Order.leaves[leaf++]=start+1;
      zz->TPL_Order.ancestor[start+1] = start + 2;
      start+=2;
      do {
        converttab[start] = father;
        if (father %2 == 0) {
          int nextoffset;
          int level;

          level = mylog2(2*zz->Num_Proc - 1 - father);
          nextoffset = (1<<(levelmax-level));
          zz->TPL_Order.ancestor[start] = start+nextoffset;
          start++;
          break;
        }
        else {
          zz->TPL_Order.ancestor[start] = start+1;
          start++;
          father = zz->Num_Proc + father/2;
        }
      } while (father < 2*zz->Num_Proc - 1);
    }

    zz->TPL_Order.start[0] = 0;
    zz->TPL_Order.ancestor [2*zz->Num_Proc - 2] = -1;
    for (numbloc = 1 ; numbloc < 2*zz->Num_Proc ; numbloc++) {
      int oldblock=converttab[numbloc-1];
      zz->TPL_Order.start[numbloc] = zz->TPL_Order.start[numbloc-1] + ord.sep_sizes[oldblock];
    }

    ZOLTAN_FREE(&converttab);
    ZOLTAN_FREE(&ord.sep_sizes);

    zz->TPL_Order.leaves[zz->Num_Proc] = -1;
    zz->TPL_Order.nbr_leaves = zz->Num_Proc;
    zz->TPL_Order.nbr_blocks = 2*zz->Num_Proc-1;
  }
  else { /* No tree */
    zz->TPL_Order.nbr_blocks = 0;
    zz->TPL_Order.start = NULL;
    zz->TPL_Order.ancestor = NULL;
    zz->TPL_Order.leaves = NULL;
  }

  /* Correct because no redistribution */
  memcpy(gids, l_gids, n*zz->Num_GID*sizeof(ZOLTAN_ID_TYPE));
  memcpy(lids, l_lids, n*zz->Num_LID*sizeof(ZOLTAN_ID_TYPE));

  ierr = Zoltan_Postprocess_Graph (zz, l_gids, l_lids, &gr, NULL, NULL, NULL, &ord, NULL);

  ZOLTAN_FREE(&l_gids); 
  ZOLTAN_FREE(&l_lids);

  /* Get a time here */
  if (get_times) times[3] = Zoltan_Time(zz->Timer);

  if (get_times) Zoltan_Third_DisplayTime(zz, times);

  if (use_timers)
    ZOLTAN_TIMER_STOP(zz->ZTime, timer_p, zz->Communicator);

  if (sizeof(indextype) == sizeof(ZOLTAN_ID_TYPE)){
    memcpy(rank, ord.rank, gr.num_obj*sizeof(indextype));
  }
  else{
    for (i=0; i < gr.num_obj; i++){
      rank[i] = (ZOLTAN_ID_TYPE)ord.rank[i];
    }
  }

  if ((ord.iperm != NULL) && (iperm != NULL)){
    if (sizeof(indextype) == sizeof(int)){
      memcpy(iperm, ord.iperm, gr.num_obj*sizeof(indextype));
    }
    else{
      for (i=0; i < gr.num_obj; i++){
        iperm[i] = (int)ord.iperm[i];
      }
    }
  }

  if (ord.iperm != NULL)  ZOLTAN_FREE(&ord.iperm);
  ZOLTAN_FREE(&ord.rank);

  /* Free all other "graph" stuff */
  Zoltan_Third_Exit(&gr, NULL, NULL, NULL, NULL, NULL);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ZOLTAN_OK);
}

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
  }
  return(status);
}

#ifdef __cplusplus
}
#endif
