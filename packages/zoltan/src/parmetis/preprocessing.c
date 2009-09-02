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
#ifdef TPL_NEW_GRAPH
#include "graph.h"
#endif

/**********  parameters structure used by both ParMetis and Jostle **********/
static PARAM_VARS Graph_params[] = {
	{ "CHECK_GRAPH", NULL, "INT", 0 },
	{ "SCATTER_GRAPH", NULL, "INT", 0 },
	{ "FINAL_OUTPUT", NULL, "INT", 0 },
	{ "USE_TIMERS", NULL, "INT", 0 },
	{ "ADD_OBJ_WEIGHT",  NULL,  "STRING", 0},
	{ NULL, NULL, NULL, 0 } };

/* Extern function prototypes. Should be in a separate header file? */
extern int Zoltan_Verify_Graph(MPI_Comm comm, indextype *vtxdist, indextype *xadj,
       indextype *adjncy, indextype *vwgt, indextype *adjwgt,
       int vwgt_dim, int ewgt_dim, int graph_type, int check_graph,
       int debug_level);
extern int Zoltan_Scatter_Graph(indextype **vtxdist, indextype **xadj,
       indextype **adjncy, indextype **vwgt, indextype **vsize, indextype **adjwgt,
       float **xyz, int ndims, int, ZZ *zz, ZOLTAN_COMM_OBJ **plan);
extern int Zoltan_Compare_Ints(const void *key, const void *arg);


/* Auxiliary function prototypes. */
static int Zoltan_Preprocess_Add_Weight (ZZ*, ZOLTAN_Third_Graph * gr,
				  ZOLTAN_Third_Part  * prt,
				  char * add_obj_weight);


static int
Zoltan_Preprocess_Scale_Weights (ZOLTAN_Third_Graph *gr, float *flt_wgt, weighttype** rnd_wgt,
				 int number, int ndim, int mode, ZZ* zz,
				 char * name, int offset);
static int
Zoltan_Preprocess_Extract_Geom (ZZ *zz,
				ZOLTAN_ID_PTR *global_ids,
				ZOLTAN_ID_PTR *local_ids,
				ZOLTAN_Third_Graph *gr,
				ZOLTAN_Third_Geom *geo);
static int
Zoltan_Preprocess_Extract_Vsize (ZZ *zz,
				 ZOLTAN_ID_PTR *global_ids,
				 ZOLTAN_ID_PTR *local_ids,
				 ZOLTAN_Third_Graph *gr,
				 ZOLTAN_Third_Vsize *vsp);
static int
Zoltan_Preprocess_Scatter_Graph (ZZ *zz,
				 ZOLTAN_Third_Graph *gr,
				 ZOLTAN_Third_Part *prt,
				 ZOLTAN_Third_Geom *geo,
				 ZOLTAN_Third_Vsize *vsp);

static int scale_round_weights(float *fwgts, indextype *iwgts, int n, int dim,
			       int mode, int max_wgt_sum, int debug_level, MPI_Comm comm);



/**
 *  This function fills data structures in order to call a third party
 *  library like ParMetis or Scotch
 **/

int Zoltan_Preprocess_Graph(
  ZZ *zz,                               /* Zoltan structure */
  ZOLTAN_ID_PTR *global_ids,
  ZOLTAN_ID_PTR *local_ids,
  ZOLTAN_Third_Graph *gr,              /* Graph for third part libs */
  ZOLTAN_Third_Geom  *geo,
  ZOLTAN_Third_Part  *prt,
  ZOLTAN_Third_Vsize *vsp
)
{
  static char *yo = "Zoltan_Preprocess_Graph";

  int ierr;
  float *float_vwgt, *float_ewgts;
  char msg[256];
#ifdef TPL_NEW_GRAPH
  ZG *graph = &(gr->graph);
  int local;
#endif /* TPL_NEW_GRAPH */

  char add_obj_weight[MAX_PARAM_STRING_LEN+1];

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  gr->vtxdist = gr->xadj = gr->adjncy = NULL;
  gr->vwgt = gr->ewgts = NULL;
  float_vwgt = float_ewgts = NULL;

  if (gr->obj_wgt_dim >= 0) {
    /* Check weight dimensions */
    if (zz->Obj_Weight_Dim<0){
      sprintf(msg, "Object weight dimension is %d, "
	      "but should be >= 0. Using Obj_Weight_Dim = 0.",
	      zz->Obj_Weight_Dim);
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
      gr->obj_wgt_dim = 0;
    }
    else {
      gr->obj_wgt_dim = zz->Obj_Weight_Dim;
    }
  }
  else
    gr->obj_wgt_dim = 0;
  if (gr->edge_wgt_dim >= 0) {
    if (zz->Edge_Weight_Dim<0){
      sprintf(msg, "Edge weight dimension is %d, "
	      "but should be >= 0. Using Edge_Weight_Dim = 0.",
	      zz->Edge_Weight_Dim);
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
      gr->edge_wgt_dim = 0;
    }
    else if (zz->Edge_Weight_Dim>1){
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "This method does not support "
			"multidimensional edge weights. Using Edge_Weight_Dim = 1.");
      gr->edge_wgt_dim = 1;
    }
    else {
      gr->edge_wgt_dim = zz->Edge_Weight_Dim;
    }
  }
  else
      gr->edge_wgt_dim = 0;

  /* Default graph type is GLOBAL. */

  /* Get parameter options shared by ParMetis and Jostle */
  gr->check_graph = 1;          /* default */
  gr->scatter = 1;              /* default */
  gr->final_output = 0;
  strcpy(add_obj_weight, "NONE");  /* default */
  Zoltan_Bind_Param(Graph_params, "CHECK_GRAPH", (void *) &gr->check_graph);
  Zoltan_Bind_Param(Graph_params, "SCATTER_GRAPH", (void *) &gr->scatter);
  Zoltan_Bind_Param(Graph_params, "FINAL_OUTPUT", (void *) &gr->final_output);
  Zoltan_Bind_Param(Graph_params, "ADD_OBJ_WEIGHT", (void *) add_obj_weight);
  Zoltan_Assign_Param_Vals(zz->Params, Graph_params, zz->Debug_Level, zz->Proc,
			   zz->Debug_Proc);


  /* If reorder is true, we already have the id lists. Ignore weights. */
  if ((*global_ids == NULL) || (!gr->id_known)){
    int *input_part = NULL;
    ierr = Zoltan_Get_Obj_List(zz, &gr->num_obj, global_ids, local_ids,
			       gr->obj_wgt_dim, &float_vwgt, &input_part);
    if (prt) {
      prt->input_part = input_part;
    }
    else if (input_part) { /* Ordering, dont need part */
      ZOLTAN_FREE(&input_part);
    }
    if (ierr){
      /* Return error */
      ZOLTAN_PARMETIS_ERROR(ierr, "Get_Obj_List returned error.");
    }
  }
  /* Build Graph for third party library data structures, or just get vtxdist. */

  if (gr->get_data) {
#ifndef TPL_NEW_GRAPH
    ierr = Zoltan_Build_Graph(zz, &gr->graph_type, gr->check_graph, gr->num_obj,
			      *global_ids, *local_ids, gr->obj_wgt_dim, &gr->edge_wgt_dim,
			      &gr->vtxdist, &gr->xadj, &gr->adjncy, &float_ewgts, &gr->adjproc);
#else
    local = IS_LOCAL_GRAPH(gr->graph_type);
    ierr = Zoltan_ZG_Build (zz, graph, local); /* Normal graph */
    ierr = Zoltan_ZG_Export (zz, graph,
			     &gr->num_obj, &gr->num_obj, &gr->obj_wgt_dim, &gr->edge_wgt_dim,
			     &gr->vtxdist, &gr->xadj, &gr->adjncy, &gr->adjproc,
			     &float_vwgt, &float_ewgts, NULL);
    /* TODO: support graph redistribution */
/*   if (prt) */
/*     ierr = Zoltan_ZG_Vertex_Info(zz, &graph, global_ids, &prt->input_part); */
/*   else */
/*     ierr = Zoltan_ZG_Vertex_Info(zz, &graph, global_ids, NULL); */

/*     /\* Just to try *\/ */
/*     if (prt) { */
/*       prt->input_part = (int *)ZOLTAN_CALLOC(gr->num_obj, sizeof(int)); */
/*     } */

#endif

    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      ZOLTAN_PARMETIS_ERROR(ierr, "Zoltan_Build_Graph returned error.");
    }
  }
  else{ /* Only geometry */
    int i;
    /* No graph but still needs vtxdist*/
    gr->vtxdist = (int*) ZOLTAN_MALLOC ((zz->Num_Proc+1)*sizeof(int));
    if (gr->vtxdist == NULL)
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");

    gr->vtxdist[0] = 0;
    MPI_Allgather(&gr->num_obj, 1, MPI_INT, gr->vtxdist+1, 1, MPI_INT, zz->Communicator);
    for (i=1 ; i <= zz->Num_Proc ; ++i) {
      gr->vtxdist[i] += gr->vtxdist[i-1];
    }
  }

  if (prt) {
    prt->part_sizes = prt->input_part_sizes;
    /* ParMETIS needs prt->part to be allocated, even when num_obj=0. */
    /* See Zoltan bug 4299. */
    prt->part = (indextype *)ZOLTAN_MALLOC((gr->num_obj+1) * sizeof(indextype));
    if (!prt->part) {
      /* Not enough memory */
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_MEMERR, "Out of memory.");
    }
    if (gr->num_obj >0) {
      memcpy (prt->part, prt->input_part, (gr->num_obj) * sizeof(indextype));
    }
    else {
      prt->input_part = NULL;
    }
  }

  /* Convert from float. */

  /* Get vertex weights if needed */
  if (gr->obj_wgt_dim){
    ierr = Zoltan_Preprocess_Scale_Weights (gr, float_vwgt, &gr->vwgt,
					    gr->num_obj, gr->obj_wgt_dim, 1, zz,
					    "vertex", gr->vtxdist[zz->Proc]);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      /* Return error code */
      ZOLTAN_PARMETIS_ERROR(ierr, "Error in scaling of weights.");
    }
    ZOLTAN_FREE(&float_vwgt);
  }

  if (strcasecmp(add_obj_weight, "NONE")){
    if (Zoltan_Preprocess_Add_Weight(zz, gr, prt, add_obj_weight) != ZOLTAN_OK)
      ZOLTAN_PARMETIS_ERROR(ierr, "Error in adding  vertex weights.");
  }

  /* Get edge weights if needed */
  if (gr->get_data)
    gr->num_edges = gr->xadj[gr->num_obj];
  else {
    gr->num_edges = 0;
    gr->edge_wgt_dim = 0;
  }

  if (gr->edge_wgt_dim){
    ierr = Zoltan_Preprocess_Scale_Weights (gr, float_ewgts, &gr->ewgts,
					    gr->num_edges, gr->edge_wgt_dim, 1, zz,
					    "edge", 0);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
      /* Return error code */
      ZOLTAN_PARMETIS_ERROR(ierr, "Error in scaling of weights.");
    }
    if (!gr->final_output) {
      ZOLTAN_FREE(&float_ewgts);
    }
    else
      gr->float_ewgts = float_ewgts;
  }
  else {
    ZOLTAN_FREE(&float_ewgts);
  }

  if (geo){
    ierr = Zoltan_Preprocess_Extract_Geom (zz, global_ids, local_ids, gr, geo);
    if (ierr) {
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL,
			    "Error returned from Zoltan_Preprocess_Extract_Geom");
    }
  }


  if (vsp) {
    ierr = Zoltan_Preprocess_Extract_Vsize (zz, global_ids, local_ids, gr, vsp);
    if (ierr) {
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL,
			  "Error returned from Zoltan_Preprocess_Extract_Vsize");
    }
  }

  /* Scatter graph?
   * If data distribution is highly imbalanced, it is better to
   * redistribute the graph data structure before calling ParMetis.
   * After partitioning, the results must be mapped back.
   */
  if (gr->scatter < gr->scatter_min) gr->scatter = gr->scatter_min;

  if (gr->scatter>0) {
    ierr = Zoltan_Preprocess_Scatter_Graph (zz, gr, prt, geo, vsp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PARMETIS_ERROR(ZOLTAN_FATAL,
			    "Error returned from Zoltan_Preprocess_Scatter_Graph");
    }
  }


  /* Verify that graph is correct */
  if (gr->get_data){
    int flag;

    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      flag = 2; /* Maximum output level */
    else
      flag = 1; /* Medium output level */
    ierr = Zoltan_Verify_Graph(zz->Communicator, gr->vtxdist, gr->xadj, gr->adjncy, gr->vwgt,
	      gr->ewgts, gr->obj_wgt_dim, gr->edge_wgt_dim, gr->graph_type, gr->check_graph, flag);

  }

 End:

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

static int
Zoltan_Preprocess_Add_Weight (ZZ *zz,
			      ZOLTAN_Third_Graph * gr,
			      ZOLTAN_Third_Part  * prt,
			      char * add_obj_weight)
{
  static char * yo = "Zoltan_Preprocess_Add_Weight";
  /* Add a vertex weight? */
  int add_type = 0;
  weighttype  *vwgt_new;
  int ierr = ZOLTAN_OK;
  int i,j;

  vwgt_new = (weighttype *)ZOLTAN_MALLOC((gr->obj_wgt_dim+1)*gr->num_obj
					    * sizeof(weighttype));
  if ((!strcasecmp(add_obj_weight, "UNIT")) || (!strcasecmp(add_obj_weight, "VERTICES"))){
    add_type = 1;
  }
  else if ((!strcasecmp(add_obj_weight, "VERTEX DEGREE")) ||
	   (!strcasecmp(add_obj_weight, "PINS")) ||
	   (!strcasecmp(add_obj_weight, "NONZEROS"))){
    add_type = 2;
  }
  else {
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Invalid parameter value for ADD_OBJ_WEIGHT!\n");
    ierr = ZOLTAN_WARN;
    add_type = 0;
  }
  if (add_type){
    if (prt != NULL) {
      /* update part_sizes array */
      ierr = Zoltan_LB_Add_Part_Sizes_Weight(zz,
					     (gr->obj_wgt_dim ? gr->obj_wgt_dim : 1),
					     gr->obj_wgt_dim+1,
					     prt->input_part_sizes, &prt->part_sizes);
    }
	/* Add implicit weights in new array */
    for (i=0; i<gr->num_obj; i++){
      /* First copy old weights */
      for (j=0; j<gr->obj_wgt_dim; j++){
	vwgt_new[i*(gr->obj_wgt_dim+1)+j] = gr->vwgt[i*gr->obj_wgt_dim+j];
      }
      if (add_type==1)
	/* unit weight */
	vwgt_new[i*(gr->obj_wgt_dim+1)+gr->obj_wgt_dim] = 1;
      else if (add_type==2)
	/* weight is vertex degree (EBEB should we add +1?) */
	/*	vwgt_new[i*(gr->obj_wgt_dim+1)+gr->obj_wgt_dim] = 10*(gr->xadj[i+1] -gr->xadj[i]) + 1; */
	vwgt_new[i*(gr->obj_wgt_dim+1)+gr->obj_wgt_dim] = gr->xadj[i+1] -gr->xadj[i];
    }
	/* Use new vwgt array */
    if (gr->vwgt) ZOLTAN_FREE(&gr->vwgt);
    gr->vwgt = vwgt_new;
    gr->obj_wgt_dim += 1;
  }

  return (ierr);
}


static int
Zoltan_Preprocess_Scale_Weights (ZOLTAN_Third_Graph *gr, float *flt_wgt, weighttype** rnd_wgt,
				 int number, int ndim, int mode, ZZ* zz,
				 char * name, int offset)
{
  static char * yo = "Zoltan_Preprocess_Scale_Weights";
  int ierr = ZOLTAN_OK;
  int i99;
  int k;
  char msg[256];

  if (ndim == 0)
    return ierr;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
    for (i99=0; i99 < (number <3 ? number : 3); i99++){
      for (k=0; k<gr->obj_wgt_dim; k++)
	sprintf(msg+10*k, "%.9f ", flt_wgt[i99*gr->obj_wgt_dim+k]);
      printf("[%1d] Debug: before scaling weights for %s %d = %s\n",
	     zz->Proc, name, offset+i99, msg);
    }
  }

  *rnd_wgt = (weighttype *)ZOLTAN_MALLOC(ndim * number * sizeof(weighttype));
  if ((number >0 ) && (*rnd_wgt == NULL)){
	/* Not enough memory */
    ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
  }

  ierr = scale_round_weights(flt_wgt, *rnd_wgt, number, ndim, mode,
			     (int) MAX_WGT_SUM, zz->Debug_Level, zz->Communicator);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
    ZOLTAN_FREE(rnd_wgt);
    /* Return error code */
    ZOLTAN_THIRD_ERROR(ierr, "Error in scaling of weights.");
  }

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
    for (i99=0; i99 < (number < 3 ? number : 3); i99++){
      for (k=0; k<gr->obj_wgt_dim; k++)
	sprintf(msg+10*k, "%9d ", (*rnd_wgt)[i99*gr->obj_wgt_dim+k]);
      printf("[%1d] Debug: scaled weights for %s %d = %s\n",
	     zz->Proc, name, offset+i99, msg);
    }
  }

  return (ierr);
}



static int
Zoltan_Preprocess_Extract_Geom (ZZ *zz,
				ZOLTAN_ID_PTR *global_ids,
				ZOLTAN_ID_PTR *local_ids,
				ZOLTAN_Third_Graph *gr,
				ZOLTAN_Third_Geom *geo)
{
  static char * yo = "Zoltan_Preprocess_Extract_Geom";
  int ierr;
  double *geom_vec;
  int i;

  geom_vec = NULL;
  /* Get coordinate information */
  ierr = Zoltan_Get_Coordinates(zz, gr->num_obj, *global_ids, *local_ids,
				&geo->ndims, &geom_vec);
  if (ierr) {
    ZOLTAN_THIRD_ERROR(ZOLTAN_FATAL,
			  "Error returned from Zoltan_Get_Coordinates");
  }
  /* Convert geometry info from double to float for ParMETIS */
  if (gr->num_obj && geo->ndims) {
    geo->xyz = (float *) ZOLTAN_MALLOC(gr->num_obj * geo->ndims * sizeof(float));
    if (geo->xyz == NULL)  {
      ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Memory error.");
    }
    for (i = 0; i < gr->num_obj * geo->ndims; i++)
      geo->xyz[i] = (float) geom_vec[i];
    ZOLTAN_FREE(&geom_vec);
  }

  return ierr;
}


static int
Zoltan_Preprocess_Extract_Vsize (ZZ *zz,
				 ZOLTAN_ID_PTR *global_ids,
				 ZOLTAN_ID_PTR *local_ids,
				 ZOLTAN_Third_Graph *gr,
				 ZOLTAN_Third_Vsize *vsp)
{
  static char * yo = "Zoltan_Preprocess_Extract_Vsize";
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int ierr= ZOLTAN_OK;
  int i;

  if (gr->obj_wgt_dim)
    vsp->vsize_malloc = 0;

  if (!vsp->vsize) {                    /* If not already allocated */
    if (vsp->vsize_malloc)
      vsp->vsize = (indextype *) malloc(gr->num_obj*sizeof(indextype));
    else
      vsp->vsize = (indextype *) ZOLTAN_MALLOC(gr->num_obj*sizeof(indextype));
  }
  vsp->vsizeBACKUP = (indextype *) ZOLTAN_MALLOC(gr->num_obj*sizeof(indextype));
  if (gr->num_obj && !vsp->vsize){
      /* Not enough space */
    ZOLTAN_THIRD_ERROR(ZOLTAN_MEMERR, "Out of memory.");
  }
  if (zz->Get_Obj_Size_Multi) {
    zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data,
			   num_gid_entries, num_lid_entries, gr->num_obj,
			   *global_ids, *local_ids, vsp->vsize, &ierr);
  }
  else if (zz->Get_Obj_Size) {
    ZOLTAN_ID_PTR lid;
    for (i=0; i<gr->num_obj; i++){
      lid = (num_lid_entries ? &((*local_ids)[i*num_lid_entries]) : NULL);
      vsp->vsize[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data,
				  num_gid_entries, num_lid_entries,
				  &((*global_ids)[i*num_gid_entries]),
				  lid, &ierr);
    }
  }
  else {
    /* Assume uniform sizes if no Obj_Size callbacks are registered. */
    for (i = 0; i < gr->num_obj; i++)
    vsp->vsize[i] = 1;
  }
  memcpy(vsp->vsizeBACKUP, vsp->vsize, sizeof(indextype)*gr->num_obj);
  return (ierr);
}


static int
Zoltan_Preprocess_Scatter_Graph (ZZ *zz,
				 ZOLTAN_Third_Graph *gr,
				 ZOLTAN_Third_Part *prt,
				 ZOLTAN_Third_Geom *geo,
				 ZOLTAN_Third_Vsize *vsp)
{
  static char * yo = "Zoltan_Preprocess_Scatter_Graph";
  int ierr = ZOLTAN_OK;
  int tmp;
  int procid;
  int i;

  if ((gr->scatter>0) && (gr->scatter<3)){
    int j;
    /* Decide if the data imbalance is so bad that we should scatter the graph. */
    /* scatter==1: Scatter if all the objects are on a single processor.        */
    /* scatter==2: Scatter if any processor has no objects.                     */
    gr->num_obj_orig = gr->num_obj; /* Save value for later. */
    if (gr->num_obj==0)
      j = 1;
    else
      j = 0;
    MPI_Allreduce(&j, &tmp, 1, MPI_INT, MPI_SUM, zz->Communicator);
    if (gr->scatter == 1){
      if (tmp < zz->Num_Proc-1)
	gr->scatter = 0;
    }
    else if (gr->scatter==2){
      if (tmp == 0)
	gr->scatter = 0;
    }
  }

  /* We need to make sure we don't scatter the graph
   * if graph_type = LOCAL_GRAPH, i.e. METIS is used.
   */
  if (gr->scatter && (IS_LOCAL_GRAPH(gr->graph_type))){
    gr->scatter = 0;
    ZOLTAN_PRINT_WARN(zz->Proc, yo, "Setting scatter_graph=0 since the graph"
      " is local on each proc");
    ierr = ZOLTAN_WARN;
  }

  if (gr->scatter){

    if (geo)
      ierr = Zoltan_Scatter_Graph(&gr->vtxdist, &gr->xadj, &gr->adjncy,
                                  &gr->vwgt, (vsp ? &vsp->vsize : NULL),
				  &gr->ewgts, &geo->xyz, geo->ndims,
                                  gr->obj_wgt_dim, zz, &gr->comm_plan);
    else {
      float* xyz = NULL;
      ierr = Zoltan_Scatter_Graph(&gr->vtxdist, &gr->xadj, &gr->adjncy,
                                  &gr->vwgt, (vsp ? &vsp->vsize : NULL),
				  &gr->ewgts, &xyz, 0,
                                  gr->obj_wgt_dim, zz, &gr->comm_plan);
    }
    if ((ierr == ZOLTAN_FATAL) || (ierr == ZOLTAN_MEMERR)){
      ZOLTAN_THIRD_ERROR(ierr, "Zoltan_Scatter_Graph returned error.");
    }
    gr->num_obj = gr->vtxdist[zz->Proc+1]-gr->vtxdist[zz->Proc];

    /* Construct the adjproc array */
    ZOLTAN_FREE(&gr->adjproc);
    gr->adjproc =  (int *) ZOLTAN_MALLOC(gr->xadj[gr->num_obj] * sizeof(int));
    for (i = 0, procid=zz->Proc ; i < gr->xadj[gr->num_obj] ; ++ i) {
      gr->adjproc[i] = give_proc (gr->adjncy[i], gr->vtxdist, zz->Num_Proc, &procid);
    }

    if (prt) {
      prt->part_orig = prt->part;
      prt->part =  (indextype *) ZOLTAN_MALLOC((gr->num_obj+1) * sizeof(indextype));
      Zoltan_Comm_Do(gr->comm_plan, TAG1, (char *) prt->part_orig, sizeof(indextype),
		     (char *) prt->part);

    }

  }
  return ierr;
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

static int scale_round_weights(float *fwgts, indextype *iwgts, int n, int dim,
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
	  else /* sum_wgt[j] != 0 */
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


int Zoltan_Preprocess_Timer(ZZ *zz, int *use_timers)
{
  static int timer_p = -1;

  *use_timers = 0;
  Zoltan_Bind_Param(Graph_params, "USE_TIMERS", (void *) use_timers);
  if (*use_timers) {
    if (timer_p < 0)
      timer_p = Zoltan_Timer_Init(zz->ZTime, 1, "ThirdLibrary");
    ZOLTAN_TIMER_START(zz->ZTime, timer_p, zz->Communicator);
  }

  return (timer_p);
}



void Zoltan_Third_DisplayTime(ZZ* zz, double* times)
{
  if (zz->Proc == zz->Debug_Proc) printf("\nZOLTAN timing statistics:\n");
  Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[1]-times[0],
		     " Partitioner Pre-processing time  ");
  Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[2]-times[1],
		     " Partitioner Library time         ");
  Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[4]-times[3],
		     " Partitioner Post-processing time ");
  Zoltan_Print_Stats(zz->Communicator, zz->Debug_Proc, times[4]-times[3] + times[2] - times[0],
		     " Partitioner Total time           ");
  if (zz->Proc==zz->Debug_Proc) printf("\n");
}


int Zoltan_Third_Set_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
  int  index;
  PARAM_UTYPE result;

  return Zoltan_Check_Param(name, val, Graph_params, &result, &index);
}


#ifdef __cplusplus
}
#endif
