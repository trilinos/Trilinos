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
#include "phg.h"
#include "zoltan_eval.h"

/************************************************************************/
static void iget_strided_stats(int *v, int stride, int offset, int len,
                             float *min, float *max, float *sum);

static void fget_strided_stats(float *v, int stride, int offset, int len,
                             float *min, float *max, float *sum);

static int get_nbor_parts( ZZ *zz, int nobj, ZOLTAN_ID_PTR global_ids, 
  ZOLTAN_ID_PTR local_ids, int *part, int nnbors, ZOLTAN_ID_PTR nbors_global,
  int *nbors_part);

static int *objects_by_partition(ZZ *zz, int num_obj, int *part,
  int *nparts, int *nonempty);

static int
object_metrics(ZZ *zz, int num_obj, int *parts, float *vwgts, int wgt_dim, 
               int *nparts, int *nonempty, float *obj_imbalance, float *imbalance, float *nobj,
               float *obj_wgt, float *xtra_imbalance, float (*xtra_obj_wgt)[EVAL_SIZE]);

static int 
add_graph_extra_weight(ZZ *zz, int num_obj, int *edges_per_obj, int *vwgt_dim, float **vwgts);

/*****************************************************************************/

int Zoltan_LB_Eval_Balance(ZZ *zz, int print_stats, BALANCE_EVAL *eval)
{
  /*****************************************************************************/
  /* Return performance metrics in BALANCE_EVAL structure.                     */ 
  /* Also print them out if print_stats is true.                               */
  /*****************************************************************************/

  char *yo = "Zoltan_LB_Eval_Balance";
  int vwgt_dim = zz->Obj_Weight_Dim;
  int i, ierr;
  int nparts, nonempty_nparts, req_nparts;
  int num_obj = 0;
  BALANCE_EVAL localEval;

  int *parts=NULL;
  float *vwgts=NULL;

  ZOLTAN_ID_PTR global_ids=NULL, local_ids=NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = ZOLTAN_OK;

  if (!eval)
    eval = &localEval;

  memset(eval, 0, sizeof(BALANCE_EVAL));

  /* Get requested number of partitions.  Actual number may differ  */

  ierr = Zoltan_LB_Build_PartDist(zz);
  if (ierr != ZOLTAN_OK){
    goto End;
  }

  req_nparts = zz->LB.Num_Global_Parts;

  /* Get object weights and partitions */

  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, vwgt_dim, &vwgts, &parts);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);

  /* Get metrics based on number of objects and object weights */

  ierr = object_metrics(zz, num_obj, parts, vwgts, vwgt_dim,
          &nparts,           /* actual number of partitions */
          &nonempty_nparts,  /* number of non-empty partitions */
          &eval->obj_imbalance,
          &eval->imbalance,
          eval->nobj,
          eval->obj_wgt,
          eval->xtra_imbalance,
          eval->xtra_obj_wgt);

  if (ierr != ZOLTAN_OK)
    goto End;
     
  /************************************************************************
   * Print results
   */

  if (print_stats && (zz->Proc == zz->Debug_Proc)){

    printf("\n%s  Part count: %1d requested, %1d actual , %1d non-empty\n", 
      yo, req_nparts, nparts, nonempty_nparts);

    printf("%s  Statistics with respect to %1d partitions: \n", yo, nparts);
    printf("%s                             Min      Max      Sum  Imbalance\n", yo);

    printf("%s  Number of objects  :  %8.3g %8.3g %8.3g     %5.3f\n", yo, 
        eval->nobj[EVAL_GLOBAL_MIN], eval->nobj[EVAL_GLOBAL_MAX], 
        eval->nobj[EVAL_GLOBAL_SUM], eval->obj_imbalance);

    if (vwgt_dim > 0){
      printf("%s  Object weight      :  %8.3g %8.3g %8.3g     %5.3f\n", yo, 
        eval->obj_wgt[EVAL_GLOBAL_MIN], eval->obj_wgt[EVAL_GLOBAL_MAX], 
        eval->obj_wgt[EVAL_GLOBAL_SUM], eval->imbalance);

      for (i=0; i < vwgt_dim-1; i++){
        if (i == EVAL_MAX_XTRA_VWGTS){
          break;
        }
        printf("%s  Object weight %d   :  %8.3g %8.3g %8.3g     %5.3f\n", yo, i+2,
          eval->xtra_obj_wgt[i][EVAL_GLOBAL_MIN], eval->xtra_obj_wgt[i][EVAL_GLOBAL_MAX], 
          eval->xtra_obj_wgt[i][EVAL_GLOBAL_SUM], eval->imbalance);
      }
      if (vwgt_dim-1 > EVAL_MAX_XTRA_VWGTS){
        printf("(We calculate up to %d extra object weights.  This can be changed.)\n",
              EVAL_MAX_XTRA_VWGTS);
      }
    }

    printf("\n\n");
  }

End:

  /* Free data */

  ZOLTAN_FREE(&vwgts);
  ZOLTAN_FREE(&parts);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}
/*****************************************************************************/

int Zoltan_LB_Eval_Graph(ZZ *zz, int print_stats, GRAPH_EVAL *graph)
{
  /*****************************************************************************/
  /* Return performance metrics in GRAPH_EVAL structure.                       */ 
  /* Also print them out if print_stats is true.                               */
  /*****************************************************************************/

  char *yo = "Zoltan_LB_Eval_Graph";
  MPI_Comm comm = zz->Communicator;
  int vwgt_dim = zz->Obj_Weight_Dim;
  int ewgt_dim = zz->Edge_Weight_Dim;

  ZOLTAN_ID_PTR global_ids=NULL, local_ids=NULL, nbors_global=NULL;

  int i, j, k, e, ierr;
  int nparts, nonempty_nparts, req_nparts;
  int num_weights, obj_part, nbor_part, ncuts;
  int num_obj = 0;
  int num_edges = 0;

  int *localCount = NULL, *globalCount = NULL;
  int *parts=NULL, *nbors_part=NULL;
  int *edges_per_obj=NULL, *nbors_proc=NULL;
  int *num_boundary=NULL, *cuts=NULL;

  float obj_edge_weights;

  float *vwgts=NULL, *ewgts=NULL, *wgt=NULL;
  float *localVals = NULL, *globalVals = NULL;
  float *cutn=NULL, *cutl=NULL, *cut_wgt=NULL;

  char **cute=NULL;

  GRAPH_EVAL localEval;

  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = ZOLTAN_OK;

  if (!graph)
    graph = &localEval;

  memset(graph, 0, sizeof(GRAPH_EVAL));

  /* Get requested number of partitions.  Actual number may differ  */

  ierr = Zoltan_LB_Build_PartDist(zz);
  if (ierr != ZOLTAN_OK){
    goto End;
  }

  req_nparts = zz->LB.Num_Global_Parts;

  /* Get object weights and partitions */

  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, vwgt_dim, &vwgts, &parts);

  if (ierr != ZOLTAN_OK)
    goto End;

  /*****************************************************************
   * Get graph from query functions
   */

  ierr = Zoltan_Graph_Queries(zz, num_obj, global_ids, local_ids,
                              &num_edges, &edges_per_obj, 
                              &nbors_global, &nbors_proc, &ewgts);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&nbors_proc);

  /*****************************************************************
   * Add a vertex weight if ADD_OBJ_WEIGHT is set
   */

  ierr = add_graph_extra_weight(zz, num_obj, edges_per_obj, &vwgt_dim, &vwgts);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  /*****************************************************************
   * Get metrics based on number of objects and object weights 
   */

  ierr = object_metrics(zz, num_obj, parts, vwgts, vwgt_dim,
          &nparts,          /* actual number of partitions */
          &nonempty_nparts,  /* number of non-empty partitions */
          &graph->obj_imbalance,
          &graph->imbalance,
          graph->nobj,
          graph->obj_wgt,
          graph->xtra_imbalance,
          graph->xtra_obj_wgt);

  if (ierr != ZOLTAN_OK)
    goto End;

  /*****************************************************************
   * Compute the partition number of neighboring objects
   */

  nbors_part = (int *)ZOLTAN_MALLOC(num_edges * sizeof(int));
  if (num_edges && !nbors_part){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  ierr = get_nbor_parts(zz, num_obj, global_ids, local_ids, parts, 
                        num_edges, nbors_global, nbors_part);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&nbors_global);
  ZOLTAN_FREE(&local_ids);

  /*****************************************************************
   * Compute cut statisticics
   */

  cuts = (int *)ZOLTAN_CALLOC(nparts, sizeof(int));
  num_boundary = (int *)ZOLTAN_CALLOC(nparts, sizeof(int));
  cutn = (float *)ZOLTAN_CALLOC(nparts, sizeof(float));
  cutl = (float *)ZOLTAN_CALLOC(nparts, sizeof(float));
  cute = (char **)ZOLTAN_MALLOC(nparts * sizeof(char *));

  if (nparts && (!cuts || !cutn || !cutl || !num_boundary || !cute)){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  for (i=0; i < nparts; i++){
    cute[i] = (char *)ZOLTAN_CALLOC(nparts , sizeof(char));
    if (!cute[i]){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  }

  cut_wgt = (float *)ZOLTAN_CALLOC(nparts * ewgt_dim, sizeof(float));

  if (nparts && ewgt_dim && !cut_wgt){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  for (i=0,k=0; i < num_obj; i++){   /* object */

    obj_edge_weights = 0;
    obj_part = parts[i];
    ncuts = 0;

    for (j=0; j < edges_per_obj[i]; j++,k++){    /* neighbor in graph */

      nbor_part = nbors_part[k];

      if (ewgt_dim > 0){
        obj_edge_weights += ewgts[k * ewgt_dim];  /* "hypergraph" weight */
      }

      if (nbor_part != obj_part){
        /* 
         * number of edges that have nbor in a different partition 
         */
        cuts[obj_part]++; 

        /*
         * save this info so that we can count for each partition,
         * the number of partitions that it has neighbors in
         */
        cute[obj_part][nbor_part] = 1;

        for (e=0; e < ewgt_dim; e++){
          /*
           * For each partition, the sum of the weights of the edges
           * whos neighbor is in a different partition
           */
          cut_wgt[obj_part * ewgt_dim + e] += ewgts[k * ewgt_dim + e];
        }

        ncuts++;
      }
    }

    if (ncuts){
      /*
       * hypergraph ConCut measure
       */
      cutn[obj_part] += (obj_edge_weights * ncuts);

      /*
       * hypergraph NetCut measure
       */
      cutl[obj_part] += obj_edge_weights;

      /*
       * for each partition, the number of objects with a neighbor outside
       * the partition
       */
      num_boundary[obj_part]++;
    }
  }

  ZOLTAN_FREE(&parts);
  ZOLTAN_FREE(&edges_per_obj);
  ZOLTAN_FREE(&ewgts);
  ZOLTAN_FREE(&nbors_part);

  /************************************************************************
   * Write cut statistics to the return structure.
   */

  k = ((ewgt_dim > 0) ? ewgt_dim : 1);

  globalVals = (float *)ZOLTAN_MALLOC(nparts * ewgt_dim * sizeof(float));
  if (nparts && ewgt_dim && !globalVals){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  globalCount = (int *)ZOLTAN_MALLOC(nparts * sizeof(int));
  if (nparts && !globalCount){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /*
   * CUTN - ConCut
   */

  MPI_Allreduce(cutn, globalVals, nparts, MPI_FLOAT, MPI_SUM, comm);

  ZOLTAN_FREE(&cutn);

  fget_strided_stats(globalVals, 1, 0, nparts,
               graph->cutn + EVAL_GLOBAL_MIN,
               graph->cutn + EVAL_GLOBAL_MAX,
               graph->cutn + EVAL_GLOBAL_SUM);

  graph->cutn[EVAL_GLOBAL_AVG] = graph->cutn[EVAL_GLOBAL_SUM] / nparts;

  /*
   * CUTL - NetCut
   */

  MPI_Allreduce(cutl, globalVals, nparts, MPI_FLOAT, MPI_SUM, comm);

  ZOLTAN_FREE(&cutl);

  fget_strided_stats(globalVals, 1, 0, nparts,
               graph->cutl + EVAL_GLOBAL_MIN,
               graph->cutl + EVAL_GLOBAL_MAX,
               graph->cutl + EVAL_GLOBAL_SUM);

  graph->cutl[EVAL_GLOBAL_AVG] = graph->cutl[EVAL_GLOBAL_SUM] / nparts;

  /*
   * CUTS - Number of cut edges in each partition
   */

  MPI_Allreduce(cuts, globalCount, nparts, MPI_INT, MPI_SUM, comm);

  ZOLTAN_FREE(&cuts);

  iget_strided_stats(globalCount, 1, 0, nparts,
               graph->cuts + EVAL_GLOBAL_MIN,
               graph->cuts + EVAL_GLOBAL_MAX,
               graph->cuts + EVAL_GLOBAL_SUM);

  graph->cuts[EVAL_GLOBAL_AVG] = graph->cuts[EVAL_GLOBAL_SUM] / nparts;

  /*
   * CUTE - The number of neighboring partitions
   */

  localCount = (int *)ZOLTAN_MALLOC(nparts * sizeof(int));

  for (i=0; i < nparts; i++){
    localCount[i] = 0;
    for (j=0; j < nparts; j++){
      if (cute[i][j]) localCount[i]++;
    }
    ZOLTAN_FREE(&cute[i]);
  }
  ZOLTAN_FREE(&cute);

  MPI_Allreduce(localCount, globalCount, nparts, MPI_INT, MPI_SUM, comm);

  ZOLTAN_FREE(&localCount);

  iget_strided_stats(globalCount, 1, 0, nparts,
               graph->cute + EVAL_GLOBAL_MIN,
               graph->cute + EVAL_GLOBAL_MAX,
               graph->cute + EVAL_GLOBAL_SUM);

  graph->cute[EVAL_GLOBAL_AVG] = graph->cute[EVAL_GLOBAL_SUM] / nparts;

  /*
   * CUT WEIGHT - The sum of the weights of the cut edges.
   */

  num_weights = nparts * ewgt_dim;

  MPI_Allreduce(cut_wgt, globalVals, num_weights, MPI_FLOAT, MPI_SUM, comm);

  ZOLTAN_FREE(&cut_wgt);

  fget_strided_stats(globalVals, ewgt_dim, 0, num_weights,
               graph->cut_wgt + EVAL_GLOBAL_MIN,
               graph->cut_wgt + EVAL_GLOBAL_MAX,
               graph->cut_wgt + EVAL_GLOBAL_SUM);

  graph->cut_wgt[EVAL_GLOBAL_AVG] = graph->cut_wgt[EVAL_GLOBAL_SUM] / nparts;

  for (i=0; i < ewgt_dim-1; i++){
    /* end of calculations for multiple edge weights */

    if (i == EVAL_MAX_XTRA_EWGTS){
      break;
    }
   
    wgt = graph->xtra_cut_wgt[i];

    fget_strided_stats(globalVals, ewgt_dim, i+1, num_weights,
                       wgt + EVAL_GLOBAL_MIN,
                       wgt + EVAL_GLOBAL_MAX,
                       wgt + EVAL_GLOBAL_SUM);

    wgt[EVAL_GLOBAL_AVG] = wgt[EVAL_GLOBAL_SUM]/(float)nparts;
  }

  ZOLTAN_FREE(&globalVals);

  /*
   * Number of objects in partition that have an off-partition neighbor.
   */

  MPI_Allreduce(num_boundary, globalCount, nparts, MPI_INT, MPI_SUM, comm);

  iget_strided_stats(globalCount, 1, 0, nparts,
               graph->num_boundary + EVAL_GLOBAL_MIN,
               graph->num_boundary + EVAL_GLOBAL_MAX,
               graph->num_boundary + EVAL_GLOBAL_SUM);

  graph->num_boundary[EVAL_GLOBAL_AVG] = graph->num_boundary[EVAL_GLOBAL_SUM] / nparts;

  ZOLTAN_FREE(&num_boundary);
  ZOLTAN_FREE(&globalCount);

  /************************************************************************
   * Print results
   */

  if (print_stats && (zz->Proc == zz->Debug_Proc)){

    printf("\n%s  Part count: %1d requested, %1d actual, %1d non-empty\n", 
      yo, req_nparts, nparts, nonempty_nparts);

    printf("%s  Statistics with respect to %1d partitions: \n", yo, nparts);
    printf("%s                             Min      Max      Sum  Imbalance\n", yo);

    printf("%s  Number of objects  :  %8.3g %8.3g %8.3g     %5.3g\n", yo, 
      graph->nobj[EVAL_GLOBAL_MIN], graph->nobj[EVAL_GLOBAL_MAX],
      graph->nobj[EVAL_GLOBAL_SUM], graph->obj_imbalance);

    if (vwgt_dim > 0){
      printf("%s  Object weight      :  %8.3g %8.3g %8.3g     %5.3f\n", yo, 
        graph->obj_wgt[EVAL_GLOBAL_MIN], graph->obj_wgt[EVAL_GLOBAL_MAX], 
        graph->obj_wgt[EVAL_GLOBAL_SUM], graph->imbalance);
  
      for (i=0; i < vwgt_dim-1; i++){
        if (i == EVAL_MAX_XTRA_VWGTS){
          break;
        }
        printf("%s  Object weight %d    :  %8.3g %8.3g %8.3g     %5.3f\n", yo, i+2,
          graph->xtra_obj_wgt[i][EVAL_GLOBAL_MIN], graph->xtra_obj_wgt[i][EVAL_GLOBAL_MAX], 
          graph->xtra_obj_wgt[i][EVAL_GLOBAL_SUM], graph->imbalance);
      }
      if (vwgt_dim-1 > EVAL_MAX_XTRA_VWGTS){
        printf("(We calculate up to %d extra object weights.  This can be changed.)\n",
              EVAL_MAX_XTRA_VWGTS);
      }
    }

    printf("\n\n");

    printf("%s  Statistics with respect to %1d partitions: \n", yo, nparts);
    printf("%s                               Min      Max    Average    Sum\n", yo);

    printf("%s  Num boundary objects :  %8.3g %8.3g %8.3g %8.3g\n", yo, 
      graph->num_boundary[EVAL_GLOBAL_MIN], graph->num_boundary[EVAL_GLOBAL_MAX], 
      graph->num_boundary[EVAL_GLOBAL_AVG], graph->num_boundary[EVAL_GLOBAL_SUM]);

    printf("%s  Number of cut edges  :  %8.3g %8.3g %8.3g %8.3g\n", yo, 
      graph->cuts[EVAL_GLOBAL_MIN], graph->cuts[EVAL_GLOBAL_MAX], graph->cuts[EVAL_GLOBAL_AVG],
      graph->cuts[EVAL_GLOBAL_SUM]);

    printf("%s  Weight of cut edges  :  %8.3g %8.3g %8.3g %8.3g\n", yo, 
      graph->cut_wgt[EVAL_GLOBAL_MIN], graph->cut_wgt[EVAL_GLOBAL_MAX], 
      graph->cut_wgt[EVAL_GLOBAL_AVG], graph->cut_wgt[EVAL_GLOBAL_SUM]);

    for (i=0; i < ewgt_dim-1; i++){
      if (i == EVAL_MAX_XTRA_EWGTS){
        break;
      }
      printf("%s  Weight %d            :  %8.3g %8.3g %8.3g %8.3g\n", yo, i+2,
        graph->xtra_cut_wgt[i][EVAL_GLOBAL_MIN], graph->xtra_cut_wgt[i][EVAL_GLOBAL_MAX], 
        graph->xtra_cut_wgt[i][EVAL_GLOBAL_AVG], graph->xtra_cut_wgt[i][EVAL_GLOBAL_SUM]);
    }
    if (ewgt_dim-1 > EVAL_MAX_XTRA_EWGTS){
      printf("(We calculate up to %d extra edge weights.  This can be changed.)\n",
            EVAL_MAX_XTRA_EWGTS);
    }

    printf("%s  Num Nbor Parts (CUTE):  %8.3g %8.3g %8.3g %8.3g\n", yo, 
      graph->cute[EVAL_GLOBAL_MIN], graph->cute[EVAL_GLOBAL_MAX], 
      graph->cute[EVAL_GLOBAL_AVG], graph->cute[EVAL_GLOBAL_SUM]);

    printf("%s  Conn Cut (CUTN)      :  %8.3g %8.3g %8.3g %8.3g\n", yo, 
      graph->cutn[EVAL_GLOBAL_MIN], graph->cutn[EVAL_GLOBAL_MAX], 
      graph->cutn[EVAL_GLOBAL_AVG], graph->cutn[EVAL_GLOBAL_SUM]);

    printf("%s  Net Cut (CUTL)       :  %8.3g %8.3g %8.3g %8.3g\n", yo, 
      graph->cutl[EVAL_GLOBAL_MIN], graph->cutl[EVAL_GLOBAL_MAX], 
      graph->cutl[EVAL_GLOBAL_AVG], graph->cutl[EVAL_GLOBAL_SUM]);
    
    printf("\n\n");
  }

End:

  /* Free data */

  ZOLTAN_FREE(&localVals);
  ZOLTAN_FREE(&localCount);
  ZOLTAN_FREE(&globalVals);
  ZOLTAN_FREE(&globalCount);
  ZOLTAN_FREE(&vwgts);
  ZOLTAN_FREE(&parts);
  ZOLTAN_FREE(&nbors_proc);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&edges_per_obj);
  ZOLTAN_FREE(&nbors_global);
  ZOLTAN_FREE(&ewgts);
  ZOLTAN_FREE(&nbors_part);
  ZOLTAN_FREE(&num_boundary);
  ZOLTAN_FREE(&cut_wgt);
  ZOLTAN_FREE(&cutn);
  ZOLTAN_FREE(&cutl);
  ZOLTAN_FREE(&cuts);

  if (cute){
    for (i=0; i < nparts; i++)
      ZOLTAN_FREE(&cute[i]);
    ZOLTAN_FREE(&cute);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}
/*****************************************************************************/

int Zoltan_LB_Eval_HG(ZZ *zz, int print_stats, HG_EVAL *hg)
{
  /*****************************************************************************/
  /* Return performance metrics in HG_EVAL structure.  Also print them out     */
  /* if print_stats is true.  Results are per partition, not per process.      */
  /*****************************************************************************/

  char *yo = "Zoltan_LB_Eval_HG";
  MPI_Comm comm = zz->Communicator;

  float *part_sizes=NULL;
  float *localVals=NULL;

  double local[2], global[2];

  int ierr, debug_level;
  int nparts, nonempty_nparts, req_nparts;
  int *localCount = NULL;
  PHGPartParams hgp;

  ZHG* zhg=NULL;

  HG_EVAL localEval;
  
  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Set default error code */
  ierr = ZOLTAN_OK;

  if (!hg)
    hg = &localEval;

  memset(hg, 0, sizeof(HG_EVAL));

  /* get the requested number of partitions */ 

  ierr = Zoltan_LB_Build_PartDist(zz);
  if (ierr != ZOLTAN_OK){
    goto End;
  }

  req_nparts = zz->LB.Num_Global_Parts;

  /*****************************************************************
   * Get the hypergraph, via the hypergraph or graph query functions
   */

  part_sizes = (float*)ZOLTAN_MALLOC(sizeof(float) * req_nparts);
  if (req_nparts && !part_sizes){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  Zoltan_LB_Get_Part_Sizes(zz, zz->LB.Num_Global_Parts, 1, part_sizes);

  debug_level = zz->Debug_Level;
  zz->Debug_Level = 0;

  ierr = Zoltan_PHG_Initialize_Params(zz, part_sizes, &hgp);
  if (ierr != ZOLTAN_OK)
    goto End;

  zz->Debug_Level = debug_level;

  zhg = (ZHG*) ZOLTAN_MALLOC (sizeof(ZHG));
  if (zhg == NULL){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  Zoltan_Input_HG_Init(zhg);

  ierr = Zoltan_Get_Hypergraph_From_Queries(zz, &hgp, zhg);
  if (ierr != ZOLTAN_OK)
    goto End;

  if (zhg->globalObj == 0){
    if (zz->Proc == zz->Debug_Proc){
      printf("%s: No objects to evaluate\n",yo);
    }
    goto End;
  }

  /* Get metrics based on number of objects and object weights */

  ierr = object_metrics(zz, zhg->nObj, zhg->Input_Parts, zhg->objWeight, zhg->objWeightDim,
          &nparts,          /* actual number of partitions */
          &nonempty_nparts,  /* number of non-empty partitions */
          &hg->obj_imbalance,
          &hg->imbalance,
          hg->nobj,
          hg->obj_wgt,
          NULL, NULL);

  if (ierr != ZOLTAN_OK)
    goto End;

  /************************************************************************
   * Compute the cutn and cutl 
   */

  if (!zhg->Output_Parts)
    zhg->Output_Parts = zhg->Input_Parts;

  ierr = Zoltan_PHG_Cuts(zz, zhg, local);

  if (zhg->Output_Parts == zhg->Input_Parts)
    zhg->Output_Parts = NULL;

  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN){
    goto End;
  }

  MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, comm);

  hg->cutn = (float)global[0];
  hg->cutl = (float)global[1];
            
  /************************************************************************
   * Print results
   */

  if (print_stats && (zz->Proc == zz->Debug_Proc)){

    printf("\n%s  Part count: %1d requested, %1d actual, %1d non-empty\n", 
      yo, req_nparts, nparts, nonempty_nparts);

    printf("%s  Statistics with respect to %1d partitions: \n", yo, nparts);
    printf("%s                            Min      Max     Sum  Imbalance\n", yo);

    printf("%s  Number of objects  :  %8.3g %8.3g %8.3g   %5.3f\n", yo, 
      hg->nobj[EVAL_GLOBAL_MIN], hg->nobj[EVAL_GLOBAL_MAX], 
      hg->nobj[EVAL_GLOBAL_SUM], hg->obj_imbalance);

    if (zhg->objWeight > 0){
      printf("%s  Object weight      :  %8.3g %8.3g %8.3g   %5.3f\n", yo, 
        hg->obj_wgt[EVAL_GLOBAL_MIN], hg->obj_wgt[EVAL_GLOBAL_MAX], 
        hg->obj_wgt[EVAL_GLOBAL_SUM], hg->imbalance);
    }
    
    printf("\n");

    printf("%s  Hyperedge (k-1)-connectivity cut:     %8.0f\n", yo, global[0]);
    printf("%s  No. cut hyperedges:                   %8.0f\n\n", yo, global[1]);
  }

End:

  /* Free data */

  ZOLTAN_FREE(&part_sizes);
  if (zhg){
    Zoltan_PHG_Free_Hypergraph_Data(zhg);
    ZOLTAN_FREE(&zhg);
  }
  ZOLTAN_FREE(&localVals);
  ZOLTAN_FREE(&localCount);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
/************************************************************************/

static void iget_strided_stats(int *v, int stride, int offset, int len,
                             float *min, float *max, float *sum)
{
  int i;
  int *v0;
  float fmin, fmax, fsum;

  if (len < 1) return;

  v0 = v + offset;

  fmin = fmax = fsum = v0[0];

  for (i=stride; i < len; i += stride){
    if (v0[i] < fmin) fmin = (float)v0[i];
    else if (v0[i] > fmax) fmax = (float)v0[i];
    fsum += (float)v0[i];
  }

  *min = fmin;
  *max = fmax;
  *sum = fsum;
}
/************************************************************************/
static void fget_strided_stats(float *v, int stride, int offset, int len,
                             float *min, float *max, float *sum)
{
  int i;
  float *v0;
  float fmin, fmax, fsum;

  if (len < 1) return;

  v0 = v + offset;

  fmin = fmax = fsum = v0[0];

  for (i=stride; i < len; i += stride){
    if (v0[i] < fmin) fmin = v0[i];
    else if (v0[i] > fmax) fmax = v0[i];
    fsum += v0[i];
  }

  *min = fmin;
  *max = fmax;
  *sum = fsum;
}

/************************************************************************/
#define TEST_DD_ERROR(ierr, yo, proc, fn) \
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) { \
    char msg[64];  \
    sprintf(msg, "Error in %s\n", fn); \
    ZOLTAN_PRINT_ERROR(proc, yo, msg); \
    goto End;\
  }


static int get_nbor_parts(
  ZZ *zz,
  int nobj,                     /* Input:  number of local objs */
  ZOLTAN_ID_PTR global_ids,     /* Input:  GIDs of local objs */
  ZOLTAN_ID_PTR local_ids,      /* Input:  LIDs of local objs */
  int *part,                    /* Input:  part assignments of local objs */
  int nnbors,                   /* Input:  number of neighboring objs */
  ZOLTAN_ID_PTR nbors_global,   /* Input:  GIDs of neighboring objs */
  int *nbors_part               /* Output:  part assignments of neighbor objs */
)
{
/* Function to retrieve the partition number for neighboring nodes. */
char *yo = "get_nbor_parts";
struct Zoltan_DD_Struct *dd = NULL;
int *owner = NULL;
int maxnobj;
int ierr;
int start, size, i_am_done, alldone;
const int MAXSIZE = 200000;

  ZOLTAN_TRACE_ENTER(zz, yo);

  MPI_Allreduce(&nobj, &maxnobj, 1, MPI_INT, MPI_MAX, zz->Communicator);
  ierr = Zoltan_DD_Create(&dd, zz->Communicator, zz->Num_GID, zz->Num_LID,
                          0, MIN(maxnobj,MAXSIZE), 0);
  TEST_DD_ERROR(ierr, yo, zz->Proc, "Zoltan_DD_Create");

  ierr = Zoltan_DD_Update(dd, global_ids, local_ids, NULL, part, nobj);
  TEST_DD_ERROR(ierr, yo, zz->Proc, "Zoltan_DD_Update");

  /* Do the find in chunks to avoid swamping memory. */
  owner = (int *) ZOLTAN_MALLOC(MIN(MAXSIZE,nnbors) * sizeof(int));
  start = 0;
  alldone = 0;
  while (!alldone) {
    size = MIN(MAXSIZE, (nnbors-start > 0 ? nnbors-start : 0));
    if (start < nnbors)
      ierr = Zoltan_DD_Find(dd, &nbors_global[start*zz->Num_GID],
                            NULL, NULL, &nbors_part[start],
                            size, owner);
    else  /* I'm done, but other processors might not be */
      ierr = Zoltan_DD_Find(dd, NULL, NULL, NULL, NULL, 0, NULL);
    start += size;
    i_am_done = (nnbors - start > 0 ? 0 : 1);
    MPI_Allreduce(&i_am_done, &alldone, 1, MPI_INT, MPI_MIN, zz->Communicator);
  }

  ZOLTAN_FREE(&owner);
  TEST_DD_ERROR(ierr, yo, zz->Proc, "Zoltan_DD_Find");

End:
  Zoltan_DD_Destroy(&dd);
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/************************************************************************/

static int *
objects_by_partition(ZZ *zz, int num_obj, int *part, int *nparts, int *nonempty)
{
  char *yo = "objects_by_partition";
  int i, num_parts, num_nonempty, max_part, gmax_part;
  int *partCounts = NULL, *totalCounts;
  int *returnBuf = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  max_part = part[0];
  for (i=1; i < num_obj; i++){
    if (part[i] > max_part) max_part = part[i];
  }

  MPI_Allreduce(&max_part, &gmax_part, 1, MPI_INT, MPI_MAX, zz->Communicator);

  max_part = gmax_part + 1;

  /* Allocate and return a buffer containing the local count of objects by partition,
     followed by the global count.  Set "nparts" to the number of partitions, and
     set "nonempty" to the count of those that have objects in them.
   */

  partCounts = (int *)ZOLTAN_CALLOC(max_part * 2, sizeof(int));
  if (max_part && !partCounts)
    return NULL;

  totalCounts = partCounts + max_part;

  for (i=0; i < num_obj; i++)
    partCounts[part[i]]++;

  MPI_Allreduce(partCounts, totalCounts, max_part, MPI_INT, MPI_SUM, zz->Communicator);

  num_parts = max_part;

  for (i=max_part-1; i > 0; i--){
    if (totalCounts[i] > 0) break;
    num_parts--;
  }

  returnBuf = (int *)ZOLTAN_MALLOC(num_parts * sizeof(int));
  if (num_parts && !returnBuf)
    return NULL;

  memcpy(returnBuf, totalCounts, sizeof(int) * num_parts);

  ZOLTAN_FREE(&partCounts);

  num_nonempty = 0;

  for (i=0; i < num_parts; i++){
    if (returnBuf[i] > 0) num_nonempty++;
  }

  *nparts = num_parts;
  *nonempty = num_nonempty;

  ZOLTAN_TRACE_EXIT(zz, yo);
  return returnBuf;
}

/************************************************************************/

static int
object_metrics(ZZ *zz, int num_obj, int *parts, float *vwgts, int vwgt_dim,
               int *nparts,       /* return actual number of partitions */
               int *nonempty,     /* return number of those that are non-empty */
               float *obj_imbalance,  /* object # imbalance wrt partitions */
               float *imbalance,      /* object wgt imbalance wrt partitions */
               float *nobj,       /* number of objects */
               float *obj_wgt,    /* object weights */
               float *xtra_imbalance,  /* return if vertex weight dim > 1 */
    float (*xtra_obj_wgt)[EVAL_SIZE])  /* return if vertex weight dim > 1 */
{
  char *yo = "object_metrics";
  MPI_Comm comm = zz->Communicator;

  int i, j, ierr;
  int num_weights, num_parts, num_nonempty_parts;

  int *globalCount = NULL;

  float *wgt=NULL;
  float *localVals = NULL, *globalVals = NULL;

  ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Get the actual number of partitions, and number of objects per partition */

  globalCount = objects_by_partition(zz, num_obj, parts,
                &num_parts,           /* actual number of partitions */
                &num_nonempty_parts); /* number of those which are non-empty */

  if (!globalCount){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  *nparts = num_parts;
  *nonempty = num_nonempty_parts;

  iget_strided_stats(globalCount, 1, 0, num_parts,
                     nobj + EVAL_GLOBAL_MIN,
                     nobj + EVAL_GLOBAL_MAX,
                     nobj + EVAL_GLOBAL_SUM);

  nobj[EVAL_GLOBAL_AVG] = nobj[EVAL_GLOBAL_SUM]/(float)num_parts;
 
  *obj_imbalance = (nobj[EVAL_GLOBAL_SUM] > 0 ? 
      nobj[EVAL_GLOBAL_MAX]*num_parts / nobj[EVAL_GLOBAL_SUM]: 1);

  ZOLTAN_FREE(&globalCount);

  if (vwgt_dim > 0){

    num_weights = num_parts * vwgt_dim;
  
    localVals = (float *)ZOLTAN_CALLOC(num_weights*2, sizeof(float));
  
    if (num_weights && !localVals){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  
    globalVals = localVals + num_weights;
  
    if (vwgt_dim>0){
      for (i=0; i<num_obj; i++){
        for (j=0; j<vwgt_dim; j++){
          localVals[parts[i]*vwgt_dim+j] += vwgts[i*vwgt_dim+j];
        }
      }
    }
  
    MPI_Allreduce(localVals, globalVals, num_weights, MPI_FLOAT, MPI_SUM, comm);
  
    fget_strided_stats(globalVals, vwgt_dim, 0, num_weights,
                       obj_wgt + EVAL_GLOBAL_MIN,
                       obj_wgt + EVAL_GLOBAL_MAX,
                       obj_wgt + EVAL_GLOBAL_SUM);
  
    obj_wgt[EVAL_GLOBAL_AVG] = obj_wgt[EVAL_GLOBAL_SUM]/(float)num_parts;
  
    *imbalance = (obj_wgt[EVAL_GLOBAL_SUM] > 0 ? 
        obj_wgt[EVAL_GLOBAL_MAX]*num_parts / obj_wgt[EVAL_GLOBAL_SUM]: 1);
  
    for (i=0; i < vwgt_dim-1; i++){
      /* calculations for multiple vertex weights */
  
      if (i == EVAL_MAX_XTRA_VWGTS){
        break;
      }
     
      wgt = xtra_obj_wgt[i];
  
      fget_strided_stats(globalVals, vwgt_dim, i+1, num_weights,
                         wgt + EVAL_GLOBAL_MIN,
                         wgt + EVAL_GLOBAL_MAX,
                         wgt + EVAL_GLOBAL_SUM);
  
      wgt[EVAL_GLOBAL_AVG] = wgt[EVAL_GLOBAL_SUM]/(float)num_parts;
  
      xtra_imbalance[i] = (wgt[EVAL_GLOBAL_SUM] > 0 ? 
            (wgt[EVAL_GLOBAL_MAX]*num_parts / wgt[EVAL_GLOBAL_SUM]) : 1);
    }
  }
  else{
    /* assume an object weight of one per object */
    for (i=0; i < EVAL_SIZE; i++){
      obj_wgt[i] = nobj[i];
    }
    *imbalance = *obj_imbalance;
  }

  globalVals = NULL;

End:
  ZOLTAN_FREE(&localVals);
  ZOLTAN_FREE(&globalCount);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/************************************************************************/

static int 
add_graph_extra_weight(ZZ *zz, int num_obj, int *edges_per_obj, int *vwgt_dim, float **vwgts)
{
  char *yo = "add_graph_extra_weight";
  PARAM_VARS params[2] = 
    {{ "ADD_OBJ_WEIGHT",  NULL,  "STRING", 0},
     { NULL, NULL, NULL, 0 } };
  char add_obj_weight[64];
  int ierr, i, j;
  int add_type = 0;
  float *tmpnew, *tmpold;
  float *weights = NULL;
  int weight_dim = 0;

  ierr = ZOLTAN_OK;

  strcpy(add_obj_weight, "NONE");
  Zoltan_Bind_Param(params, "ADD_OBJ_WEIGHT", (void *) add_obj_weight);
  Zoltan_Assign_Param_Vals(zz->Params, params, 0, zz->Proc, zz->Debug_Proc);

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

  if (add_type > 0){

    weight_dim = *vwgt_dim + 1;
    weights = (float *)ZOLTAN_CALLOC(num_obj * weight_dim , sizeof(float));
    if (num_obj && weight_dim && !weights){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
   
    tmpnew = weights;
    tmpold = *vwgts;

    for (i=0; i < num_obj; i++){
      for (j=0; j < *vwgt_dim; j++){
        *tmpnew++ = *tmpold++;
      }
      if (add_type == 1){
        *tmpnew++ = 1.0;
      }
      else{
        *tmpnew++ = edges_per_obj[i];
      }
    }

    tmpold = *vwgts;
    ZOLTAN_FREE(&tmpold);
    *vwgt_dim = weight_dim;
    *vwgts = weights;
  }

End:

  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
