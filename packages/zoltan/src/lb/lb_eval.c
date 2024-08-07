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
#include "zz_util_const.h"
#include "phg.h"
#include "zoltan_eval.h"

/************************************************************************/
static void iget_strided_stats(int *v, int stride, int offset, int len,
                               float *min, float *max, float *sum);

static void fget_strided_stats(float *v, int stride, int offset, int len,
                               float *min, float *max, float *sum);

static int get_nbor_parts(ZZ *zz, int nobj, ZOLTAN_ID_PTR global_ids, 
                          ZOLTAN_ID_PTR local_ids, int *part, int nnbors,
                          ZOLTAN_ID_PTR nbors_global, int *nbors_part);

static int *objects_by_part(ZZ *zz, int num_obj, int *part,
                            int *nparts, int *nonempty);

static int object_metrics(ZZ *zz, int num_obj, int *parts, float *part_sizes,
                          int req_nparts, float *vwgts, int wgt_dim, 
                          int *nparts, int *nonempty, float *obj_imbalance, 
                          float *imbalance, float *nobj, float *obj_wgt, 
                          float *xtra_imbalance, 
                          float (*xtra_obj_wgt)[EVAL_SIZE]);

static int add_graph_extra_weight(ZZ *zz, int num_obj, int *edges_per_obj, 
                                  int *vwgt_dim, float **vwgts);

extern int zoltan_lb_eval_sort_increasing(const void *a, const void *b);

/*****************************************************************************/

int Zoltan_LB_Eval_Balance(ZZ *zzin, int print_stats, ZOLTAN_BALANCE_EVAL *eval)
{
  /***************************************************************************/
  /* Return performance metrics in ZOLTAN_BALANCE_EVAL structure.            */ 
  /* Also print them out if print_stats is true.                             */
  /***************************************************************************/

  char *yo = "Zoltan_LB_Eval_Balance";
  ZZ *zz = Zoltan_Copy(zzin);  /* Some operations have side-effects in zz;
                                  Don't change the user's settings */
  int vwgt_dim = zz->Obj_Weight_Dim;
  int part_dim = 0;
  int i, j, ierr;
  int nparts, nonempty_nparts, req_nparts;
  int num_obj = 0;
  ZOLTAN_BALANCE_EVAL localEval;

  int *parts=NULL;
  float *vwgts=NULL;
  float *part_sizes = NULL;

  ZOLTAN_ID_PTR global_ids=NULL, local_ids=NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = ZOLTAN_OK;

  if (!eval)
    eval = &localEval;

  memset(eval, 0, sizeof(ZOLTAN_BALANCE_EVAL));

  /* Get requested number of parts.  Actual number may differ  */

  ierr = Zoltan_LB_Build_PartDist(zz);
  if (ierr != ZOLTAN_OK){
    goto End;
  }

  req_nparts = zz->LB.Num_Global_Parts;

  /* Get object weights and parts */

  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, 
                             vwgt_dim, &vwgts, &parts);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);

  /* Local stats */

  eval->nobj[EVAL_LOCAL_SUM] = num_obj;

  if (vwgt_dim > 0){
    for (i=0; i < num_obj; i++){
      eval->obj_wgt[EVAL_LOCAL_SUM]  += vwgts[i*vwgt_dim];
      for (j=1; j <= EVAL_MAX_XTRA_VWGTS; j++){
        if (j == vwgt_dim) break; 
        eval->xtra_obj_wgt[j-1][EVAL_LOCAL_SUM]  += vwgts[i*vwgt_dim + j];
      }
    }
  }

  /* Get the relative size of each partition */

  part_dim = (vwgt_dim > 0) ? vwgt_dim : 1;

  part_sizes = (float*)ZOLTAN_MALLOC(sizeof(float) * part_dim * req_nparts);
  if (req_nparts && !part_sizes){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  Zoltan_LB_Get_Part_Sizes(zz, part_dim, part_sizes);

  /* Get metrics based on number of objects and object weights */

  ierr = object_metrics(zz, num_obj, parts, part_sizes, req_nparts,
          vwgts, vwgt_dim,
          &nparts,           /* actual number of parts */
          &nonempty_nparts,  /* number of non-empty parts */
          &eval->obj_imbalance,
          &eval->imbalance,
          eval->nobj,
          eval->obj_wgt,
          eval->xtra_imbalance,
          eval->xtra_obj_wgt);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&part_sizes);
     
  /************************************************************************
   * Print results
   */

  if (print_stats && (zz->Proc == zz->Debug_Proc)){

    printf("\n%s  Part count: %1d requested, %1d actual , %1d non-empty\n", 
           yo, req_nparts, nparts, nonempty_nparts);

    printf("%s  Statistics with respect to %1d parts: \n", yo, nparts);
    printf("%s                             Min      Max      Sum  Imbalance\n",
           yo);

    printf("%s  Number of objects  :  %8.3g %8.3g %8.3g", yo, 
           eval->nobj[EVAL_GLOBAL_MIN], eval->nobj[EVAL_GLOBAL_MAX], 
           eval->nobj[EVAL_GLOBAL_SUM]);

    if (eval->obj_imbalance >= 0){
      printf("     %5.3f\n", eval->obj_imbalance);
    }
    else{
      printf("     ----\n");
    }

    if (vwgt_dim > 0){
      printf("%s  Object weight      :  %8.3g %8.3g %8.3g", yo, 
             eval->obj_wgt[EVAL_GLOBAL_MIN], eval->obj_wgt[EVAL_GLOBAL_MAX], 
             eval->obj_wgt[EVAL_GLOBAL_SUM]);

      if (eval->imbalance >= 0){
        printf("     %5.3f\n", eval->imbalance);
      }
      else{
        printf("     ----\n");
      }

      for (i=0; i < vwgt_dim-1; i++){
        if (i == EVAL_MAX_XTRA_VWGTS){
          break;
        }
        printf("%s  Object weight %d    :  %8.3g %8.3g %8.3g", yo, i+2,
               eval->xtra_obj_wgt[i][EVAL_GLOBAL_MIN], 
               eval->xtra_obj_wgt[i][EVAL_GLOBAL_MAX], 
               eval->xtra_obj_wgt[i][EVAL_GLOBAL_SUM]);

        if (eval->xtra_imbalance[i] >= 0){
          printf("     %5.3f\n", eval->xtra_imbalance[i] );
        }
        else{
          printf("     ----\n");
        }
      }
      if (vwgt_dim-1 > EVAL_MAX_XTRA_VWGTS){
        printf("(We calculate up to %d extra object weights.  "
               "This can be changed.)\n",
              EVAL_MAX_XTRA_VWGTS);
      }
    }

    printf("\n\n");
  }

End:

  /* Free data */

  ZOLTAN_FREE(&vwgts);
  ZOLTAN_FREE(&parts);
  ZOLTAN_FREE(&part_sizes);

  ZOLTAN_TRACE_EXIT(zz, yo);
  Zoltan_Destroy(&zz);
  return ierr;
}
/*****************************************************************************/

int Zoltan_LB_Eval_Graph(ZZ *zzin, int print_stats, ZOLTAN_GRAPH_EVAL *graph)
{
  /****************************************************************************/
  /* Return performance metrics in ZOLTAN_GRAPH_EVAL structure.               */
  /* Also print them out if print_stats is true.                              */
  /****************************************************************************/

  char *yo = "Zoltan_LB_Eval_Graph";
  ZZ *zz = Zoltan_Copy(zzin);  /* Some operations have side-effects in zz;
                                  Don't change the user's settings */
  MPI_Comm comm = zz->Communicator;
  int vwgt_dim = zz->Obj_Weight_Dim;
  int ewgt_dim = zz->Edge_Weight_Dim;
  int orig_vwgt_dim, part_dim, eval_vwgt_dim;
  int fromidx, toidx;

  ZOLTAN_ID_PTR global_ids=NULL, local_ids=NULL, nbors_global=NULL;

  int i, j, k, e, ierr, count;
  int nparts, nonempty_nparts, req_nparts;
  int num_weights=0, obj_part, nbor_part, nother_parts;
  int num_pairs, num_parts;
  ZOLTAN_MAP *map = NULL;
  int num_obj = 0;
  int num_edges = 0;
  int hashTableSize = 0;

  int *localCount = NULL, *globalCount = NULL;
  int *parts=NULL, *nbors_part=NULL, *part_check=NULL;
  int *edges_per_obj=NULL, *nbors_proc=NULL;
  int *num_boundary=NULL, *cuts=NULL;
  int *partNbors= NULL, *partCount=NULL;
  int *key=NULL;

  float obj_edge_weights;

  float *vwgts=NULL, *ewgts=NULL, *wgt=NULL;
  float *globalVals = NULL;
  float *cut_wgt=NULL;
  float *part_sizes = NULL;

  int partPair[2];
  intptr_t keyValue, dummyValue = 0;

  ZOLTAN_GRAPH_EVAL localEval;

  ZOLTAN_TRACE_ENTER(zz, yo);
  ierr = ZOLTAN_OK;

  if (!graph)
    graph = &localEval;

  memset(graph, 0, sizeof(ZOLTAN_GRAPH_EVAL));

  if ((zz->Get_Num_Edges == NULL && zz->Get_Num_Edges_Multi == NULL) ||
           (zz->Get_Edge_List == NULL && zz->Get_Edge_List_Multi == NULL)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "This function requires caller-defined graph query functions.\n");
    return ZOLTAN_FATAL;
  }

  /* Get requested number of parts.  Actual number may differ  */

  ierr = Zoltan_LB_Build_PartDist(zz);
  if (ierr != ZOLTAN_OK){
    goto End;
  }

  req_nparts = zz->LB.Num_Global_Parts;

  /* Get object weights and parts */

  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids,
                             vwgt_dim, &vwgts, &parts);

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

  orig_vwgt_dim = vwgt_dim;

  ierr = add_graph_extra_weight(zz, num_obj, edges_per_obj, &vwgt_dim, &vwgts);

  if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
    goto End;
  }

  /*****************************************************************
   * Get the user defined partition sizes for each weight.  Create
   * partition sizes for the additional (ADD_OBJ_WEIGHT) weight, if any.  
   */

  part_dim = ((orig_vwgt_dim > 0) ? orig_vwgt_dim : 1);

  eval_vwgt_dim = ((vwgt_dim > 0) ? vwgt_dim : 1);

  part_sizes = (float*)ZOLTAN_MALLOC(sizeof(float)*req_nparts*eval_vwgt_dim);
  if (req_nparts && !part_sizes){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  Zoltan_LB_Get_Part_Sizes(zz, part_dim, part_sizes);

  if (eval_vwgt_dim > part_dim){
    for (i=req_nparts-1; i >= 0; i--){
      toidx = i * vwgt_dim;
      fromidx = i * part_dim;
      for (j=0; j < eval_vwgt_dim -1; j++){
        part_sizes[toidx + j] = part_sizes[fromidx + j];
      }

      /* for the library-added weight, set the partition sizes equal to the
       * partition sizes associated with the first user-supplied weight.
       */

      part_sizes[toidx + j] = part_sizes[fromidx];
    }
  }

  /*****************************************************************
   * Local stats 
   */

  graph->nobj[EVAL_LOCAL_SUM] = num_obj;

  if (vwgt_dim > 0){
    for (i=0; i < num_obj; i++){
      graph->obj_wgt[EVAL_LOCAL_SUM]  += vwgts[i*vwgt_dim];
      for (j=1; j <= EVAL_MAX_XTRA_VWGTS; j++){
        if (j == vwgt_dim) break;
        graph->xtra_obj_wgt[j-1][EVAL_LOCAL_SUM]  += vwgts[i*vwgt_dim + j];
      }
    }
  }

  /*****************************************************************
   * Get metrics based on number of objects and object weights 
   */

  ierr = object_metrics(zz, num_obj, parts, part_sizes, req_nparts, 
          vwgts, vwgt_dim,
          &nparts,          /* actual number of parts */
          &nonempty_nparts,  /* number of non-empty parts */
          &graph->obj_imbalance,
          &graph->imbalance,
          graph->nobj,
          graph->obj_wgt,
          graph->xtra_imbalance,
          graph->xtra_obj_wgt);

  if (ierr != ZOLTAN_OK)
    goto End;

  ZOLTAN_FREE(&part_sizes);

  /*****************************************************************
   * Compute the part number of neighboring objects
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

  if (num_edges){

    hashTableSize = Zoltan_Recommended_Hash_Size(num_edges);

    /*  
     * For calculation of each part's number of neighbors,
     * initialize a set of part number pairs using Zoltan_Hash.
     * Alternative is a nparts*nparts array, which uses too much memory.
     */

    map = Zoltan_Map_Create(zz,
                   hashTableSize,     /* size of hash table            */
                   2 * sizeof(int),   /* size of key */
                   1,             /* yes, store a copy of the key */
                   0);            /* dynamically allocate hash table entries */

    if (map == NULL){
      ierr = ZOLTAN_FATAL;
      goto End;
    }
  }

  /*****************************************************************
   * Compute cut statisticics
   */

  cuts = (int *)ZOLTAN_CALLOC(nparts, sizeof(int));
  num_boundary = (int *)ZOLTAN_CALLOC(nparts, sizeof(int));
  part_check = (int *) ZOLTAN_CALLOC(nparts, sizeof(int));

  if (nparts && (!cuts || !num_boundary || !part_check)){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  cut_wgt = (float *)ZOLTAN_CALLOC(nparts * ewgt_dim, sizeof(float));

  if (nparts && ewgt_dim && !cut_wgt){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  num_parts = 0;

  for (i=0,k=0; i < num_obj; i++){   /* object */

    obj_edge_weights = 0;
    obj_part = parts[i];
    nother_parts= 0;

    for (j=0; j < edges_per_obj[i]; j++,k++){    /* neighbor in graph */

      nbor_part = nbors_part[k];

      if (ewgt_dim > 0){
        /* "hypergraph" edge weight is max of the relevant graph edge weights */
        if (ewgts[k * ewgt_dim] > obj_edge_weights)
           obj_edge_weights = ewgts[k * ewgt_dim];  
      }
      else{
        obj_edge_weights = 1.0;
      }

      if (nbor_part != obj_part){

        /* 
         * number of edges that have nbor in a different part 
         */

        cuts[obj_part]++; 

        for (e=0; e < ewgt_dim; e++){
          /*
           * For each part, the sum of the weights of the edges
           * whos neighbor is in a different part
           */
          cut_wgt[obj_part * ewgt_dim + e] += ewgts[k * ewgt_dim + e];
        }

        if (part_check[nbor_part] < i+1){

          nother_parts++;
          part_check[nbor_part] = i + 1;

          partPair[0] = obj_part;
          partPair[1] = nbor_part;

          ierr = Zoltan_Map_Add(zz, map, (char *)partPair, dummyValue);
          if (ierr != ZOLTAN_OK){
            goto End;
          }
        } /* first time neighbor part is seen for this object */
      } /* if neighbor's part is different than object's part */
    } /* next neighbor */

    if (nother_parts){
      /*
       * for each part, the number of objects with a neighbor outside
       * the part
       */
      num_boundary[obj_part]++;
    }
  } /* next object */

  ZOLTAN_FREE(&part_check);
  ZOLTAN_FREE(&parts);
  ZOLTAN_FREE(&edges_per_obj);
  ZOLTAN_FREE(&ewgts);
  ZOLTAN_FREE(&nbors_part);

  /*
   * List of each local part, followed by the number of neighbors of that part
   */
  partCount = (int *)ZOLTAN_CALLOC(2 * nparts, sizeof(int));
  if (!partCount){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  count = 0;   /* number of parts listed in partCount */

  if (num_edges){

    num_pairs = Zoltan_Map_Size(zz, map);
    num_parts = 0;

    if (num_pairs > 0){

      partNbors = (int *)ZOLTAN_MALLOC(sizeof(int) * num_pairs * 2);
      if (!partNbors){
        ierr = ZOLTAN_MEMERR; 
        goto End;
      }

      /* Zoltan_Map "iterator */

      ierr = Zoltan_Map_First(zz, map, (char **) (void *) &key, &keyValue);

      if ( ((ierr == ZOLTAN_OK) && !key) ||  /* must be at least one pair */
           (ierr != ZOLTAN_OK)){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_Map_First\n");
        goto End;
      }

      while (key) {
        partNbors[num_parts++] = key[0];
        partNbors[num_parts++] = key[1];

        ierr = Zoltan_Map_Next(zz, map, (char **) (void *) &key, &keyValue);
  
        if (ierr != ZOLTAN_OK){
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_Map_Next\n");
          goto End;
        }
      } 

      if (num_parts != num_pairs * 2){
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in Zoltan_Map_First/Next\n");
        ierr = ZOLTAN_FATAL;
        goto End;
      }
    }

    Zoltan_Map_Destroy(zz, &map);

    if (num_pairs > 0){
      qsort(partNbors, num_pairs, sizeof(int) * 2,
            zoltan_lb_eval_sort_increasing);

      obj_part = -1;
      count = -1;
      for (i=0; i < num_pairs; i++){

        if (partNbors[2*i] != obj_part){
          obj_part = partNbors[2*i];
          partCount[++count] = obj_part;  /* object part */
          partCount[++count] = 0;       /* number of neighbor parts */
        }
        partCount[count]++;
      }

      ZOLTAN_FREE(&partNbors);
      count++;
      count >>= 1;
    } 
  }

  /************************************************************************
   * Write cut statistics to the return structure.
   */

  k = ((ewgt_dim > 0) ? ewgt_dim : 1);

  globalVals = (float *)ZOLTAN_MALLOC(nparts * k * sizeof(float));
  if (nparts && !globalVals){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  globalCount = (int *)ZOLTAN_MALLOC(nparts * sizeof(int));
  if (nparts && !globalCount){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  /*
   * NNBORPARTS - The number of neighboring parts
   */

  localCount = (int *)ZOLTAN_CALLOC(nparts , sizeof(int));

  for (i=0; i < count; i++){
    localCount[partCount[2*i]] = partCount[2*i + 1];
  }

  ZOLTAN_FREE(&partCount);

  /* Note: this is incorrect if parts are split across processes, as
   * they could be if migration has not yet occured.
   */

  MPI_Allreduce(localCount, globalCount, nparts, MPI_INT, MPI_SUM, comm);

  ZOLTAN_FREE(&localCount);

  iget_strided_stats(globalCount, 1, 0, nparts,
               graph->nnborparts + EVAL_GLOBAL_MIN,
               graph->nnborparts + EVAL_GLOBAL_MAX,
               graph->nnborparts + EVAL_GLOBAL_SUM);

  graph->nnborparts[EVAL_GLOBAL_AVG] = 
                    graph->nnborparts[EVAL_GLOBAL_SUM] / nparts;

  /*
   * CUTS - Number of cut edges in each part
   */

  MPI_Allreduce(cuts, globalCount, nparts, MPI_INT, MPI_SUM, comm);

  ZOLTAN_FREE(&cuts);

  iget_strided_stats(globalCount, 1, 0, nparts,
               graph->cuts + EVAL_GLOBAL_MIN,
               graph->cuts + EVAL_GLOBAL_MAX,
               graph->cuts + EVAL_GLOBAL_SUM);

  graph->cuts[EVAL_GLOBAL_AVG] = graph->cuts[EVAL_GLOBAL_SUM] / nparts;

  /*
   * CUT WEIGHT - The sum of the weights of the cut edges.
   */

  if (ewgt_dim) {
    num_weights = nparts * ewgt_dim;

    MPI_Allreduce(cut_wgt, globalVals, num_weights, MPI_FLOAT, MPI_SUM, comm);

    ZOLTAN_FREE(&cut_wgt);

    fget_strided_stats(globalVals, ewgt_dim, 0, num_weights,
                 graph->cut_wgt + EVAL_GLOBAL_MIN,
                 graph->cut_wgt + EVAL_GLOBAL_MAX,
                 graph->cut_wgt + EVAL_GLOBAL_SUM);

    graph->cut_wgt[EVAL_GLOBAL_AVG] = graph->cut_wgt[EVAL_GLOBAL_SUM] / nparts;
  }

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
   * Number of objects in part that have an off-part neighbor.
   */

  MPI_Allreduce(num_boundary, globalCount, nparts, MPI_INT, MPI_SUM, comm);

  iget_strided_stats(globalCount, 1, 0, nparts,
               graph->num_boundary + EVAL_GLOBAL_MIN,
               graph->num_boundary + EVAL_GLOBAL_MAX,
               graph->num_boundary + EVAL_GLOBAL_SUM);

  graph->num_boundary[EVAL_GLOBAL_AVG] = 
                      graph->num_boundary[EVAL_GLOBAL_SUM] / nparts;

  ZOLTAN_FREE(&num_boundary);
  ZOLTAN_FREE(&globalCount);

  /************************************************************************
   * Print results
   */

  if (print_stats && (zz->Proc == zz->Debug_Proc)){

    printf("\n%s  Part count: %1d requested, %1d actual, %1d non-empty\n", 
      yo, req_nparts, nparts, nonempty_nparts);

    printf("%s  Statistics with respect to %1d parts: \n", yo, nparts);
    printf("%s                             Min      Max      Sum  Imbalance\n",
           yo);

    printf("%s  Number of objects  :  %8.3g %8.3g %8.3g", yo, 
           graph->nobj[EVAL_GLOBAL_MIN], graph->nobj[EVAL_GLOBAL_MAX],
           graph->nobj[EVAL_GLOBAL_SUM]);

    if (graph->obj_imbalance >= 0){
      printf("    %5.3g\n", graph->obj_imbalance);
    }
    else{
      printf("    ----\n");
    }

    if (vwgt_dim > 0){
      printf("%s  Object weight      :  %8.3g %8.3g %8.3g", yo, 
             graph->obj_wgt[EVAL_GLOBAL_MIN], graph->obj_wgt[EVAL_GLOBAL_MAX], 
             graph->obj_wgt[EVAL_GLOBAL_SUM]);

      if (graph->imbalance >= 0){
        printf("     %5.3f\n", graph->imbalance);
      }
      else{
        printf("     ----\n");
      }
  
      for (i=0; i < vwgt_dim-1; i++){
        if (i == EVAL_MAX_XTRA_VWGTS){
          break;
        }
        printf("%s  Object weight %d    :  %8.3g %8.3g %8.3g", yo, i+2,
               graph->xtra_obj_wgt[i][EVAL_GLOBAL_MIN], 
               graph->xtra_obj_wgt[i][EVAL_GLOBAL_MAX], 
               graph->xtra_obj_wgt[i][EVAL_GLOBAL_SUM]);

        if (graph->xtra_imbalance[i] >= 0){
          printf("     %5.3f\n", graph->xtra_imbalance[i]);
        }
        else{
          printf("     ----\n");
        }
      }
      if (vwgt_dim-1 > EVAL_MAX_XTRA_VWGTS){
        printf("(We calculate up to %d extra object weights.  "
               "This can be changed.)\n",
              EVAL_MAX_XTRA_VWGTS);
      }
    }

    printf("\n");

    printf("%s  Statistics with respect to %1d parts: \n", yo, nparts);
    printf("%s                                    "
           "Min      Max    Average    Sum\n", yo);

    printf("%s  Num boundary objects      :  %8.3g %8.3g %8.3g %8.3g\n", yo, 
           graph->num_boundary[EVAL_GLOBAL_MIN], 
           graph->num_boundary[EVAL_GLOBAL_MAX], 
           graph->num_boundary[EVAL_GLOBAL_AVG], 
           graph->num_boundary[EVAL_GLOBAL_SUM]);

    printf("%s  Number of cut edges       :  %8.3g %8.3g %8.3g %8.3g\n", yo, 
           graph->cuts[EVAL_GLOBAL_MIN], graph->cuts[EVAL_GLOBAL_MAX], 
           graph->cuts[EVAL_GLOBAL_AVG],
           graph->cuts[EVAL_GLOBAL_SUM]);

    if (ewgt_dim)
      printf("%s  Weight of cut edges (CUTE):  %8.3g %8.3g %8.3g %8.3g\n", yo, 
             graph->cut_wgt[EVAL_GLOBAL_MIN], graph->cut_wgt[EVAL_GLOBAL_MAX], 
             graph->cut_wgt[EVAL_GLOBAL_AVG], graph->cut_wgt[EVAL_GLOBAL_SUM]);

    for (i=0; i < ewgt_dim-1; i++){
      if (i == EVAL_MAX_XTRA_EWGTS){
        break;
      }
      printf("%s  Weight %d                 :  %8.3g %8.3g %8.3g %8.3g\n", 
             yo, i+2,
             graph->xtra_cut_wgt[i][EVAL_GLOBAL_MIN], 
             graph->xtra_cut_wgt[i][EVAL_GLOBAL_MAX], 
             graph->xtra_cut_wgt[i][EVAL_GLOBAL_AVG], 
             graph->xtra_cut_wgt[i][EVAL_GLOBAL_SUM]);
    }
    if (ewgt_dim-1 > EVAL_MAX_XTRA_EWGTS){
      printf("(We calculate up to %d extra edge weights.  "
             "This can be changed.)\n",
             EVAL_MAX_XTRA_EWGTS);
    }

    printf("%s  Num Nbor Parts            :  %8.3g %8.3g %8.3g %8.3g\n", yo, 
           graph->nnborparts[EVAL_GLOBAL_MIN], 
           graph->nnborparts[EVAL_GLOBAL_MAX], 
           graph->nnborparts[EVAL_GLOBAL_AVG], 
           graph->nnborparts[EVAL_GLOBAL_SUM]);

    printf("\n\n");
  }

End:

  /* Free data */

  Zoltan_Map_Destroy(zz, &map);

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
  ZOLTAN_FREE(&cuts);
  ZOLTAN_FREE(&partCount);
  ZOLTAN_FREE(&part_check);
  ZOLTAN_FREE(&part_sizes);

  ZOLTAN_TRACE_EXIT(zz, yo);
  Zoltan_Destroy(&zz);
  return ierr;
}
/*****************************************************************************/

int Zoltan_LB_Eval_HG(ZZ *zzin, int print_stats, ZOLTAN_HG_EVAL *hg)
{
  /****************************************************************************/
  /* Return performance metrics in ZOLTAN_HG_EVAL structure.  Also print them */
  /* if print_stats is true.  Results are per part, not per process.      */
  /****************************************************************************/

  char *yo = "Zoltan_LB_Eval_HG";
  ZZ *zz = Zoltan_Copy(zzin);  /* Some operations have side-effects in zz;
                                  Don't change the user's settings */
  MPI_Comm comm = zz->Communicator;

  float *part_sizes=NULL;
  float *localVals=NULL;

  double local[2], global[2];

  int ierr, debug_level, i;
  int nparts, nonempty_nparts, req_nparts;
  int vwgt_dim = zz->Obj_Weight_Dim;
  int part_dim = (vwgt_dim > 0 ? vwgt_dim : 1);

  int *localCount = NULL;
  PHGPartParams hgp;

  ZHG* zhg=NULL;

  ZOLTAN_HG_EVAL localEval;
  
  ZOLTAN_TRACE_ENTER(zz, yo);
  /* Set default error code */
  ierr = ZOLTAN_OK;

  if (!hg)
    hg = &localEval;

  memset(hg, 0, sizeof(ZOLTAN_HG_EVAL));

  if ((zz->Get_HG_Size_CS==NULL) && (zz->Get_HG_CS==NULL) &&
      (zz->Get_Num_Edges==NULL) && (zz->Get_Num_Edges_Multi==NULL) &&
      (zz->Get_Edge_List==NULL) && (zz->Get_Edge_List_Multi==NULL)) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
      "This function requires graph or hypergraph query functions.\n");
    return ZOLTAN_FATAL;
  }

  /* get the requested number of parts */ 

  ierr = Zoltan_LB_Build_PartDist(zz);
  if (ierr != ZOLTAN_OK){
    goto End;
  }

  req_nparts = zz->LB.Num_Global_Parts;

  /*****************************************************************
   * Get the hypergraph, via the hypergraph or graph query functions
   */

  part_sizes = (float*)ZOLTAN_MALLOC(sizeof(float) * part_dim * req_nparts);
  if (req_nparts && !part_sizes){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  Zoltan_LB_Get_Part_Sizes(zz, part_dim, part_sizes);

  debug_level = zz->Debug_Level;
  zz->Debug_Level = 0;

  ierr = Zoltan_PHG_Initialize_Params(zz, part_sizes, HYPERGRAPH, &hgp);
  if (ierr != ZOLTAN_OK)
    goto End;

  zz->Debug_Level = debug_level;

  zhg = (ZHG*) ZOLTAN_MALLOC (sizeof(ZHG));
  if (zhg == NULL){
    ierr = ZOLTAN_MEMERR;
    goto End;
  }

  Zoltan_Input_HG_Init(zhg);

  ierr = Zoltan_Get_Hypergraph_From_Queries(zz, &hgp, HYPERGRAPH, zhg);
  if (ierr != ZOLTAN_OK)
    goto End;

  if (zhg->globalObj == 0){
    if (zz->Proc == zz->Debug_Proc){
      printf("%s: No objects to evaluate\n",yo);
    }
    goto End;
  }

  /* Get metrics based on number of objects and object weights */

  ierr = object_metrics(zz, zhg->nObj, zhg->Input_Parts, part_sizes, req_nparts,
          zhg->objWeight, zhg->objWeightDim,
          &nparts,          /* actual number of parts */
          &nonempty_nparts,  /* number of non-empty parts */
          &hg->obj_imbalance,
          &hg->imbalance,
          hg->nobj,
          hg->obj_wgt,
          hg->xtra_imbalance,
          hg->xtra_obj_wgt);

  if (ierr != ZOLTAN_OK)
    goto End;

  /*****************************************************************
   * Local stats 
   */

  hg->nobj[EVAL_LOCAL_SUM] = zhg->nObj;

  if (zhg->objWeightDim > 0){
    for (i=0; i < zhg->nObj; i++){
      hg->obj_wgt[EVAL_LOCAL_SUM]  += zhg->objWeight[i];
    }
  }

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
  hg->cutl[EVAL_GLOBAL_SUM] = (float)global[0];
  hg->cutn[EVAL_GLOBAL_SUM] = (float)global[1];

  MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MIN, comm);
  hg->cutl[EVAL_GLOBAL_MIN] = (float)global[0];
  hg->cutn[EVAL_GLOBAL_MIN] = (float)global[1];

  MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, comm);
  hg->cutl[EVAL_GLOBAL_MAX] = (float)global[0];
  hg->cutn[EVAL_GLOBAL_MAX] = (float)global[1];

  hg->cutl[EVAL_GLOBAL_AVG] = hg->cutl[EVAL_GLOBAL_SUM] / nparts;
  hg->cutn[EVAL_GLOBAL_AVG] = hg->cutn[EVAL_GLOBAL_SUM] / nparts;
            
  /************************************************************************
   * Print results
   */

  if (print_stats && (zz->Proc == zz->Debug_Proc)){

    printf("\n%s  Part count: %1d requested, %1d actual, %1d non-empty\n", 
      yo, req_nparts, nparts, nonempty_nparts);

    printf("%s  Statistics with respect to %1d parts: \n", yo, nparts);
    printf("%s                            Min      Max     Sum  Imbalance\n", 
           yo);

    printf("%s  Number of objects :  %8.3g %8.3g %8.3g", yo, 
      hg->nobj[EVAL_GLOBAL_MIN], hg->nobj[EVAL_GLOBAL_MAX], 
      hg->nobj[EVAL_GLOBAL_SUM]);

    if (hg->obj_imbalance >= 0){
      printf("   %5.3f\n", hg->obj_imbalance);
    }
    else{
      printf("   ----\n");
    }

    if (zhg->objWeightDim > 0){
      printf("%s  Object weight     :  %8.3g %8.3g %8.3g", yo, 
        hg->obj_wgt[EVAL_GLOBAL_MIN], hg->obj_wgt[EVAL_GLOBAL_MAX], 
        hg->obj_wgt[EVAL_GLOBAL_SUM]);

      if (hg->imbalance >= 0){
        printf("   %5.3f\n", hg->imbalance);
      }
      else{
        printf("   ----\n");
      }

      for (i=0; i < zhg->objWeightDim-1; i++){
        if (i == EVAL_MAX_XTRA_VWGTS){
          break;
        }
        printf("%s  Object weight %d   :  %8.3g %8.3g %8.3g", 
               yo, i+2,
               hg->xtra_obj_wgt[i][EVAL_GLOBAL_MIN],
               hg->xtra_obj_wgt[i][EVAL_GLOBAL_MAX],
               hg->xtra_obj_wgt[i][EVAL_GLOBAL_SUM]);

        if (hg->xtra_imbalance[i] >= 0){
          printf("   %5.3f\n", hg->xtra_imbalance[i]);
        }
        else{
          printf("   ----\n");
        }
      }
      if (zhg->objWeightDim-1 > EVAL_MAX_XTRA_VWGTS){
        printf("(We calculate up to %d extra object weights.  "
               "This can be changed.)\n",
               EVAL_MAX_XTRA_VWGTS);
      }
    }
    printf("\n");

    printf("%s  CUTN (Sum_edges( (#parts(edge)>1)*ewgt )): %8.3f\n", 
           yo, hg->cutn[EVAL_GLOBAL_SUM]);
    printf("%s  CUTL (Sum_edges( (#parts(edge)-1)*ewgt )): %8.3f\n", 
           yo, hg->cutl[EVAL_GLOBAL_SUM]);
    printf("%s  CUTL-MAX (Max_procs( comm. volume ):       %8.3f\n", 
           yo, hg->cutl[EVAL_GLOBAL_MAX]);


    printf("\n\n");
  }

End:

  /* Free data */

  if (hgp.globalcomm.row_comm != MPI_COMM_NULL)
    MPI_Comm_free(&(hgp.globalcomm.row_comm));
  if (hgp.globalcomm.col_comm != MPI_COMM_NULL)
    MPI_Comm_free(&(hgp.globalcomm.col_comm));
  if (hgp.globalcomm.Communicator != MPI_COMM_NULL)
    MPI_Comm_free(&(hgp.globalcomm.Communicator));

  ZOLTAN_FREE(&part_sizes);
  if (zhg){
    Zoltan_PHG_Free_Hypergraph_Data(zhg);
    ZOLTAN_FREE(&zhg);
  }
  ZOLTAN_FREE(&localVals);
  ZOLTAN_FREE(&localCount);

  ZOLTAN_TRACE_EXIT(zz, yo);
  Zoltan_Destroy(&zz);
  return ierr;
}

/************************************************************************/
/************************************************************************/
/************************************************************************/
int Zoltan_LB_Eval(ZZ *zz, int print_stats, 
                   ZOLTAN_BALANCE_EVAL *obj, 
                   ZOLTAN_GRAPH_EVAL *graph, 
                   ZOLTAN_HG_EVAL *hg
)
/*
 * Input:
 *   zzin        - pointer to Zoltan structure
 *   print_stats - if > 0, compute and print partition quality metrics
 *                 if == 0, stay silent but update EVAL structures with metrics
 *
 * Input/Output:
 *   obj         - pointer to a ZOLTAN_BALANCE_EVAL structure,
 *                 if non-NULL partitioning
 *                 metrics will be written to the structure
 *   graph       - pointer to a ZOLTAN_GRAPH_EVAL structure, 
 *                 if non_null and graph query
 *                 functions are defined by the application, graph partitioning
 *                 quality metrics will be written to the structure
 *   hg          - pointer to an ZOLTAN_HG_EVAL structure, 
 *                 if non_null and graph or
 *                 hypergraph query functions are defined by the application,
 *                 hypergraph partitioning quality metrics will be written to 
 *                 the structure
 *
 */
{
  char *yo = "Zoltan_LB_Eval";
  int ierr = ZOLTAN_OK;
  int hypergraph_callbacks = 0;
  int graph_callbacks = 0;

  if (!print_stats && !obj && !graph && !hg){
    return ierr;
  }

  if (zz->Get_HG_Size_CS && zz->Get_HG_CS){
    hypergraph_callbacks = 1;
  }
  if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
           (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) {
    graph_callbacks = 1;
  }

  if (print_stats || obj){
    ierr = Zoltan_LB_Eval_Balance(zz, print_stats, obj);
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Zoltan_LB_Eval_Balance");
      goto End;
    }
  }

  if ((print_stats || graph) && graph_callbacks){
    ierr = Zoltan_LB_Eval_Graph(zz, print_stats, graph);
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from Zoltan_LB_Eval_Graph");
      goto End;
    }
  }

  if ((print_stats || hg) && (graph_callbacks || hypergraph_callbacks)){
    ierr = Zoltan_LB_Eval_HG(zz, print_stats, hg);
    if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                         "Error returned from Zoltan_LB_Eval_HG");
      goto End;
    }
  }

End:
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
/* Function to retrieve the part number for neighboring nodes. */
char *yo = "get_nbor_parts";
struct Zoltan_DD_Struct *dd = NULL;
int *owner = NULL;
int maxnobj;
int ierr;
int i, start, size, i_am_done, alldone;
const int MAXSIZE = 200000;

  for (i=0; i < nnbors; i++){
    /* 
     * a check on validity of data supplied by query functions
     */
    nbors_part[i] = -1;
  }

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

  for (i=0, size=0; i < nnbors; i++){
    if (nbors_part[i] < 0){
      if (size == 10){
        fprintf(stderr,
                "%s (%d) more uninitialized entries omitted from print out\n",
                yo, zz->Proc);
        break;
      }
      else{
        fprintf(stderr,
                "%s (%d) ERROR part array index %d is uninitialized\n",
                yo, zz->Proc,i);
        size++;
      }
    }
  }
  if (size){
    fprintf(stderr,
            "%s (%d) Most likely cause is incorrect edge data from app\n",
            yo, zz->Proc);
    ierr = ZOLTAN_FATAL;
  }

End:
  Zoltan_DD_Destroy(&dd);
  return ierr;
}

/************************************************************************/

static int *
objects_by_part(ZZ *zz, int num_obj, int *part, int *nparts, int *nonempty)
{
  int i, num_parts, num_nonempty, max_part, gmax_part;
  int *partCounts = NULL, *totalCounts;
  int *returnBuf = NULL;

  max_part = 0;
  for (i=0; i < num_obj; i++){
    if (part[i] > max_part) max_part = part[i];
  }

  MPI_Allreduce(&max_part, &gmax_part, 1, MPI_INT, MPI_MAX, zz->Communicator);

  max_part = gmax_part + 1;

  /* Allocate and return a buffer containing the local count of objects by part,
     followed by the global count.  Set "nparts" to the number of parts, and
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

  return returnBuf;
}

/************************************************************************/

static int
object_metrics(ZZ *zz, int num_obj, int *parts, float *part_sizes, int req_nparts, 
               float *vwgts, int vwgt_dim,
               int *nparts,       /* return actual number of parts */
               int *nonempty,     /* return number of those that are non-empty*/
               float *obj_imbalance,  /* object # imbalance wrt parts */
               float *imbalance,      /* object wgt imbalance wrt parts */
               float *nobj,       /* number of objects */
               float *obj_wgt,    /* object weights */
               float *xtra_imbalance,  /* return if vertex weight dim > 1 */
    float (*xtra_obj_wgt)[EVAL_SIZE])  /* return if vertex weight dim > 1 */
{
  MPI_Comm comm = zz->Communicator;

  int i, j, idx, ierr, part_dim;
  int num_weights; 
  int num_parts = zz->LB.Num_Global_Parts; 
  int num_nonempty_parts = zz->LB.Num_Global_Parts;

  int *globalCount = NULL;

  float *wgt=NULL;
  float *localVals = NULL, *globalVals = NULL;

  float imbal, tmp;

  ierr = ZOLTAN_OK;

  part_dim = (vwgt_dim > 0) ? vwgt_dim : 1;

  /* Get the actual number of parts, and number of objects per part */

  globalCount = objects_by_part(zz, num_obj, parts,
                &num_parts,           /* actual number of parts */
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

  /* Weight object count imbalance by desired partition sizes (for first weight) */

  if (req_nparts >= num_parts){
    imbal = 0.0;
  
    if (nobj[EVAL_GLOBAL_SUM] > 0) {
  
      for (i=0, idx=0; i < num_parts; i++, idx+=part_dim){
        if (part_sizes[idx] > 0){
          tmp = globalCount[i] / (nobj[EVAL_GLOBAL_SUM] * part_sizes[idx]);
          if (tmp > imbal) imbal = tmp;
        }
      }
    }
    *obj_imbalance = (imbal > 0 ? imbal : 1.0);
  }
  else{
    /* 
     * flag that imbalance is infinite (part_size is 0 for non-empty part) 
     */
    *obj_imbalance = -1;
  }

  ZOLTAN_FREE(&globalCount);

  if (vwgt_dim > 0){

    /* if vwgt_dim > 0, then part_dim is the same as vwgt_dim */

    num_weights = num_parts * vwgt_dim;
  
    localVals = (float *)ZOLTAN_CALLOC(num_weights*2, sizeof(float));
  
    if (num_weights && !localVals){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
  
    globalVals = localVals + num_weights;
  
    for (i=0; i<num_obj; i++){
      for (j=0; j<vwgt_dim; j++){
        localVals[parts[i]*vwgt_dim+j] += vwgts[i*vwgt_dim+j];
      }
    }
  
    MPI_Allreduce(localVals, globalVals, num_weights, MPI_FLOAT, MPI_SUM, comm);
  
    fget_strided_stats(globalVals, vwgt_dim, 0, num_weights,
                       obj_wgt + EVAL_GLOBAL_MIN,
                       obj_wgt + EVAL_GLOBAL_MAX,
                       obj_wgt + EVAL_GLOBAL_SUM);
  
    obj_wgt[EVAL_GLOBAL_AVG] = obj_wgt[EVAL_GLOBAL_SUM]/(float)num_parts;

    /* imbalance scaled by user specified part_sizes */

    if (req_nparts >= num_parts){
      imbal = 0.0;
    
      if (obj_wgt[EVAL_GLOBAL_SUM] > 0){
    
        for (i=0, idx = 0; i < num_parts; i++, idx += vwgt_dim){
          if (part_sizes[idx] > 0){
            tmp = globalVals[idx] / (obj_wgt[EVAL_GLOBAL_SUM] * part_sizes[idx]);
            if (tmp > imbal) imbal = tmp;
          }
        }
      }
     
      *imbalance = (imbal > 0 ? imbal : 1.0);
    }
    else{
      *imbalance = -1;  /* flag some part_sizes are zero */
    }

    /* calculations for multiple vertex weights */
  
    for (i=0; i < vwgt_dim-1; i++){

      if (i == EVAL_MAX_XTRA_VWGTS){
        break;
      }
     
      wgt = xtra_obj_wgt[i];
  
      fget_strided_stats(globalVals, vwgt_dim, i+1, num_weights,
                         wgt + EVAL_GLOBAL_MIN,
                         wgt + EVAL_GLOBAL_MAX,
                         wgt + EVAL_GLOBAL_SUM);
  
      wgt[EVAL_GLOBAL_AVG] = wgt[EVAL_GLOBAL_SUM]/(float)num_parts;
 
      /* imbalance scaled by user specified part_sizes */

      if (req_nparts >= num_parts){
        imbal = 0.0;
      
        if (wgt[EVAL_GLOBAL_SUM] > 0){
      
          for (j=0; j < num_parts; j++){
            idx = (j * vwgt_dim) + i + 1;
            if (part_sizes[idx] > 0){
              tmp = globalVals[idx] / (wgt[EVAL_GLOBAL_SUM] * part_sizes[idx]);
              if (tmp > imbal) imbal = tmp;
            }
          }
        }
        xtra_imbalance[i] = (imbal > 0 ? imbal : 1.0);
      }
      else{
        xtra_imbalance[i] = -1;  /* flag some part_sizes are 0 */
      }
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

  if (!strcasecmp(add_obj_weight, "NONE")){
    return ierr;
  }
  else if ((!strcasecmp(add_obj_weight, "UNIT")) || (!strcasecmp(add_obj_weight, "VERTICES"))){
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

/************************************************************************/

int zoltan_lb_eval_sort_increasing(const void *a, const void *b)
{
  const int *val_a = (const int *)a;
  const int *val_b = (const int *)b;

  if (*val_a < *val_b){
    return -1;
  }
  else if (*val_a > *val_b){
    return 1;
  }
  else{
    return 0;
  }
}

/************************************************************************/

void Zoltan_LB_Eval_Print_Balance(ZOLTAN_BALANCE_EVAL *lb)
{
  int i;

  printf("               Minimum     Maximum      Average      Sum         Sum\n");
  printf("                across      across        of          of       on local\n");
  printf("                parts       parts        parts       parts     process\n");
  printf("               ========    ========    ========    ========    ========\n");

  printf("num objects %11.4f %11.4f %11.4f %11.4f %11.4f\n",
    lb->nobj[EVAL_GLOBAL_MIN], lb->nobj[EVAL_GLOBAL_MAX], lb->nobj[EVAL_GLOBAL_AVG],
    lb->nobj[EVAL_GLOBAL_SUM], lb->nobj[EVAL_LOCAL_SUM]);

  printf("weight objs %11.4f %11.4f %11.4f %11.4f %11.4f\n",
    lb->obj_wgt[EVAL_GLOBAL_MIN], lb->obj_wgt[EVAL_GLOBAL_MAX], lb->obj_wgt[EVAL_GLOBAL_AVG],
    lb->obj_wgt[EVAL_GLOBAL_SUM], lb->obj_wgt[EVAL_LOCAL_SUM]);

  for (i=0; i < EVAL_MAX_XTRA_VWGTS; i++){
    if (lb->xtra_obj_wgt[i][EVAL_GLOBAL_SUM] == 0) break;
    printf(" weight #%d %11.4f %11.4f %11.4f %11.4f %11.4f\n", i+2,
      lb->xtra_obj_wgt[i][EVAL_GLOBAL_MIN], lb->xtra_obj_wgt[i][EVAL_GLOBAL_MAX], lb->xtra_obj_wgt[i][EVAL_GLOBAL_AVG],
      lb->xtra_obj_wgt[i][EVAL_GLOBAL_SUM], lb->xtra_obj_wgt[i][EVAL_LOCAL_SUM]);
  }

  printf("object count imbalance     %11.4f\n",lb->obj_imbalance);
  printf("object weight imbalance    %11.4f\n",lb->imbalance);

  for (i=0; i < EVAL_MAX_XTRA_VWGTS; i++){
    if (lb->xtra_imbalance[i] == 0) break;
    printf("  object weight %d         %11.4f\n",i+2,lb->xtra_imbalance[i]);
  }
  printf("\n");
}
/************************************************************************/

void Zoltan_LB_Eval_Print_HG(ZOLTAN_HG_EVAL *hg)
{
  int i;

  printf("               Minimum     Maximum      Average      Sum         Sum\n");
  printf("                across      across        of          of       on local\n");
  printf("                parts       parts        parts       parts     process\n");
  printf("               ========    ========    ========    ========    ========\n");

  printf("num vtxs    %11.4f %11.4f %11.4f %11.4f %11.4f\n",
    hg->nobj[EVAL_GLOBAL_MIN], hg->nobj[EVAL_GLOBAL_MAX], hg->nobj[EVAL_GLOBAL_AVG],
    hg->nobj[EVAL_GLOBAL_SUM], hg->nobj[EVAL_LOCAL_SUM]);

  printf("weight vtxs %11.4f %11.4f %11.4f %11.4f %11.4f\n",
    hg->obj_wgt[EVAL_GLOBAL_MIN], hg->obj_wgt[EVAL_GLOBAL_MAX], hg->obj_wgt[EVAL_GLOBAL_AVG],
    hg->obj_wgt[EVAL_GLOBAL_SUM], hg->obj_wgt[EVAL_LOCAL_SUM]);

  for (i=0; i < EVAL_MAX_XTRA_VWGTS; i++){
    if (hg->xtra_obj_wgt[i][EVAL_GLOBAL_SUM] == 0) break;
      printf("  weight %d  %11.4f %11.4f %11.4f %11.4f %11.4f\n", i+2,
      hg->xtra_obj_wgt[i][EVAL_GLOBAL_MIN], hg->xtra_obj_wgt[i][EVAL_GLOBAL_MAX], hg->xtra_obj_wgt[i][EVAL_GLOBAL_AVG],
      hg->xtra_obj_wgt[i][EVAL_GLOBAL_SUM], hg->xtra_obj_wgt[i][EVAL_LOCAL_SUM]);
  }

  printf(" cutn       %11.4f %11.4f %11.4f %11.4f\n",
    hg->cutn[EVAL_GLOBAL_MIN], hg->cutn[EVAL_GLOBAL_MAX], hg->cutn[EVAL_GLOBAL_AVG],
    hg->cutn[EVAL_GLOBAL_SUM]);

  printf(" cutl       %11.4f %11.4f %11.4f %11.4f\n",
    hg->cutl[EVAL_GLOBAL_MIN], hg->cutl[EVAL_GLOBAL_MAX], hg->cutl[EVAL_GLOBAL_AVG],
    hg->cutl[EVAL_GLOBAL_SUM]);

  printf("vertex number imbalance    %11.4f\n",hg->obj_imbalance);
  printf("vertex weight imbalance    %11.4f\n",hg->imbalance);

  for (i=0; i < EVAL_MAX_XTRA_VWGTS; i++){
    if (hg->xtra_imbalance[i] == 0) break;
    printf("  weight %d               %11.4f\n",i+2,hg->xtra_imbalance[i]);
  }
  printf("\n");
}

/************************************************************************/

void Zoltan_LB_Eval_Print_Graph(ZOLTAN_GRAPH_EVAL *graph)
{
  int i;

  printf("               Minimum     Maximum      Average      Sum         Sum\n");
  printf("                across      across        of          of       on local\n");
  printf("                parts       parts        parts       parts     process\n");
  printf("               ========    ========    ========    ========    ========\n");

  printf("num vtxs    %11.4f %11.4f %11.4f %11.4f %11.4f\n",
    graph->nobj[EVAL_GLOBAL_MIN], graph->nobj[EVAL_GLOBAL_MAX], graph->nobj[EVAL_GLOBAL_AVG],
    graph->nobj[EVAL_GLOBAL_SUM], graph->nobj[EVAL_LOCAL_SUM]);

  printf("weight vtxs %11.4f %11.4f %11.4f %11.4f %11.4f\n",
    graph->obj_wgt[EVAL_GLOBAL_MIN], graph->obj_wgt[EVAL_GLOBAL_MAX], graph->obj_wgt[EVAL_GLOBAL_AVG],
    graph->obj_wgt[EVAL_GLOBAL_SUM], graph->obj_wgt[EVAL_LOCAL_SUM]);

  for (i=0; i < EVAL_MAX_XTRA_VWGTS; i++){
    if (graph->xtra_obj_wgt[i][EVAL_GLOBAL_SUM] == 0) break;
    printf("  weight %d  %11.4f %11.4f %11.4f %11.4f %11.4f\n", i+2,
      graph->xtra_obj_wgt[i][EVAL_GLOBAL_MIN], graph->xtra_obj_wgt[i][EVAL_GLOBAL_MAX], graph->xtra_obj_wgt[i][EVAL_GLOBAL_AVG],
      graph->xtra_obj_wgt[i][EVAL_GLOBAL_SUM], graph->xtra_obj_wgt[i][EVAL_LOCAL_SUM]);
  }

  printf("# bdry vtxs %11.4f %11.4f %11.4f %11.4f %11.4f\n",
    graph->num_boundary[EVAL_GLOBAL_MIN], graph->num_boundary[EVAL_GLOBAL_MAX], graph->num_boundary[EVAL_GLOBAL_AVG],
    graph->num_boundary[EVAL_GLOBAL_SUM], graph->num_boundary[EVAL_LOCAL_SUM]);

  printf(" cuts       %11.4f %11.4f %11.4f %11.4f\n",
    graph->cuts[EVAL_GLOBAL_MIN], graph->cuts[EVAL_GLOBAL_MAX], graph->cuts[EVAL_GLOBAL_AVG],
    graph->cuts[EVAL_GLOBAL_SUM]);

  printf("cut wgt     %11.4f %11.4f %11.4f %11.4f\n",
    graph->cut_wgt[EVAL_GLOBAL_MIN], graph->cut_wgt[EVAL_GLOBAL_MAX], graph->cut_wgt[EVAL_GLOBAL_AVG],
    graph->cut_wgt[EVAL_GLOBAL_SUM]);

  for (i=0; i < EVAL_MAX_XTRA_EWGTS; i++){
    if (graph->xtra_cut_wgt[i][EVAL_GLOBAL_SUM] == 0) break;
      printf("  weight %d  %11.4f %11.4f %11.4f %11.4f\n", i+2,
      graph->xtra_cut_wgt[i][EVAL_GLOBAL_MIN], graph->xtra_cut_wgt[i][EVAL_GLOBAL_MAX], graph->xtra_cut_wgt[i][EVAL_GLOBAL_AVG],
      graph->xtra_cut_wgt[i][EVAL_GLOBAL_SUM]);
  }

  printf("#nbor parts %11.4f %11.4f %11.4f %11.4f\n",
    graph->nnborparts[EVAL_GLOBAL_MIN], graph->nnborparts[EVAL_GLOBAL_MAX], graph->nnborparts[EVAL_GLOBAL_AVG],
    graph->nnborparts[EVAL_GLOBAL_SUM]);

  printf("vertex number imbalance    %11.4f\n",graph->obj_imbalance);
  printf("vertex weight imbalance    %11.4f\n",graph->imbalance);

  for (i=0; i < EVAL_MAX_XTRA_VWGTS; i++){
    if (graph->xtra_imbalance[i] == 0) break;
    printf("  weight %d                 %11.4f\n",i+2,graph->xtra_imbalance[i]);
  }
  printf("\n");
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
