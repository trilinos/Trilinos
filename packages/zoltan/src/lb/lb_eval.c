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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for evaluating the quality of the
 *  current partitioning/balance.
 *  These functions are all callable by the application. 
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#define NUM_GSTATS 4 /* Number of graph statistics */

void LB_Eval (LB *lb, int print_stats, 
     int *nobj, float *obj_wgt, int *cut_wgt, int *nboundary,
     int *nadj, int *ierr)
/* 
 * Input:
 *   lb          - pointer to lb structure
 *   print_stats - if > 0, compute and print max, min and sum of the metrics
 *
 * Output:
 *   nobj      - number of objects (for each proc)
 *   obj_wgt   - obj_wgt[0:lb->Obj_Weight_Dim-1] are the object weights (on 
 *               each proc)
 *   cut_wgt   - cut size/weight (for each proc)
 *   nboundary - number of boundary objects (for each proc)
 *   nadj      - the number of adjacent procs (for each proc)
 *   ierr      - error code
 *
 * Ouput parameters will only be returned if they are 
 * not NULL on entry (except for the error code ierr).
 */

{
  char *yo = "LB_Eval";
  char *yo2 = "ZOLTAN LB_Eval";
  int i, j, num_obj, max_edges, flag, nedges;
  int num_adj, num_boundary, cut_weight;
  int stats[4*NUM_GSTATS], *ewgts;
  int *proc, *nbors_proc;
  float *tmp_wgt, *vwgts, nproc;
  LB_LID * local_ids;
  LB_GID *global_ids, *nbors_global;
  
  LB_TRACE_ENTER(lb, yo);
  /* Set default error code */
  *ierr = LB_OK;

  /* Set all pointers to NULL */
  global_ids = NULL;
  local_ids = NULL;
  tmp_wgt = NULL;
  vwgts = NULL;
  ewgts = NULL;
  nbors_global = NULL;
  nbors_proc = NULL;
  proc = NULL;

  /* First compute number of objs and object weight on each proc */
  num_obj = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, ierr);

  if (num_obj>0){

    /* Allocate space for object data */
    global_ids = (LB_GID *) LB_MALLOC(num_obj * sizeof(LB_GID));
    local_ids  = (LB_LID *) LB_MALLOC(num_obj * sizeof(LB_LID));
      
    if ((!global_ids) || (!local_ids)){
      *ierr = LB_MEMERR;
      LB_FREE(&global_ids);
      LB_FREE(&local_ids);
      LB_TRACE_EXIT(lb, yo);
      return;
    }
  }
  

  /* Allocate space for weights if needed */
  if (lb->Obj_Weight_Dim>0){
    vwgts   = (float  *) LB_MALLOC(lb->Obj_Weight_Dim*num_obj * sizeof(float));
    tmp_wgt = (float *) LB_MALLOC(4*lb->Obj_Weight_Dim * sizeof(float));
    if ((num_obj && !vwgts) || (!tmp_wgt)){
      *ierr = LB_MEMERR;
      LB_FREE(&global_ids);
      LB_FREE(&local_ids);
      LB_FREE(&vwgts);
      LB_FREE(&tmp_wgt);
      LB_TRACE_EXIT(lb, yo);
      return;
    }
  } 
  
  LB_Get_Obj_List(lb, global_ids, local_ids, lb->Obj_Weight_Dim, vwgts, ierr);
  if (*ierr == LB_FATAL){
    LB_FREE(&global_ids);
    LB_FREE(&local_ids);
    LB_TRACE_EXIT(lb, yo);
    return;
  }
  

  /* Compute object weight sums */
  if (lb->Obj_Weight_Dim>0){
    for (j=0; j<lb->Obj_Weight_Dim; j++)
      tmp_wgt[j] = 0;
    for (i=0; i<num_obj; i++){
      for (j=0; j<lb->Obj_Weight_Dim; j++){
        tmp_wgt[j] += vwgts[i*lb->Obj_Weight_Dim+j];
      }
    }
  }


  /* Compute (weighted) edge cuts, #boundary vertices,
     and # adjacent procs if possible */

  cut_weight = 0;
  num_boundary = 0;
  num_adj = 0;

  if (lb->Get_Num_Edges != NULL) {
    /* Use the basic graph query functions */

    /* First compute max no. of edges so we can allocate the right
       amount of space */
    max_edges = 0;
    for (i=0; i< num_obj; i++){
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], ierr);
      if (*ierr){
        fprintf(stderr, "Error in %s: Get_Num_Edges returned error code %d\n", 
          yo, *ierr);
        LB_FREE(&global_ids);
        LB_FREE(&local_ids);
        LB_FREE(&vwgts);
        LB_FREE(&tmp_wgt);
        LB_TRACE_EXIT(lb, yo);
        return;
      }
      if (nedges>max_edges) max_edges = nedges;
    }

    /* Allocate edge list space */
    nbors_global = (LB_GID *)LB_MALLOC(max_edges * sizeof(LB_GID));
    nbors_proc = (int *)LB_MALLOC(max_edges * sizeof(int));
    ewgts = (int *)LB_MALLOC(lb->Comm_Weight_Dim*max_edges * sizeof(int));
    /* Allocate a proc list for computing nadjacent */
    proc = (int *)LB_MALLOC((lb->Num_Proc)* sizeof(int));

    if (max_edges && ((!nbors_global) || (!nbors_proc) || 
        (lb->Comm_Weight_Dim && (!ewgts)) || (!proc))){
      *ierr = LB_MEMERR;
      LB_FREE(&global_ids);
      LB_FREE(&local_ids);
      LB_FREE(&vwgts);
      LB_FREE(&tmp_wgt);
      LB_FREE(&nbors_global);
      LB_FREE(&nbors_proc);
      LB_FREE(&ewgts);
      LB_FREE(&proc);
      LB_TRACE_EXIT(lb, yo);
      return;
    }

    for (i=0; i<lb->Num_Proc; i++)
      proc[i] = 0;

    for (i=0; i<num_obj; i++){
      flag = 0;
      nedges = lb->Get_Num_Edges(lb->Get_Edge_List_Data, global_ids[i], 
               local_ids[i], ierr);
      if (*ierr == LB_FATAL){
        LB_FREE(&global_ids);
        LB_FREE(&local_ids);
        LB_FREE(&vwgts);
        LB_FREE(&tmp_wgt);
        LB_FREE(&nbors_global);
        LB_FREE(&nbors_proc);
        LB_FREE(&ewgts);
        LB_FREE(&proc);
        LB_TRACE_EXIT(lb, yo);
        return;
      }
      lb->Get_Edge_List(lb->Get_Edge_List_Data, global_ids[i], local_ids[i],
          nbors_global, nbors_proc, lb->Comm_Weight_Dim, ewgts, ierr);
      if (*ierr == LB_FATAL){
        LB_FREE(&global_ids);
        LB_FREE(&local_ids);
        LB_FREE(&vwgts);
        LB_FREE(&tmp_wgt);
        LB_FREE(&nbors_global);
        LB_FREE(&nbors_proc);
        LB_FREE(&ewgts);
        LB_FREE(&proc);
        LB_TRACE_EXIT(lb, yo);
        return;
      }
      /* Check for cut edges */
      for (j=0; j<nedges; j++){
        if (nbors_proc[j] != lb->Proc){
          if (lb->Comm_Weight_Dim == 0)
            cut_weight++;
          else if (lb->Comm_Weight_Dim == 1)
            cut_weight += ewgts[j];
          else{
            fprintf(stderr, "Error in %s: Comm_Weight_Dim=%d not supported\n", 
                    yo, lb->Comm_Weight_Dim);
            *ierr = LB_WARN;
          }
          if (flag==0){
            num_boundary++;
            flag = 1;
          }
          proc[nbors_proc[j]]++;
        }
      }
    }
    /* Compute the number of adjacent procs */
    for (j=0; j<lb->Num_Proc; j++)
      if (proc[j]>0) num_adj++;
  }
  else{
    /* No graph query functions available */
  }
  
  if (print_stats){
    /* Global reduction for object weights. */
    if (lb->Obj_Weight_Dim>0){
      MPI_Allreduce(tmp_wgt, &tmp_wgt[lb->Obj_Weight_Dim], lb->Obj_Weight_Dim, 
                    MPI_FLOAT, MPI_MAX, lb->Communicator);
      MPI_Allreduce(tmp_wgt, &tmp_wgt[2*lb->Obj_Weight_Dim], lb->Obj_Weight_Dim,
                    MPI_FLOAT, MPI_SUM, lb->Communicator);
      MPI_Allreduce(tmp_wgt, &tmp_wgt[3*lb->Obj_Weight_Dim], lb->Obj_Weight_Dim,
                    MPI_FLOAT, MPI_MIN, lb->Communicator);
    }
    stats[0] = num_obj;
    stats[1] = cut_weight;
    stats[2] = num_boundary;
    stats[3] = num_adj;

    /* Compute max, min and sum in the upper portions of the stats array. */
    MPI_Allreduce(stats, &stats[NUM_GSTATS], NUM_GSTATS, MPI_INT, MPI_MAX, 
                  lb->Communicator);
    MPI_Allreduce(stats, &stats[2*NUM_GSTATS], NUM_GSTATS, MPI_INT, MPI_SUM, 
                  lb->Communicator);
    MPI_Allreduce(stats, &stats[3*NUM_GSTATS], NUM_GSTATS, MPI_INT, MPI_MIN, 
                  lb->Communicator);

    /* Print max-sum of results */

    /* Print max-sum of results */
    nproc = lb->Num_Proc; /* convert to float */
    if (lb->Proc == lb->Debug_Proc){
      printf("\n%s  Statistics for current partitioning/balance:\n", yo2);
      for (i=0; i<lb->Obj_Weight_Dim; i++)
        printf("%s  Object weight #%1d :  Max = %6.1f, Min = %6.1f, "
          "Sum = %7.1f, Imbal. = %5.3f\n",
          yo2, i+1, tmp_wgt[lb->Obj_Weight_Dim+i], 
          tmp_wgt[3*lb->Obj_Weight_Dim+i], 
          tmp_wgt[2*lb->Obj_Weight_Dim+i], 
          (tmp_wgt[2*lb->Obj_Weight_Dim+i] > 0 
           ? tmp_wgt[lb->Obj_Weight_Dim+i]*nproc/tmp_wgt[2*lb->Obj_Weight_Dim+i]
           : 1.));
      printf("%s  No. of objects   :  Max = %6d, Min = %6d, "
             "Sum = %7d, Imbal. = %5.3f\n",
        yo2, stats[NUM_GSTATS], 
        stats[3*NUM_GSTATS],
        stats[2*NUM_GSTATS], 
        (stats[2*NUM_GSTATS] > 0
              ? stats[NUM_GSTATS]*nproc/stats[2*NUM_GSTATS]
              : 1.));
      if (lb->Get_Num_Edges != NULL){
        printf("%s  Cut weight       :  Max = %6d, Min = %6d, "
               "Sum = %7d, Imbal. = %5.3f\n",
          yo2, stats[NUM_GSTATS+1], 
          stats[3*NUM_GSTATS+1],
          stats[2*NUM_GSTATS+1], 
          (stats[2*NUM_GSTATS+1] > 0 
                ? stats[NUM_GSTATS+1]*nproc/stats[2*NUM_GSTATS+1]
                : 1.));
        printf("%s  Boundary objects :  Max = %6d, Min = %6d, "
               "Sum = %7d, Imbal. = %5.3f\n",
          yo2, stats[NUM_GSTATS+2], 
          stats[3*NUM_GSTATS+2],
          stats[2*NUM_GSTATS+2], 
          (stats[2*NUM_GSTATS+2] > 0
                ? stats[NUM_GSTATS+2]*nproc/stats[2*NUM_GSTATS+2]
                : 1.));
        printf("%s  Adjacent procs   :  Max = %6d, Min = %6d, "
               "Sum = %7d, Imbal. = %5.3f\n",
          yo2, stats[NUM_GSTATS+3], 
          stats[3*NUM_GSTATS+3],
          stats[2*NUM_GSTATS+3], 
          (stats[2*NUM_GSTATS+3] > 0
                ? stats[NUM_GSTATS+3]*nproc/stats[2*NUM_GSTATS+3]
                : 1.));
      }
      printf("\n");
    }
  }

  /* Copy results to output parameters if desired */
  if (nobj) *nobj = num_obj;
  if (nadj) *nadj = num_adj;
  if (nboundary) *nboundary = num_boundary;
  if (cut_wgt) *cut_wgt = cut_weight;
  if (obj_wgt){
    for (i=0; i<lb->Obj_Weight_Dim; i++) 
      obj_wgt[i] = tmp_wgt[i];
  }

  /* Free data */
  LB_FREE(&global_ids);
  LB_FREE(&local_ids);
  LB_FREE(&tmp_wgt);
  LB_FREE(&vwgts);
  LB_FREE(&ewgts);
  LB_FREE(&nbors_global);
  LB_FREE(&nbors_proc);
  LB_FREE(&proc);
  LB_TRACE_EXIT(lb, yo);
}
