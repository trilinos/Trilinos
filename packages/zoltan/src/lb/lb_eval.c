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
#include <limits.h>
#include <float.h>
#include "parmetis_jostle.h"
#ifndef FLT_MAX
#define FLT_MAX (1e38)
#endif

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

#define NUM_STATS 4 /* Number of graph stats other than vertex/edge weights */
#define NUM_STATS_PART 3 /* Number of graph stats for partitions. */

static int eval_edge_list(ZZ *, int, int, ZOLTAN_ID_PTR, int *, ZOLTAN_ID_PTR,
  int *, float *, int, int, int, int *, int *, int *, int *, int *, float *,
  int *, float *, int *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_LB_Eval (ZZ *zz, int print_stats, 
     int *nobj, float *obj_wgt, 
     int *ncuts, float *cut_wgt, 
     int *nboundary, int *nadj)
/* 
 * Input:
 *   zz          - pointer to Zoltan structure
 *   print_stats - if > 0, compute and print max, min and sum of the metrics
 *
 * Output:
 *   nobj      - number of objects (for each proc)
 *   obj_wgt   - obj_wgt[0:zz->Obj_Weight_Dim-1] are the object weights
 *               (for each proc)
 *   ncuts     - number of cuts (for each proc)
 *   cut_wgt   - cut_wgt[0:zz->Obj_Weight_Dim-1] are the cut weights (for each proc)
 *   nboundary - number of boundary objects (for each proc)
 *   nadj      - the number of adjacent procs (for each proc)
 *
 * Output parameters will only be returned if they are 
 * not NULL on entry.
 *
 * EBEB: We should change the function interface to return min/max/sum stats. 
 * EBEB: These data can be computed both w.r.t. processors and partitions.
 * EBEB: Not all partition stats have been implemented yet.
 *
 */

{
  char *yo = "Zoltan_LB_Eval";
  int i, j, k, max_edges, num_edges;
  int cuts, proc_flag, part_flag;
  int num_obj = 0, num_adj, num_boundary, ierr, compute_part;
  int nproc = zz->Num_Proc, nparts, maxpart, obj_wgt_dim;
  int stats[4*NUM_STATS];
  int imin, imax, isum, iimbal;
  int *proc_count, *nbors_proc;
  float imbal[NUM_STATS];
  float *tmp_vwgt, *vwgts, *ewgts, *tmp_cutwgt, *part_sizes, temp;
  int *edges_per_obj;
  ZOLTAN_ID_PTR local_ids; 
  ZOLTAN_ID_PTR global_ids, nbors_global;
  int *part;
  ZOLTAN_ID_PTR lid;   /* Temporary pointer to a local id; used to pass NULL
                      pointers to query functions when NUM_LID_ENTRIES = 0. */
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int gid_off, lid_off;
  int have_graph_callbacks;
  int edge_list_size;
  int sum;
  char msg[256];
  /* Arrays for partition data. */
  int *nobj_arr, *cut_arr, *bndry_arr, *all_arr, *all_arr_glob;
  float *vwgt_arr, *vwgt_arr_glob, *cutwgt_arr, *cutwgt_arr_glob;
  
  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Set default error code */
  ierr = ZOLTAN_OK;

  /* Set all pointers to NULL */
  global_ids = NULL;
  local_ids = NULL;
  part = NULL;
  tmp_vwgt = tmp_cutwgt = NULL;
  vwgts = ewgts = NULL;
  nbors_global = NULL;
  nbors_proc = proc_count = NULL;
  vwgt_arr = vwgt_arr_glob = cutwgt_arr = cutwgt_arr_glob = NULL;
  nobj_arr = cut_arr = bndry_arr = all_arr = all_arr_glob = NULL;
  edges_per_obj = NULL;

  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, 
                             zz->Obj_Weight_Dim, &vwgts, &part);

  if (ierr == ZOLTAN_FATAL){
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ierr;
  }

  /* Have graph query functions? */
  have_graph_callbacks = 
       ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
        (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL));

  /* Compute statistics w.r.t. partitions? */
  compute_part = (zz->Get_Partition != NULL || zz->Get_Partition_Multi != NULL);

  if (compute_part){
    /* Get partition information. We need LB.Num_Global_Parts */
    /* We don't know if the application uses a Zoltan-style
       mapping of partitions to processsors, so avoid Part2Proc! */
    ierr = Zoltan_LB_Build_PartDist(zz);
    if (ierr == ZOLTAN_FATAL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_LB_Build_PartDist returned error");
      goto End;
    }
    /* Requested number of partitions. */
    nparts = zz->LB.Num_Global_Parts;

    /* Compute actual number of partitions. */
    maxpart = 0;
    for (i=0; i<num_obj; i++)
      if (part[i]>maxpart) 
        maxpart = part[i];
    /* Find max over all procs. */
    i = maxpart;
    MPI_Allreduce(&i, &maxpart, 1, MPI_INT, MPI_MAX,
                  zz->Communicator);

    if (maxpart+1 > nparts){
      ierr = ZOLTAN_WARN;
      sprintf(msg, "Actual number of partitions (%1d) is greater than requested # partitions (%1d)", maxpart+1, nparts);
      ZOLTAN_PRINT_WARN(zz->Proc, yo, msg);
    }

    /* Allocate space. */
    all_arr = (int *)  ZOLTAN_CALLOC(2*NUM_STATS*nparts, sizeof(int));
    vwgt_arr = (float *) ZOLTAN_CALLOC(2*nparts*(zz->Obj_Weight_Dim +
                           zz->Edge_Weight_Dim), sizeof(float));
    if (!all_arr || (zz->Obj_Weight_Dim+zz->Edge_Weight_Dim && !vwgt_arr)){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    nobj_arr = all_arr;
    cut_arr  = nobj_arr + nparts;
    bndry_arr  = cut_arr + nparts;
    all_arr_glob = all_arr + NUM_STATS*nparts;
    vwgt_arr_glob = vwgt_arr + nparts*(zz->Obj_Weight_Dim);
    cutwgt_arr = vwgt_arr_glob + nparts*(zz->Obj_Weight_Dim);
    cutwgt_arr_glob = cutwgt_arr + nparts*(zz->Edge_Weight_Dim);

    /* Count number of objects. */
    for (i=0; i<num_obj; i++){
      nobj_arr[part[i]]++;
    }
  }

  /* Compute object weight sums */
  if (zz->Obj_Weight_Dim>0){
    tmp_vwgt = (float *) ZOLTAN_CALLOC(4*zz->Obj_Weight_Dim,  sizeof(float));
    if (!tmp_vwgt){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    for (i=0; i<num_obj; i++){
      for (j=0; j<zz->Obj_Weight_Dim; j++){
        tmp_vwgt[j] += vwgts[i*zz->Obj_Weight_Dim+j];
        if (compute_part)
          vwgt_arr[part[i]*zz->Obj_Weight_Dim+j] += 
                   vwgts[i*zz->Obj_Weight_Dim+j];
      }
    }
  }


  /* Compute (weighted) edge cuts, #boundary vertices,
     and # adjacent procs if possible */

  num_boundary = 0;
  num_adj = 0;
  cuts = 0;

  if (have_graph_callbacks) {
    /* Use the basic graph query functions */

    /* First compute max no. of edges so we can allocate the right
       amount of space */

    Zoltan_Get_Num_Edges_Per_Obj(zz, num_obj, global_ids, local_ids, 
                                 &edges_per_obj, &max_edges, &num_edges);

    if (zz->Get_Edge_List_Multi) 
      edge_list_size = num_edges;  /* Get all edge info at once */
    else
      edge_list_size = max_edges;  /* Get edge info one obj at a time */

    /* Allocate edge list space */
    nbors_global = ZOLTAN_MALLOC_GID_ARRAY(zz, edge_list_size);
    nbors_proc = (int *)ZOLTAN_MALLOC(edge_list_size * sizeof(int));
    /* Allocate a proc list for computing nadjacent */
    proc_count = (int *)ZOLTAN_MALLOC((zz->Num_Proc)* sizeof(int));
    /* Allocate space for edge weights if needed */
    if (zz->Edge_Weight_Dim){
      ewgts = (float *)ZOLTAN_MALLOC((zz->Edge_Weight_Dim) * edge_list_size 
                                                           * sizeof(float));
      tmp_cutwgt = (float *) ZOLTAN_MALLOC(4 * (zz->Edge_Weight_Dim) 
                                             * sizeof(float));
    }

    if ((edge_list_size && ((!nbors_global) || (!nbors_proc) || (zz->Edge_Weight_Dim && !ewgts))) || 
        (zz->Edge_Weight_Dim && (!tmp_cutwgt)) || (!proc_count)){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i=0; i<zz->Num_Proc; i++)
      proc_count[i] = 0;
    for (i=0; i<zz->Edge_Weight_Dim; i++)
      tmp_cutwgt[i] = 0;

    if (zz->Get_Edge_List_Multi) {
      zz->Get_Edge_List_Multi(zz->Get_Edge_List_Multi_Data, 
                              num_gid_entries, num_lid_entries, num_obj,
                              global_ids, local_ids, edges_per_obj,
                              nbors_global, nbors_proc,
                              zz->Edge_Weight_Dim, ewgts, &ierr);
      if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
        goto End;
      }
      sum = 0;
      for (k = 0; k < num_obj; k++) {
        proc_flag = part_flag = 0;
        ierr = eval_edge_list(zz, k, num_gid_entries, global_ids, part, 
                              &(nbors_global[sum]), &(nbors_proc[sum]), 
                              &(ewgts[sum*zz->Edge_Weight_Dim]), 
                              edges_per_obj[k],
                              num_obj, compute_part, proc_count,
                              &proc_flag, &part_flag, &cuts, &num_boundary, 
                              tmp_cutwgt, cut_arr, cutwgt_arr, bndry_arr);
        if (ierr) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in eval_edge_list");
          goto End;
        }
        sum += edges_per_obj[k];
      }
    }
    else {
      for (k=0; k<num_obj; k++){
        proc_flag = part_flag = 0;
        gid_off = k * num_gid_entries;
        lid_off = k * num_lid_entries;
        lid = (num_lid_entries > 0 ? &(local_ids[lid_off]) : NULL);
        if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
          goto End;
        }
        zz->Get_Edge_List(zz->Get_Edge_List_Data, 
                          num_gid_entries, num_lid_entries,
                          &(global_ids[gid_off]), lid,
                          nbors_global, nbors_proc, 
                          zz->Edge_Weight_Dim, ewgts, &ierr);
        if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)){
          goto End;
        }
  
        ierr = eval_edge_list(zz, k, num_gid_entries, global_ids, part, 
                              nbors_global, nbors_proc, ewgts, edges_per_obj[k],
                              num_obj, compute_part, proc_count,
                              &proc_flag, &part_flag, &cuts, &num_boundary, 
                              tmp_cutwgt, cut_arr, cutwgt_arr, bndry_arr);
        if (ierr) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in eval_edge_list");
          goto End;
        }
      }
    }
    /* Compute the number of adjacent procs */
    for (j=0; j<zz->Num_Proc; j++)
      if (proc_count[j]>0) num_adj++;
  }
  else{
    /* No graph query functions available */
  }
  
  if (print_stats){
    imin = 1;
    imax = 2;
    isum = 3;

    /* Global reduction for object weights. */
    if (zz->Obj_Weight_Dim>0){
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[imin*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_MIN, zz->Communicator);
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[imax*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[isum*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_SUM, zz->Communicator);
    }
    /* Global reduction for comm weights. */
    if (zz->Edge_Weight_Dim>0 && have_graph_callbacks) {
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[imin*(zz->Edge_Weight_Dim)], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_MIN, zz->Communicator);
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[imax*zz->Edge_Weight_Dim], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[isum*(zz->Edge_Weight_Dim)], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_SUM, zz->Communicator);
    }
     
    stats[0] = num_obj;
    stats[1] = cuts;
    stats[2] = num_boundary;
    stats[3] = num_adj;

    /* Compute min, max, sum in the upper portions of the stats array. */
    MPI_Allreduce(stats, &stats[imin*NUM_STATS], NUM_STATS, MPI_INT, MPI_MIN, 
                  zz->Communicator);
    MPI_Allreduce(stats, &stats[imax*NUM_STATS], NUM_STATS, MPI_INT, MPI_MAX, 
                  zz->Communicator);
    MPI_Allreduce(stats, &stats[isum*NUM_STATS], NUM_STATS, MPI_INT, MPI_SUM, 
                  zz->Communicator);

    /* Print min-max-sum of results */
    if (zz->Proc == zz->Debug_Proc){
      printf("\n%s  Statistics with respect to %1d processors:\n", yo, nproc);
      printf("%s                         Min.     Max.      Sum  Imbalance\n", yo);
      printf("%s  No. of objects   :  %8d %8d %8d   %5.3f\n",
        yo, stats[imin*NUM_STATS], 
        stats[imax*NUM_STATS],
        stats[isum*NUM_STATS], 
        (stats[isum*NUM_STATS] > 0
              ? stats[imax*NUM_STATS]*nproc/ (float) stats[isum*NUM_STATS]
              : 1.));
 
      for (i=0; i<zz->Obj_Weight_Dim; i++){
        printf("%s  Object weight #%1d :  %8.3g %8.3g %8.3g   %5.3f\n",
          yo, i, tmp_vwgt[imin*zz->Obj_Weight_Dim+i], 
          tmp_vwgt[imax*zz->Obj_Weight_Dim+i], 
          tmp_vwgt[isum*zz->Obj_Weight_Dim+i], 
          (tmp_vwgt[isum*zz->Obj_Weight_Dim+i] > 0 
           ? tmp_vwgt[imax*zz->Obj_Weight_Dim+i]*nproc/tmp_vwgt[isum*zz->Obj_Weight_Dim+i]
           : 1.));
      }

      if (have_graph_callbacks) {
        for (i=0; i<zz->Edge_Weight_Dim; i++){
          printf("%s  Cut weight  #%1d   :  %8.3g %8.3g %8.3g   %5.3f\n",
            yo, i, tmp_cutwgt[imin*(zz->Edge_Weight_Dim)+i], 
            tmp_cutwgt[imax*(zz->Edge_Weight_Dim)+i], 
            tmp_cutwgt[isum*(zz->Edge_Weight_Dim)+i], 
            (tmp_cutwgt[isum*(zz->Edge_Weight_Dim)+i] > 0 
             ? tmp_cutwgt[imax*(zz->Edge_Weight_Dim)+i]*nproc/tmp_cutwgt[isum*(zz->Edge_Weight_Dim)+i]
             : 1.));
        }

        for (i=1; i<NUM_STATS; i++){
          switch(i){
          case 1:  sprintf(msg, "%s", "No. of cuts      : ");
                   break;
          case 2:  sprintf(msg, "%s", "Boundary objects : ");
                   break;
          case 3:  sprintf(msg, "%s", "Adjacent procs   : ");
                   break;
          default: sprintf(msg, "%s", "                   ");
                   break;
          }

          printf("%s  %s %8d %8d %8d   %5.3f\n",
            yo, msg, stats[imin*NUM_STATS+i], 
            stats[imax*NUM_STATS+i],
            stats[isum*NUM_STATS+i], 
            (stats[isum*NUM_STATS+i] > 0 
                  ? stats[imax*NUM_STATS+i]*nproc/ (float) stats[isum*NUM_STATS+i]
                  : 1.));
        }
      }
      printf("\n");
    }
  }

  /* Copy results to output parameters if desired */
  if (nobj) *nobj = num_obj;
  if (ncuts) *ncuts = cuts;
  if (nadj) *nadj = num_adj;
  if (nboundary) *nboundary = num_boundary;
  if (obj_wgt){
    for (i=0; i<zz->Obj_Weight_Dim; i++) 
      obj_wgt[i] = tmp_vwgt[i];
  }
  if (cut_wgt){
    for (i=0; i<zz->Edge_Weight_Dim; i++) 
      cut_wgt[i] = tmp_cutwgt[i];
  }


  if (compute_part && print_stats){
    /* Compute statistics w.r.t. partitions. */
    imin = 0;
    imax = 1;
    isum = 2;
    iimbal = 3;

    /* Get partition sizes. */
    obj_wgt_dim = (zz->Obj_Weight_Dim > 0 ? zz->Obj_Weight_Dim : 1);
    part_sizes = (float *) ZOLTAN_MALLOC(nparts * obj_wgt_dim * sizeof(float));
    if (!part_sizes){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    Zoltan_LB_Get_Part_Sizes(zz, nparts, obj_wgt_dim, part_sizes);

    /* Allreduce data w.r.t. partitions onto all procs. */
    MPI_Allreduce(all_arr, all_arr_glob, NUM_STATS_PART*nparts, 
                  MPI_INT, MPI_SUM, zz->Communicator);
    if (zz->Obj_Weight_Dim > 0)
      MPI_Allreduce(vwgt_arr, vwgt_arr_glob, nparts*(zz->Obj_Weight_Dim), 
                    MPI_FLOAT, MPI_SUM, zz->Communicator);
    if (zz->Edge_Weight_Dim > 0)
      MPI_Allreduce(cutwgt_arr, cutwgt_arr_glob, nparts*(zz->Edge_Weight_Dim), 
                    MPI_FLOAT, MPI_SUM, zz->Communicator);

    /* Find min, max, sum. */
    for (i=0; i<NUM_STATS; i++)
      stats[i] = INT_MAX; /* set min to very large number */
    for (i=NUM_STATS; i<3*NUM_STATS; i++)
      stats[i] = 0;       /* max and sum */
 
    for (i=0; i<nparts; i++){
      for (j=0; j<NUM_STATS_PART; j++){
        if (all_arr_glob[j*nparts+i] < stats[j])           /* min */
          stats[j] = all_arr_glob[j*nparts+i];
        if (all_arr_glob[j*nparts+i] > stats[imax*NUM_STATS+j]) /* max */
          stats[imax*NUM_STATS+j] = all_arr_glob[j*nparts+i];
        stats[isum*NUM_STATS+j] += all_arr_glob[j*nparts+i];  /* sum */
      }
    }

    /* Compute scaled imbalance. */
    for (j=0; j<NUM_STATS_PART; j++)
      imbal[j] = 0.0;
    k = (zz->Obj_Weight_Dim>0 ? zz->Obj_Weight_Dim : 1); 
    for (i=0; i<nparts; i++){
      for (j=0; j<NUM_STATS_PART; j++){
        if (all_arr_glob[j*nparts+i]*part_sizes[i*k] != 0.0)
          temp = all_arr_glob[j*nparts+i] / (stats[isum*NUM_STATS+j]
                 * part_sizes[i*k]);
        else
          temp = 1.0;
        if (temp > imbal[j])
          imbal[j] = temp;
      }
    }

    /* Min, max, sum for object weights. */
    if (zz->Obj_Weight_Dim>0){
      for (i=0; i<zz->Obj_Weight_Dim; i++)
        tmp_vwgt[i] = FLT_MAX; /*min */
      for (i=zz->Obj_Weight_Dim; i<4*zz->Obj_Weight_Dim; i++)
        tmp_vwgt[i] = 0.0; /* max and sum and imbalance */

      for (j=0; j<nparts; j++){
        for (i=0; i<zz->Obj_Weight_Dim; i++){
          if (vwgt_arr_glob[j*zz->Obj_Weight_Dim+i] < tmp_vwgt[i]) 
            tmp_vwgt[i] = vwgt_arr_glob[j*zz->Obj_Weight_Dim+i];
          if (vwgt_arr_glob[j*zz->Obj_Weight_Dim+i] > tmp_vwgt[imax*zz->Obj_Weight_Dim+i]) 
            tmp_vwgt[imax*zz->Obj_Weight_Dim+i] = vwgt_arr_glob[j*zz->Obj_Weight_Dim+i];
          tmp_vwgt[isum*zz->Obj_Weight_Dim+i] += vwgt_arr_glob[j*zz->Obj_Weight_Dim+i]; 
        }
      }

      /* Compute scaled imbalance. */
      for (j=0; j<nparts; j++){
        for (i=0; i<zz->Obj_Weight_Dim; i++){
          if (tmp_vwgt[isum*zz->Obj_Weight_Dim+i]*part_sizes[j*zz->Obj_Weight_Dim+i] != 0.0)
            temp = vwgt_arr_glob[j*zz->Obj_Weight_Dim+i]/(tmp_vwgt[isum*zz->Obj_Weight_Dim+i]*part_sizes[j*zz->Obj_Weight_Dim+i]);
          else
            temp = 1.0;
          if (temp > tmp_vwgt[iimbal*zz->Obj_Weight_Dim+i])
            tmp_vwgt[iimbal*zz->Obj_Weight_Dim+i] = temp;
        }
      }
    }

    /* Min, max, sum for cut weights. */
    if (zz->Edge_Weight_Dim>0){
      for (i=0; i<zz->Edge_Weight_Dim; i++)
        tmp_cutwgt[i] = FLT_MAX; /*min */
      for (i=zz->Edge_Weight_Dim; i<3*zz->Edge_Weight_Dim; i++)
        tmp_cutwgt[i] = 0; /* max and sum */

      for (j=0; j<nparts; j++){
        for (i=0; i<zz->Edge_Weight_Dim; i++){
          if (cutwgt_arr_glob[j*zz->Edge_Weight_Dim+i] < tmp_cutwgt[i]) 
            tmp_cutwgt[i] = cutwgt_arr_glob[j*zz->Edge_Weight_Dim+i];
          if (cutwgt_arr_glob[j*zz->Edge_Weight_Dim+i] > tmp_cutwgt[imax*zz->Edge_Weight_Dim+i]) 
            tmp_cutwgt[imax*zz->Edge_Weight_Dim+i] = cutwgt_arr_glob[j*zz->Edge_Weight_Dim+i];
          tmp_cutwgt[isum*zz->Edge_Weight_Dim+i] += cutwgt_arr_glob[j*zz->Edge_Weight_Dim+i]; 
        }
      }

    }

    /* Print min-max-sum of results */
    if (zz->Proc == zz->Debug_Proc){
      printf("\n%s  Requested # partitions = %1d, actual # partitions = %1d\n", yo, nparts, maxpart+1);
      printf("%s  Statistics with respect to %1d partitions: \n", yo, nparts);
      printf("%s                         Min.     Max.      Sum  Imbalance\n", yo);
      printf("%s  No. of objects   :  %8d %8d %8d   %5.3f\n",
        yo, stats[imin*NUM_STATS], 
        stats[imax*NUM_STATS],
        stats[isum*NUM_STATS], 
        imbal[0]);
 
      for (i=0; i<zz->Obj_Weight_Dim; i++){
        printf("%s  Object weight #%1d :  %8.3g %8.3g %8.3g   %5.3f\n",
          yo, i, tmp_vwgt[imin*zz->Obj_Weight_Dim+i], 
        tmp_vwgt[imax*zz->Obj_Weight_Dim+i], 
        tmp_vwgt[isum*zz->Obj_Weight_Dim+i], 
        tmp_vwgt[iimbal*zz->Obj_Weight_Dim+i]);
      }

      if (have_graph_callbacks) {
        for (i=0; i<zz->Edge_Weight_Dim; i++){
          printf("%s  Cut weight  #%1d   :  %8.3g %8.3g %8.3g   %5.3f\n",
            yo, i, tmp_cutwgt[imin*(zz->Edge_Weight_Dim)+i], 
          tmp_cutwgt[imax*(zz->Edge_Weight_Dim)+i], 
          tmp_cutwgt[isum*(zz->Edge_Weight_Dim)+i], 
          (tmp_cutwgt[isum*(zz->Edge_Weight_Dim)+i] > 0 
             ? tmp_cutwgt[imax*(zz->Edge_Weight_Dim)+i]*nparts/tmp_cutwgt[isum*(zz->Edge_Weight_Dim)+i]
             : 1.));
        }
      }

      for (i=1; i<NUM_STATS_PART; i++){
        switch(i){
        case 1:  sprintf(msg, "%s", "No. of cuts      : ");
                 break;
        case 2:  sprintf(msg, "%s", "Boundary objects : ");
                 break;
        case 3:  sprintf(msg, "%s", "Adj. partitions  : ");
                 break;
        default: sprintf(msg, "%s", "                   ");
                 break;
        }

        printf("%s  %s %8d %8d %8d   %5.3f\n",
          yo, msg, stats[imin*NUM_STATS+i], 
          stats[imax*NUM_STATS+i],
          stats[isum*NUM_STATS+i], 
          imbal[i]);
      }
      printf("\n");
    }

  }

End:
  /* Free data */
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&part); 
  ZOLTAN_FREE(&tmp_vwgt);
  ZOLTAN_FREE(&tmp_cutwgt);
  ZOLTAN_FREE(&vwgts);
  ZOLTAN_FREE(&ewgts);
  ZOLTAN_FREE(&nbors_global);
  ZOLTAN_FREE(&nbors_proc);
  ZOLTAN_FREE(&proc_count);
  ZOLTAN_FREE(&edges_per_obj);
  if (compute_part){
    ZOLTAN_FREE(&all_arr);
    ZOLTAN_FREE(&vwgt_arr);
    ZOLTAN_FREE(&part_sizes);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}
/***************************************************************************/
static int eval_edge_list(
  ZZ *zz, 
  int k,                      /* Index of object being processed */
  int num_gid_entries,        /* # of unsigned ints in a GID */
  ZOLTAN_ID_PTR global_ids,   /* Array of GIDs on this proc */
  int *part,                  /* Partition numbers of GIDs on this proc */
  ZOLTAN_ID_PTR nbors_global, /* GIDs of this obj's neighboring objs. */
  int *nbors_proc,            /* Procs owning the neighboring objs. */
  float *ewgts,               /* Edge weights for each of obj's edges */
  int nedges,                 /* # of edges for this object */
  int num_obj,                /* # of objs on this proc */
  int compute_part,           /* Flag:  compute partition stats */
  int *proc_count,            /* # of nbors on each other proc */
  int *proc_flag,             /* Flag:  indicates whether this obj has been
                                 counted as a boundary object yet. */
  int *part_flag,             /* Flag:  indicates whether this obj has been
                                 counted as a partition boundary obj yet. */
  int *cuts,                  /* # of cut edges */
  int *num_boundary,          /* # of boundary objs */
  float *tmp_cutwgt,          /* total weight of cut edges (and other stats) */
  int *cut_arr,               /* # of partition cut edges. */
  float *cutwgt_arr,          /* weights of partition cut edges */
  int *bndry_arr              /* # of partition boundary objs */
)
{
/* Function to evaluate edge cuts, etc., for an object's edge list. */
char *yo = "eval_edge_list";
int i, j, p, found;
int ierr = ZOLTAN_OK;

  /* Check for cut edges */
  for (j=0; j<nedges; j++){
    if (nbors_proc[j] != zz->Proc){
      (*cuts)++;
      for (i=0; i<zz->Edge_Weight_Dim; i++)
        tmp_cutwgt[i] += ewgts[j*(zz->Edge_Weight_Dim)+i];
      if ((*proc_flag)==0){
        (*num_boundary)++;
        (*proc_flag) = 1;
      }
      proc_count[nbors_proc[j]]++;
    }
    if (compute_part){
      if (nbors_proc[j] == zz->Proc){
        /* Need to find the nbors_global ID in global_ids
         * in order to determine the partition numbers.
         * For now, look through the global_id list and compare
         * every time. This is quite slow, O(n1*n2) where
         * n1, n2 are the length of the lists, so
         * in the future we should use a hash table
         * or the Ddirectory.
         */
        found = -1;
        for (i=0; i<num_obj; i++){
          if (ZOLTAN_EQ_GID(zz, &(global_ids[i*num_gid_entries]),
                            &(nbors_global[j*num_gid_entries]))){
            found = i;
            break;
          }
        }
        if (found<0){
          /* This should never happen. */
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Sanity check failed!");
          ierr = ZOLTAN_FATAL;
          goto End;
        }
        p = part[found];
      }
      else
        p = -1; /* Assume remote data belong to different partition. */
                /* EBEB: This is not necessarily true! */
                /* Answers may be incorrect. Fix: Use Ddirectory. */

      if (p != part[k]){
        cut_arr[part[k]]++;
        for (i=0; i<zz->Edge_Weight_Dim; i++)
          cutwgt_arr[part[k]*zz->Edge_Weight_Dim+i] 
            += ewgts[j*zz->Edge_Weight_Dim+i];
        if ((*part_flag)==0){
          bndry_arr[part[k]]++;
          (*part_flag) = 1;
        }
      }
    }
  }
End:
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
