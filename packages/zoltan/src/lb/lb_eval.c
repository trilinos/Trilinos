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

#define NUM_STATS 4 /* Number of graph statistics other than vertex/edge weights */

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
 * EBEB: Current version only works if every partition is wholly
 *       contained within a processor.
 */

{
  char *yo = "Zoltan_LB_Eval";
  int i, j, k, p, num_obj = 0, max_edges, nedges, cuts, flag;
  int num_adj, num_boundary, ierr, compute_part;
  int imin, imax, isum, nproc = zz->Num_Proc, nparts = zz->LB.Num_Global_Parts;
  int stats[6*NUM_STATS];
  int *proc_count, *nbors_proc;
  float *tmp_vwgt, *vwgts, *ewgts, *tmp_cutwgt;
  ZOLTAN_ID_PTR local_ids; 
  ZOLTAN_ID_PTR global_ids, nbors_global;
  int *part;
  ZOLTAN_ID_PTR lid;   /* Temporary pointer to a local id; used to pass NULL
                      pointers to query functions when NUM_LID_ENTRIES = 0. */
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int gid_off, lid_off;
  char msg[256];
  /* Arrays for partition data. */
  int *nobj_arr, *cut_arr, *bndry_arr, *nadj_arr, *all_arr;
  float *vwgt_arr, *cutwgt_arr;
  
  
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
  vwgt_arr = cutwgt_arr = NULL;
  nobj_arr = cut_arr = bndry_arr = nadj_arr = all_arr = NULL;

  ierr = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids, 
                             zz->Obj_Weight_Dim, &vwgts, &part);

  if (ierr == ZOLTAN_FATAL){
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ierr;
  }

  /* Compute statistics w.r.t. partitions? */
  if (zz->LB.Num_Global_Parts != zz->Num_Proc)
    compute_part = 1;
  else
    compute_part = 0;

  if (compute_part){
    /* Allocate space. */
    all_arr = (int *)  ZOLTAN_CALLOC(NUM_STATS*nparts, sizeof(int));
    vwgt_arr = (float *) ZOLTAN_CALLOC(nparts*(zz->Obj_Weight_Dim +
                           zz->Edge_Weight_Dim), sizeof(float));
    if (!all_arr || (zz->Obj_Weight_Dim+zz->Edge_Weight_Dim && !vwgt_arr)){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    nobj_arr = all_arr;
    cut_arr  = nobj_arr + nparts;
    bndry_arr  = cut_arr + nparts;
    nadj_arr = bndry_arr + nparts;
    cutwgt_arr = vwgt_arr + nparts*(zz->Obj_Weight_Dim);

    /* Count number of objects. */
    for (i=0; i<num_obj; i++){
      nobj_arr[part[i]]++;
    }
  }

  /* Compute object weight sums */
  if (zz->Obj_Weight_Dim>0){
    tmp_vwgt = (float *) ZOLTAN_CALLOC(6*zz->Obj_Weight_Dim,  sizeof(float));
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

  if (zz->Get_Num_Edges && zz->Get_Edge_List) {
    /* Use the basic graph query functions */

    /* First compute max no. of edges so we can allocate the right
       amount of space */
    max_edges = 0;
    for (i=0; i< num_obj; i++){
      lid = (num_lid_entries > 0 ? &(local_ids[i*num_lid_entries]) : NULL);
      nedges = zz->Get_Num_Edges(zz->Get_Num_Edges_Data, 
                                 num_gid_entries, num_lid_entries,
                                 &(global_ids[i*num_gid_entries]), 
                                 lid, &ierr);
      if (ierr){
        sprintf(msg, "Get_Num_Edges returned error code %d.", ierr);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
        goto End;
      }
      if (nedges>max_edges) max_edges = nedges;
    }

    /* Allocate edge list space */
    nbors_global = ZOLTAN_MALLOC_GID_ARRAY(zz, max_edges);
    nbors_proc = (int *)ZOLTAN_MALLOC(max_edges * sizeof(int));
    /* Allocate a proc list for computing nadjacent */
    proc_count = (int *)ZOLTAN_MALLOC((zz->Num_Proc)* sizeof(int));
    /* Allocate space for edge weights if needed */
    if (zz->Edge_Weight_Dim){
      ewgts = (float *)ZOLTAN_MALLOC((zz->Edge_Weight_Dim)*max_edges * sizeof(float));
      tmp_cutwgt = (float *) ZOLTAN_MALLOC(6*(zz->Edge_Weight_Dim) * sizeof(float));
    }

    if ((max_edges && ((!nbors_global) || (!nbors_proc) || (zz->Edge_Weight_Dim && !ewgts))) || 
        (zz->Edge_Weight_Dim && (!tmp_cutwgt)) || (!proc_count)){
      ierr = ZOLTAN_MEMERR;
      goto End;
    }

    for (i=0; i<zz->Num_Proc; i++)
      proc_count[i] = 0;
    for (i=0; i<zz->Edge_Weight_Dim; i++)
      tmp_cutwgt[i] = 0;

    for (k=0; k<num_obj; k++){
      flag = 0;
      gid_off = k * num_gid_entries;
      lid_off = k * num_lid_entries;
      lid = (num_lid_entries > 0 ? &(local_ids[lid_off]) : NULL);
      nedges = zz->Get_Num_Edges(zz->Get_Num_Edges_Data, 
                                 num_gid_entries, num_lid_entries,
                                 &(global_ids[gid_off]),
                                 lid, &ierr);
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
      /* Check for cut edges */
      for (j=0; j<nedges; j++){
        if (nbors_proc[j] != zz->Proc){
          cuts++;
          for (i=0; i<zz->Edge_Weight_Dim; i++)
            tmp_cutwgt[i] += ewgts[j*(zz->Edge_Weight_Dim)+i];
          if (flag==0){
            num_boundary++;
            flag = 1;
          }
          proc_count[nbors_proc[j]]++;
        }
        if (compute_part){
          if (nbors_proc[j] == zz->Proc)
            p = zz->Get_Partition(zz->Get_Partition_Data,
                    num_gid_entries, 0,
                    &(nbors_global[j*num_gid_entries]), NULL, &ierr);
          else
            p = -1;
          if ((ierr != ZOLTAN_OK) && (ierr != ZOLTAN_WARN)) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                               "Error returned from ZOLTAN_PARTITION_FN");
            goto End;
          }
          if (p != part[k]){
            cut_arr[part[k]]++;
            for (i=0; i<zz->Edge_Weight_Dim; i++)
              cutwgt_arr[part[k]*zz->Edge_Weight_Dim+i] 
                += ewgts[j*zz->Edge_Weight_Dim+i];
          }
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
    /* Global reduction for object weights. */
    if (zz->Obj_Weight_Dim>0){
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[1*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_MIN, zz->Communicator);
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[2*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[3*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_SUM, zz->Communicator);
    }
    /* Global reduction for comm weights. */
    if (zz->Edge_Weight_Dim>0 && zz->Get_Num_Edges && zz->Get_Edge_List){
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[1*(zz->Edge_Weight_Dim)], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_MIN, zz->Communicator);
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[2*zz->Edge_Weight_Dim], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[3*(zz->Edge_Weight_Dim)], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_SUM, zz->Communicator);
    }
    /* fflush(stdout); */
     
    stats[0] = num_obj;
    stats[1] = cuts;
    stats[2] = num_boundary;
    stats[3] = num_adj;

    /* Compute min, max, sum in the upper portions of the stats array. */
    imin = 1;
    imax = 2;
    isum = 3;
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

      if (zz->Get_Num_Edges && zz->Get_Edge_List){
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
    /* Assume that each partition is WHOLLY CONTAINED WITHIN A PROC. */

    /* Find local min, max, sum. */
    for (i=0; i<6*NUM_STATS; i++)
      stats[i] = 0;
    for (i=0; i<NUM_STATS; i++)
      stats[i] = INT_MAX; /* set min to very large number */
 
    Zoltan_LB_Proc_To_Part(zz, zz->Proc, &p, &k);
    for (i=k; i<k+p; i++){
      for (j=0; j<NUM_STATS; j++){
        if (all_arr[j*nparts+i] < stats[j])           /* min */
          stats[j] = all_arr[j*nparts+i];
        if (all_arr[j*nparts+i] > stats[NUM_STATS+j]) /* max */
          stats[NUM_STATS+j] = all_arr[j*nparts+i];
        stats[2*NUM_STATS+j] += all_arr[j*nparts+i];  /* sum */
      }
    }

    /* Compute global min, max, sum in the upper part of the stats array. */
    imin = 3;
    imax = 4;
    isum = 5;
    MPI_Allreduce(stats, &stats[imin*NUM_STATS], NUM_STATS, MPI_INT, MPI_MIN, 
                  zz->Communicator);
    MPI_Allreduce(&stats[NUM_STATS], &stats[imax*NUM_STATS], NUM_STATS, 
                  MPI_INT, MPI_MAX, zz->Communicator);
    MPI_Allreduce(&stats[2*NUM_STATS], &stats[isum*NUM_STATS], NUM_STATS, 
                  MPI_INT, MPI_SUM, zz->Communicator);

    /* Local and global reduction for object weights. */
    if (zz->Obj_Weight_Dim>0){
      for (i=0; i<zz->Obj_Weight_Dim; i++)
        tmp_vwgt[i] = FLT_MAX; /*min */
      for (i=zz->Obj_Weight_Dim; i<6*zz->Obj_Weight_Dim; i++)
        tmp_vwgt[i] = 0; /* max and sum */
      Zoltan_LB_Proc_To_Part(zz, zz->Proc, &p, &k);
      for (j=k; j<k+p; j++){
        for (i=0; i<zz->Obj_Weight_Dim; i++){
          if (vwgt_arr[j*zz->Obj_Weight_Dim+i] < tmp_vwgt[i]) 
            tmp_vwgt[i] = vwgt_arr[j*zz->Obj_Weight_Dim+i];
          if (vwgt_arr[j*zz->Obj_Weight_Dim+i] > tmp_vwgt[zz->Obj_Weight_Dim+i]) 
            tmp_vwgt[zz->Obj_Weight_Dim+i] = vwgt_arr[j*zz->Obj_Weight_Dim+i];
          tmp_vwgt[2*zz->Obj_Weight_Dim+i] += vwgt_arr[j*zz->Obj_Weight_Dim+i]; 
        }
      }

      
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[imin*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_MIN, zz->Communicator);
      MPI_Allreduce(&tmp_vwgt[zz->Obj_Weight_Dim], 
                    &tmp_vwgt[imax*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(&tmp_vwgt[2*zz->Obj_Weight_Dim], 
                    &tmp_vwgt[isum*zz->Obj_Weight_Dim], 
                    zz->Obj_Weight_Dim, MPI_FLOAT, MPI_SUM, zz->Communicator);
    }

    /* Local and global reduction for cut weights. */
    if (zz->Edge_Weight_Dim>0){
      for (i=0; i<zz->Edge_Weight_Dim; i++)
        tmp_cutwgt[i] = FLT_MAX; /*min */
      for (i=zz->Edge_Weight_Dim; i<6*zz->Edge_Weight_Dim; i++)
        tmp_cutwgt[i] = 0; /* max and sum */
      Zoltan_LB_Proc_To_Part(zz, zz->Proc, &p, &k);
      for (j=k; j<k+p; j++){
        for (i=0; i<zz->Edge_Weight_Dim; i++){
          if (cutwgt_arr[j*zz->Edge_Weight_Dim+i] < tmp_cutwgt[i]) 
            tmp_cutwgt[i] = cutwgt_arr[j*zz->Edge_Weight_Dim+i];
          if (cutwgt_arr[j*zz->Edge_Weight_Dim+i] > tmp_cutwgt[zz->Edge_Weight_Dim+i]) 
            tmp_cutwgt[zz->Edge_Weight_Dim+i] = cutwgt_arr[j*zz->Edge_Weight_Dim+i];
          tmp_cutwgt[2*zz->Edge_Weight_Dim+i] += cutwgt_arr[j*zz->Edge_Weight_Dim+i]; 
        }
      }

      
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[imin*zz->Edge_Weight_Dim], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_MIN, zz->Communicator);
      MPI_Allreduce(&tmp_cutwgt[zz->Edge_Weight_Dim], 
                    &tmp_cutwgt[imax*zz->Edge_Weight_Dim], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(&tmp_cutwgt[2*zz->Edge_Weight_Dim], 
                    &tmp_cutwgt[isum*zz->Edge_Weight_Dim], 
                    zz->Edge_Weight_Dim, MPI_FLOAT, MPI_SUM, zz->Communicator);
    }

    /* Print min-max-sum of results */
    if (zz->Proc == zz->Debug_Proc){
      printf("\n%s  Statistics with respect to %1d partitions: \n", yo, nparts);
      printf("%s                         Min.     Max.      Sum  Imbalance\n", yo);
      printf("%s  No. of objects   :  %8d %8d %8d   %5.3f\n",
        yo, stats[imin*NUM_STATS], 
        stats[imax*NUM_STATS],
        stats[isum*NUM_STATS], 
        (stats[isum*NUM_STATS] > 0
              ? stats[imax*NUM_STATS]*nparts/ (float) stats[isum*NUM_STATS]
              : 1.));
 
      for (i=0; i<zz->Obj_Weight_Dim; i++){
        printf("%s  Object weight #%1d :  %8.3g %8.3g %8.3g   %5.3f\n",
          yo, i, tmp_vwgt[imin*zz->Obj_Weight_Dim+i], 
        tmp_vwgt[imax*zz->Obj_Weight_Dim+i], 
        tmp_vwgt[isum*zz->Obj_Weight_Dim+i], 
        (tmp_vwgt[isum*zz->Obj_Weight_Dim+i] > 0 
           ? tmp_vwgt[imax*zz->Obj_Weight_Dim+i]*nparts/
             tmp_vwgt[isum*zz->Obj_Weight_Dim+i]
           : 1.));
      }

      if (zz->Get_Num_Edges && zz->Get_Edge_List){
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

      for (i=1; i<2 /* NUM_STATS */; i++){
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
          (stats[isum*NUM_STATS+i] > 0 
                ? stats[imax*NUM_STATS+i]*nparts/ (float) stats[isum*NUM_STATS+i]
                : 1.));
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
  if (compute_part){
    ZOLTAN_FREE(&all_arr);
    ZOLTAN_FREE(&vwgt_arr);
  }
  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
