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
 */

{
  char *yo = "Zoltan_LB_Eval";
  int i, j, k, num_obj, max_edges, nedges, cuts, flag;
  int num_adj, num_boundary, ierr;
  int stats[4*NUM_STATS];
  int *proc, *nbors_proc;
  float *tmp_vwgt, *vwgts, *ewgts, *tmp_cutwgt, nproc;
  ZOLTAN_ID_PTR local_ids; 
  ZOLTAN_ID_PTR global_ids, nbors_global;
  ZOLTAN_ID_PTR lid;   /* Temporary pointer to a local id; used to pass NULL
                      pointers to query functions when NUM_LID_ENTRIES = 0. */
  int num_gid_entries = zz->Num_GID;
  int num_lid_entries = zz->Num_LID;
  int gid_off, lid_off;
  char msg[256];
  
  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Set default error code */
  ierr = ZOLTAN_OK;

  /* Set all pointers to NULL */
  global_ids = NULL;
  local_ids = NULL;
  tmp_vwgt = NULL;
  tmp_cutwgt = NULL;
  vwgts = NULL;
  ewgts = NULL;
  nbors_global = NULL;
  nbors_proc = NULL;
  proc = NULL;

  /* First compute number of objs and object weight on each proc */
  num_obj = zz->Get_Num_Obj(zz->Get_Num_Obj_Data, &ierr);

  if (num_obj>0){

    /* Allocate space for object data */
    global_ids = ZOLTAN_MALLOC_GID_ARRAY(zz, num_obj);
    local_ids  = ZOLTAN_MALLOC_LID_ARRAY(zz, num_obj);
      
    if ((!global_ids) || (num_lid_entries && !local_ids)){
      ZOLTAN_FREE(&global_ids);
      ZOLTAN_FREE(&local_ids);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
  }

  /* Allocate space for weights if needed */
  if (zz->Obj_Weight_Dim>0){
    vwgts   = (float  *) ZOLTAN_MALLOC(zz->Obj_Weight_Dim*num_obj * sizeof(float));
    tmp_vwgt = (float *) ZOLTAN_MALLOC(4*zz->Obj_Weight_Dim * sizeof(float));
    if ((num_obj && !vwgts) || (!tmp_vwgt)){
      ZOLTAN_FREE(&global_ids);
      ZOLTAN_FREE(&local_ids);
      ZOLTAN_FREE(&vwgts);
      ZOLTAN_FREE(&tmp_vwgt);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }
  } 
  
  Zoltan_Get_Obj_List(zz, global_ids, local_ids, zz->Obj_Weight_Dim, vwgts, 
    &ierr);
  if (ierr == ZOLTAN_FATAL){
    ZOLTAN_FREE(&global_ids);
    ZOLTAN_FREE(&local_ids);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ierr;
  }
  

  /* Compute object weight sums */
  if (zz->Obj_Weight_Dim>0){
    for (j=0; j<zz->Obj_Weight_Dim; j++)
      tmp_vwgt[j] = 0;
    for (i=0; i<num_obj; i++){
      for (j=0; j<zz->Obj_Weight_Dim; j++){
        tmp_vwgt[j] += vwgts[i*zz->Obj_Weight_Dim+j];
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
        ZOLTAN_FREE(&global_ids);
        ZOLTAN_FREE(&local_ids);
        ZOLTAN_FREE(&vwgts);
        ZOLTAN_FREE(&tmp_vwgt);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ierr;
      }
      if (nedges>max_edges) max_edges = nedges;
    }

    /* Allocate edge list space */
    nbors_global = ZOLTAN_MALLOC_GID_ARRAY(zz, max_edges);
    nbors_proc = (int *)ZOLTAN_MALLOC(max_edges * sizeof(int));
    /* Allocate a proc list for computing nadjacent */
    proc = (int *)ZOLTAN_MALLOC((zz->Num_Proc)* sizeof(int));
    /* Allocate space for edge weights if needed */
    if (zz->Edge_Weight_Dim){
      ewgts = (float *)ZOLTAN_MALLOC((zz->Edge_Weight_Dim)*max_edges * sizeof(float));
      tmp_cutwgt = (float *) ZOLTAN_MALLOC(4*(zz->Edge_Weight_Dim) * sizeof(float));
    }

    if ((max_edges && ((!nbors_global) || (!nbors_proc) || (zz->Edge_Weight_Dim && !ewgts))) || 
        (zz->Edge_Weight_Dim && (!tmp_cutwgt)) || (!proc)){
      ZOLTAN_FREE(&global_ids);
      ZOLTAN_FREE(&local_ids);
      ZOLTAN_FREE(&vwgts);
      ZOLTAN_FREE(&tmp_vwgt);
      ZOLTAN_FREE(&nbors_global);
      ZOLTAN_FREE(&nbors_proc);
      ZOLTAN_FREE(&ewgts);
      ZOLTAN_FREE(&tmp_cutwgt);
      ZOLTAN_FREE(&proc);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return ZOLTAN_MEMERR;
    }

    for (i=0; i<zz->Num_Proc; i++)
      proc[i] = 0;
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
      if (ierr == ZOLTAN_FATAL){
        ZOLTAN_FREE(&global_ids);
        ZOLTAN_FREE(&local_ids);
        ZOLTAN_FREE(&vwgts);
        ZOLTAN_FREE(&tmp_vwgt);
        ZOLTAN_FREE(&nbors_global);
        ZOLTAN_FREE(&nbors_proc);
        ZOLTAN_FREE(&ewgts);
        ZOLTAN_FREE(&tmp_cutwgt);
        ZOLTAN_FREE(&proc);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ierr;
      }
      zz->Get_Edge_List(zz->Get_Edge_List_Data, 
                        num_gid_entries, num_lid_entries,
                        &(global_ids[gid_off]), lid,
                        nbors_global, nbors_proc, 
                        zz->Edge_Weight_Dim, ewgts, &ierr);
      if (ierr == ZOLTAN_FATAL){
        ZOLTAN_FREE(&global_ids);
        ZOLTAN_FREE(&local_ids);
        ZOLTAN_FREE(&vwgts);
        ZOLTAN_FREE(&tmp_vwgt);
        ZOLTAN_FREE(&nbors_global);
        ZOLTAN_FREE(&nbors_proc);
        ZOLTAN_FREE(&ewgts);
        ZOLTAN_FREE(&tmp_cutwgt);
        ZOLTAN_FREE(&proc);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return ierr;
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
          proc[nbors_proc[j]]++;
        }
      }
    }
    /* Compute the number of adjacent procs */
    for (j=0; j<zz->Num_Proc; j++)
      if (proc[j]>0) num_adj++;
  }
  else{
    /* No graph query functions available */
  }
  
  if (print_stats){
    /* Global reduction for object weights. */
    if (zz->Obj_Weight_Dim>0){
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[zz->Obj_Weight_Dim], zz->Obj_Weight_Dim, 
                    MPI_FLOAT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[2*zz->Obj_Weight_Dim], zz->Obj_Weight_Dim,
                    MPI_FLOAT, MPI_SUM, zz->Communicator);
      MPI_Allreduce(tmp_vwgt, &tmp_vwgt[3*zz->Obj_Weight_Dim], zz->Obj_Weight_Dim,
                    MPI_FLOAT, MPI_MIN, zz->Communicator);
    }
    /* Global reduction for comm weights. */
    if (zz->Edge_Weight_Dim>0 && zz->Get_Num_Edges && zz->Get_Edge_List){
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[zz->Edge_Weight_Dim], zz->Edge_Weight_Dim, 
                    MPI_FLOAT, MPI_MAX, zz->Communicator);
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[2*(zz->Edge_Weight_Dim)], zz->Edge_Weight_Dim,
                    MPI_FLOAT, MPI_SUM, zz->Communicator);
      MPI_Allreduce(tmp_cutwgt, &tmp_cutwgt[3*(zz->Edge_Weight_Dim)], zz->Edge_Weight_Dim,
                    MPI_FLOAT, MPI_MIN, zz->Communicator);
    }
    fflush(stdout);
     
    stats[0] = num_obj;
    stats[1] = cuts;
    stats[2] = num_boundary;
    stats[3] = num_adj;

    /* Compute max, min and sum in the upper portions of the stats array. */
    MPI_Allreduce(stats, &stats[NUM_STATS], NUM_STATS, MPI_INT, MPI_MAX, 
                  zz->Communicator);
    MPI_Allreduce(stats, &stats[2*NUM_STATS], NUM_STATS, MPI_INT, MPI_SUM, 
                  zz->Communicator);
    MPI_Allreduce(stats, &stats[3*NUM_STATS], NUM_STATS, MPI_INT, MPI_MIN, 
                  zz->Communicator);

    /* Print max-sum of results */
    nproc = zz->Num_Proc; /* convert to float */
    if (zz->Proc == zz->Debug_Proc){
      printf("\n%s  Statistics for current partitioning/balance:\n", yo);
      printf("%s  No. of objects   :  Max = %7d, Min = %7d, "
             "Sum = %7d, Imbal. = %5.3f\n",
        yo, stats[NUM_STATS], 
        stats[3*NUM_STATS],
        stats[2*NUM_STATS], 
        (stats[2*NUM_STATS] > 0
              ? stats[NUM_STATS]*nproc/stats[2*NUM_STATS]
              : 1.));
 
      for (i=0; i<zz->Obj_Weight_Dim; i++){
        printf("%s  Object weight #%1d :  Max =%8.3g, Min =%8.3g, "
          "Sum =%8.3g, Imbal. = %5.3f\n",
          yo, i+1, tmp_vwgt[zz->Obj_Weight_Dim+i], 
          tmp_vwgt[3*zz->Obj_Weight_Dim+i], 
          tmp_vwgt[2*zz->Obj_Weight_Dim+i], 
          (tmp_vwgt[2*zz->Obj_Weight_Dim+i] > 0 
           ? tmp_vwgt[zz->Obj_Weight_Dim+i]*nproc/tmp_vwgt[2*zz->Obj_Weight_Dim+i]
           : 1.));
      }

      if (zz->Get_Num_Edges && zz->Get_Edge_List){
        for (i=0; i<zz->Edge_Weight_Dim; i++){
          printf("%s  Comm. weight  #%1d :  Max =%8.3g, Min =%8.3g, "
            "Sum =%8.3g, Imbal. = %5.3f\n",
            yo, i+1, tmp_cutwgt[(zz->Edge_Weight_Dim)+i], 
            tmp_cutwgt[3*(zz->Edge_Weight_Dim)+i], 
            tmp_cutwgt[2*(zz->Edge_Weight_Dim)+i], 
            (tmp_cutwgt[2*(zz->Edge_Weight_Dim)+i] > 0 
             ? tmp_cutwgt[(zz->Edge_Weight_Dim)+i]*nproc/tmp_cutwgt[2*(zz->Edge_Weight_Dim)+i]
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

          printf("%s  %s Max = %7d, Min = %7d, Sum = %7d, Imbal. = %5.3f\n",
            yo, msg, stats[NUM_STATS+i], 
            stats[3*NUM_STATS+i],
            stats[2*NUM_STATS+i], 
            (stats[2*NUM_STATS+i] > 0 
                  ? stats[NUM_STATS+i]*nproc/stats[2*NUM_STATS+i]
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

  /* Free data */
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&tmp_vwgt);
  ZOLTAN_FREE(&tmp_cutwgt);
  ZOLTAN_FREE(&vwgts);
  ZOLTAN_FREE(&ewgts);
  ZOLTAN_FREE(&nbors_global);
  ZOLTAN_FREE(&nbors_proc);
  ZOLTAN_FREE(&proc);
  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
