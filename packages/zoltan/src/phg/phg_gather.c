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
#include "phg.h"
#include <limits.h>


/****************************************************************************/
int Zoltan_PHG_Gather_To_All_Procs(
  ZZ *zz, 
  HGraph *phg,           /* Input:   Local part of distributed hypergraph */
  PHGPartParams *hgp,        /* Input:   Hypergraph parameters */
  PHGComm *scomm,        /* Input:   Serial PHGComm for use by shg. */
  HGraph **gathered_hg   /* Output:  combined hypergraph combined to proc */
)
{
/* 
 * Function to gather distributed hypergraph onto each processor for
 * coarsest partitioning.
 * First hypergraph arrays for the hypergraph on a column of processors
 * are built using MPI_Allgathers down the processor columns.
 * These hypergraph arrays contain complete info about a subset of vertices.
 * Second the column hypergraphs are gathered along processor rows.
 * Each processor then has a complete description of the hypergraph.
 */
char *yo = "Zoltan_PHG_Gather_To_All_Procs";
int ierr = ZOLTAN_OK;
int i, tmp, sum;
int *each = NULL,
    *disp = NULL;      /* Size and displacement arrays for MPI_Allgatherv */
int *send_buf = NULL;    /* Buffer of values to be sent */
int send_size;           /* Size of buffer send_buf */
int *col_vedge = NULL;   /* vedge array for the proc-column hypergraph */
int *col_vindex = NULL;  /* vindex array for the proc-column hypergraph */
int *col_hvertex = NULL; /* hvertex array for the proc-column hypergraph */
int *col_hindex = NULL;  /* hindex array for the proc-column hypergraph */
int col_nVtx;            /* Number of vertices in processor column */
int col_nEdge;           /* Number of edges in processor column */
int col_nPin;            /* Number of pins in processor column */

int *recv_size = NULL;   /* nPins for each proc in col or row */

HGraph *shg;             /* Pointer to the serial hypergraph to be
                            returned by this function. */

int myProc_x = phg->comm->myProc_x;
int nProc_x = phg->comm->nProc_x;
int nProc_y = phg->comm->nProc_y;
int max_nProc_xy = MAX(nProc_x, nProc_y);

  if (phg->comm->nProc == 1) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Do not call this routine on one proc.");
    return ZOLTAN_FATAL;
  }

#ifdef KDDKDD_CHECK
  Zoltan_HG_Print(zz, phg, NULL, stdout, "GatherBefore");/* NULL parts for now;
                                                           add non-NULL later */
#endif

  /******************************************************************
   *  0. Allocate the hypergraph to be returned. 
   *  Set values that we already know. 
   ******************************************************************/

  shg = *gathered_hg = (HGraph *) ZOLTAN_MALLOC(sizeof(HGraph));
  if (!shg) MEMORY_ERROR;

  Zoltan_HG_HGraph_Init(shg);
  shg->nVtx = phg->dist_x[nProc_x];    /* TODO64 - can this exceed 2B? */
  shg->nEdge = phg->dist_y[nProc_y];

  shg->dist_x = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(2 * sizeof(ZOLTAN_GNO_TYPE));
  shg->dist_y = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(2 * sizeof(ZOLTAN_GNO_TYPE));
  if (!shg->dist_x || !shg->dist_y) MEMORY_ERROR;

  shg->dist_x[0] = shg->dist_y[0] = 0;
  shg->dist_x[1] = shg->nVtx;
  shg->dist_y[1] = shg->nEdge;

  shg->comm = scomm;

  shg->EdgeWeightDim = phg->EdgeWeightDim;
  shg->VtxWeightDim = phg->VtxWeightDim;
  if (shg->VtxWeightDim && shg->nVtx)
    shg->vwgt = (float *) ZOLTAN_MALLOC(shg->nVtx * shg->VtxWeightDim 
                                                  * sizeof(float));
  if (shg->EdgeWeightDim && shg->nEdge)
    shg->ewgt = (float *) ZOLTAN_MALLOC(shg->nEdge * shg->EdgeWeightDim 
                                                  * sizeof(float));
  /* Fixed vertices */
  shg->bisec_split = phg->bisec_split;
  if (hgp->UseFixedVtx)
    shg->fixed_part = (int *) ZOLTAN_MALLOC(shg->nVtx * sizeof(int));
  if (hgp->UsePrefPart)
    shg->pref_part = (int *) ZOLTAN_MALLOC(shg->nVtx * sizeof(int));
  
  /* Allocate arrays for use in gather operations */
  recv_size = (int *) ZOLTAN_MALLOC(3 * max_nProc_xy * sizeof(int));
  each = recv_size + max_nProc_xy;
  disp = each + max_nProc_xy;
 
  /* TODO64 - phg->dist_y[nProc_y] could exceed 2 Billion, NO? */
  send_size = MAX(phg->dist_x[myProc_x+1] - phg->dist_x[myProc_x], 
                  phg->dist_y[nProc_y]);
  send_buf = (int *) ZOLTAN_MALLOC(send_size * sizeof(int));
  

  if ((shg->VtxWeightDim && shg->nVtx && !shg->vwgt) ||
      (shg->EdgeWeightDim && shg->nEdge && !shg->ewgt) || !recv_size ||
      (send_size && !send_buf)) 
    MEMORY_ERROR;
  

  /*************************************************************
   *  1. Gather all non-zeros for vertices in processor column *
   *************************************************************/
  
  if (nProc_y == 1) {
    /* 
     * Don't need a gather; just set pointers appropriately for row-gather
     * in Step 2 below.
     */

    col_nVtx = phg->nVtx;
    col_nEdge = phg->nEdge;
    col_nPin = phg->nPins;
    col_vindex = phg->vindex;
    col_vedge = phg->vedge;
    col_hindex = phg->hindex;
    col_hvertex = phg->hvertex;

    for (i = 0; i < shg->EdgeWeightDim * shg->nEdge; i++)
      shg->ewgt[i] = phg->ewgt[i];
  }

  else {

    /* Gather local size info for each proc in column */

    MPI_Allgather(&(phg->nPins), 1, MPI_INT, recv_size, 1, MPI_INT, 
                  phg->comm->col_comm);
  
    /* Compute number of vtx, edge, and nnz in column */
    col_nVtx = (int)(phg->dist_x[myProc_x+1] - phg->dist_x[myProc_x]);
    col_nEdge = phg->dist_y[nProc_y];   /* SCHEMEA */
    col_nPin = 0;
    for (i = 0; i < nProc_y; i++) {
      col_nPin += recv_size[i];
    }
    
    /* Allocate arrays for column hypergraph */
    col_hindex = (int *) ZOLTAN_CALLOC((col_nEdge+1), sizeof(int));
    col_hvertex = (int *) ZOLTAN_MALLOC(col_nPin * sizeof(int));
  
    col_vindex = (int *) ZOLTAN_CALLOC((col_nVtx+1), sizeof(int));
    col_vedge = (int *) ZOLTAN_MALLOC(col_nPin * sizeof(int));
  
    if (!col_vindex || !col_hindex || 
        (col_nPin && (!col_vedge || !col_hvertex)))
      MEMORY_ERROR;
    
    /* Gather hvertex data for all procs in column */
  
    /* SCHEMEA uses same vertex LNO on each proc in column. */
    /* SCHEMEB would require conversion from vertex LNO to GNO here. */
  
    disp[0] = 0;
    for (i = 1; i < nProc_y; i++)
      disp[i] = disp[i-1] + recv_size[i-1];
  
    MPI_Allgatherv(phg->hvertex, phg->nPins, MPI_INT,
                   col_hvertex, recv_size, disp, MPI_INT, phg->comm->col_comm);
  
    /* SCHEMEA uses same vertex LNO on each proc in column. */
    /* SCHEMEB would require conversion from vertex GNO to LNO here */
  
    /* Gather hindex data for all procs in column */
  
    for (i = 0; i < phg->nEdge; i++)
      send_buf[i] = phg->hindex[i+1] - phg->hindex[i];
  
    /* SCHEMEA can assume a recv for each edge;
     * SCHEMEB needs to gather the number of edges recv'd from each proc. */
  
    for (i = 0; i < nProc_y; i++) 
      each[i] = phg->dist_y[i+1] - phg->dist_y[i];

    disp[0] = 0;  /* Can't use dist_y because it may not be sizeof(int) */
    for (i=1; i < nProc_y; i++){
      disp[i] = disp[i-1] + each[i-1];
    }
  
    /* SCHEMEA can use phg->dist_y for displacement array.
     * SCHEMEB requires separate displacement array. */

    MPI_Allgatherv(send_buf, phg->nEdge, MPI_INT, 
                   col_hindex, each, disp, MPI_INT, phg->comm->col_comm);
  
    /* Perform prefix sum on col_hindex */
    sum = 0;
    for (i = 0; i < col_nEdge; i++) {
      tmp = col_hindex[i];
      col_hindex[i] = sum;
      sum += tmp;
    }
    col_hindex[col_nEdge] = sum;

    /* Sanity check */
    if (col_hindex[col_nEdge] != col_nPin) {
      printf("%d Sanity check failed:  "
             "col_hindex[col_nEdge] %d != col_nPin %d\n", 
              zz->Proc, col_hindex[col_nEdge], col_nPin);
      exit(-1);
    }
  
    /* Gather edge weights, if any. */
    if (shg->EdgeWeightDim) {
  
      /* Can use nearly the same each array. */
      /* Need to compute new disp array. */
  
      disp[0] = 0;
      each[0] *= phg->EdgeWeightDim;
      for (i = 1; i < nProc_y; i++) {
        each[i] *= phg->EdgeWeightDim;
        disp[i] = disp[i-1] + each[i-1];
      }
      
      MPI_Allgatherv(phg->ewgt, phg->nEdge*phg->EdgeWeightDim, MPI_FLOAT, 
                     shg->ewgt, each, disp, MPI_FLOAT, phg->comm->col_comm);
    }
   
  
    Zoltan_HG_Mirror(col_nEdge, col_hindex, col_hvertex, 
                     col_nVtx, col_vindex, col_vedge);
  
  }  /* End column-gather */
  
  /*************************************************************
   *  2. Gather all non-zeros for edges in processor rows      *
   *  All processors in a processor column now have the same   *
   *  hypergraph; we now gather it across rows.                *
   *************************************************************/

  if (nProc_x == 1) {
    /* 
     * Don't need a gather across the row; just set pointers appropriately
     * in shg.
     */
    shg->vindex = col_vindex;
    shg->vedge = col_vedge;
    shg->hindex = col_hindex;
    shg->hvertex = col_hvertex;

    /* Copy vwgt and fixed arrays so shg owns this memory */
    for (i = 0; i < shg->VtxWeightDim*shg->nVtx; i++)
      shg->vwgt[i] = phg->vwgt[i];
    if (hgp->UseFixedVtx)
      for (i = 0; i < shg->nVtx; i++)
        shg->fixed_part[i] = phg->fixed_part[i];
    if (hgp->UsePrefPart)
      for (i = 0; i < shg->nVtx; i++)
        shg->pref_part[i] = phg->pref_part[i];
  }

  else {

    /* Gather info about size within the row */

    MPI_Allgather(&col_nPin, 1, MPI_INT, recv_size, 1, MPI_INT, 
                  phg->comm->row_comm);
  
    tmp = 0;
    for (i = 0; i < nProc_x; i++) 
      tmp += recv_size[i];

    shg->nPins = tmp;
  
    shg->vindex = (int *) ZOLTAN_CALLOC((shg->nVtx+1), sizeof(int));
    shg->vedge = (int *) ZOLTAN_MALLOC(shg->nPins * sizeof(int));
    shg->hindex = (int *) ZOLTAN_CALLOC((shg->nEdge+1), sizeof(int));
    shg->hvertex = (int *) ZOLTAN_MALLOC(shg->nPins * sizeof(int));
   
    if (!shg->vindex || !shg->hindex ||
        (shg->nPins && (!shg->vedge || !shg->hvertex)))
      MEMORY_ERROR;
    
    /* Gather vedge data for all procs in row */
  
    /* SCHEMEA can send local edge numbers; 
       SCHEMEB requires edge LNO to GNO conversion. */
  
    disp[0] = 0;
    for (i = 1; i < nProc_x; i++)
      disp[i] = disp[i-1] + recv_size[i-1];
  
    MPI_Allgatherv(col_vedge, col_nPin, MPI_INT,
                   shg->vedge, recv_size, disp, MPI_INT, phg->comm->row_comm);
  
    /* Gather vindex data for all procs in row */
  
    for (i = 0; i < col_nVtx; i++)
      send_buf[i] = col_vindex[i+1] - col_vindex[i];
  
    /* SCHEMEA can assume a recv for each vertex;
     * SCHEMEB would need to gather the number of vtxs recv'd from each proc. */
  
    for (i = 0; i < nProc_x; i++) 
      each[i] = (int)(phg->dist_x[i+1] - phg->dist_x[i]);

    disp[0] = 0;  /* Can't use dist_x, may not be sizeof(int) */
    for (i = 1; i < nProc_x; i++) 
      disp[i] = disp[i-1] + each[i-1];

    /* SCHEMEA can use phg->dist_x as displacement array;
     * SCHEMEB requires separate displacement array. */

    MPI_Allgatherv(send_buf, col_nVtx, MPI_INT, 
                   shg->vindex, each, disp,
                   MPI_INT, phg->comm->row_comm);

    /* Perform prefix sum on shg->vindex */
    sum = 0;
    for (i = 0; i < shg->nVtx; i++) {
      tmp = shg->vindex[i];
      shg->vindex[i] = sum;
      sum += tmp;
    }
    shg->vindex[shg->nVtx] = sum;
  
    /* Sanity check */
    if (shg->vindex[shg->nVtx] != shg->nPins) {
      printf("%d Sanity check failed:  "
             "shg->vindex %d != nPins %d\n", 
              zz->Proc, shg->vindex[shg->nVtx], shg->nPins);
      exit(-1);
    }
  
    /* Gather fixed array, if any  */
    if (hgp->UseFixedVtx){
  
#ifdef DEBUG_
      uprintf(phg->comm, "Debug in PHG_gather before gather. phg->fixed =");
      for (i=0; i<phg->nVtx; i++){
        printf(" %d ", phg->fixed_part[i]);
      }
      printf("\n");
#endif

      /* Can use the same each array. */
      /* Need to compute new disp array. */
  
      disp[0] = 0;
      for (i = 1; i < nProc_x; i++) {
        disp[i] = disp[i-1] + each[i-1];
      }
      
      MPI_Allgatherv(phg->fixed_part, phg->nVtx, MPI_FLOAT, 
                     shg->fixed_part, each, disp, MPI_FLOAT, phg->comm->row_comm);

#ifdef DEBUG_
      uprintf(phg->comm, "Debug in PHG_gather after gather. shg->fixed =");
      for (i=0; i<shg->nVtx; i++){
        printf(" %d ", shg->fixed_part[i]);
      }
      printf("\n");
#endif
    }
    /* Gather pref part array, if any  */
    if (hgp->UsePrefPart){
      /* Can use the same each array. */
      /* Need to compute new disp array. */
      disp[0] = 0;
      for (i = 1; i < nProc_x; i++) {
        disp[i] = disp[i-1] + each[i-1];
      }
      
      MPI_Allgatherv(phg->pref_part, phg->nVtx, MPI_FLOAT, 
                     shg->pref_part, each, disp, MPI_FLOAT, phg->comm->row_comm);
    }
    
    /* Gather vertex weights, if any. */
    if (shg->VtxWeightDim) {
  
      /* Can use nearly the same each array. */
      /* Need to compute new disp array. */
  
      disp[0] = 0;
      each[0] *= phg->VtxWeightDim;
      for (i = 1; i < nProc_x; i++) {
        each[i] *= phg->VtxWeightDim;
        disp[i] = disp[i-1] + each[i-1];
      }
      
      MPI_Allgatherv(phg->vwgt, phg->nVtx*phg->VtxWeightDim, MPI_FLOAT, 
                     shg->vwgt, each, disp, MPI_FLOAT, phg->comm->row_comm);
    }
  
    Zoltan_HG_Mirror(shg->nVtx, shg->vindex, shg->vedge, 
                     shg->nEdge, shg->hindex, shg->hvertex);

  }  /* End row gather */
  
#ifdef KDDKDD_CHECK
  Zoltan_HG_Print(zz, shg, NULL, stdout, "GatherAfter");/* NULL parts for now;
                                                           add non-NULL later */
  Zoltan_PHG_Plot_2D_Distrib(zz, phg);
  Zoltan_PHG_Plot_2D_Distrib(zz, shg);
#endif

End:

  if (ierr < 0) {
    Zoltan_HG_HGraph_Free(*gathered_hg);
    ZOLTAN_FREE(gathered_hg);
  }

  Zoltan_Multifree(__FILE__, __LINE__, 2, &send_buf, 
                                          &recv_size);

  if (nProc_x > 1 && nProc_y > 1) 
    Zoltan_Multifree(__FILE__, __LINE__, 4, &col_vedge,
                                            &col_vindex,
                                            &col_hvertex,
                                            &col_hindex);
  return ierr;
}

/******************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
