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
#include <limits.h>

/****************************************************************************/

#define NIDX 3
#define VIDX 0
#define EIDX 1
#define PIDX 2

#define MEMORY_ERROR {ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
                      ierr = ZOLTAN_MEMERR;                              \
                      goto End;                                          \
                     }

/****************************************************************************/
int Zoltan_PHG_Gather_To_All_Procs(
  ZZ *zz, 
  PHGraph *phg,           /* Input:   Local part of distributed hypergraph */
  PHGraph **gathered_hg   /* Output:  combined hypergraph combined to proc */
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
int i, tmpv, tmpp;
int *each = NULL, 
    *disp = NULL;        /* Size and displacement arrays for MPI_Allgatherv */
int *send_buf = NULL;    /* Buffer of values to be sent */
int send_size;           /* Size of buffer send_buf */
int *col_vedge = NULL;   /* vedge array for the proc-column hypergraph */
int *col_vindex = NULL;  /* vindex array for the proc-column hypergraph */
int *col_hvertex = NULL; /* hvertex array for the proc-column hypergraph */
int *col_hindex = NULL;  /* hindex array for the proc-column hypergraph */
int col_size[NIDX];      /* array of nvtx, nedge & nnz for the proc-column
                            hypergraph.  */
int size[NIDX];          /* array of nvtx, nedge & nnz for the input
                            parallel hypergraph.  */

int *recv_size = NULL;   /* nvtx, nedge, & nnz for each proc in col or row */

PHGraph *shg;            /* Pointer to the serial hypergraph to be
                            returned by this function. */
PHGComm scomm;           /* Serial PHGComm for use by shg. */

int myProc_x = phg->comm->myProc_x;
int nProc_x = phg->comm->nProc_x;
int nProc_y = phg->comm->nProc_y;
int max_nProc_xy;        /* Max of nProc_x, nProc_y */
int first_vtx = phg->dist_x[myProc_x];


  /******************************************************************
   *  0. Allocate the hypergraph to be returned. 
   *  Set values that we already know. 
   ******************************************************************/

  shg = *gathered_hg = (PHGraph *) ZOLTAN_MALLOC(sizeof(PHGraph));
  if (!shg) MEMORY_ERROR;

  Zoltan_PHG_PHGraph_Init(shg);
  shg->nVtx = phg->dist_x[nProc_x];
  shg->nEdge = phg->dist_y[nProc_y];
  MPI_Allreduce(&(phg->nNonZero), &(shg->nNonZero), 1, MPI_INT, MPI_SUM,
                zz->Communicator);

  shg->vindex = (int *) ZOLTAN_CALLOC((shg->nVtx+1), sizeof(int));
  shg->vedge = (int *) ZOLTAN_MALLOC(shg->nNonZero * sizeof(int));
  shg->hindex = (int *) ZOLTAN_CALLOC((shg->nEdge+1), sizeof(int));
  shg->hvertex = (int *) ZOLTAN_MALLOC(shg->nNonZero * sizeof(int));

  shg->dist_x = (int *) ZOLTAN_MALLOC(2 * sizeof(int));
  shg->dist_y = (int *) ZOLTAN_MALLOC(2 * sizeof(int));
  shg->dist_x[0] = shg->dist_y[0] = 0;
  shg->dist_x[1] = shg->nVtx;
  shg->dist_x[1] = shg->nEdge;

  /* KDDKDD -- Need allocation error check here. */

  scomm.Communicator = MPI_COMM_SELF;
  scomm.Proc = 0;
  scomm.Num_Proc = 1;
  scomm.nProc_x = 1;
  scomm.nProc_y = 1;
  scomm.myProc_x = 0;
  scomm.myProc_y = 0;
  scomm.row_comm = MPI_COMM_SELF;
  scomm.col_comm = MPI_COMM_SELF;
  shg->comm = &scomm;

  shg->EdgeWeightDim = phg->EdgeWeightDim;
  shg->VtxWeightDim = phg->VtxWeightDim;
  
/* KDDKDD WEIGHTS ARE NOT YET DONE!!!! */
/* KDDKDD
  shg->vwgt = vwgt;
  shg->ewgt = ewgt;
*/

  if (!shg->vindex || !shg->hindex || 
      (shg->nNonZero && (!shg->vedge || !shg->hvertex))) 
    MEMORY_ERROR;
  

  /*************************************************************
   *  1. Gather all non-zeros for vertices in processor column *
   *************************************************************/

  /* Gather local size info for each proc in column */

  max_nProc_xy = MAX(nProc_x, nProc_y); 
            /* Use MAX so can reuse recv_size in step 2 */
  recv_size = (int *) ZOLTAN_MALLOC(NIDX * max_nProc_xy * sizeof(int));
  if (!recv_size) MEMORY_ERROR;

  size[VIDX] = phg->nVtx;
  size[EIDX] = phg->nEdge;
  size[PIDX] = phg->nNonZero;

  MPI_Allgather(size, 3, MPI_INT, recv_size, 3, MPI_INT, 
                phg->comm->col_comm);
  
  /* Compute number of vtx, edge, and nnz in column */
  col_size[VIDX] = phg->dist_x[myProc_x+1] - phg->dist_x[myProc_x];
  col_size[EIDX] = col_size[PIDX] = 0;
  for (i = 0; i < nProc_y; i++) {
    col_size[EIDX] += recv_size[NIDX*i + EIDX];
    col_size[PIDX] += recv_size[NIDX*i + PIDX];
  }
  
  /* Allocate arrays for column hypergraph */
  col_hindex = (int *) ZOLTAN_CALLOC((col_size[EIDX]+1), sizeof(int));
  col_hvertex = (int *) ZOLTAN_MALLOC(col_size[PIDX] * sizeof(int));

  col_vindex = (int *) ZOLTAN_CALLOC((col_size[VIDX]+1), sizeof(int));
  col_vedge = (int *) ZOLTAN_MALLOC(col_size[PIDX] * sizeof(int));

  /* Allocate arrays for use in gather operations */

  send_size = MAX(col_size[PIDX],col_size[EIDX]+1);
  send_buf = (int *) ZOLTAN_MALLOC(send_size * sizeof(int));

  each = (int *) ZOLTAN_MALLOC(2 * max_nProc_xy * sizeof(int));
  disp = each + max_nProc_xy;

  if (!each || !send_buf || !col_vindex || !col_hindex || 
      (col_size[PIDX] && (!col_vedge || !col_hvertex)))
    MEMORY_ERROR;
  
  /* Gather hvertex data for all procs in column */

  for (i = 0; i < phg->nNonZero; i++)
    send_buf[i] = VTX_LNO_TO_GNO(phg, phg->hvertex[i]);

  for (i = 0; i < nProc_y; i++)
    each[i] = recv_size[NIDX * i + PIDX];

  disp[0] = 0;
  for (i = 1; i < nProc_y; i++)
    disp[i] = disp[i-1] + each[i-1];

  MPI_Allgatherv(send_buf, phg->nNonZero, MPI_INT,
                 col_hvertex, each, disp, MPI_INT, phg->comm->col_comm);

  /* 
   * Convert vertex GNOs to LNOs with respect to the column hypergraph
   * Assume LNOs in densely within range 
   * phg->dist_x[myProc_x]:(phg->distx[myProc_x+1]-1).
   * This dense distribution is correct regardless of whether we use
   * Scheme A or Scheme B for parallel distribution.
   * Conversion is needed for Zoltan_PHG_Mirror to work correctly;
   * it assumes a local numbering system.
   */

  for (i = 0; i < col_size[PIDX]; i++)
    col_hvertex[i] = col_hvertex[i] - first_vtx;

  /* Gather hindex data for all procs in column */

  for (i = 0; i < phg->nEdge; i++)
    send_buf[i] = phg->hindex[i+1] - phg->hindex[i];

  for (i = 0; i < nProc_y; i++) 
    each[i] = recv_size[NIDX * i + EIDX] + 1;

  disp[0] = 0;
  for (i = 1; i < nProc_y; i++)
    disp[i] = disp[i-1] + each[i-1];

  MPI_Allgatherv(send_buf, phg->nEdge, MPI_INT, 
                 col_hindex, each, disp, MPI_INT, phg->comm->col_comm);

  /* Perform prefix sum on col_hindex */

  for (i = 1; i <= col_size[EIDX]; i++)
    col_hindex[i] += col_hindex[i-1];
  col_hindex[0] = 0;

  /* Sanity check */
  if (col_hindex[col_size[EIDX]] != col_size[PIDX]) {
    printf("%d Sanity check failed:  "
           "col_hindex[col_size[EIDX]] %d != col_size[PIDX] %d\n", 
            zz->Proc, col_hindex[col_size[EIDX]], col_size[PIDX]);
    exit(-1);
  }

  Zoltan_PHG_Mirror(col_size[EIDX], col_hindex, col_hvertex, 
                    col_size[VIDX], col_vindex, col_vedge);


  /*************************************************************
   *  2. Gather all non-zeros for edges in processor rows      *
   *  All processors in a processor column now have the same   *
   *  hypergraph; we now gather it across rows.                *
   *************************************************************/

  /* Gather info about size within the row */

  MPI_Allgather(col_size, 3, MPI_INT, recv_size, 3, MPI_INT, 
                phg->comm->row_comm);

  /* Sanity checks  KDDKDD remove later */
  tmpv = tmpp = 0;
  for (i = 0; i < nProc_x; i++) {
    tmpv += recv_size[NIDX*i + VIDX];
    tmpp += recv_size[NIDX*i + PIDX];
  }

  if (tmpp != shg->nNonZero) {
    printf("%d sum_size[PIDX] %d != shg->nNonZero %d\n", 
            zz->Proc, tmpp, shg->nNonZero);
    exit(-1);
  }

  if (tmpv != shg->nVtx) {
    printf("%d sum_size[VIDX] %d != shg->nVtx %d\n", 
            zz->Proc, tmpv, shg->nVtx);
    exit(-1);
  }
  /* End sanity checks */
  
  /* Gather vedge data for all procs in row */

  for (i = 0; i < col_size[PIDX]; i++)
    send_buf[i] = EDGE_LNO_TO_GNO(phg, col_vedge[i]);

  for (i = 0; i < nProc_x; i++)
    each[i] = recv_size[NIDX * i + PIDX];

  disp[0] = 0;
  for (i = 1; i < nProc_x; i++)
    disp[i] = disp[i-1] + each[i-1];

  MPI_Allgatherv(send_buf, col_size[PIDX], MPI_INT,
                 shg->vedge, each, disp, MPI_INT, phg->comm->row_comm);

  /* Gather vindex data for all procs in row */

  for (i = 0; i < col_size[VIDX]; i++)
    send_buf[i] = col_vindex[i+1] - col_vindex[i];

  for (i = 0; i < nProc_x; i++) 
    each[i] = recv_size[NIDX * i + VIDX] + 1;

  disp[0] = 0;
  for (i = 1; i < nProc_x; i++)
    disp[i] = disp[i-1] + each[i-1];

  MPI_Allgatherv(send_buf, col_size[VIDX], MPI_INT, 
                 shg->vindex, each, disp, MPI_INT, phg->comm->row_comm);

  /* Perform prefix sum on shg->vindex */
  for (i = 1; i <= shg->nVtx; i++)
    shg->vindex[i] += shg->vindex[i-1];
  shg->vindex[0] = 0;

  /* Sanity check */
  if (shg->vindex[shg->nVtx] != shg->nNonZero) {
    printf("%d Sanity check failed:  "
           "shg->vindex %d != nNonZero %d\n", 
            zz->Proc, shg->vindex[shg->nVtx], shg->nNonZero);
    exit(-1);
  }

  Zoltan_PHG_Mirror(shg->nVtx, shg->vindex, shg->vedge, 
                    shg->nEdge, shg->hindex, shg->hvertex);

  Zoltan_PHG_Plot_2D_Distrib(zz, phg);
  Zoltan_PHG_Plot_2D_Distrib(zz, shg);

End:

  if (ierr < 0) {
    Zoltan_PHG_HGraph_Free(*gathered_hg);
    ZOLTAN_FREE(gathered_hg);
  }

  Zoltan_Multifree(__FILE__, __LINE__, 7, &each,
                                          &send_buf, 
                                          &col_vedge,
                                          &col_vindex,
                                          &col_hvertex,
                                          &col_hindex,
                                          &recv_size);
  return ierr;
}

/******************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
