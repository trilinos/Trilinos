/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2009 Sandia National Laboratories.                          *
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

#include <math.h>
#include "zz_const.h"
#include "zz_util_const.h"
#include "zoltan_dd.h"
#include "phg.h"
#include "matrix.h"

static int
compar_couple (const int* e1, const int* e2)
{
  if (e1[0] == e2[0])
    return (e1[1] - e1[1]);
  return (e1[0] - e1[0]);
}


int Zoltan_Distribute_Square (ZZ * zz, PHGComm *layout)
{
  return Zoltan_Distribute_layout(zz, NULL, 0, zz->Num_Proc-1, -1, -1, layout);
}

int Zoltan_Distribute_LinearY (ZZ * zz, PHGComm *layout)
{
  return Zoltan_Distribute_layout(zz, NULL, 0, zz->Num_Proc-1, 1, zz->Num_Proc, layout);
}

int
Zoltan_Distribute_layout (ZZ *zz, const PHGComm * const inlayout,
			  int loRank, int hiRank,
			  int reqx, int reqy,
			  PHGComm *outlayout)
{
  MPI_Group allgrp, newgrp;
  int *ranks;
  MPI_Comm  nmpicomm;
  MPI_Comm  ompicomm;
  int myProc;
  int i;
  int nProc;

  ompicomm = (inlayout != NULL)?inlayout->Communicator:zz->Communicator;
  myProc = (inlayout != NULL)?inlayout->myProc:zz->Proc;
  nProc= (inlayout != NULL)?inlayout->nProc:zz->Num_Proc;
  if (((reqx != 1) && (reqy != 1) && (nProc > 3)) && Zoltan_PHG_isPrime(nProc)) nProc--;

  Zoltan_PHGComm_Init(outlayout);

  /* create a new communicator for procs[lo..hi] */

  MPI_Comm_group(ompicomm, &allgrp);
  ranks = (int *) ZOLTAN_MALLOC(nProc * sizeof(int));
  for (i=loRank; i<=hiRank; ++i)
    ranks[i-loRank] = i;

  MPI_Group_incl(allgrp, nProc, ranks, &newgrp);
  MPI_Comm_create(ompicomm, newgrp, &nmpicomm);
  MPI_Group_free(&newgrp);
  MPI_Group_free(&allgrp);
  ZOLTAN_FREE(&ranks);

  return (Zoltan_PHG_Set_2D_Proc_Distrib(zz, nmpicomm,
					myProc-loRank, nProc,
					reqx, reqy, outlayout));
}


/* if !copy, inmat is not usable after this call */
/* for pin wgt, we may do a "savage cast" to avoid padding problem */

int
Zoltan_Matrix2d_Distribute (ZZ* zz, const Zoltan_matrix inmat,
			    Zoltan_matrix_2d *outmat, int copy)
{
  static char *yo = "Zoltan_Matrix_Build2d";
  int ierr = ZOLTAN_OK;
  int nProc_x, nProc_y;
  int myProc_x, myProc_y;
  int *dist_x=NULL, *dist_y=NULL;
  int frac_x, frac_y;
  int i, j, cnt;
  int *proclist = NULL, *sendbuf = NULL;
  int nPins, nEdge, nVtx;
  int myProcId;
  int *procptr = NULL, *tmparray = NULL;
  int msg_tag = 1021982;
  int *nonzeros=NULL;
  int offset;
  int prev_x, prev_y;
  ZOLTAN_COMM_OBJ *plan;
  int elem_size = 2;

  ZOLTAN_TRACE_ENTER(zz, yo);

  memcpy(&outmat->mtx, &inmat, sizeof(Zoltan_matrix));
  if(copy) {
    Zoltan_Matrix_Reset (&outmat->mtx);
    /* Copy also directories */
    outmat->mtx.ddX = Zoltan_DD_Copy (inmat.ddX);
    if (inmat.ddY == inmat.ddX)
      outmat->mtx.ddY = outmat->mtx.ddX;
    else
      outmat->mtx.ddY = Zoltan_DD_Copy (inmat.ddY);
  }

  /****************************************************************************************
   * Compute the distribution of vertices and edges to the 2D data distribution's processor
   * columns and rows. For now, these distributions are described by arrays dist_x  and dist_y;
   * in the future, we may prefer a hashing function mapping GIDs to processor columns and rows.
   ****************************************************************************************/

  nProc_x = outmat->comm->nProc_x;
  nProc_y = outmat->comm->nProc_y;
  myProc_x = outmat->comm->myProc_x;
  myProc_y = outmat->comm->myProc_y;

  outmat->dist_x = dist_x = (int *) ZOLTAN_CALLOC((nProc_x+1), sizeof(int));
  outmat->dist_y = dist_y = (int *) ZOLTAN_CALLOC((nProc_y+1), sizeof(int));

  if (!dist_x || !dist_y) MEMORY_ERROR;

  frac_x = (float) inmat.globalX / (float) nProc_x;
  for (i = 1; i < nProc_x; i++)
    dist_x[i] = (int) (i * frac_x);
  dist_x[nProc_x] = inmat.globalX;

  frac_y = (float) inmat.globalY / (float) nProc_y;
  for (i = 1; i < nProc_y; i++)
    dist_y[i] = (int) (i * frac_y);
  dist_y[nProc_y] = inmat.globalY;

  /* myProc_y and myProc_x can be -1 when we use a 2D decomposition.
   * Some processor may be excluded from the 2D communicator; for it,
   * myProc_y and myProc_x == -1. */
  nEdge = (myProc_y >= 0 ? dist_y[myProc_y+1] - dist_y[myProc_y] : 0);
  nVtx  = (myProc_x >= 0 ? dist_x[myProc_x+1] - dist_x[myProc_x] : 0);

  /* Construct procptr, the array to know what are the correct procnumber in
   * the big communicator.
   */
  tmparray = (int*)ZOLTAN_CALLOC (zz->Num_Proc, sizeof(int));
  /* Trick : +1 to avoid a test in the loop */
  procptr = (int*)ZOLTAN_CALLOC (nProc_x*nProc_y+1, sizeof(int));
  if (tmparray == NULL || procptr == NULL) MEMORY_ERROR;
  procptr ++;

  myProcId = (myProc_x >= 0)?(myProc_y*nProc_x+myProc_x):-1;
  MPI_Allgather (&myProcId, 1, MPI_INT, tmparray, 1, MPI_INT, zz->Communicator);

  for (i=0 ; i < zz->Num_Proc ; ++i)
    procptr[tmparray[i]] = i;            /* Don't worry about the -1 */
  ZOLTAN_FREE(&tmparray);

  offset = dist_y[myProc_y];
  nPins = inmat.nPins;
  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution.
   */
  if (inmat.pinwgtdim > 0)
    elem_size = 3;
  proclist = (int *)ZOLTAN_MALLOC(nPins *sizeof(int));
  sendbuf = (int *) ZOLTAN_MALLOC(nPins * elem_size * sizeof(int));

  if ((nPins >0) && (proclist == NULL || sendbuf == NULL)) MEMORY_ERROR;

  cnt = 0;
  for (i = 0; i < inmat.nY; i++) {
    int edge_gno=-1, edge_Proc_y=-1;
    /* processor row for the edge */
    edge_gno = inmat.yGNO[i];
    edge_Proc_y = EDGE_TO_PROC_Y(outmat, edge_gno);

    for (j = inmat.ystart[i]; j < inmat.yend[i]; j++) {
      int vtx_gno=-1, vtx_Proc_x=-1;
      /* processor column for the vertex */
      vtx_gno = inmat.pinGNO[j];
      vtx_Proc_x = VTX_TO_PROC_X(outmat, vtx_gno);

      proclist[cnt] = procptr[edge_Proc_y * nProc_x + vtx_Proc_x];
      sendbuf[elem_size*cnt] = edge_gno;
      sendbuf[elem_size*cnt+1] = vtx_gno;
      if (inmat.pinwgtdim > 0)
	sendbuf[elem_size*cnt+2] = j; /* Where is the corresponding pinwgt */
      cnt++;
    }
  }
  procptr -= 1;
  ZOLTAN_FREE(&procptr);

  /* May be interesting to remove local duplicates now ! */
  /* Will require more memory as we have to copy "sendbuf" */
  /* However we can add the weights only now */

  /*
   * Send pins to their target processors.
   * They become non-zeros in the 2D data distribution.
   */

  msg_tag--;
  ierr = Zoltan_Comm_Create(&plan, cnt, proclist, zz->Communicator, msg_tag, &outmat->mtx.nPins);

  ZOLTAN_FREE(&proclist);

  if (outmat->mtx.nPins) {
    nonzeros = (int *) ZOLTAN_MALLOC(outmat->mtx.nPins * elem_size * sizeof(int));
    if (!nonzeros) MEMORY_ERROR;
  }

  msg_tag--;
  Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, elem_size*sizeof(int),
		 (char *) nonzeros);

  ZOLTAN_FREE(&sendbuf);

  /* Not yet destruction, I need to transfert pinWgt */
  Zoltan_Comm_Destroy(&plan);

  /* Unpack the non-zeros received. */

  outmat->mtx.ystart = (int *) ZOLTAN_REALLOC(outmat->mtx.ystart, (nEdge + 1)*sizeof(int));
  outmat->mtx.pinGNO = (int *) ZOLTAN_REALLOC(outmat->mtx.pinGNO, (outmat->mtx.nPins) * sizeof(int));

  if ( outmat->mtx.ystart == NULL || ((outmat->mtx.nPins && outmat->mtx.pinGNO == NULL)))
     MEMORY_ERROR;

  /* Sort Edges: Allow to easily remove duplicates */
  /* inconvenient: make pin weight association harder */
  qsort ((void*)nonzeros, outmat->mtx.nPins, sizeof(int)*2,
	 (int (*)(const void*,const void*))compar_couple);

  if (!copy && inmat.yend != inmat.ystart + 1)
    ZOLTAN_FREE(&outmat->mtx.yend);

  outmat->mtx.yend= outmat->mtx.ystart + 1; /* Keep in compact mode, no new edges */
  outmat->mtx.ystart[0] = 0;

  for (i = 0, prev_x=-1, prev_y=-1; i < outmat->mtx.nPins; i++) {
    int x = VTX_GNO_TO_LNO(outmat, nonzeros[elem_size*i+1]);
    int y = EDGE_GNO_TO_LNO(outmat, nonzeros[elem_size*i]);
    /* Why do not deal with GNO instead of LNO ? << We own only on the row */

    if (y == prev_y && x == prev_x) /* If it is a duplicate edge, skip ! */
      continue;
    for (j=prev_y+1 ; j <= y ; j++)               /* Some y may be empty */
      outmat->mtx.yend[j] = outmat->mtx.ystart[j]; /* Compact mode */

    prev_y = y; prev_x = x;
    outmat->mtx.pinGNO[outmat->mtx.yend[y]++] = x;
  }
  for (j=prev_y+1 ; j < nEdge ; j++)               /* Some y may be empty */
    outmat->mtx.yend[j] = outmat->mtx.ystart[j]; /* Compact mode */

  ZOLTAN_FREE(&nonzeros);

  outmat->mtx.nPins = outmat->mtx.yend[nEdge - 1];
  /* Try to minimize memory */
  outmat->mtx.pinGNO = (int *) ZOLTAN_REALLOC(outmat->mtx.pinGNO,
			       (outmat->mtx.nPins) * sizeof(int));

  outmat->mtx.nY = nEdge;

  /* Now construct yGNO array */
  outmat->mtx.yGNO = (int *)ZOLTAN_REALLOC(outmat->mtx.yGNO, nEdge*sizeof(int));
  if (nEdge && outmat->mtx.yGNO == NULL) MEMORY_ERROR;
  for (i = 0 ; i < nEdge; ++i) {
    outmat->mtx.yGNO[i] = offset + i;
  }

 End:
  if (procptr != NULL) {
    procptr -= 1;
    ZOLTAN_FREE(&procptr);
  }
  ZOLTAN_FREE(&proclist);
  ZOLTAN_FREE(&sendbuf);
  ZOLTAN_FREE(&nonzeros);
  ZOLTAN_FREE(&tmparray);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

#ifdef __cplusplus
}
#endif
