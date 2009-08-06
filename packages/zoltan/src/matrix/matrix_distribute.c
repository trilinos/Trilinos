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


int
Zoltan_Matrix2d_Distribute (ZZ* zz, const Zoltan_matrix inmat,
			    Zoltan_matrix_2d *outmat)
{
  static char *yo = "Zoltan_Matrix_Build2d";
  int ierr = ZOLTAN_OK;
  int nProc_x, nProc_y;
  int myProc_x, myProc_y;
  int *dist_x, *dist_y;
  int frac_x, frac_y;
  int i, j, cnt;
  int *proclist = NULL, *sendbuf = NULL;
  int nPins, nEdge, nVtx;
  int myProcId;
  int *procptr = NULL, *tmparray = NULL;
  int msg_tag = 1021982;
  int *nonzeros=NULL;
  ZOLTAN_COMM_OBJ *plan;
  int final_output = 0;

  ZOLTAN_TRACE_ENTER(zz, yo);

  memset(&outmat->mtx, 0, sizeof(Zoltan_matrix));

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

  nPins = inmat.nPins;
  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution.
   */
  proclist = (int *)ZOLTAN_MALLOC(nPins *sizeof(int));
  sendbuf = (int *) ZOLTAN_MALLOC(nPins * 2 * sizeof(int));

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
      sendbuf[2*cnt] = edge_gno;
      sendbuf[2*cnt+1] = vtx_gno;
      cnt++;
    }
  }
  procptr -= 1;
  ZOLTAN_FREE(&procptr);

  /*
   * Send pins to their target processors.
   * They become non-zeros in the 2D data distribution.
   */

  msg_tag--;
  ierr = Zoltan_Comm_Create(&plan, cnt, proclist, zz->Communicator, msg_tag, &outmat->mtx.nPins);

  ZOLTAN_FREE(&proclist);

  if (outmat->mtx.nPins) {
    nonzeros = (int *) ZOLTAN_MALLOC(outmat->mtx.nPins * 2 * sizeof(int));
    if (!nonzeros) MEMORY_ERROR;
  }

  msg_tag--;
  Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(int),
		 (char *) nonzeros);

  ZOLTAN_FREE(&sendbuf);

  /* Not yet destruction, I need to transfert pinWgt */
  Zoltan_Comm_Destroy(&plan);

  /* Unpack the non-zeros received. */

  tmparray = (int *) ZOLTAN_CALLOC(nEdge + 1, sizeof(int));
  outmat->mtx.ystart = (int *) ZOLTAN_CALLOC(nEdge + 1 , sizeof(int));
  outmat->mtx.pinGNO = (int *) ZOLTAN_MALLOC((outmat->mtx.nPins) * sizeof(int));

  if (!tmparray || !outmat->mtx.ystart || ((outmat->mtx.nPins && !outmat->mtx.pinGNO)))
    MEMORY_ERROR;

  /* Count the number of nonzeros per hyperedge */
  for (i = 0; i < outmat->mtx.nPins; i++) {
    j = EDGE_GNO_TO_LNO(outmat, nonzeros[2*i]);
    tmparray[j]++;
  }

  outmat->mtx.yend= outmat->mtx.ystart + 1; /* Keep in compact mode, no new edges */
  /* Compute prefix sum to represent hindex correctly. */
  for (i = 0; i < nEdge; i++)  {
    outmat->mtx.ystart[i+1] = outmat->mtx.ystart[i] + tmparray[i];
    tmparray[i] = 0;
  }

  for (i = 0; i < outmat->mtx.nPins; i++) {
    j = EDGE_GNO_TO_LNO(outmat, nonzeros[2*i]);
    outmat->mtx.pinGNO[outmat->mtx.ystart[j]+tmparray[j]] = VTX_GNO_TO_LNO(outmat, nonzeros[2*i+1]);
    tmparray[j]++;
  }

  ZOLTAN_FREE(&nonzeros);
  ZOLTAN_FREE(&tmparray);

  outmat->mtx.nY = nEdge;

/*   /\* Create ObjGNO array from the dist informations *\/ */
/*   if (myProc_x >= 0) { */
/*     outmat->mtx.objGNO = (int*) ZOLTAN_MALLOC(nVtx*sizeof(int)); */
/*     if (nVtx && !outmat->mtx.objGNO) MEMORY_ERROR; */

/*     for (i = 0 ; i < nVtx ; ++i) */
/*       outmat->mtx.xGNO[i] = dist_x[myProc_x] + i; */
/*   } */

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
