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

/* Layout related functions */

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
Zoltan_Matrix2d_Distribute (ZZ* zz, Zoltan_matrix inmat, /* Cannot be const as we can share it inside outmat */
			    Zoltan_matrix_2d *outmat, int copy)
{
  static char *yo = "Zoltan_Matrix_Build2d";
  int ierr = ZOLTAN_OK;
  int nProc_x, nProc_y;
  int myProc_x, myProc_y;
  int *dist_x=NULL, *dist_y=NULL;
  int frac_x, frac_y;
  int i, j, cnt;
  int *proclist = NULL;
  int nPins, nEdge, nVtx;
  int myProcId;
  int *procptr = NULL;
  Zoltan_Arc *nonzeros= NULL, *sendbuf= NULL;
  int *tmparray=NULL;
  float *tmpwgtarray = NULL;
  int msg_tag = 1021982;
  int offset;
  ZOLTAN_COMM_OBJ *plan;
  MPI_Comm communicator = MPI_COMM_NULL;
  int nProc;
  int *yGNO = NULL;
  int *pinGNO = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  memcpy(&outmat->mtx, &inmat, sizeof(Zoltan_matrix));
  if(copy) {
    /* TODO: We need to copy the arrays also */
    Zoltan_Matrix_Reset (&outmat->mtx);
    /* Copy also directories */
    outmat->mtx.ddX = Zoltan_DD_Copy (inmat.ddX);
    if (inmat.ddY == inmat.ddX)
      outmat->mtx.ddY = outmat->mtx.ddX;
    else
      outmat->mtx.ddY = Zoltan_DD_Copy (inmat.ddY);
  }

  communicator = outmat->comm->Communicator;
  nProc = outmat->comm->nProc;
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

  /* Do we need to redistribute the graph ? *
   * redistribute if the matrix is defined using hypergraph queries or if there is
   * more than one proc on x axis.
   */
  if (inmat.redist || (nProc_x> 1) || !inmat.opts.keep_distribution) {
    /* Do a redistribution by "slice" on X and Y */
    frac_x = (float) inmat.globalX / (float) nProc_x;
    for (i = 1; i < nProc_x; i++)
      dist_x[i] = (int) (i * frac_x);
    dist_x[nProc_x] = inmat.globalX;

    frac_y = (float) inmat.globalY / (float) nProc_y;
    for (i = 1; i < nProc_y; i++)
      dist_y[i] = (int) (i * frac_y);
    dist_y[nProc_y] = inmat.globalY;
  }
  else {
    /* No redistribution, code is only here to finish "symmetrization" */
    int nY;
/*     int flag; */
/*     int general_flag; */

    nY = inmat.nY;
    /* This code only works for a linear distribution on Y */
    MPI_Allgather (&inmat.nY_ori, 1, MPI_INT, dist_y + 1, 1, MPI_INT, communicator);
    dist_y[0] = 0;
    for (i = 1 ; i <= nProc_y ; i++) {
      dist_y[i] += dist_y[i-1];
    }
    dist_x[0] = 0; dist_x[1] = inmat.globalX;

    /* TODO: we need more to keep the same order in the permutation ! */
/*     /\* Perhaps we have to insure that the permutation is correct ? *\/ */
/*     /\* I will check to avoid a call to really permute things *\/ */
/*     for (i = 0, flag=0 ; i < nY ; ++i) */
/*       flag |= ((inmat.yGNO[i] < dist_y[myProc_y]) || (inmat.yGNO[i] >= dist_y[myProc_y+1])); */
/*     MPI_Allreduce(&flag, &general_flag, 1, MPI_INT, MPI_MAX, communicator); */
/*     if (general_flag) { /\* We have to compute the "permutation" *\/ */
/*       int offset; */

/*       tmparray = (int*)ZOLTAN_MALLOC(nY *sizeof(int)); */
/*       if (tmparray == NULL) MEMORY_ERROR; */
/*       offset = dist_y[myProc_y]; */
/*       for (i=0 ; i < nY ; ++i) */
/* 	tmparray[i] = offset + i; */
/*       Zoltan_Matrix_Permute(zz, &outmat->mtx, tmparray); */
/*       ZOLTAN_FREE(&tmparray); */
/*     } */
  }

  /* myProc_y and myProc_x can be -1 when we use a 2D decomposition.
   * Some processor may be excluded from the 2D communicator; for it,
   * myProc_y and myProc_x == -1. */
  nEdge = (myProc_y >= 0 ? dist_y[myProc_y+1] - dist_y[myProc_y] : 0);
  nVtx  = (myProc_x >= 0 ? dist_x[myProc_x+1] - dist_x[myProc_x] : 0);

  /* Construct procptr, the array to know what are the correct procnumber in
   * the big communicator.
   */
  tmparray = (int*)ZOLTAN_CALLOC (nProc, sizeof(int));
  /* Trick : +1 to avoid a test in the loop */
  procptr = (int*)ZOLTAN_CALLOC (nProc_x*nProc_y+1, sizeof(int));
  if (tmparray == NULL || procptr == NULL) MEMORY_ERROR;
  procptr ++;

  myProcId = (myProc_x >= 0)?(myProc_y*nProc_x+myProc_x):-1;
  MPI_Allgather (&myProcId, 1, MPI_INT, tmparray, 1, MPI_INT, communicator);

  for (i=0 ; i < nProc ; ++i)
    procptr[tmparray[i]] = i;            /* Don't worry about the -1 */
  ZOLTAN_FREE(&tmparray);

  offset = dist_y[myProc_y];

  ierr = Zoltan_Matrix_Remove_Duplicates(zz, outmat->mtx, &outmat->mtx);
  nPins = outmat->mtx.nPins;

  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution.
   */
  proclist = (int *)ZOLTAN_MALLOC(nPins *sizeof(int));
  sendbuf = (Zoltan_Arc*) ZOLTAN_MALLOC(nPins * sizeof(Zoltan_Arc));

  if ((nPins >0) && (proclist == NULL || sendbuf == NULL)) MEMORY_ERROR;

  yGNO = outmat->mtx.yGNO;
  pinGNO = outmat->mtx.pinGNO;

  cnt = 0;
  for (i = 0; i < outmat->mtx.nY; i++) {
    int edge_gno=-1, edge_Proc_y=-1;
    /* processor row for the edge */
    edge_gno = yGNO[i];
    edge_Proc_y = EDGE_TO_PROC_Y(outmat, edge_gno);

    for (j = outmat->mtx.ystart[i]; j < outmat->mtx.yend[i]; j++) {
      int vtx_gno=-1, vtx_Proc_x=-1;
      /* processor column for the vertex */
      vtx_gno = pinGNO[j];
      vtx_Proc_x = VTX_TO_PROC_X(outmat, vtx_gno);

      proclist[cnt] = procptr[edge_Proc_y * nProc_x + vtx_Proc_x];
      sendbuf[cnt].yGNO = edge_gno;
      sendbuf[cnt].pinGNO = vtx_gno;
      sendbuf[cnt].offset = j;
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
  ierr = Zoltan_Comm_Create(&plan, cnt, proclist, communicator, msg_tag, &outmat->mtx.nPins);
  ZOLTAN_FREE(&proclist);

  nonzeros = (Zoltan_Arc *) ZOLTAN_MALLOC((outmat->mtx.nPins+nEdge) * sizeof(Zoltan_Arc));
  if (outmat->mtx.nPins && nonzeros == NULL) MEMORY_ERROR;

  msg_tag--;
  Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, sizeof(Zoltan_Arc),
		 (char *) nonzeros);
  ZOLTAN_FREE(&sendbuf);

  if (outmat->mtx.pinwgtdim) { /* We have to take care about weights */
    tmpwgtarray = (float*) ZOLTAN_MALLOC(outmat->mtx.nPins*outmat->mtx.pinwgtdim*sizeof(float));
    if (outmat->mtx.nPins && tmpwgtarray == NULL) MEMORY_ERROR;

    msg_tag--;
    Zoltan_Comm_Do_Post(plan, msg_tag, (char *) outmat->mtx.pinwgt, outmat->mtx.pinwgtdim*sizeof(float),
		   (char *) tmpwgtarray);
    for (i = 0 ; i < outmat->mtx.nPins ; ++i) { /* Fill indirection structure */
      nonzeros[i].offset = i;
    }
    Zoltan_Comm_Do_Wait(plan, msg_tag, (char *) outmat->mtx.pinwgt, outmat->mtx.pinwgtdim*sizeof(float),
		   (char *) tmpwgtarray);
    ZOLTAN_FREE(&outmat->mtx.pinwgt);
  }
  Zoltan_Comm_Destroy(&plan);

  /* Unpack the non-zeros received. */

  /* First, add fake edges for empty vertices */
  for (i = 0 ; i < nEdge ; ++i) {
    int offset = dist_y[myProc_y];
    j=outmat->mtx.nPins + i;
    nonzeros[j].yGNO = offset + i;
    nonzeros[j].pinGNO = -1;
    nonzeros[j].offset = 0;
  }
  outmat->mtx.ystart = (int *) ZOLTAN_REALLOC(outmat->mtx.ystart, (nEdge + 1)*sizeof(int));
  if (outmat->mtx.ystart == NULL) MEMORY_ERROR;
  outmat->mtx.yend = outmat->mtx.ystart + 1;
  outmat->mtx.pinGNO = (int *) ZOLTAN_REALLOC(outmat->mtx.pinGNO, (outmat->mtx.nPins) * sizeof(int));
  if (outmat->mtx.nPins && (outmat->mtx.pinGNO == NULL)) MEMORY_ERROR;
  outmat->mtx.pinwgt = (float *) ZOLTAN_REALLOC(outmat->mtx.pinwgt,
	     (outmat->mtx.nPins * outmat->mtx.pinwgtdim) * sizeof(float));
  if (outmat->mtx.nPins && outmat->mtx.pinwgtdim && outmat->mtx.pinwgt == NULL)
    MEMORY_ERROR;

  /* TODO: Ensure that we don't have empty vertex */
  Zoltan_Matrix_Remove_DupArcs(zz, outmat->mtx.nPins + nEdge, (Zoltan_Arc*)nonzeros, tmpwgtarray,
			       &outmat->mtx);

 End:
  if (procptr != NULL) {
    procptr -= 1;
    ZOLTAN_FREE(&procptr);
  }
  ZOLTAN_FREE(&proclist);
  ZOLTAN_FREE(&sendbuf);
  ZOLTAN_FREE(&nonzeros);
  ZOLTAN_FREE(&tmparray);
  ZOLTAN_FREE(&tmpwgtarray);

  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}


#ifdef __cplusplus
}
#endif
