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


int Zoltan_Distribute_Origin(int edge_gno, int vtx_gno, void* data, int *part_y)
{
  Zoltan_matrix_2d *mat;

  mat = (Zoltan_matrix_2d*) data;
  *part_y = (int)floor((double)edge_gno/((double)mat->mtx.globalY/(double)mat->comm->nProc));

  return (*part_y);
}

int Zoltan_Distribute_Linear(int edge_gno, int vtx_gno, void* data, int *part_y)
{
  Zoltan_matrix_2d *mat;

  mat = (Zoltan_matrix_2d*) data;
  *part_y = (int)floor((double)edge_gno/((double)mat->mtx.globalY/(double)mat->comm->nProc));

  return (*part_y);
}

int Zoltan_Distribute_Cyclic(int edge_gno, int vtx_gno, void* data, int *part_y)
{
  Zoltan_matrix_2d *mat;

  mat = (Zoltan_matrix_2d*) data;
  *part_y = (int)floor((double)edge_gno/((double)mat->mtx.globalY/(double)mat->comm->nProc));

  return (*part_y);
}

int Zoltan_Distribute_Partition(int edge_gno, int vtx_gno, void* data, int *part_y)
{
  Zoltan_matrix_2d *mat;

  mat = (Zoltan_matrix_2d*) data;
  *part_y = (int)floor((double)edge_gno/((double)mat->mtx.globalY/(double)mat->comm->nProc));

  return (*part_y);
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
  int i, j, cnt;
  int *proclist = NULL;
  Zoltan_Arc *nonzeros= NULL, *sendbuf= NULL;
  int *perm_y = NULL;
  float *wgtarray = NULL;
  float *tmpwgtarray = NULL;
  int msg_tag = 1021982;
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

  nProc_x = outmat->comm->nProc_x;
  nProc_y = outmat->comm->nProc_y;
  myProc_x = outmat->comm->myProc_x;
  myProc_y = outmat->comm->myProc_y;


  ierr = Zoltan_Matrix_Remove_Duplicates(zz, outmat->mtx, &outmat->mtx);

  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution.
   */
  proclist = (int *)ZOLTAN_MALLOC(outmat->mtx.nPins *sizeof(int));
  sendbuf = (Zoltan_Arc*) ZOLTAN_MALLOC(outmat->mtx.nPins * sizeof(Zoltan_Arc));

  if ((outmat->mtx.nPins >0) && (proclist == NULL || sendbuf == NULL)) MEMORY_ERROR;

  wgtarray = (float*) ZOLTAN_MALLOC(outmat->mtx.nPins*outmat->mtx.pinwgtdim*sizeof(float));

  yGNO = outmat->mtx.yGNO;
  pinGNO = outmat->mtx.pinGNO;

  cnt = 0;
  for (i = 0; i < outmat->mtx.nY; i++) {
    int edge_gno=-1;
    /* processor row for the edge */
    edge_gno = yGNO[i];

    for (j = outmat->mtx.ystart[i]; j < outmat->mtx.yend[i]; j++) {
      int vtx_gno=-1;
      /* processor column for the vertex */
      vtx_gno = pinGNO[j];

      proclist[cnt] = (outmat->hashDistFct)(edge_gno, vtx_gno, outmat->hashDistData,
					  &sendbuf[cnt].part_y);
      if (proclist[cnt] < 0) /* Discard this nnz */
        continue;
      sendbuf[cnt].GNO[0] = edge_gno;
      sendbuf[cnt].GNO[1] = vtx_gno;
      memcpy(wgtarray+cnt*outmat->mtx.pinwgtdim, outmat->mtx.pinwgt+j*outmat->mtx.pinwgtdim,
               outmat->mtx.pinwgtdim*sizeof(float));
      cnt++;
    }
  }

  if (outmat->mtx.yend != outmat->mtx.ystart + 1)
    ZOLTAN_FREE(&outmat->mtx.yend);
  outmat->mtx.yend = NULL;
  ZOLTAN_FREE(&outmat->mtx.ystart);
  ZOLTAN_FREE(&outmat->mtx.yGNO);
  ZOLTAN_FREE(&outmat->mtx.pinGNO);
  ZOLTAN_FREE(&outmat->mtx.pinwgt);
  ZOLTAN_FREE(&outmat->mtx.yGID);
  ZOLTAN_FREE(&outmat->mtx.ywgt);


  /*
   * Send pins to their target processors.
   * They become non-zeros in the 2D data distribution.
   */

  msg_tag--;
  ierr = Zoltan_Comm_Create(&plan, cnt, proclist, communicator, msg_tag, &outmat->mtx.nPins);
  ZOLTAN_FREE(&proclist);

  nonzeros = (Zoltan_Arc *) ZOLTAN_MALLOC((outmat->mtx.nPins) * sizeof(Zoltan_Arc));
  if (outmat->mtx.nPins && nonzeros == NULL) MEMORY_ERROR;

  msg_tag--;
  Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, sizeof(Zoltan_Arc),
		 (char *) nonzeros);
  ZOLTAN_FREE(&sendbuf);

  if (outmat->mtx.pinwgtdim) { /* We have to take care about weights */
    tmpwgtarray = (float*) ZOLTAN_MALLOC(outmat->mtx.nPins*outmat->mtx.pinwgtdim*sizeof(float));
    if (outmat->mtx.nPins && tmpwgtarray == NULL) MEMORY_ERROR;

    msg_tag--;
    Zoltan_Comm_Do(plan, msg_tag, (char *) wgtarray, outmat->mtx.pinwgtdim*sizeof(float),
		   (char *) tmpwgtarray);
    ZOLTAN_FREE(&wgtarray);
  }
  Zoltan_Comm_Destroy(&plan);

  /* Unpack the non-zeros received. */

  outmat->mtx.pinGNO = (int *) ZOLTAN_MALLOC(outmat->mtx.nPins * sizeof(int));
  if (outmat->mtx.nPins && (outmat->mtx.pinGNO == NULL)) MEMORY_ERROR;

  outmat->mtx.pinwgt = (float *) ZOLTAN_MALLOC(outmat->mtx.nPins * outmat->mtx.pinwgtdim
					       * sizeof(float));
  if (outmat->mtx.nPins && outmat->mtx.pinwgtdim && outmat->mtx.pinwgt == NULL)
    MEMORY_ERROR;

  /* TODO: do take care about singletons */

  Zoltan_Matrix_Remove_DupArcs(zz, outmat->mtx.nPins, (Zoltan_Arc*)nonzeros, tmpwgtarray,
			       &outmat->mtx);

  /* Now we just have to change numbering */
  outmat->dist_y = (int *) ZOLTAN_CALLOC((nProc_y+1), sizeof(int));
  outmat->dist_x = (int *) ZOLTAN_CALLOC((nProc_x+1), sizeof(int));
  if (outmat->dist_y == NULL || outmat->dist_x == NULL) MEMORY_ERROR;

  /* FIXME: Work only in 1D */
  MPI_Allgather(&outmat->mtx.nY, 1, MPI_INT,
		outmat->dist_y+1, 1, MPI_INT,
		communicator);
  for (i = 1 ; i <= nProc_y ; i ++) {
    outmat->dist_y[i] += outmat->dist_y[i-1];
  }
  outmat->dist_x[1] = outmat->mtx.globalX;

  perm_y = (int *) ZOLTAN_MALLOC(outmat->mtx.nY * sizeof(int));
  if (outmat->mtx.nY > 0 && perm_y == NULL) MEMORY_ERROR;
  for (i = 0 ; i < outmat->mtx.nY ; ++i)
    perm_y[i] = i + outmat->dist_y[myProc_y];

  Zoltan_Matrix_Permute(zz, &outmat->mtx, perm_y);

 End:
  ZOLTAN_FREE(&perm_y);
  ZOLTAN_FREE(&proclist);
  ZOLTAN_FREE(&sendbuf);
  ZOLTAN_FREE(&nonzeros);
  ZOLTAN_FREE(&tmpwgtarray);


  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}


#ifdef __cplusplus
}
#endif
