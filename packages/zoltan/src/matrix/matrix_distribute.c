/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <math.h>

#include "zz_const.h"
#include "zz_util_const.h"
#include "zoltan_dd.h"
#include "phg.h"
#include "zoltan_matrix.h"


typedef struct ZOLTAN_DIST_PART_ {
  ZZ* zz;
  ZOLTAN_MAP* map;
  int nProc;
  int nPart;
} ZOLTAN_DIST_PART;

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
  if (!ranks) return ZOLTAN_MEMERR;

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


int Zoltan_Distribute_Origin(ZOLTAN_GNO_TYPE edge_gno, ZOLTAN_GNO_TYPE vtx_gno, void* data, int *part_y)
{
  ZOLTAN_DIST_PART* part;
  intptr_t answer = 0;

  part = (ZOLTAN_DIST_PART*) data;
  Zoltan_Map_Find(part->zz, part->map, (char *)&edge_gno, &answer);
  *part_y = (int)answer;

  return (*part_y);
}

int Zoltan_Distribute_Linear(ZOLTAN_GNO_TYPE edge_gno, ZOLTAN_GNO_TYPE vtx_gno, void* data, int *part_y)
{
  Zoltan_matrix_2d *mat;

  mat = (Zoltan_matrix_2d*) data;
  *part_y = (int)floor((double)edge_gno/((double)mat->mtx.globalY/(double)mat->comm->nProc));

  return (*part_y);
}

int Zoltan_Distribute_Cyclic(ZOLTAN_GNO_TYPE edge_gno, ZOLTAN_GNO_TYPE vtx_gno, void* data, int *part_y)
{
  Zoltan_matrix_2d *mat;

  mat = (Zoltan_matrix_2d*) data;
  *part_y = edge_gno%mat->comm->nProc;

  return (*part_y);
}

int Zoltan_Distribute_Partition(ZOLTAN_GNO_TYPE edge_gno, ZOLTAN_GNO_TYPE vtx_gno, void* data, int *part_y)
{
  ZOLTAN_DIST_PART* part;
  intptr_t answer;

  part = (ZOLTAN_DIST_PART*) data;
  Zoltan_Map_Find(part->zz, part->map, (char *)&edge_gno, &answer);
  *part_y = (int)answer;

  return ((int)floor((double)*part_y/((double)part->nPart/(double)part->nProc)));
}


void* Zoltan_Distribute_Partition_Register(ZZ* zz, int size, ZOLTAN_GNO_TYPE *yGNO, int *part, int nProc, int nPart)
{
  ZOLTAN_DIST_PART* dist;
  int i;

  dist = (ZOLTAN_DIST_PART*) ZOLTAN_MALLOC(sizeof(ZOLTAN_DIST_PART));
  if (dist == NULL)
    return (NULL);

  dist->zz = zz;

  dist->map = Zoltan_Map_Create(zz, 0, sizeof(ZOLTAN_GNO_TYPE), 1, size);
  if (dist->map == NULL) {
    ZOLTAN_FREE(&dist);
    return (NULL);
  }

  for (i = 0 ; i < size ; ++i ) {
    Zoltan_Map_Add(dist->zz, dist->map, (char *)&yGNO[i], (intptr_t)part[i]);
  }

  dist->nProc = nProc;
  dist->nPart = nPart;

  return ((void*)dist);
}

void
Zoltan_Distribute_Partition_Free(void** dist)
{
  ZOLTAN_DIST_PART* part;

  if (dist == NULL || *dist == NULL)
    return;

  part = (ZOLTAN_DIST_PART*) (*dist);
  Zoltan_Map_Destroy(part->zz, &part->map);
  ZOLTAN_FREE(dist);
}

/* if !copy, inmat is not usable after this call */
int
Zoltan_Matrix2d_Distribute (ZZ* zz, Zoltan_matrix inmat, /* Cannot be const as we can share it inside outmat */
			    Zoltan_matrix_2d *outmat, int copy)
{
  static char *yo = "Zoltan_Matrix2d_Distribute";
  int ierr = ZOLTAN_OK;
  int nProc_x, nProc_y;
  int myProc_x, myProc_y;
  int i, j, cnt;
  int *proclist = NULL;
  Zoltan_Arc *nonzeros= NULL, *sendbuf= NULL;
  ZOLTAN_GNO_TYPE *perm_y = NULL;
  float *wgtarray = NULL;
  float *tmpwgtarray = NULL;
  int msg_tag = 1021982;
  ZOLTAN_COMM_OBJ *plan;
  MPI_Comm communicator = MPI_COMM_NULL;
  int nProc;
  ZOLTAN_GNO_TYPE *yGNO = NULL;
  ZOLTAN_GNO_TYPE *pinGNO = NULL;
  ZOLTAN_GNO_TYPE tmp_gno;
  void *partdata = NULL;
  MPI_Datatype zoltan_gno_mpi_type;

  ZOLTAN_TRACE_ENTER(zz, yo);

  zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

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

KDDKDDKDD(zz->Proc, "    Zoltan_Matrix_Remove_Duplicates");
  ierr = Zoltan_Matrix_Remove_Duplicates(zz, outmat->mtx, &outmat->mtx);

/* KDDKDDKDD  FIX INDENTATION OF THIS BLOCK */
if (inmat.opts.speed != MATRIX_NO_REDIST) {
  if (outmat->hashDistFct == (distFnct *)&Zoltan_Distribute_Origin) {
    /* I need to know the original distribution */
    if (outmat->mtx.ddX != outmat->mtx.ddY) { /* No initial distribution */
      outmat->hashDistFct = (distFnct *)&Zoltan_Distribute_Linear;
    }
    else {
      int *cmember = NULL;

      cmember = (int*)ZOLTAN_MALLOC(outmat->mtx.nY*sizeof(int));
      if (outmat->mtx.nY > 0 && cmember == NULL) MEMORY_ERROR;
      Zoltan_DD_Find (outmat->mtx.ddY, (ZOLTAN_ID_PTR)outmat->mtx.yGNO, NULL, (char *)cmember, NULL,
		      outmat->mtx.nY, NULL);
KDDKDDKDD(zz->Proc, "    Zoltan_Distribute_Partition_Register");
      partdata = Zoltan_Distribute_Partition_Register(zz, outmat->mtx.nY, outmat->mtx.yGNO,
						      cmember, zz->Num_Proc, zz->Num_Proc);
      ZOLTAN_FREE(&cmember);
      Zoltan_Distribute_Set(outmat, (distFnct *)&Zoltan_Distribute_Origin, partdata);
    }
  }

  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution.
   */
  /* TRICK: create fake arc (edgeno, -1) for empty Y. Upper bound for size might be nPins + nY */
  proclist = (int *)ZOLTAN_MALLOC((outmat->mtx.nPins+outmat->mtx.nY) *sizeof(int));
  sendbuf = (Zoltan_Arc*) ZOLTAN_MALLOC((outmat->mtx.nPins +outmat->mtx.nY)* sizeof(Zoltan_Arc));

  if ((outmat->mtx.nPins + outmat->mtx.nY >0) && (proclist == NULL || sendbuf == NULL)) MEMORY_ERROR;

  wgtarray = (float*) ZOLTAN_MALLOC((outmat->mtx.nPins+outmat->mtx.nY)*outmat->mtx.pinwgtdim*sizeof(float));

  if (outmat->mtx.nPins*outmat->mtx.pinwgtdim && !wgtarray) MEMORY_ERROR;

  yGNO = outmat->mtx.yGNO;
  pinGNO = outmat->mtx.pinGNO;

KDDKDDKDD(zz->Proc, "    CommPlan Hash");
  cnt = 0;
  for (i = 0; i < outmat->mtx.nY; i++) {
    ZOLTAN_GNO_TYPE edge_gno=-1;
    /* processor row for the edge */
    edge_gno = yGNO[i];

    for (j = outmat->mtx.ystart[i]; j < outmat->mtx.yend[i]; j++) {
      ZOLTAN_GNO_TYPE vtx_gno=-1;
      /* processor column for the vertex */
      vtx_gno = pinGNO[j];

      proclist[cnt] = (*outmat->hashDistFct)(edge_gno, vtx_gno, outmat->hashDistData,
					  &sendbuf[cnt].part_y);
      if (proclist[cnt] < 0) /* Discard this nnz */
        continue;
      sendbuf[cnt].GNO[0] = edge_gno;
      sendbuf[cnt].GNO[1] = vtx_gno;
      memcpy(wgtarray+cnt*outmat->mtx.pinwgtdim, outmat->mtx.pinwgt+j*outmat->mtx.pinwgtdim,
               outmat->mtx.pinwgtdim*sizeof(float));
      cnt++;
    }
    if(outmat->mtx.ystart[i] == outmat->mtx.yend[i]) {
      proclist[cnt] = (*outmat->hashDistFct)(edge_gno, -1, outmat->hashDistData,
					  &sendbuf[cnt].part_y);
      if (proclist[cnt] < 0) /* Discard this nnz */
        continue;
      sendbuf[cnt].GNO[0] = edge_gno;
      sendbuf[cnt].GNO[1] = -1;
      memset(wgtarray+cnt*outmat->mtx.pinwgtdim, 0,outmat->mtx.pinwgtdim*sizeof(float));
      cnt++;
    }
  }

  if (outmat->hashDistFct == (distFnct *)&Zoltan_Distribute_Origin)
    Zoltan_Distribute_Partition_Free(&outmat->hashDistData);

  if (outmat->mtx.yend != outmat->mtx.ystart + 1)
    ZOLTAN_FREE(&outmat->mtx.yend);
  outmat->mtx.yend = NULL;
  ZOLTAN_FREE(&outmat->mtx.ystart);
  ZOLTAN_FREE(&outmat->mtx.yGNO);
  ZOLTAN_FREE(&outmat->mtx.pinGNO);
  ZOLTAN_FREE(&outmat->mtx.pinwgt);
  ZOLTAN_FREE(&outmat->mtx.yGID);


  /*
   * Send pins to their target processors.
   * They become non-zeros in the 2D data distribution.
   */

KDDKDDKDD(zz->Proc, "    CommPlan Create");
  msg_tag--;
  ierr = Zoltan_Comm_Create(&plan, cnt, proclist, communicator, msg_tag, &outmat->mtx.nPins);
  ZOLTAN_FREE(&proclist);

  nonzeros = (Zoltan_Arc *) ZOLTAN_MALLOC((outmat->mtx.nPins) * sizeof(Zoltan_Arc));
  if (outmat->mtx.nPins && nonzeros == NULL) MEMORY_ERROR;

  msg_tag--;
  Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, sizeof(Zoltan_Arc), (char *) nonzeros);
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

KDDKDDKDD(zz->Proc, "    Zoltan_Matrix_Remove_DupArcs");
  /* TODO: do take care about singletons */
  Zoltan_Matrix_Remove_DupArcs(zz, outmat->mtx.nPins, (Zoltan_Arc*)nonzeros, tmpwgtarray,
			       &outmat->mtx);
}

  /* Now we just have to change numbering */
  outmat->dist_y = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC((nProc_y+1), sizeof(ZOLTAN_GNO_TYPE));
  outmat->dist_x = (ZOLTAN_GNO_TYPE *) ZOLTAN_CALLOC((nProc_x+1), sizeof(ZOLTAN_GNO_TYPE));
  if (outmat->dist_y == NULL || outmat->dist_x == NULL) MEMORY_ERROR;

  /* FIXME: Work only in 1D */
  tmp_gno = (ZOLTAN_GNO_TYPE)outmat->mtx.nY;
  MPI_Allgather(&tmp_gno, 1, zoltan_gno_mpi_type, outmat->dist_y+1, 1, zoltan_gno_mpi_type, communicator);
  for (i = 1 ; i <= nProc_y ; i ++) {
    outmat->dist_y[i] += outmat->dist_y[i-1];
  }
  outmat->dist_x[1] = outmat->mtx.globalX;

  perm_y = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(outmat->mtx.nY * sizeof(ZOLTAN_GNO_TYPE));
  if (outmat->mtx.nY > 0 && perm_y == NULL) MEMORY_ERROR;
  for (i = 0 ; i < outmat->mtx.nY ; ++i){
    perm_y[i] = i + outmat->dist_y[myProc_y];
  }

KDDKDDKDD(zz->Proc, "    Zoltan_Matrix_Permute");
  Zoltan_Matrix_Permute(zz, &outmat->mtx, perm_y);

KDDKDDKDD(zz->Proc, "    Zoltan_Matrix_Permute done");
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
