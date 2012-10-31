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

int
Zoltan_Matrix_Sym(ZZ* zz, Zoltan_matrix *matrix, int bipartite)
{
  static char *yo = "Zoltan_Matrix_Sym";
  int ierr = ZOLTAN_OK;
  Zoltan_Arc *tr_tab = NULL;
  int i, j, cnt;
  ZOLTAN_ID_PTR yGID = NULL;
  int *ypid=NULL;
  float *pinwgt=NULL;
  int * ybipart = NULL;
  int gno_size_for_dd;

  gno_size_for_dd = sizeof(ZOLTAN_GNO_TYPE) / sizeof(ZOLTAN_ID_TYPE);

  ZOLTAN_TRACE_ENTER(zz, yo);
  if (bipartite || !matrix->opts.enforceSquare) {
    matrix->bipartite = 1;
    matrix->redist = 1;
  }

  if (matrix->ywgtdim != zz->Obj_Weight_Dim)
      FATAL_ERROR("Cannot form bipartite graph: vertex and edge weights are not consistant");

  matrix->opts.symmetrize = 1;

  /* Update the data directories */
  tr_tab = (Zoltan_Arc*) ZOLTAN_MALLOC(sizeof(Zoltan_Arc)*(matrix->nPins*2+matrix->nY));
  if (matrix->nPins && tr_tab == NULL) MEMORY_ERROR;

  pinwgt = (float*)ZOLTAN_MALLOC((matrix->nPins*2+matrix->nY)*matrix->pinwgtdim*sizeof(float));

  if (((matrix->nPins*2+matrix->nY)*matrix->pinwgtdim) && !pinwgt) MEMORY_ERROR;

  for (i=0, cnt = 0 ; i < matrix->nY ; ++i) {
    for (j = matrix->ystart[i] ; j < matrix->yend[i] ; ++j) {
      tr_tab[cnt].GNO[0] = matrix->yGNO[i] + bipartite*matrix->globalX;   /* Normal arc */
      tr_tab[cnt].GNO[1] = matrix->pinGNO[j];
      memcpy(pinwgt + cnt*matrix->pinwgtdim, matrix->pinwgt+j*matrix->pinwgtdim,
	     matrix->pinwgtdim*sizeof(float));
      cnt ++;

      tr_tab[cnt].GNO[0] = matrix->pinGNO[j];          /* Symmetric arc */
      tr_tab[cnt].GNO[1] = matrix->yGNO[i] + bipartite*matrix->globalX; /* new ordering */
      memcpy(pinwgt + cnt*matrix->pinwgtdim, matrix->pinwgt+j*matrix->pinwgtdim,
	     matrix->pinwgtdim*sizeof(float));
      cnt ++;
    }
    if (matrix->ystart[i] == matrix->yend[i]) { /* Singleton */
      tr_tab[cnt].GNO[0] = matrix->yGNO[i] + bipartite*matrix->globalX;   /* Normal arc */
      tr_tab[cnt].GNO[1] = -1;
      cnt ++;
    }
  }
  ZOLTAN_FREE(&matrix->pinwgt);

  Zoltan_Matrix_Remove_DupArcs(zz, cnt, tr_tab, pinwgt, matrix); /* Compute new nY too */
  ZOLTAN_FREE(&tr_tab);
  ZOLTAN_FREE(&pinwgt);

  if (bipartite) {
    int endX;
    ZOLTAN_GNO_TYPE * yGNO = NULL;

    /* Update data directories */
    yGID = ZOLTAN_MALLOC_GID_ARRAY(zz, matrix->nY);
    ypid = (int*) ZOLTAN_MALLOC(matrix->nY*sizeof(int));

    ybipart = (int*) ZOLTAN_MALLOC(matrix->nY*sizeof(int));

    if (matrix->nY && (!yGID || !ypid || !ybipart)) MEMORY_ERROR;

    for (endX = 0 ; endX < matrix->nY ; ++endX) {
      if (matrix->yGNO[endX] >= matrix->globalX)
	break;
      ybipart[endX] = 0;
    }
    /* Get Informations about X */
    Zoltan_DD_Find (matrix->ddX, (ZOLTAN_ID_PTR)matrix->yGNO, yGID, (char *)ypid, NULL,
		    endX, NULL);

    yGNO = (ZOLTAN_GNO_TYPE*)ZOLTAN_MALLOC(endX*sizeof(ZOLTAN_GNO_TYPE));
    if (endX && !yGNO) MEMORY_ERROR;

    for (i = endX ; i < matrix->nY ; ++i) {
      yGNO[i-endX] = matrix->yGNO[i] - matrix->globalX;
      /* TODO: add a something to have the correct ypid */
      ybipart[i] = 1;
      ypid[i] = ypid[i-endX]; /* Be sure to update ypid */
    }
    /* Get Informations about Y */
    Zoltan_DD_Find (matrix->ddY, (ZOLTAN_ID_PTR)yGNO,
		    yGID + endX*zz->Num_GID,
		    NULL, NULL,
		    matrix->nY-endX, NULL);

    if (matrix->ddY != matrix->ddX)
      Zoltan_DD_Destroy (&matrix->ddY);
    Zoltan_DD_Destroy (&matrix->ddX);

    matrix->globalX += matrix->globalY;
    matrix->globalY = matrix->globalX;

    /* I store : xGNO, xGID, xpid, bipart */
    ierr = Zoltan_DD_Create (&matrix->ddX, zz->Communicator, gno_size_for_dd, zz->Num_GID,
			     sizeof(int), matrix->globalX/zz->Num_Proc, 0);
    matrix->ddY = matrix->ddX;
    /* Hope a linear assignment will help a little */
    if (matrix->globalX/zz->Num_Proc)
      Zoltan_DD_Set_Neighbor_Hash_Fn1(matrix->ddX, matrix->globalX/zz->Num_Proc);
    /* Associate all the data with our xyGNO */
    Zoltan_DD_Update (matrix->ddX, (ZOLTAN_ID_PTR)matrix->yGNO, yGID, (char *)ypid, ybipart,
		      matrix->nY);
    ZOLTAN_FREE(&yGNO);
  }

 End:
  ZOLTAN_FREE(&ybipart);
  ZOLTAN_FREE(&ypid);
  ZOLTAN_FREE(&pinwgt);
  ZOLTAN_FREE(&yGID);
  ZOLTAN_FREE(&tr_tab);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


#ifdef __cplusplus
}
#endif
