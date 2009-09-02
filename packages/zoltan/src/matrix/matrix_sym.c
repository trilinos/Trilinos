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

int
Zoltan_Matrix_Sym(ZZ* zz, Zoltan_matrix *matrix, int bipartite)
{
  static char *yo = "Zoltan_Matrix_Sym";
  int ierr = ZOLTAN_OK;
  Zoltan_Arc *tr_tab = NULL;
  int i, j, cnt;
  ZOLTAN_ID_PTR yGID = NULL;
  float *ywgt = NULL;
  int *Input_Parts = NULL;
  float *pinwgt=NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);
  if (bipartite || !matrix->opts.enforceSquare) {
    bipartite = 1;
    matrix->redist = 1;
  }

  if (matrix->ywgtdim != zz->Obj_Weight_Dim)
      FATAL_ERROR("Cannot form bipartite graph: vertex and edge weights are not consistant");

  /* Update the data directories */
  tr_tab = (Zoltan_Arc*) ZOLTAN_MALLOC(sizeof(Zoltan_Arc)*matrix->nPins*2);
  if (matrix->nPins && tr_tab == NULL) MEMORY_ERROR;

  for (i=0, cnt = 0 ; i < matrix->nY ; ++i) {
    for (j = matrix->ystart[i] ; j < matrix->yend[i] ; ++j) {
      tr_tab[cnt].yGNO = matrix->yGNO[i] + bipartite*matrix->globalX;   /* Normal arc */
      tr_tab[cnt].pinGNO = matrix->pinGNO[j];
      tr_tab[cnt].offset = j;
      cnt ++;

      tr_tab[cnt].yGNO = matrix->pinGNO[j];                        /* Symmetric arc */
      tr_tab[cnt].pinGNO = matrix->yGNO[i] + bipartite*matrix->globalX; /* new ordering */
      tr_tab[cnt].offset = j;
      cnt ++;
    }
  }

  matrix->nY += MIN(matrix->globalX, matrix->nPins);
  matrix->nPins*=2;
  if (matrix->yend != matrix->ystart + 1)
    ZOLTAN_FREE(&matrix->yend);
  matrix->yend=NULL;
  ZOLTAN_FREE(&matrix->ystart);
  ZOLTAN_FREE(&matrix->yGNO);
  ZOLTAN_FREE(&matrix->pinGNO);
  pinwgt = matrix->pinwgt;

  matrix->ystart = (int*) ZOLTAN_MALLOC((matrix->nY+1)*sizeof(int));
  if (matrix->ystart == NULL) MEMORY_ERROR;
  matrix->yGNO = (int*) ZOLTAN_MALLOC( matrix->nY*sizeof(int));
  if (matrix->nY && matrix->yGNO == NULL) MEMORY_ERROR;
  matrix->pinGNO = (int*) ZOLTAN_MALLOC(matrix->nPins*sizeof(int));
  if (matrix->nPins && matrix->pinGNO == NULL) MEMORY_ERROR;
  matrix->pinwgt = (float*) ZOLTAN_MALLOC(matrix->pinwgtdim*matrix->nPins*sizeof(float));
  if (matrix->nPins && matrix->pinwgtdim && matrix->pinwgt == NULL)
    MEMORY_ERROR;

  Zoltan_Matrix_Remove_DupArcs(zz, cnt, tr_tab, pinwgt, matrix);
  ZOLTAN_FREE(&tr_tab);
  ZOLTAN_FREE(&pinwgt);

  if (bipartite) {
    /* Update data directories */
    Input_Parts = (int*) ZOLTAN_MALLOC(matrix->nY * sizeof(int));
    yGID = ZOLTAN_MALLOC_GID_ARRAY(zz, matrix->nY);
    ywgt = (float*) ZOLTAN_MALLOC(matrix->nY * sizeof(float) * matrix->ywgtdim);
    if (matrix->nY && (Input_Parts == NULL || yGID == NULL
		       || (matrix->ywgtdim && ywgt == NULL)))
      MEMORY_ERROR;

    /* Get Informations about Y */
    Zoltan_DD_Find (matrix->ddY, (ZOLTAN_ID_PTR)matrix->yGNO, yGID, (ZOLTAN_ID_PTR)ywgt, NULL,
		    matrix->nY, NULL);
    /* Get Informations about X */
    Zoltan_DD_Find (matrix->ddX, (ZOLTAN_ID_PTR)matrix->yGNO + matrix->nY, yGID + matrix->nY*zz->Num_GID,
		    (ZOLTAN_ID_PTR)ywgt + matrix->nY*sizeof(float)*matrix->ywgtdim/sizeof(int),
		    Input_Parts + matrix->nY,
		    matrix->nY - matrix->nY, NULL);

    if (matrix->ddY != matrix->ddX)
      Zoltan_DD_Destroy (&matrix->ddY);
    Zoltan_DD_Destroy (&matrix->ddX);

  /* Update yGNO: new yGNO = prev yGNO + matrix->globalX */
    /* for (i=0 ; i < matrix->nY ; ++i) { */
    /*   matrix->yGNO[i] += matrix->globalX; */
    /*   Input_Parts[i] = -1; */
    /* } */

    matrix->offsetY = matrix->globalX;
    matrix->globalX += matrix->globalY;
    matrix->globalY += matrix->offsetY;

    /* I store : xGNO, xGID, xwgt, Input_Part */
    ierr = Zoltan_DD_Create (&matrix->ddX, zz->Communicator, 1, zz->Num_GID,
			     matrix->ywgtdim*sizeof(float)/sizeof(int), matrix->globalX/zz->Num_Proc, 0);
    matrix->ddY = matrix->ddX;
    /* Hope a linear assignment will help a little */
    Zoltan_DD_Set_Neighbor_Hash_Fn1(matrix->ddX, matrix->globalX/zz->Num_Proc);
    /* Associate all the data with our xyGNO */
    Zoltan_DD_Update (matrix->ddX, (ZOLTAN_ID_PTR)matrix->yGNO, yGID, (ZOLTAN_ID_PTR) ywgt,
		      Input_Parts, matrix->nY);
  }

 End:
  ZOLTAN_FREE(&pinwgt);
  ZOLTAN_FREE(&ywgt);
  ZOLTAN_FREE(&yGID);
  ZOLTAN_FREE(&Input_Parts);
  ZOLTAN_FREE(&tr_tab);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


#ifdef __cplusplus
}
#endif
