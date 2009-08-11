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

struct Transpose_Elem {
  int xGNO;
  int pinGNO;
};

int
compar_elems(const struct Transpose_Elem *elem1, const struct Transpose_Elem* elem2)
{
  if (elem1->xGNO == elem2->xGNO)            /* useful only to be determinist in our case */
    return (elem1->pinGNO - elem2->pinGNO);
  return (elem1->xGNO - elem2->xGNO);
}


int
Zoltan_Matrix_Bipart(ZZ* zz, Zoltan_matrix *matrix, int nProc, int myProc)
{
  static char *yo = "Zoltan_Matrix_Bipartite";
  int ierr = ZOLTAN_OK;
  int *tmparray = NULL;
  int *ystart = NULL;
  struct Transpose_Elem *tr_tab = NULL;
  int i, j, cnt;
  int newGNOsize = 0;
  ZOLTAN_ID_PTR yGID = NULL;
  float *ywgt = NULL;
  int *Input_Parts = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Update the data directories */

  tr_tab = (struct Transpose_Elem*) ZOLTAN_MALLOC(sizeof(struct Transpose_Elem)*matrix->nPins);
  if (matrix->nPins && tr_tab == NULL) MEMORY_ERROR;

  for (i=0, cnt = 0 ; i < matrix->nY ; ++i) {
    for (j = matrix->ystart[i] ; j < matrix->yend[i] ; ++j) {
      tr_tab[cnt].xGNO = matrix->pinGNO[j];
      tr_tab[cnt].pinGNO = matrix->yGNO[i] + matrix->globalX; /* new ordering */
      cnt ++;
    }
  }

  /* Sort the transpose part in order to be able to construct easily the CSR
     part */
  qsort ((void*)tr_tab, matrix->nPins, sizeof(struct Transpose_Elem),
	 (int (*)(const void*,const void*))compar_elems);

  /* Now we will construct the new indirection arrays */
  newGNOsize = matrix->nY + MIN(matrix->nPins, matrix->globalX); /* Maximum number of vertices we can have */

  /* I do array one by one to avoid to allocate all the arrays at the same time */
  matrix->yGNO = (int*) ZOLTAN_REALLOC(matrix->yGNO, newGNOsize*sizeof(int));
  if (newGNOsize && matrix->yGNO == NULL) MEMORY_ERROR;

  /* Update the pinGNO and ystart */
  /* The output use a compact array */
  tmparray = (int*)ZOLTAN_MALLOC(2*matrix->nPins*sizeof(int));
  if (matrix->nPins && tmparray == NULL) MEMORY_ERROR;
  ystart = (int*) ZOLTAN_MALLOC((newGNOsize + 1)*sizeof(int));
  if (ystart == NULL) MEMORY_ERROR;
  for (i=0, cnt=0 ; i < matrix->nY ; ++i) {
    ystart[i] = cnt;
    for (j = matrix->ystart[i] ; j < matrix->yend[i] ; ++j)
      tmparray[cnt++] = matrix->pinGNO[j];
  }
  ystart[matrix->nY] = cnt;

  if (matrix->yend != matrix->ystart + 1)
    ZOLTAN_FREE(&matrix->yend);
  ZOLTAN_FREE(&matrix->ystart);
  matrix->ystart = ystart;
  matrix->yend = matrix->ystart + 1;
  ystart = NULL;

  ZOLTAN_FREE(&matrix->pinGNO);
  matrix->pinGNO = tmparray;
  tmparray = NULL;

  /* Now deal with the new nnz */
  for (cnt = 0,j =matrix->nPins, i=matrix->nY-1 ;
       cnt < matrix->nPins ; ++cnt) {
    int xGNO;
    xGNO = tr_tab[cnt].xGNO;
    if (matrix->yGNO[i] != xGNO) {
      matrix->yend[i] = j;
      i++;
      /* matrix->ystart[i] = j */ /* No need because compact view ! */
      matrix->yGNO[i] = xGNO;
    }
    matrix->pinGNO[j] = tr_tab[cnt].pinGNO;
    j++;
  }
  matrix->yend[i] = j;
  newGNOsize = i + 1; /* i is 0 based */
  ZOLTAN_FREE(&tr_tab);


  /* Update data directories */
  if (matrix->ywgtdim != zz->Obj_Weight_Dim)
      FATAL_ERROR("Cannot form bipartite graph: vertex and edge weights are not consistant");
  Input_Parts = (int*) ZOLTAN_MALLOC(newGNOsize * sizeof(int));
  yGID = ZOLTAN_MALLOC_GID_ARRAY(zz, newGNOsize);
  ywgt = (float*) ZOLTAN_MALLOC(newGNOsize * sizeof(float) * matrix->ywgtdim);
  if (newGNOsize && (Input_Parts == NULL || yGID == NULL
		     || (matrix->ywgtdim && ywgt == NULL)))
    MEMORY_ERROR;

  /* Get Informations about Y */
  Zoltan_DD_Find (matrix->ddY, (ZOLTAN_ID_PTR)matrix->yGNO, yGID, (ZOLTAN_ID_PTR)ywgt, NULL,
		  matrix->nY, NULL);
  /* Get Informations about X */
  Zoltan_DD_Find (matrix->ddX, (ZOLTAN_ID_PTR)matrix->yGNO + matrix->nY, yGID + matrix->nY*zz->Num_GID,
		  (ZOLTAN_ID_PTR)ywgt + matrix->nY*sizeof(float)*matrix->ywgtdim/sizeof(int),
		  Input_Parts + matrix->nY,
		  newGNOsize - matrix->nY, NULL);

  if (matrix->ddY != matrix->ddX)
    Zoltan_DD_Destroy (&matrix->ddY);
  Zoltan_DD_Destroy (&matrix->ddX);

  /* Update yGNO: new yGNO = prev yGNO + matrix->globalX */
  for (i=0 ; i < matrix->nY ; ++i) {
    matrix->yGNO[i] += matrix->globalX;
    Input_Parts[i] = -1;
  }

  matrix->offsetY = matrix->globalX;
  matrix->nY = newGNOsize;
  matrix->nPins *= 2;

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

 End:
  ZOLTAN_FREE(&ywgt);
  ZOLTAN_FREE(&yGID);
  ZOLTAN_FREE(&Input_Parts);
  ZOLTAN_FREE(&ystart);
  ZOLTAN_FREE(&tmparray);
  ZOLTAN_FREE(&tr_tab);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


#ifdef __cplusplus
}
#endif
