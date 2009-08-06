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
  int slice;
  int offset, nX;

  ZOLTAN_TRACE_ENTER(zz, yo);

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
  /* First update yGNO: new yGNO = prev yGNO + matrix->globalX */
  matrix->yGNO = (int*) ZOLTAN_REALLOC(matrix->yGNO, newGNOsize*sizeof(int));
  if (newGNOsize && matrix->yGNO == NULL) MEMORY_ERROR;
  for (i=0 ; i < matrix->nY ; ++i)
    matrix->yGNO[i] += matrix->globalX;

  /* Then update the pinGNO and ystart */
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

  matrix->offsetY = matrix->globalX;
  matrix->nY = i + 1; /* i is 0 based */
  matrix->nPins *= 2;

  /* Now update the xGNO to be coherent ! */
  /* Ownership by slices */
  slice = matrix->globalY/nProc;
  offset = slice * myProc;
/*   nX = matrix->nX + MIN(slice*(myProc+1), matrix->globalY)-offset; */
/*   matrix->xGNO = (int*) ZOLTAN_REALLOC(matrix->xGNO, nX*sizeof(int)); */
/*   if (nX && matrix->xGNO == NULL) MEMORY_ERROR; */
/*   for (i = matrix->nX ; i < nX ; ++i) */
/*     matrix->xGNO[i] = matrix->globalX + offset + i; */

/*   matrix->nX= nX; */
  matrix->globalX += matrix->globalY;
  matrix->globalY += matrix->offsetY;

 End:
  ZOLTAN_FREE(&ystart);
  ZOLTAN_FREE(&tmparray);
  ZOLTAN_FREE(&tr_tab);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


#ifdef __cplusplus
}
#endif
