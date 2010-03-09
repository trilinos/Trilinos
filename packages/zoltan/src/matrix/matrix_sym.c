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
  float *pinwgt=NULL;
  int *ypid = NULL;
  int *yoffset = NULL;
  ZOLTAN_MAP* obj_map = NULL;

  ZOLTAN_TRACE_ENTER(zz, yo);
  if (bipartite || !matrix->opts.enforceSquare) {
    bipartite = 1;
    matrix->redist = 1;
  }

  if (matrix->ywgtdim != zz->Obj_Weight_Dim)
      FATAL_ERROR("Cannot form bipartite graph: vertex and edge weights are not consistant");

  matrix->opts.symmetrize = 1;

  /* Save objs informations */
  ypid = matrix->ypid;
  yoffset = matrix->yoffset;
/*   obj_map = Zoltan_Matrix_Save_Obj(zz, matrix->nY, ypid, yoffset); */
/*   if (obj_map == NULL) MEMORY_ERROR; */
  matrix->ypid = NULL;
  matrix->yoffset = NULL;
  matrix->xpid = NULL;
  matrix->xoffset = NULL;

  /* Construct arcs and their symmetric */
  tr_tab = (Zoltan_Arc*) ZOLTAN_MALLOC(sizeof(Zoltan_Arc)*matrix->nPins*2);
  if (matrix->nPins && tr_tab == NULL) MEMORY_ERROR;

  pinwgt = (float*)ZOLTAN_MALLOC(matrix->nPins*2*matrix->pinwgtdim*sizeof(float));
  for (i = 0 ; i < 2 ; ++i) /* Copy pin weights */
    memcpy(pinwgt + i*matrix->nPins*matrix->pinwgtdim*sizeof(float),
	   matrix->pinwgt, matrix->nPins*matrix->pinwgtdim*sizeof(float));
  ZOLTAN_FREE(&matrix->pinwgt);

  for (i=0, cnt = 0 ; i < matrix->nY ; ++i) {
    for (j = matrix->ystart[i] ; j < matrix->yend[i] ; ++j) {
      tr_tab[cnt].GNO[0] = matrix->yGNO[i] + bipartite*matrix->globalX;   /* Normal arc */
      tr_tab[cnt].GNO[1] = matrix->pinGNO[j];
      cnt ++;

      tr_tab[cnt].GNO[0] = matrix->pinGNO[j];                        /* Symmetric arc */
      tr_tab[cnt].GNO[1] = matrix->yGNO[i] + bipartite*matrix->globalX; /* new ordering */
      cnt ++;
   }
  }

  Zoltan_Matrix_Remove_DupArcs(zz, cnt, tr_tab, pinwgt, matrix);
  ZOLTAN_FREE(&tr_tab);
  ZOLTAN_FREE(&pinwgt);

  /* Complete objs informations */
  matrix->ypid = matrix->xpid = (int*) ZOLTAN_MALLOC(matrix->nY*sizeof(int));
  if (matrix->nY > 0 && matrix->ypid == NULL) MEMORY_ERROR;
  matrix->yoffset = matrix->xoffset = (int*) ZOLTAN_MALLOC(matrix->nY*sizeof(int));
  if (matrix->nY > 0 && matrix->yoffset == NULL) MEMORY_ERROR;
  matrix->ybipart = (int*) ZOLTAN_CALLOC(matrix->nY,sizeof(int));
  if (matrix->nY > 0 && matrix->ybipart == NULL) MEMORY_ERROR;

  /* TODO: code works only for square matrices */
  for (i = 0 ; i < matrix->nY ; ++i) {
    void *y_ptr;
    int yGNO;

    yGNO = matrix->yGNO[i];
    if (yGNO > matrix->globalX) { /* Bipartite vertex */
      yGNO -= matrix->globalX; /* information about ancestor */
      matrix->ybipart[i] = 1;
    }
    Zoltan_Map_Find(zz, obj_map, &yGNO, (void**)&y_ptr);

    if (y_ptr == NULL) { /* Don't know this information */
      matrix->ypid[i] = -1;
      matrix->yoffset[i] = -1;
    }
    else {
      int y_index;
      y_index = (int)(long)y_ptr;
      matrix->ypid[i] = ypid[y_index];
      matrix->yoffset[i] = yoffset[y_index];
    }
  }

 End:
  Zoltan_Map_Destroy(zz, &obj_map);
  ZOLTAN_FREE(&pinwgt);
  ZOLTAN_FREE(&tr_tab);
  ZOLTAN_FREE(&ypid);
  ZOLTAN_FREE(&yoffset);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}


#ifdef __cplusplus
}
#endif
