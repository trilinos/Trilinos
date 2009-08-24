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

void
Zoltan_Matrix_Free(ZZ *zz, Zoltan_matrix *m)
{
  if (m->yend != m->ystart + 1)
    ZOLTAN_FREE(&m->yend);
  ZOLTAN_FREE(&m->ystart);
  ZOLTAN_FREE(&m->yGNO);
  ZOLTAN_FREE(&m->pinGNO);
  ZOLTAN_FREE(&m->yGID);

  if (m->ddY != m->ddX && m->ddY != NULL)
    Zoltan_DD_Destroy(&m->ddY);
  if (m->ddX != NULL)
    Zoltan_DD_Destroy(&m->ddX);

  memset (m, 0, sizeof(Zoltan_matrix));
}

void
Zoltan_Matrix2d_Free(ZZ *zz, Zoltan_matrix_2d *m)
{
  Zoltan_Matrix_Free (zz, &m->mtx);

  ZOLTAN_FREE(&m->dist_x);
  ZOLTAN_FREE(&m->dist_y);

  Zoltan_PHGComm_Destroy(m->comm);
}

void
Zoltan_Matrix_Reset(Zoltan_matrix* m)
{
  m->yGNO = NULL;
  m->ystart = NULL;
  m->yend = NULL;
  m->pinGNO = NULL;
  m->pinwgt = NULL;
  m->ywgt = NULL;
  m->yGID = NULL;
}


/* This function compute the indices of the diagonal terms.
   This function needs that diagonal terms are declared at most
   1 time locally.
 */
int
Zoltan_Matrix_Mark_Diag(ZZ* zz, const Zoltan_matrix* const m,
			int *n_nnz, int **nnz)
{
  static char *yo = "Zoltan_Matrix_Mark_Diag";
  int ierr = ZOLTAN_OK;
  int y;

  ZOLTAN_TRACE_ENTER(zz, yo);

  (*nnz) = (int*)ZOLTAN_MALLOC(m->nY*sizeof(int));
  if (m->nY && (*nnz) == NULL)
    MEMORY_ERROR;

  (*n_nnz) = 0;
  for (y = 0 ; y < m->nY ; ++y) {
    int pin;
    for (pin = m->ystart[y] ; pin < m->yend[y] ; ++pin) {
      if (m->pinGNO[pin] == m->yGNO[y]) {
	(*nnz)[(*n_nnz)] = pin;
	(*n_nnz)++;
      }
    }
  }

  if (*n_nnz == 0) ZOLTAN_FREE(nnz); /* Avoid memory leaks */

 End:
  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}


  /* This function removes nnz which are listed as arguments (list of indexes in
   pin* arrays.
   nnz array have to be sorted.
  */
int
Zoltan_Matrix_Delete_nnz(ZZ* zz, Zoltan_matrix* m,
			 const int n_nnz, const int* const nnz)
{
  static char *yo = "Zoltan_Matrix_Delete_nnz";
  int ierr = ZOLTAN_OK;
  int i;
  int y;

  if (n_nnz == 0)
    return (ZOLTAN_OK);

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (m->yend == m->ystart + 1) { /* Cannot do this "efficiently" in compact mode */
    m->yend = (int*)ZOLTAN_MALLOC(m->nY*sizeof(int));
    if (m->nY && m->yend == NULL)
      MEMORY_ERROR;
    memcpy(m->yend, m->ystart+1, m->nY*sizeof(int));
  }

  /* Loop over elements we have to remove */
  for (i = 0, y=0; i < n_nnz ; ) {
    int lenght=0;
    int n_removed = 0;

    while (y < m->nY && !(m->ystart[y] <= nnz[i] && m->yend[y] > nnz[i] )) {
      y++;
    }
    if (y >= m->nY){
      ierr = ZOLTAN_WARN;
      break;
    }

    while (i<n_nnz && nnz[i] < m->yend[y]) {
      if (i+1 < n_nnz) lenght = MIN(nnz[i+1], m->yend[y]);
      else lenght = m->yend[y];

      lenght -= nnz[i]+1; /* We remove at least nnz[i] */
      memmove(m->pinGNO+nnz[i], m->pinGNO+nnz[i]+1, lenght*sizeof(int));
      memmove(m->pinwgt+nnz[i]*m->pinwgtdim, m->pinwgt+(nnz[i]+1)*m->pinwgtdim,
	     lenght*sizeof(float)*m->pinwgtdim);
      n_removed ++;
      i++;
    }
    m->yend[y] -= n_removed;
  }
  m->nPins -= n_nnz;

 End:
  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}



int
Zoltan_Matrix_Complete(ZZ* zz,Zoltan_matrix* m)
{
  static char *yo = "Zoltan_Matrix_Complete";
  int ierr = ZOLTAN_OK;

  if(m->completed)
    return (ZOLTAN_OK);

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (m->yend != m->ystart +1 ) {/* Not in compact mode yet */
    int y;
    /* I have to also rewrites all the pinarrays */

    for (y = 1 ; y <= m->nY ; ++y) {
      int length;
      if (m->ystart[y] == m->yend[y-1]) /* No hole */
	continue;
      length = m->yend[y]-m->ystart[y];
      memcpy(m->pinGNO+m->yend[y-1], m->pinGNO+m->ystart[y], length*sizeof(int));
      memcpy(m->pinwgt+m->yend[y-1]*m->pinwgtdim,
	     m->pinGNO+m->ystart[y]*m->pinwgtdim, length*sizeof(float)*m->pinwgtdim);
      m->ystart[y] = m->yend[y-1];
      m->yend[y] = m->ystart[y] + length;
    }

    ZOLTAN_FREE(&m->yend);
    m->yend = m->ystart + 1;
  }

  /* Update data directories */
  m->yGID = ZOLTAN_MALLOC_GID_ARRAY(zz, m->nY);
  m->ywgt = (float*) ZOLTAN_MALLOC(m->nY * sizeof(float) * m->ywgtdim);
  if (m->nY && (m->yGID == NULL || (m->ywgtdim && m->ywgt == NULL)))
    MEMORY_ERROR;

  /* Get Informations about Y */
  Zoltan_DD_Find (m->ddY, (ZOLTAN_ID_PTR)m->yGNO, m->yGID, (ZOLTAN_ID_PTR)m->ywgt, NULL,
		  m->nY, NULL);

  if (m->ddY != m->ddX) {
    Zoltan_DD_Destroy(&m->ddY);
    m->ddY = NULL;
  }

  m->completed = 1;
 End:
  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}

/* Performs a permutation of the matrix, perm_y A perm_y^t.
 * TODO: at this time we only do symmetric permutations (don't know xGNO !).
 */
int
Zoltan_Matrix_Permute(ZZ* zz, Zoltan_matrix *m, const int* const perm_y)
{
  static char *yo = "Zoltan_Matrix_Permute";
  int ierr = ZOLTAN_OK;
  int *pinGNO = NULL;
  ZOLTAN_ID_PTR yGID=NULL;
  float *ywgt=NULL;
  struct Zoltan_DD_Struct *dd;
  
  ZOLTAN_TRACE_ENTER(zz, yo);

  /* First apply y permutation */
  if (m->completed) { /* We directly know the good arrays */
    yGID = m->yGID;
    ywgt = m->ywgt;

    if (m->ddY == NULL || m->ddY != m->ddX) { /* We have to create again the DD */
      /* We have to define ddY : yGNO, yGID, ywgt */
      ierr = Zoltan_DD_Create (&m->ddY, zz->Communicator, 1, zz->Num_GID,
			       m->ywgtdim*sizeof(float)/sizeof(int), m->globalY/zz->Num_Proc, 0);
      /* Hope a linear assignment will help a little */
      Zoltan_DD_Set_Neighbor_Hash_Fn1(m->ddY, m->globalY/zz->Num_Proc);
    }
  }
  else { /* We have to get these fields */
    /* Update data directories */
    yGID = ZOLTAN_MALLOC_GID_ARRAY(zz, m->nY);
    ywgt = (float*) ZOLTAN_MALLOC(m->nY * sizeof(float) * m->ywgtdim);
    if (m->nY && (yGID == NULL || (m->ywgtdim && ywgt == NULL)))
      MEMORY_ERROR;
    /* Get Informations about Y */
    Zoltan_DD_Find (m->ddY, (ZOLTAN_ID_PTR)m->yGNO, yGID, (ZOLTAN_ID_PTR)ywgt, NULL,
		    m->nY, NULL);
  }
  
  memcpy (m->yGNO, perm_y, m->nY*sizeof(int));

  /* Get Informations about Y */
  Zoltan_DD_Update (m->ddY, (ZOLTAN_ID_PTR)m->yGNO, yGID, (ZOLTAN_ID_PTR)ywgt, NULL,
		    m->nY);
  ZOLTAN_FREE (&yGID);
  ZOLTAN_FREE (&ywgt);

  /* We have to define dd : old_yGNO, new_yGNO */
  ierr = Zoltan_DD_Create (&dd, zz->Communicator, 1, 1, 0, m->globalY/zz->Num_Proc, 0);
  /* Hope a linear assignment will help a little */
  Zoltan_DD_Set_Neighbor_Hash_Fn1(dd, m->globalY/zz->Num_Proc);

  Zoltan_DD_Update (dd, (ZOLTAN_ID_PTR)m->yGNO, perm_y, NULL, NULL, m->nY);

  pinGNO = (int*)ZOLTAN_MALLOC(m->nPins*sizeof(int));
  if (m->nPins && pinGNO == NULL)
    MEMORY_ERROR;

  Zoltan_DD_Find (dd, (ZOLTAN_ID_PTR)m->pinGNO, (ZOLTAN_ID_PTR)pinGNO, NULL, NULL,
		  m->nY, NULL);

  ZOLTAN_FREE(&m->pinGNO);
  m->pinGNO = pinGNO;
  pinGNO = NULL;
  
 End:
  ZOLTAN_FREE (&pinGNO);
  ZOLTAN_FREE (&yGID);
  ZOLTAN_FREE (&ywgt);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

int
Zoltan_Matrix2d_adjproc (ZZ* zz, const Zoltan_matrix_2d * const mat, int **adjproc)
{
  static char *yo = "Zoltan_Matrix2d_adjproc";
  int ierr = ZOLTAN_OK;
  int i;
  ZOLTAN_TRACE_ENTER(zz, yo);

  *adjproc = (int*) ZOLTAN_MALLOC(mat->mtx.nPins*sizeof(int));
  if (mat->mtx.nPins && (*adjproc == NULL))
    MEMORY_ERROR;

  for (i = 0 ; i < mat->mtx.nPins ; ++i ) {
    (*adjproc)[i] = EDGE_TO_PROC_Y(mat, mat->mtx.pinGNO[i]);
  }

 End:
  ZOLTAN_TRACE_EXIT(zz, yo);

  return (ierr);
}



#ifdef __cplusplus
}
#endif
