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


/************************************/
/* Auxiliary functions declarations */
/************************************/

static int
compar_arcs (const Zoltan_Arc* e1, const Zoltan_Arc* e2);

/* static int */
/* compar_int (const int *e1, const int *e2) { */
/*   return (e1-e2); */
/* } */

/* Functions used when we merge duplicate arcs */

/* Definition of the merge wgt operator : arguments are accumulation wgt, partial wgt, wgt dim */
/* Return !0 when an error is found */
typedef int(*WgtFctPtr)(float*, float*, int);

/* This function compare if the wgt are the same for both arcs*/
static int
wgtFctCmp(float* current, float* new, int dim);

/* This function adds the weights */
static int
wgtFctAdd(float* current, float* new, int dim);

/* This function chooses the maximum weight */
static int
wgtFctMax(float* current, float* new, int dim);

/****************************************/
/* Function definitions are here        */
/****************************************/

static int
compar_arcs (const Zoltan_Arc* e1, const Zoltan_Arc* e2)
{
  if (e1->yGNO == e2->yGNO)
    return (e1->pinGNO - e2->pinGNO);
  return (e1->yGNO - e2->yGNO);
}

/* This function compare if the wgt are the same for both arcs*/
static int
wgtFctCmp(float* current, float* new, int dim)
{
  int i; int diff;

  for (i = 0, diff=0 ; (i <dim)&&(!diff) ; ++i) {
    diff |= (new[i] != current[i]);
  }
  return (diff);
}

/* This function adds the weights */
static int
wgtFctAdd(float* current, float* new, int dim)
{
  int i;
  for (i = 0 ; i <dim ; ++i) {
    current[i] += new[i];
  }
  return (0);
}

/* This function chooses the maximum weight */
static int
wgtFctMax(float* current, float* new, int dim)
{
  int i;
  for (i = 0 ; i <dim ; ++i) {
    current[i] = MAX(current[i],new[i]);
  }
  return (0);
}

/* Function that removes locale duplicated nnz */
/* TODO: Add an option to deal with disconnected vertices */
int
Zoltan_Matrix_Remove_DupArcs(ZZ *zz, int size, Zoltan_Arc *arcs, float* pinwgt,
			     Zoltan_matrix *outmat)
{
  static char *yo = "Zoltan_Matrix_Remove_DupArcs";
  int ierr = ZOLTAN_OK;
  WgtFctPtr wgtfct;
  int nY, nPin;
  int i;
  int prev_pinGNO;

  ZOLTAN_TRACE_ENTER(zz, yo);
/*   if (addweight) */

  switch (outmat->opts.pinwgtop) {
  case MAX_WEIGHT:
    wgtfct = &wgtFctMax;
    break;
  case CMP_WEIGHT:
    wgtfct = &wgtFctCmp;
    break;
  case ADD_WEIGHT:
  default:
    wgtfct = &wgtFctAdd;
  }

  qsort ((void*)arcs, size, sizeof(Zoltan_Arc),
	 (int (*)(const void*,const void*))compar_arcs);

  prev_pinGNO = -2;
  for (i = 0, nY=-1, nPin=-1; i < size ; ++i) {
    int new = 0;
    int yGNO = arcs[i].yGNO;
    int pinGNO = arcs[i].pinGNO;

    new = ((nY < 0) || (outmat->yGNO[nY] != yGNO));
    if (new) {
      nY++;
      outmat->ystart[nY] = nPin + 1;
      outmat->yGNO[nY] = yGNO;
      prev_pinGNO = -1;
      if (pinGNO < 0)
	continue;
    }
    new = new ||(pinGNO != prev_pinGNO);
    if (new) { /* New edge */
      nPin ++;
      outmat->pinGNO[nPin] = pinGNO;
      memcpy(outmat->pinwgt + nPin*outmat->pinwgtdim, pinwgt + arcs[i].offset*outmat->pinwgtdim,
	     outmat->pinwgtdim*sizeof(float));
    }
    else { /* Duplicate */
      wgtfct(outmat->pinwgt + nPin* outmat->pinwgtdim,
	     pinwgt + arcs[i].offset*outmat->pinwgtdim, outmat->pinwgtdim);
    }
    prev_pinGNO = outmat->pinGNO[nPin];
  }
  nY ++;
  outmat->ystart[nY] = nPin+1; /* compact mode */
  outmat->nPins = nPin+1;
  outmat->nY = nY;

  /* Try to minimize memory */
  /* We reduce memory, thus I don't think these realloc can fail */
  if (outmat->yend != outmat->ystart + 1)
    ZOLTAN_FREE(&outmat->yend);

  outmat->pinGNO = (int *) ZOLTAN_REALLOC(outmat->pinGNO, outmat->nPins * sizeof(int));
  outmat->pinwgt = (float *) ZOLTAN_REALLOC(outmat->pinwgt,
			       outmat->nPins*outmat->pinwgtdim*sizeof(float));
  outmat->yend = outmat->ystart + 1;

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
}

/* Function that removes locale duplicated nnz */
int
Zoltan_Matrix_Remove_Duplicates(ZZ *zz, Zoltan_matrix inmat, Zoltan_matrix *outmat)
{
  static char *yo = "Zoltan_Matrix_Remove_Duplicates";
  int ierr = ZOLTAN_OK;
  Zoltan_Arc *arcs = NULL;
  float *pinwgt = NULL;
  int freeflag = 0;
  int size;
  int i, j, cnt;

  ZOLTAN_TRACE_ENTER(zz, yo);

  size = inmat.nPins + inmat.nY; /* We add fake arcs for non connected vertices */
  arcs = (Zoltan_Arc*) ZOLTAN_MALLOC(size*sizeof(Zoltan_Arc));
  if (inmat.nPins && arcs == NULL) MEMORY_ERROR;

  for (i = 0, cnt=0 ; i < inmat.nY; i++) {
    /* Fake arc in order to be sure to keep this vertex */
    arcs[cnt].yGNO = inmat.yGNO[i];
    arcs[cnt].pinGNO = -1;
    arcs[cnt].offset = -1;
    cnt++;
    for (j = inmat.ystart[i]; j < inmat.yend[i]; j++) {
      arcs[cnt].yGNO = inmat.yGNO[i];
      arcs[cnt].pinGNO = inmat.pinGNO[j];
      arcs[cnt].offset = j;
      cnt ++;
    }
  }

  pinwgt = inmat.pinwgt;
  if (pinwgt == outmat->pinwgt) {
    freeflag = 1;
    outmat->pinwgt = (float*) ZOLTAN_MALLOC(inmat.pinwgtdim*inmat.nPins*sizeof(float));
    if (inmat.pinwgtdim && inmat.nPins && outmat->pinwgt == NULL) MEMORY_ERROR;
  }

  ierr = Zoltan_Matrix_Remove_DupArcs(zz, size, arcs, pinwgt,outmat);

  if (freeflag)
    ZOLTAN_FREE(&pinwgt);

 End:
  ZOLTAN_FREE(&arcs);

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ierr);
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




/* Performs a permutation of the matrix, perm_y A perm_y^t.
 * TODO: at this time we only do symmetric permutations (don't know xGNO !).
 */
int
Zoltan_Matrix_Permute(ZZ* zz, Zoltan_matrix *m, int* perm_y)
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

  Zoltan_DD_Update (dd, (ZOLTAN_ID_PTR)m->yGNO, (ZOLTAN_ID_PTR)perm_y, NULL, NULL, m->nY);

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


#ifdef __cplusplus
}
#endif
