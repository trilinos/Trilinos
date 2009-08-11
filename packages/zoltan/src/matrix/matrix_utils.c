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
  /* TODO : Update to free the whole PHGComm struct */
/*   MPI_Comm_free(&m->comm->row_comm); */
/*   MPI_Comm_free(&m->comm->col_comm); */
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

int
Zoltan_Matrix_Complete(ZZ* zz,Zoltan_matrix* m)
{
  static char *yo = "Zoltan_Matrix_Complete";
  int ierr = ZOLTAN_OK;

  if(m->completed)
    return (ZOLTAN_OK);

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (m->yend != m->ystart +1 ) {/* Not in compact mode yet */
    /* Not properly handled yet */
    /* I have to also rewrites all the pinarrays */

/*     ystart = (int*) ZOLTAN_MALLOC((m->nY +1)*sizeof(int)); */
/*     if (ystart == NULL) MEMORY_ERROR; */
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


#ifdef __cplusplus
}
#endif
