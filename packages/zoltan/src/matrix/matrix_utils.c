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

  if (m->ddY != m->ddX)
    Zoltan_DD_Destroy(&m->ddY);
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
  MPI_Comm_free(&m->comm->row_comm);
  MPI_Comm_free(&m->comm->col_comm);
}

#ifdef __cplusplus
}
#endif
