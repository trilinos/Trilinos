/*BHEADER**********************************************************************
 * (c) 1999   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/
/******************************************************************************
 *
 * PrunedRows.h header file.
 *
 *****************************************************************************/

#include <stdio.h>
#include "Mem.h"
#include "DiagScale.h"

#ifndef _PRUNEDROWS_H
#define _PRUNEDROWS_H

typedef struct
{
    Mem      *mem;   /* storage for arrays, indices, and values */
    int      size;

    int     *len;
    int    **ind;
}
PrunedRows;

PrunedRows *PrunedRowsCreate(Matrix *mat, int size, DiagScale *diag_scale,
  double thresh);
void PrunedRowsDestroy(PrunedRows *p);
int *PrunedRowsAlloc(PrunedRows *p, int len);
void PrunedRowsPut(PrunedRows *p, int index, int len, int *ind);
void PrunedRowsGet(PrunedRows *p, int index, int *lenp, int **indp);

#endif /* _PRUNEDROWS_H */
