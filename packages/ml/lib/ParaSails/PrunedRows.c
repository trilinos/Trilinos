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
 * PrunedRows - Collection of pruned rows that are cached on the local 
 * processor.  Direct access to these rows is available, via the local
 * index number.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include "Common.h"
#include "Mem.h"
#include "Matrix.h"
#include "DiagScale.h"
#include "PrunedRows.h"

/*--------------------------------------------------------------------------
 * PrunedRowsCreate - Return (a pointer to) a pruned rows object.
 *
 * mat        - matrix used to construct the local pruned rows (input)
 *              assumes the matrix uses local indexing
 * size       - number of unique local indices on this processor;
 *              an array of this size will be allocated to access the
 *              pruned rows (input) - includes the number of local nodes
 * diag_scale - diagonal scale object used to scale the thresholding (input)
 * thresh     - threshold for pruning the matrix (input)
 *
 * The local pruned rows are stored in the first part of the len and ind 
 * arrays.
 *--------------------------------------------------------------------------*/

PrunedRows *PrunedRowsCreate(Matrix *mat, int size, DiagScale *diag_scale, 
  double thresh)
{
    int row, len, *ind, count, j, *data;
    double *val, temp;

    PrunedRows *p = (PrunedRows *) malloc(sizeof(PrunedRows));

    p->mem  = MemCreate();
    p->size = MAX(size, mat->end_row - mat->beg_row + 1);

    p->len = (int *)  malloc(p->size * sizeof(int));
    p->ind = (int **) malloc(p->size * sizeof(int *));

    /* Prune and store the rows on the local processor */

    for (row=0; row<=mat->end_row - mat->beg_row; row++)
    {
        MatrixGetRow(mat, row, &len, &ind, &val);

        count = 1; /* automatically include the diagonal */
        for (j=0; j<len; j++)
        {
            temp = DiagScaleGet(diag_scale, row);
            if (temp*ABS(val[j])*DiagScaleGet(diag_scale, ind[j]) 
              >= thresh && ind[j] != row)
                count++;
        }

        p->ind[row] = (int *) MemAlloc(p->mem, count*sizeof(int));
        p->len[row] = count;

        data = p->ind[row];
        *data++ = row; /* the diagonal entry */
        for (j=0; j<len; j++)
        {
            temp = DiagScaleGet(diag_scale, row);
            if (temp*ABS(val[j])*DiagScaleGet(diag_scale, ind[j]) 
              >= thresh && ind[j] != row)
                *data++ = ind[j];
        }
    }

    return p;
}

/*--------------------------------------------------------------------------
 * PrunedRowsDestroy - Destroy a pruned rows object "p".
 *--------------------------------------------------------------------------*/

void PrunedRowsDestroy(PrunedRows *p)
{
    MemDestroy(p->mem);
    free(p->len);
    free(p->ind);
    free(p);
}

/*--------------------------------------------------------------------------
 * PrunedRowsAllocInd - Return space allocated for "len" indices in the
 * pruned rows object "p".  The indices may span several rows.
 *--------------------------------------------------------------------------*/

int *PrunedRowsAlloc(PrunedRows *p, int len)
{
    return (int *) MemAlloc(p->mem, len*sizeof(int));
}

/*--------------------------------------------------------------------------
 * PrunedRowsPut - Given a pruned row (len, ind), store it as row "index" in
 * the pruned rows object "p".  Only nonlocal pruned rows should be put using
 * this interface; the local pruned rows are put using the create function.
 *--------------------------------------------------------------------------*/

void PrunedRowsPut(PrunedRows *p, int index, int len, int *ind)
{
    if (index >= p->size)
    {
	p->size = index*2;
#ifdef PARASAILS_DEBUG
	printf("StoredRows resize %d\n", p->size);
#endif
	p->len = (int *)  realloc(p->len, p->size * sizeof(int));
	p->ind = (int **) realloc(p->ind, p->size * sizeof(int *));
    }

    p->len[index] = len;
    p->ind[index] = ind;
}

/*--------------------------------------------------------------------------
 * PrunedRowsGet - Return the row with index "index" through the pointers 
 * "lenp" and "indp" in the pruned rows object "p".
 *--------------------------------------------------------------------------*/

void PrunedRowsGet(PrunedRows *p, int index, int *lenp, int **indp)
{
    *lenp = p->len[index];
    *indp = p->ind[index];
}
