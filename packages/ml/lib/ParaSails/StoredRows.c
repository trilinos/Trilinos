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
 * StoredRows - Local storage of rows from other processors.  Although only
 * off-processor rows are stored, if an on-processor row is requested, it
 * is returned by referring to the local matrix.  Local indexing is used to
 * access the stored rows.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include "Common.h"
#include "Mem.h"
#include "Matrix.h"
#include "StoredRows.h"

/*--------------------------------------------------------------------------
 * StoredRowsCreate - Return (a pointer to) a stored rows object.
 *
 * mat  - matrix used for returning on-processor rows (input)
 * size - the maximum number of (off-processor) rows that can be stored 
 *        (input).  See below for more a precise description.
 *
 * A slot is available for "size" off-processor rows.  The slot for the
 * row with local index i is (i - num_loc).  Therefore, if max_i is the
 * largest local index expected, then size should be set to 
 * (max_i - num_loc + 1).  StoredRows will automatically increase its 
 * size if a row with a larger local index needs to be put in StoredRows.
 *--------------------------------------------------------------------------*/

StoredRows *StoredRowsCreate(Matrix *mat, int size)
{
    StoredRows *p = (StoredRows *) malloc(sizeof(StoredRows));

    p->mat  = mat;
    p->mem  = MemCreate();

    p->size = size;
    p->num_loc = mat->end_row - mat->beg_row + 1;

    p->len = (int *)     calloc(size,  sizeof(int));
    p->ind = (int **)    malloc(size * sizeof(int *));
    p->val = (double **) malloc(size * sizeof(double *));

    p->count = 0;

    return p;
}

/*--------------------------------------------------------------------------
 * StoredRowsDestroy - Destroy a stored rows object "p".
 *--------------------------------------------------------------------------*/

void StoredRowsDestroy(StoredRows *p)
{
    MemDestroy(p->mem);
    free(p->len);
    free(p->ind);
    free(p->val);
    free(p);
}

/*--------------------------------------------------------------------------
 * StoredRowsAllocInd - Return space allocated for "len" indices in the
 * stored rows object "p".  The indices may span several rows.
 *--------------------------------------------------------------------------*/

int *StoredRowsAllocInd(StoredRows *p, int len)
{
    return (int *) MemAlloc(p->mem, len*sizeof(int));
}

/*--------------------------------------------------------------------------
 * StoredRowsAllocVal - Return space allocated for "len" values in the
 * stored rows object "p".  The values may span several rows.
 *--------------------------------------------------------------------------*/

double *StoredRowsAllocVal(StoredRows *p, int len)
{
    return (double *) MemAlloc(p->mem, len*sizeof(double));
}

/*--------------------------------------------------------------------------
 * StoredRowsPut - Given a row (len, ind, val), store it as row "index" in
 * the stored rows object "p".  Only nonlocal stored rows should be put using
 * this interface; the local stored rows are put using the create function.
 *--------------------------------------------------------------------------*/

void StoredRowsPut(StoredRows *p, int index, int len, int *ind, double *val)
{
    int i = index - p->num_loc;

    /* Reallocate if necessary */
    if (i >= p->size)
    {
        int j;
        int newsize;

	newsize = i*2;
#ifdef PARASAILS_DEBUG
		    printf("StoredRows resize %d\n", newsize);
#endif
        p->len = (int *)     realloc(p->len, newsize * sizeof(int));
        p->ind = (int **)    realloc(p->ind, newsize * sizeof(int *));
        p->val = (double **) realloc(p->val, newsize * sizeof(double *));

	/* set lengths to zero */
        for (j=p->size; j<newsize; j++)
	    p->len[j] = 0;

        p->size = newsize;
    }

    /* check that row has not been put already */
    assert(p->len[i] == 0);

    p->len[i] = len;
    p->ind[i] = ind;
    p->val[i] = val;

    p->count++;
}

/*--------------------------------------------------------------------------
 * StoredRowsGet - Return the row with index "index" through the pointers 
 * "lenp", "indp" and "valp" in the stored rows object "p".
 *--------------------------------------------------------------------------*/

void StoredRowsGet(StoredRows *p, int index, int *lenp, int **indp, 
  double **valp)
{
    if (index < p->num_loc)
    {
        MatrixGetRow(p->mat, index, lenp, indp, valp);
    }
    else
    {
	index = index - p->num_loc;

        *lenp = p->len[index];
        *indp = p->ind[index];
        *valp = p->val[index];
    }
}
