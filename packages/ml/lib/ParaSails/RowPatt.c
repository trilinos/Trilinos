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
 * RowPatt - Pattern of a row, and functions to manipulate the pattern of
 * a row, particularly merging a pattern with a set of nonzero indices.
 *
 * Implementation and Notes: a full-length array is used to mark nonzeros
 * in the pattern.  Indices must not equal -1, which is the "empty" marker
 * used in the full length array.  It is expected that RowPatt will only be 
 * presented with local indices, otherwise the full length array may be very 
 * large.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include "Common.h"
#include "RowPatt.h"

/*--------------------------------------------------------------------------
 * resize - local function for automatically increasing the size of RowPatt
 *--------------------------------------------------------------------------*/

static void resize(RowPatt *p, int newlen)
{
    int oldlen, i;

#ifdef PARASAILS_DEBUG
    printf("RowPatt resize %d\n", newlen);
#endif

    oldlen = p->maxlen;
    p->maxlen = newlen;

    p->ind  = (int *) realloc(p->ind,  p->maxlen * sizeof(int));
    p->mark = (int *) realloc(p->mark, p->maxlen * sizeof(int));

    /* initialize the new portion of the mark array */
    for (i=oldlen; i<p->maxlen; i++)
	p->mark[i] = -1;
}

/*--------------------------------------------------------------------------
 * RowPattCreate - Return (a pointer to) a pattern of a row with a maximum
 * of "maxlen" nonzeros.
 *--------------------------------------------------------------------------*/

RowPatt *RowPattCreate(int maxlen)
{
    int i;
    RowPatt *p = (RowPatt *) malloc(sizeof(RowPatt));

    p->maxlen   = maxlen;
    p->len      = 0;
    p->prev_len = 0;
    p->ind      = (int *) malloc(maxlen * sizeof(int));
    p->mark     = (int *) malloc(maxlen * sizeof(int));
    p->buffer   = NULL;
    p->buflen   = 0;

    for (i=0; i<maxlen; i++)
        p->mark[i] = -1;

    return p;
}

/*--------------------------------------------------------------------------
 * RowPattDestroy - Destroy a row pattern object "p".
 *--------------------------------------------------------------------------*/

void RowPattDestroy(RowPatt *p)
{
    free(p->ind);
    free(p->mark);
    free(p);
}

/*--------------------------------------------------------------------------
 * RowPattReset - Empty the pattern of row pattern object "p".
 *--------------------------------------------------------------------------*/

void RowPattReset(RowPatt *p)
{
    int i;

    for (i=0; i<p->len; i++)
        p->mark[p->ind[i]] = -1;

    p->len      = 0;
    p->prev_len = 0;
}

/*--------------------------------------------------------------------------
 * RowPattMerge - Merge the "len" nonzeros in array "ind" with pattern "p".
 *--------------------------------------------------------------------------*/

void RowPattMerge(RowPatt *p, int len, int *ind)
{
    int i;

    for (i=0; i<len; i++)
    {
	if (ind[i] >= p->maxlen)
	    resize(p, ind[i]*2);

	if (p->mark[ind[i]] == -1)
	{
	    assert(p->len < p->maxlen);

	    p->mark[ind[i]] = p->len;
            p->ind[p->len] = ind[i];
            p->len++;
	}
    }
}

/*--------------------------------------------------------------------------
 * RowPattMergeExt - Merge the external nonzeros in the array "ind" of 
 * length "len" with the pattern "p".  The external indices are those
 * that are less than "beg" or greater than "end".
 *--------------------------------------------------------------------------*/

void RowPattMergeExt(RowPatt *p, int len, int *ind, int num_loc)
{
    int i;

    for (i=0; i<len; i++)
    {
        if (ind[i] < num_loc)
	    continue;

	if (ind[i] >= p->maxlen)
	    resize(p, ind[i]*2);

	if (p->mark[ind[i]] == -1)
	{
	    assert(p->len < p->maxlen);

	    p->mark[ind[i]] = p->len;
            p->ind[p->len] = ind[i];
            p->len++;
	}
    }
}

/*--------------------------------------------------------------------------
 * RowPattGet - Return the pattern of "p".  The length and pointer to the
 * pattern indices are returned through the parameters "lenp" and "indp".
 * A copy of the indices is returned; this copy is destroyed on the next
 * call to RowPattGet or RowPattPrevLevel.
 *--------------------------------------------------------------------------*/

void RowPattGet(RowPatt *p, int *lenp, int **indp)
{
    int len;

    len = p->len;

    if (len > p->buflen)
    {
	free(p->buffer);
	p->buflen = len + 100;
	p->buffer = (int *) malloc(p->buflen * sizeof(int));
    }

    memcpy(p->buffer, p->ind, len*sizeof(int));

    *lenp = len;
    *indp = p->buffer;
}

/*--------------------------------------------------------------------------
 * RowPattPrevLevel - Return the new indices added to the pattern of "p"
 * since the last call to RowPattPrevLevel (or all the indices if never
 * called).  The length and pointer to the pattern indices are returned 
 * through the parameters "lenp" and "indp".
 * A copy of the indices is returned; this copy is destroyed on the next
 * call to RowPattGet or RowPattPrevLevel.
 *--------------------------------------------------------------------------*/

void RowPattPrevLevel(RowPatt *p, int *lenp, int **indp)
{
    int len;

    len = p->len - p->prev_len;

    if (len > p->buflen)
    {
	free(p->buffer);
	p->buflen = len + 100;
	p->buffer = (int *) malloc(p->buflen * sizeof(int));
    }

    memcpy(p->buffer, &p->ind[p->prev_len], len*sizeof(int));

    *lenp = len;
    *indp = p->buffer;

    p->prev_len = p->len;
}
