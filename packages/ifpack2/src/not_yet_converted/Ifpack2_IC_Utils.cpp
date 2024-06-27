// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_IC_Utils.hpp"
#include <cmath>

#define SYMSTR 1 /* structurally symmetric version */
#include <stdio.hpp>

#define SHORTCUT(p, a, ja, ia) \
        (a)  = (p)->val; \
        (ja) = (p)->col; \
        (ia) = (p)->ptr;

#define MATNULL(p) \
        (p).val = NULL; \
        (p).col = NULL; \
        (p).ptr = NULL;

void Ifpack2_AIJMatrix_alloc(Ifpack2_AIJMatrix *a, int n, int nnz)
{
    a->val = new double[nnz];
    a->col = new int[nnz];
    a->ptr = new int[n+1];
}
 
void Ifpack2_AIJMatrix_dealloc(Ifpack2_AIJMatrix *a)
{
    delete [] (a->val);
    delete [] (a->col);
    delete [] (a->ptr);
    a->val = 0;
    a->col = 0;
    a->ptr = 0;
}

static void qsplit(double *a, int *ind, int n, int ncut)
{
    double tmp, abskey;
    int itmp, first, last, mid;
    int j;
 
    ncut--;
    first = 0;
    last = n-1;
    if (ncut < first || ncut > last) 
        return;
 
    /* outer loop while mid != ncut */
    while (1)
    {
        mid = first;
        abskey = std::abs(a[mid]);
        for (j=first+1; j<=last; j++)
        {
            if (std::abs(a[j]) > abskey)
            {
                mid = mid+1;
                /* interchange */
                tmp = a[mid];
                itmp = ind[mid];
                a[mid] = a[j];
                ind[mid] = ind[j];
                a[j]  = tmp;
                ind[j] = itmp;
            }
        }
       
        /* interchange */
        tmp = a[mid];
        a[mid] = a[first];
        a[first]  = tmp;
        itmp = ind[mid];
        ind[mid] = ind[first];
        ind[first] = itmp;
       
        /* test for while loop */
        if (mid == ncut) 
            return;

        if (mid > ncut)
            last = mid-1;
        else
            first = mid+1;
    }
}

/* update column k using previous columns */
/* assumes that column of A resides in the work vector */
/* this function can also be used to update rows */
 
static void update_column(
    int k,
    const int *ia, const int *ja, const double *a,
    const int *ifirst,
    const int *ifirst2,
    const int *list2,
    const double *multipliers, /* the val array of the other factor */
    const double *d, /* diagonal of factorization */
    int *marker,
    double *ta,
    int *itcol,
    int *ptalen)
{
    int j, i, isj, iej, row;
    int talen, pos;
    double mult;
 
    talen = *ptalen;

    j = list2[k];
    while (j != -1)
    {
        /* update column k using column j */
 
        isj = ifirst[j];

        /* skip first nonzero in column, since it does not contribute */
        /* if symmetric structure */
        /* isj++; */
 
        /* now do the actual update */
       if (isj != -1)
       {
        /* multiplier */
        mult = multipliers[ifirst2[j]] * d[j];

        /* section of column used for update */
        iej = ia[j+1]-1;

        for (i=isj; i<=iej; i++)
        {
            row = ja[i];

            /* if nonsymmetric structure */
            if (row <= k)
                continue;
 
            if ((pos = marker[row]) != -1)
            {
                /* already in pattern of working row */
                ta[pos] -= mult*a[i];
            }
            else
            {
                /* not yet in pattern of working row */
                itcol[talen] = row;
                ta[talen] = -mult*a[i];
                marker[row] = talen++;
            }
        }
       }
 
        j = list2[j];
    }

    *ptalen = talen;
}

/* update ifirst and list */
 
static void update_lists(
    int k,
    const int *ia,
    const int *ja,
    int *ifirst,
    int *list)
{
    int j, isj, iej, iptr;

    j = list[k];
    while (j != -1)
    {
        isj = ifirst[j];
        iej = ia[j+1]-1;
 
        isj++;
 
        if (isj != 0 && isj <= iej) /* nonsymmetric structure */
        {
            /* increment ifirst for column j */
            ifirst[j] = isj;

            /* add j to head of list for list[ja[isj]] */
            iptr = j;
            j = list[j];
            list[iptr] = list[ja[isj]];
            list[ja[isj]] = iptr;
        }
        else
        {
            j = list[j];
        }
    }
}

static void update_lists_newcol(
    int k,
    int isk,
    int iptr,
    int *ifirst,
    int *list)
{
    ifirst[k] = isk;
    list[k] = list[iptr];
    list[iptr] = k;
}

/*
 * crout_ict - Crout version of ICT - Incomplete Cholesky by Threshold
 *
 * The incomplete factorization L D L^T is computed,
 * where L is unit triangular, and D is diagonal
 *
 * INPUTS
 * n = number of rows or columns of the matrices
 * AL = unit lower triangular part of A stored by columns
 *            the unit diagonal is implied (not stored)
 * Adiag = diagonal part of A
 * droptol = drop tolerance
 * lfil  = max nonzeros per col in L factor or per row in U factor
 *
 * OUTPUTS
 * L     = lower triangular factor stored by columns
 * pdiag = diagonal factor stored in an array
 *
 * NOTE: calling function must free the memory allocated by this
 * function, i.e., L, pdiag.
 */

void crout_ict(
    int n,
#ifdef TIFPACK
    void * A,
    int maxentries;
    int (*getcol)( void * A, int col, int ** nentries, double * val, int * ind),
    int (*getdiag)( void *A, double * diag),
#else
    const Ifpack2_AIJMatrix *AL,
    const double *Adiag,
#endif
    double droptol,
    int lfil,
    Ifpack2_AIJMatrix *L,
    double **pdiag)
{
    int k, j, i, index;
    int count_l;
    double norm_l;

    /* work arrays; work_l is list of nonzero values */
    double *work_l = new double[n];
    int    *ind_l  = new int[n];
    int     len_l;

    /* list and ifirst data structures */
    int *list_l    = new int[n];
    int *ifirst_l  = new int[n];
    
    /* aliases */
    int *marker_l  = ifirst_l;

    /* matrix arrays */
    double *al; int *jal, *ial;
    double *l;  int *jl,  *il;

    int kl = 0;

    double *diag = new double[n];
    *pdiag = diag;

    Ifpack2_AIJMatrix_alloc(L,  n, lfil*n);

    SHORTCUT(AL, al, jal, ial);
    SHORTCUT(L,  l,  jl,  il);

    /* initialize work arrays */
    for (i=0; i<n; i++)
    {
        list_l[i]    = -1;
        ifirst_l[i]  = -1;
        marker_l[i]  = -1;
    }

    /* copy the diagonal */
#ifdef TIFPACK
    getdiag(A, diag);
#else
    for (i=0; i<n; i++)
        diag[i] = Adiag[i];
#endif

    /* start off the matrices right */
    il[0]  = 0;

    /*
     * Main loop over columns and rows
     */

    for (k=0; k<n; k++)
    {
        /*
         * lower triangular factor update by columns
         * (need ifirst for L and lists for U)
         */
 
        /* copy column of A into work vector */
        norm_l = 0.;
#ifdef TIFPACK
      getcol(A, k, len_l, work_l, ind_l);
      for (j=0; j<len_l; j++)
	{
	  norm_l += std::abs(work_l[j]);
	  marker_l[ind_l[j]] = j;
	}
#else
        len_l = 0;
        for (j=ial[k]; j<ial[k+1]; j++)
        {
            index = jal[j];
            work_l[len_l] = al[j];
            norm_l += std::abs(al[j]);
            ind_l[len_l] = index;
            marker_l[index] = len_l++; /* points into work array */
        }
#endif
        norm_l = (len_l == 0) ? 0.0 : norm_l/((double) len_l);
 
        /* update and scale */

        update_column(k, il, jl, l, ifirst_l, ifirst_l, list_l, l,
            diag, marker_l, work_l, ind_l, &len_l);

        for (j=0; j<len_l; j++)
            work_l[j] /= diag[k];

	i = 0;
        for (j=0; j<len_l; j++)
	{
	    if (std::abs(work_l[j]) < droptol * norm_l)
	    {
	        /* zero out marker array for this */
		marker_l[ind_l[j]] = -1;
	    }
	    else
	    {
	        work_l[i] = work_l[j];
		ind_l[i]  = ind_l[j];
		i++;
	    }
	}
	len_l = i;

	/*
	 * dropping:  for each work vector, perform qsplit, and then
	 * sort each part by increasing index; then copy each sorted
	 * part into the factors
	 */

	count_l = std::min(len_l, lfil);
	qsplit(work_l, ind_l, len_l, count_l);
	quicksort(ind_l, work_l, count_l);

	for (j=0; j<count_l; j++)
	{
	    l[kl] = work_l[j];
	    jl[kl++] = ind_l[j];
	}
	il[k+1] = kl;

	/*
	 * update lists
	 */
	
        update_lists(k, il, jl, ifirst_l, list_l);

	if (kl - il[k] > SYMSTR)
	    update_lists_newcol(k, il[k], jl[il[k]], ifirst_l, list_l);

        /* zero out the marker arrays */
        for (j=0; j<len_l; j++)
            marker_l[ind_l[j]] = -1;

        /*
         * update the diagonal (after dropping)
         */

        for (j=0; j<count_l; j++)
        {
            index = ind_l[j];
            diag[index] -= work_l[j] * work_l[j] * diag[k];
        }
    }

    delete [] work_l;
    delete [] ind_l;
    delete [] list_l;
    delete [] ifirst_l;
}
