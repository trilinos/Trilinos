/*BHEADER**********************************************************************
 * (c) 2002   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
*********************************************************************EHEADER*/

/*
 * Crout-form incomplete factorization
 *
 * To create a Matlab interface for icrout_cholesky_mex, define the symbol
 * MATLAB and compile this file as a Matlab mex program.
 * 
 * Contact Info for this code
 * --------------------------
 * Edmond Chow, Lawrence Livermore National Laboratory, echow@llnl.gov
 *
 * Revision History
 * ----------------
 * 06/04/02 Cleaned up for Boeing
 * 11/06/02 Added icrout_cholesky_mex
 */

#undef IFPACK

#define SYMSTR 1 /* structurally symmetric version */

#ifdef MATLAB

#include "matrix.h"
#include "mex.h"
#define malloc mxMalloc
#define calloc mxCalloc
#define realloc mxRealloc
#define free mxFree
#define printf mexPrintf

#else

#include <stdio.h>

#endif

#include "Ifpack_config.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define MIN(a,b) ((a)<=(b) ? (a) : (b))
#define MAX(a,b) ((a)>=(b) ? (a) : (b))
#define ABS(a) ((a)>=0 ? (a) : -(a))

/*
 * Data structure for sparse matrices is CSR, 0-based indexing.
 */
typedef struct {
    double *val;  /* also known as A  */
    int    *col;  /* also known as JA; first column is column 0 */
    int    *ptr;  /* also known as IA; with ptr[0] = 0 */
} Matrix;

IFPACK_DEPRECATED void ifpack_quicksort (int *const pbase, double *const daux, size_t total_elems);

#define SHORTCUT(p, a, ja, ia) \
        (a)  = (p)->val; \
        (ja) = (p)->col; \
        (ia) = (p)->ptr;

#define MATNULL(p) \
        (p).val = NULL; \
        (p).col = NULL; \
        (p).ptr = NULL;

IFPACK_DEPRECATED void Matrix_alloc(Matrix *a, int n, int nnz)
{
    a->val = (double *) malloc(nnz * sizeof(double));
    a->col = (int *) malloc(nnz * sizeof(int));
    a->ptr = (int *) malloc((n+1) * sizeof(int));
}
 
IFPACK_DEPRECATED void Matrix_dealloc(Matrix *a)
{
    free(a->val);
    free(a->col);
    free(a->ptr);
    a->val = NULL;
    a->col = NULL;
    a->ptr = NULL;
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
        abskey = ABS(a[mid]);
        for (j=first+1; j<=last; j++)
        {
            if (ABS(a[j]) > abskey)
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

IFPACK_DEPRECATED void crout_ict(
    int n,
#ifdef IFPACK
    void * A,
    int maxentries;
    int (*getcol)( void * A, int col, int ** nentries, double * val, int * ind),
    int (*getdiag)( void *A, double * diag),
#else
    const Matrix *AL,
    const double *Adiag,
#endif
    double droptol,
    int lfil,
    Matrix *L,
    double **pdiag)
{
    int k, j, i, index;
    int count_l;
    double norm_l;

    /* work arrays; work_l is list of nonzero values */
    double *work_l = (double *) malloc(n * sizeof(double));
    int    *ind_l  = (int *)    malloc(n * sizeof(int));
    int     len_l;

    /* list and ifirst data structures */
    int *list_l    = (int *) malloc(n * sizeof(int));
    int *ifirst_l  = (int *) malloc(n * sizeof(int));
    
    /* aliases */
    int *marker_l  = ifirst_l;

    /* matrix arrays */
    double *al; int *jal, *ial;
    double *l;  int *jl,  *il;

    int kl = 0;

    double *diag = (double *) malloc(n * sizeof(double));
    *pdiag = diag;

    Matrix_alloc(L,  n, lfil*n);

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
#ifdef IFPACK
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
#ifdef IFPACK
      getcol(A, k, len_l, work_l, ind_l);
      for (j=0; j<len_l; j++)
	{
	  norm_l += ABS(work_l[j]);
	  marker_l[ind_l[j]] = j;
	}
#else
        len_l = 0;
        for (j=ial[k]; j<ial[k+1]; j++)
        {
            index = jal[j];
            work_l[len_l] = al[j];
            norm_l += ABS(al[j]);
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
	    if (ABS(work_l[j]) < droptol * norm_l)
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

	count_l = MIN(len_l, lfil);
	qsplit(work_l, ind_l, len_l, count_l);
	ifpack_quicksort(ind_l, work_l, count_l);

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

    free(work_l);
    free(ind_l);
    free(list_l);
    free(ifirst_l);
}

#ifdef MATLAB
/*
 * [l, d] = icrout_cholesky_mex(AL, Adiag, droptol, lfil)
 * 
 * where AL is a sparse matrix (strictly lower triangular part)
 * Adiag is a vector
 * droptol is a scalar real
 * lfil is a scalar integer
 *
 * d must be converted to a matrix (from a vector) by icrout_cholesky.m
 */
IFPACK_DEPRECATED void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n;
    Matrix AL;
    double *Adiag;
    double droptol;
    int lfil;
    Matrix L;
    double *D;

    if (nrhs != 4)
        mexErrMsgTxt("mex function called with bad number of arguments");

    n = mxGetN(prhs[0]);
    AL.ptr = mxGetJc(prhs[0]);
    AL.col = mxGetIr(prhs[0]);
    AL.val = mxGetPr(prhs[0]);
    Adiag  = mxGetPr(prhs[1]);

    droptol = (double) *mxGetPr(prhs[2]);
    lfil   = (int) *mxGetPr(prhs[3]);

    crout_ict(n, &AL, Adiag, droptol, lfil, &L, &D);

    /* create output matrices */

    /* L */
    plhs[0] = mxCreateSparse(n, n, 1, mxREAL);
    mxFree(mxGetJc(plhs[0]));
    mxFree(mxGetIr(plhs[0]));
    mxFree(mxGetPr(plhs[0]));
    mxSetJc(plhs[0], L.ptr);
    mxSetIr(plhs[0], L.col);
    mxSetPr(plhs[0], L.val);
    mxSetNzmax(plhs[0], n*lfil); /* must agree */

    /* D */
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    mxFree(mxGetPr(plhs[1]));
    mxSetPr(plhs[1], D);
}
#endif
