/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include <math.h>
#include "superlu_ddefs.h"


static void create_msr_matrix(SuperMatrix *, int_t [], int_t,
			      double **, int_t **);
static void PrintMSRmatrix(int, double [], int_t [], gridinfo_t *);


int pdgsmv_AXglobal_setup
(
 SuperMatrix *A,       /* Matrix A permuted by columns (input).
			  The type of A can be:
			  Stype = NCP; Dtype = D; Mtype = GE. */
 Glu_persist_t *Glu_persist, /* input */
 gridinfo_t *grid,     /* input */
 int_t *m,             /* output */
 int_t *update[],      /* output */
 double *val[],        /* output */
 int_t *bindx[],       /* output */
 int_t *mv_sup_to_proc /* output */
 )
{
    int n;
    int input_option;
    int N_update;    /* Number of variables updated on this process (output) */
    int iam = grid->iam;
    int nprocs = grid->nprow * grid->npcol;
    int_t *xsup = Glu_persist->xsup;
    int_t *supno = Glu_persist->supno;
    int_t nsupers;
    int i, nsup, p, t1, t2, t3;


    /* Initialize the list of global indices.
     * NOTE: the list of global indices must be in ascending order.
     */
    n = A->nrow;
    input_option = SUPER_LINEAR;
    nsupers = supno[n-1] + 1;

#if ( DEBUGlevel>=1 )
    if ( !iam ) {
	PrintInt10("xsup", supno[n-1]+1, xsup);
	PrintInt10("supno", n, supno);
    }
#endif

    if ( input_option == SUPER_LINEAR ) { /* Block partitioning based on
					     individual rows.  */
	/* Figure out mv_sup_to_proc[] on all processes. */
	for (p = 0; p < nprocs; ++p) {
	    t1 = n / nprocs;       /* Number of rows */
	    t2 = n - t1 * nprocs;  /* left-over, which will be assigned
				      to the first t2 processes.  */
	    if ( p >= t2 ) t2 += (p * t1); /* Starting row number */
	    else { /* First t2 processes will get one more row. */
 	        ++t1;              /* Number of rows. */
		t2 = p * t1;       /* Starting row. */
	    }
	    /* Make sure the starting and ending rows are at the
	       supernode boundaries. */
	    t3 = t2 + t1;      /* Ending row. */
	    nsup = supno[t2];
	    if ( t2 > xsup[nsup] ) { /* Round up the starting row. */
		t1 -= xsup[nsup+1] - t2;
		t2 = xsup[nsup+1];
	    }
	    nsup = supno[t3];
	    if ( t3 > xsup[nsup] ) /* Round up the ending row. */
		t1 += xsup[nsup+1] - t3;
	    t3 = t2 + t1 - 1;
	    if ( t1 ) {
		for (i = supno[t2]; i <= supno[t3]; ++i) {
		    mv_sup_to_proc[i] = p;
#if ( DEBUGlevel>=1 )
		    if ( mv_sup_to_proc[i] == p-1 ) {
			fprintf(stderr, 
				"mv_sup_to_proc conflicts at supno %d\n", i);
			exit(-1);
		    }
#endif
		}
	    }
	    
	    if ( iam == p ) {
		N_update = t1;
		if ( N_update ) {
		    if ( !(*update = intMalloc_dist(N_update)) )
			ABORT("Malloc fails for update[]");
		}
		for (i = 0; i < N_update; ++i) (*update)[i] = t2 + i;
#if ( DEBUGlevel>=1 )
		printf("(%2d) N_update = %4d\t"
		       "supers %4d to %4d\trows %4d to %4d\n",
		       iam, N_update, supno[t2], supno[t3], t2, t3);
#endif
	    }
	} /* for p ... */
    } else if ( input_option == SUPER_BLOCK ) { /* Block partitioning based on
						   individual supernodes.  */
	/* This may cause bad load balance, because the blocks are usually
	   small in the beginning and large toward the end.   */
	t1 = nsupers / nprocs;
	t2 = nsupers - t1 * nprocs; /* left-over */
	if ( iam >= t2 ) t2 += (iam * t1);
	else {
	    ++t1;          /* Number of blocks. */
	    t2 = iam * t1; /* Starting block. */
	}
	N_update = xsup[t2+t1] - xsup[t2];
	if ( !(*update = intMalloc_dist(N_update)) )
	    ABORT("Malloc fails for update[]");
	for (i = 0; i < N_update; ++i) (*update)[i] = xsup[t2] + i;
    }


    /* Create an MSR matrix in val/bindx to be used by pdgsmv(). */
    create_msr_matrix(A, *update, N_update, val, bindx);

#if ( DEBUGlevel>=1 )
    PrintInt10("mv_sup_to_proc", nsupers, mv_sup_to_proc);
    PrintMSRmatrix(N_update, *val, *bindx, grid);
#endif

    *m = N_update;
    return 0;
} /* PDGSMV_Aglobal_SETUP */


/* Create the distributed modified sparse row (MSR) matrix: bindx/val.
 * For a submatrix of size m-by-n, the MSR arrays are as follows:
 *    bindx[0]      = m + 1
 *    bindx[0..m]   = pointer to start of each row
 *    bindx[ks..ke] = column indices of the off-diagonal nonzeros in row k,
 *                    where, ks = bindx[k], ke = bindx[k+1]-1
 *    val[k]        = A(k,k), k < m, diagonal elements
 *    val[m]        = not used
 *    val[ki]       = A(k, bindx[ki]), where ks <= ki <= ke
 * Both arrays are of length nnz + 1.
 */
static void create_msr_matrix
(
 SuperMatrix *A,       /* Matrix A permuted by columns (input).
			  The type of A can be:
			  Stype = NCP; Dtype = D; Mtype = GE. */
 int_t update[],       /* input (local) */
 int_t N_update,       /* input (local) */
 double **val,         /* output */
 int_t **bindx           /* output */
)
{
    int hi, i, irow, j, lo, n, nnz_local;
    NCPformat *Astore;
    double *nzval;
    int_t *rowcnt;
    
    if ( !N_update ) return;

    n = A->ncol;
    Astore = A->Store;
    nzval = Astore->nzval;

    /* One pass of original matrix A to count nonzeros of each row. */
    if ( !(rowcnt = (int_t *) intCalloc_dist(N_update)) )
	ABORT("Malloc fails for rowcnt[]");
    lo = update[0];
    hi = update[N_update-1];
    nnz_local = 0;
    for (j = 0; j < n; ++j)
	for (i = Astore->colbeg[j]; i < Astore->colend[j]; ++i) {
	    irow = Astore->rowind[i];
	    if ( irow >= lo && irow <= hi ) {
		if ( irow != j ) /* Exclude diagonal */
		    ++rowcnt[irow - lo];
		++nnz_local;
	    }
	}

    /* Allocate storage for bindx[] and val[]. */
    if ( !(*val = (double *) doubleMalloc_dist(nnz_local+1)) )
	ABORT("Malloc fails for val[]");
    if ( !(*bindx = (int_t *) SUPERLU_MALLOC((nnz_local+1) * sizeof(int_t))) )
	ABORT("Malloc fails for bindx[]");

    /* Set up row pointers. */
    (*bindx)[0] = N_update + 1;
    for (j = 1; j <= N_update; ++j) {
	(*bindx)[j] = (*bindx)[j-1] + rowcnt[j-1];
	rowcnt[j-1] = (*bindx)[j-1];
    }

    /* One pass of original matrix A to fill in matrix entries. */
    for (j = 0; j < n; ++j)
	for (i = Astore->colbeg[j]; i < Astore->colend[j]; ++i) {
	    irow = Astore->rowind[i];
	    if ( irow >= lo && irow <= hi ) {
		if ( irow == j ) /* Diagonal */
		    (*val)[irow - lo] = nzval[i];
		else {
		    irow -= lo;
		    (*bindx)[rowcnt[irow]] = j;
		    (*val)[rowcnt[irow]] = nzval[i];
		    ++rowcnt[irow];
		}
	    }
	}

    SUPERLU_FREE(rowcnt);
}

/*
 * Performs sparse matrix-vector multiplication.
 *   - val/bindx stores the distributed MSR matrix A
 *   - X is global
 *   - ax product is distributed the same way as A
 */
int
pdgsmv_AXglobal(int_t m, int_t update[], double val[], int_t bindx[],
	       double X[], double ax[])
{
    int_t i, j, k;

    if ( m <= 0 ) return; /* number of rows (local) */

    for (i = 0; i < m; ++i) {
	ax[i] = 0.0;
	for (k = bindx[i]; k < bindx[i+1]; ++k) {
	    j = bindx[k];       /* column index */
	    ax[i] += val[k] * X[j];
	}
	ax[i] += val[i] * X[update[i]]; /* diagonal */
    }
} /* PDGSMV_AXglobal */
 
/*
 * Performs sparse matrix-vector multiplication.
 *   - val/bindx stores the distributed MSR matrix A
 *   - X is global
 *   - ax product is distributed the same way as A
 */
int
pdgsmv_AXglobal_abs(int_t m, int_t update[], double val[], int_t bindx[],
		   double X[], double ax[])
{
    int_t i, j, k;

    if ( m <= 0 ) return; /* number of rows (local) */

    for (i = 0; i < m; ++i) {
	ax[i] = 0.0;
	for (k = bindx[i]; k < bindx[i+1]; ++k) {
	    j = bindx[k];       /* column index */
	    ax[i] += fabs(val[k]) * fabs(X[j]);
	}
	ax[i] += fabs(val[i]) * fabs(X[update[i]]); /* diagonal */
    }
} /* PDGSMV_AXglobal_ABS */

/*
 * Print the local MSR matrix
 */
static void PrintMSRmatrix
(
 int m,       /* Number of rows of the submatrix. */
 double val[],
 int_t bindx[],
 gridinfo_t *grid
)
{
    int iam;
    int_t nnzp1;

    if ( !m ) return;

    iam = grid->iam;
    nnzp1 = bindx[m];
    printf("(%2d) MSR submatrix has %d rows -->\n", iam, m);
    PrintDouble5("val", nnzp1, val);
    PrintInt10("bindx", nnzp1, bindx);
}
