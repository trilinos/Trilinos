/*
 * -- SuperLU routine (version 1.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
#include <stdio.h>
#include "mex.h"
#include "supermatrix.h"
#include "util.h"
#define  MatlabMatrix Matrix

/* Aliases for input and output arguments */
#define A_in		prhs[0]
#define Pc_in		prhs[1]
#define L_out    	plhs[0]
#define U_out          	plhs[1]
#define Pr_out     	plhs[2]
#define Pc_out   	plhs[3]

void LUextract(SuperMatrix *, SuperMatrix *, double *, int *, int *, 
	       double *, int *, int *, int *, int*);

#define verbose (SPUMONI>0)
#define babble  (SPUMONI>1)
#define burble  (SPUMONI>2)

void mexFunction(
    int          nlhs,           /* number of expected outputs */
    MatlabMatrix *plhs[],        /* matrix pointer array returning outputs */
    int          nrhs,           /* number of inputs */
    MatlabMatrix *prhs[]         /* matrix pointer array for inputs */
    )
{
    int SPUMONI;             /* ... as should the sparse monitor flag */
    Real FlopsInSuperLU;     /* ... as should the flop counter */
    extern flops_t LUFactFlops();
    
    /* Arguments to C dgstrf(). */
    SuperMatrix A;
    SuperMatrix Ac;        /* Matrix postmultiplied by Pc */
    SuperMatrix L, U;
    int	   	m, n, nnz;
    double      *val;
    int       	*rowind;
    int		*colptr;
    int    	*etree, *perm_r, *perm_c;
    char        refact[1], trans[1];
    int         panel_size, relax;
    double      thresh = 1.0;       /* diagonal pivoting threshold */
    double      drop_tol = 0.0;     /* drop tolerance parameter */
    int		info;
    MatlabMatrix *X, *Y;            /* args to calls back to Matlab */
    int         i, mexerr;
    double      *dp;
    double      *Lval, *Uval;
    int         *Lrow, *Urow;
    int         *Lcol, *Ucol;
    int         nnzL, nnzU, snnzL, snnzU;

    /* Check number of arguments passed from Matlab. */
    if (nrhs != 2) {
	mexErrMsgTxt("SUPERLU requires 2 input arguments.");
    } else if (nlhs != 4) {
      	mexErrMsgTxt("SUPERLU requires 4 output arguments.");
    }   

    /* Read the Sparse Monitor Flag */
    X = mxCreateString("spumoni");
    mexerr = mexCallMATLAB(1, &Y, 1, &X, "sparsfun");
    SPUMONI = mxGetScalar(Y);
    mxFreeMatrix(Y);
    mxFreeMatrix(X);

    m = mxGetM(A_in);
    n = mxGetN(A_in);
    etree = (int *) mxCalloc(n, sizeof(int));
    perm_r = (int *) mxCalloc(m, sizeof(int));
    perm_c = mxGetIr(Pc_in); 
    val = mxGetPr(A_in);
    rowind = mxGetIr(A_in);
    colptr = mxGetJc(A_in);
    nnz = colptr[n];
    dCreate_CompCol_Matrix(&A, m, n, nnz, val, rowind, colptr, NC, _D, GE);
    *refact    = 'N';
    *trans     = 'N';
    panel_size = sp_ienv(1);
    relax      = sp_ienv(2);
    thresh     = 1.0;
    drop_tol   = 0.0;
    FlopsInSuperLU      = 0;

    StatInit(panel_size, relax);

    if ( verbose ) mexPrintf("Apply column perm to A and compute etree...\n");
    sp_preorder(refact, &A, perm_c, etree, &Ac);

    if ( verbose ) {
	mexPrintf("LU factorization...\n");
	mexPrintf("\tpanel_size %d, relax %d, diag_pivot_thresh %.2g\n",
		  panel_size, relax, thresh);
    }
    dgstrf(refact, &Ac, thresh, drop_tol, relax, panel_size, etree,
	   NULL, 0, perm_r, perm_c, &L, &U, &info);

    if ( verbose ) mexPrintf("INFO from dgstrf %d\n", info);

    /* Tell Matlab how many flops we did. */
    FlopsInSuperLU += LUFactFlops();
    if (verbose) mexPrintf("SUPERLU flops: %.f\n", FlopsInSuperLU);
    mexerr = mexCallMATLAB(1, &X, 0, NULL, "flops");
    *(mxGetPr(X)) += FlopsInSuperLU;
    mexerr = mexCallMATLAB(1, &Y, 1, &X, "flops");
    mxFreeMatrix(Y);
    mxFreeMatrix(X);
	
    /* Construct output arguments for Matlab. */
    if ( info >= 0 && info <= n ) {
	Pr_out = mxCreateFull(m, 1, REAL);
	dp = mxGetPr(Pr_out);
	for (i = 0; i < m; *dp++ = (double) perm_r[i++]+1);
	Pc_out = mxCreateFull(n, 1, REAL);
	dp = mxGetPr(Pc_out);
	for (i = 0; i < n; *dp++ = (double) perm_c[i++]+1);
	
	/* Now for L and U */
	nnzL = ((SCformat*)L.Store)->nnz; /* count diagonals */
   	nnzU = ((NCformat*)U.Store)->nnz;

	L_out = mxCreateSparse(m, n, nnzL, REAL);
	Lval = mxGetPr(L_out);
	Lrow = mxGetIr(L_out);
	Lcol = mxGetJc(L_out);

	U_out = mxCreateSparse(m, n, nnzU, REAL);
	Uval = mxGetPr(U_out);
	Urow = mxGetIr(U_out);
	Ucol = mxGetJc(U_out);

	LUextract(&L, &U, Lval, Lrow, Lcol, Uval, Urow, Ucol, &snnzL, &snnzU);
	
        Destroy_CompCol_Permuted(&Ac);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
        StatFree();

	if (babble) mexPrintf("factor nonzeros: %d unsqueezed, %d squeezed.\n",
			      nnzL + nnzU, snnzL + snnzU);
    } else {
	mexErrMsgTxt("Error returned from C dgstrf().");
    }

    mxFree(etree);
    mxFree(perm_r);
    StatFree();
    return;
}

void
LUextract(SuperMatrix *L, SuperMatrix *U, double *Lval, int *Lrow,
	  int *Lcol, double *Uval, int *Urow, int *Ucol, int *snnzL,
	  int *snnzU)
{
    int         i, j, k;
    int         upper;
    int         fsupc, istart, nsupr;
    int         lastl = 0, lastu = 0;
    SCformat    *Lstore;
    NCformat    *Ustore;
    double      *SNptr;

    Lstore = L->Store;
    Ustore = U->Store;
    Lcol[0] = 0;
    Ucol[0] = 0;
    
    /* for each supernode */
    for (k = 0; k <= Lstore->nsuper; ++k) {
	
	fsupc = L_FST_SUPC(k);
	istart = L_SUB_START(fsupc);
	nsupr = L_SUB_START(fsupc+1) - istart;
	upper = 1;
	
	/* for each column in the supernode */
	for (j = fsupc; j < L_FST_SUPC(k+1); ++j) {
	    SNptr = &((double*)Lstore->nzval)[L_NZ_START(j)];

	    /* Extract U */
	    for (i = U_NZ_START(j); i < U_NZ_START(j+1); ++i) {
		Uval[lastu] = ((double*)Ustore->nzval)[i];
 		/* Matlab doesn't like explicit zero. */
		if (Uval[lastu] != 0.0) Urow[lastu++] = U_SUB(i);
	    }
	    for (i = 0; i < upper; ++i) { /* upper triangle in the supernode */
		Uval[lastu] = SNptr[i];
 		/* Matlab doesn't like explicit zero. */
		if (Uval[lastu] != 0.0) Urow[lastu++] = L_SUB(istart+i);
	    }
	    Ucol[j+1] = lastu;

	    /* Extract L */
	    Lval[lastl] = 1.0; /* unit diagonal */
	    Lrow[lastl++] = L_SUB(istart + upper - 1);
	    for (i = upper; i < nsupr; ++i) {
		Lval[lastl] = SNptr[i];
 		/* Matlab doesn't like explicit zero. */
		if (Lval[lastl] != 0.0) Lrow[lastl++] = L_SUB(istart+i);
	    }
	    Lcol[j+1] = lastl;

	    ++upper;
	    
	} /* for j ... */
	
    } /* for k ... */

    *snnzL = lastl;
    *snnzU = lastu;
}
