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
#define b_in    	prhs[1]
#define Pc_in		prhs[2]
#define x_out		plhs[0]

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
    extern flops_t LUSolveFlops();
    
    /* Arguments to dgssv(). */
    SuperMatrix A;
    SuperMatrix B;
    SuperMatrix L, U;
    int	   	m, n, nnz;
    int         numrhs;
    double 	*vb, *x;
    double      *val;
    int       	*rowind;
    int		*colptr;
    int    	*perm_r, *perm_c;
    int		info;
    MatlabMatrix *X, *Y;            /* args to calls back to Matlab */
    int         i, mexerr;

    /* Check number of arguments passed from Matlab. */
    if (nrhs != 3) {
	mexErrMsgTxt("LUSOLVE requires 3 input arguments.");
    } else if (nlhs != 1) {
      	mexErrMsgTxt("LUSOLVE requires 1 output argument.");
    }   

    /* Read the Sparse Monitor Flag */
    X = mxCreateString("spumoni");
    mexerr = mexCallMATLAB(1, &Y, 1, &X, "sparsfun");
    SPUMONI = mxGetScalar(Y);
    mxFreeMatrix(Y);
    mxFreeMatrix(X);

    m = mxGetM(A_in);
    n = mxGetN(A_in);
    numrhs = mxGetN(b_in);
    if ( babble ) printf("m=%d, n=%d, numrhs=%d\n", m, n, numrhs);
    vb = mxGetPr(b_in);
    x_out = mxCreateFull(m, numrhs, REAL);
    x = mxGetPr(x_out);
    perm_r = (int *) mxCalloc(m, sizeof(int));
    perm_c = mxGetIr(Pc_in); 

    val = mxGetPr(A_in);
    rowind = mxGetIr(A_in);
    colptr = mxGetJc(A_in);
    nnz = colptr[n];
    dCreate_CompCol_Matrix(&A, m, n, nnz, val, rowind, colptr, NC, _D, GE);
    dCopy_Dense_Matrix(m, numrhs, vb, m, x, m);
    dCreate_Dense_Matrix(&B, m, numrhs, x, m, DN, _D, GE);

    FlopsInSuperLU = 0;
    
    /* Call simple driver */
    if ( verbose ) mexPrintf("Call LUSOLVE, use SUPERLU to factor first...\n");
    dgssv(&A, perm_c, perm_r, &L, &U, &B, &info);
	
    /* Tell Matlab how many flops we did. */
    FlopsInSuperLU += LUFactFlops() + LUSolveFlops();
    if ( verbose ) mexPrintf("LUSOLVE flops: %.f\n", FlopsInSuperLU);
    mexerr = mexCallMATLAB(1, &X, 0, NULL, "flops");
    *(mxGetPr(X)) += FlopsInSuperLU;
    mexerr = mexCallMATLAB(1, &Y, 1, &X, "flops");
    mxFreeMatrix(Y);
    mxFreeMatrix(X);

    /* Construct Matlab solution matrix. */
    if ( !info ) {
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);
        if ( babble ) printf("Destroy L & U from SuperLU...\n");
    } else {
	mexErrMsgTxt("Error returned from C dgssv().");
    }

    mxFree(perm_r);

    return;
 
}

