/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

/*
  Copyright (c) 1994 by Xerox Corporation.  All rights reserved.
 
  THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
  EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 
  Permission is hereby granted to use or copy this program for any
  purpose, provided the above notices are retained on all copies.
  Permission to modify the code and to distribute modified code is
  granted, provided the above notices are retained, and a notice that
  the code was modified is included with the above copyright notice.
*/

#include <math.h>
#include "superlu_zdefs.h"


void
zCreate_CompCol_Matrix(SuperMatrix *A, int_t m, int_t n, int_t nnz, 
			doublecomplex *nzval, int_t *rowind, int_t *colptr,
		      	Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NCformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NCformat) );
    if ( !(A->Store) ) ABORT("SUPERLU_MALLOC fails for A->Store");
    Astore = (NCformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colptr = colptr;
}

/* Copy matrix A into matrix B. */
void
zCopy_CompCol_Matrix(SuperMatrix *A, SuperMatrix *B)
{
    NCformat *Astore, *Bstore;
    int      ncol, nnz, i;

    B->Stype = A->Stype;
    B->Dtype = A->Dtype;
    B->Mtype = A->Mtype;
    B->nrow  = A->nrow;;
    B->ncol  = ncol = A->ncol;
    Astore   = (NCformat *) A->Store;
    Bstore   = (NCformat *) B->Store;
    Bstore->nnz = nnz = Astore->nnz;
    for (i = 0; i < nnz; ++i)
	((doublecomplex *)Bstore->nzval)[i] = ((doublecomplex *)Astore->nzval)[i];
    for (i = 0; i < nnz; ++i) Bstore->rowind[i] = Astore->rowind[i];
    for (i = 0; i <= ncol; ++i) Bstore->colptr[i] = Astore->colptr[i];
}


void zPrint_CompCol_Matrix(SuperMatrix *A)
{
    NCformat     *Astore;
    register int i;
    doublecomplex       *dp;
    
    printf("\nCompCol matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (NCformat *) A->Store;
    dp = (doublecomplex *) Astore->nzval;
    printf("nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
    printf("\nnzval: ");
    for (i = 0; i < Astore->nnz; ++i) printf("%f  ", dp[i]);
    printf("\nrowind: ");
    for (i = 0; i < Astore->nnz; ++i) printf("%d  ", Astore->rowind[i]);
    printf("\ncolptr: ");
    for (i = 0; i <= A->ncol; ++i) printf("%d  ", Astore->colptr[i]);
    printf("\nend CompCol matrix.\n");
}

void zPrint_Dense_Matrix(SuperMatrix *A)
{
    DNformat     *Astore;
    register int i;
    doublecomplex       *dp;
    
    printf("\nDense matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (DNformat *) A->Store;
    dp = (doublecomplex *) Astore->nzval;
    printf("nrow %d, ncol %d, lda %d\n", A->nrow,A->ncol,Astore->lda);
    printf("\nnzval: ");
    for (i = 0; i < A->nrow; ++i) printf("%f  ", dp[i]);
    printf("\nend Dense matrix.\n");
}

void
zCreate_Dense_Matrix(SuperMatrix *X, int_t m, int_t n, doublecomplex *x,
		     int_t ldx, Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    DNformat    *Xstore;
    
    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    X->Store = (void *) SUPERLU_MALLOC( sizeof(DNformat) );
    if ( !(X->Store) ) ABORT("SUPERLU_MALLOC fails for X->Store");
    Xstore = (DNformat *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = (doublecomplex *) x;
}

void
zCopy_Dense_Matrix(int_t M, int_t N, doublecomplex *X, int_t ldx,
		   doublecomplex *Y, int_t ldy)
{
/*
 *
 *  Purpose
 *  =======
 *
 *  Copies a two-dimensional matrix X to another matrix Y.
 */
    int    i, j;
    
    for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i)
            Y[i + j*ldy] = X[i + j*ldx];
}

void
zCreate_SuperNode_Matrix(SuperMatrix *L, int_t m, int_t n, int_t nnz, 
			 doublecomplex *nzval, int_t *nzval_colptr,
			 int_t *rowind, int_t *rowind_colptr,
			 int_t *col_to_sup, int_t *sup_to_col,
			 Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    SCformat *Lstore;

    L->Stype = stype;
    L->Dtype = dtype;
    L->Mtype = mtype;
    L->nrow = m;
    L->ncol = n;
    L->Store = (void *) SUPERLU_MALLOC( sizeof(SCformat) );
    if ( !(L->Store) ) ABORT("SUPERLU_MALLOC fails for L->Store");
    Lstore = L->Store;
    Lstore->nnz = nnz;
    Lstore->nsuper = col_to_sup[n];
    Lstore->nzval = nzval;
    Lstore->nzval_colptr = nzval_colptr;
    Lstore->rowind = rowind;
    Lstore->rowind_colptr = rowind_colptr;
    Lstore->col_to_sup = col_to_sup;
    Lstore->sup_to_col = sup_to_col;

}


/*
 * Check whether tempv[] == 0. This should be true before and after 
 * calling any numeric routines, i.e., "panel_bmod" and "column_bmod". 
 */
void zcheck_tempv(int n, doublecomplex *tempv)
{
    int i;
	
    for (i = 0; i < n; i++) {
	if ((tempv[i].r != 0.0) || (tempv[i].i != 0.0))
	{
	    fprintf(stderr,"tempv[%d] = %f\n", i,tempv[i]);
	    ABORT("zcheck_tempv");
	}
    }
}


void
zGenXtrue(int_t n, int_t nrhs, doublecomplex *x, int_t ldx)
{
    int  i, j;
    for (j = 0; j < nrhs; ++j)
	for (i = 0; i < n; ++i) {
	    if ( i % 2 ) x[i + j*ldx].r = 1.0;
	    else x[i + j*ldx].r = 2.0;
	    x[i + j*ldx].i = 0.0;
	}
}

/*
 * Let rhs[i] = sum of i-th row of A, so the solution vector is all 1's
 */
void
zFillRHS(char *trans, int_t nrhs, doublecomplex *x, int_t ldx,
	 SuperMatrix *A, doublecomplex *rhs, int_t ldb)
{
    doublecomplex one = {1.0, 0.0};
    doublecomplex zero = {0.0, 0.0};

    sp_zgemm(trans, "N", A->nrow, nrhs, A->ncol, one, A,
	     x, ldx, zero, rhs, ldb);

}

/* 
 * Fills a doublecomplex precision array with a given value.
 */
void 
zfill(doublecomplex *a, int_t alen, doublecomplex dval)
{
    register int_t i;
    for (i = 0; i < alen; i++) a[i] = dval;
}



/* 
 * Check the inf-norm of the error vector 
 */
void zinf_norm_error(int_t n, int_t nrhs, doublecomplex *x, int_t ldx,
		     doublecomplex *xtrue, int_t ldxtrue)
{
    double err, xnorm;
    doublecomplex *x_work, *xtrue_work;
    doublecomplex temp;
    int i, j;

    for (j = 0; j < nrhs; j++) {
      x_work = &x[j*ldx];
      xtrue_work = &xtrue[j*ldxtrue];
      err = xnorm = 0.0;
      for (i = 0; i < n; i++) {
        z_sub(&temp, &x_work[i], &xtrue_work[i]);
	err = MAX(err, z_abs(&temp));
	xnorm = MAX(xnorm, z_abs(&x_work[i]));
      }
      err = err / xnorm;
      printf(".. ||X - Xtrue||/||X|| = %e\n", err);
    }
}



/* Print performance of the code. */
void
zPrintPerf(SuperMatrix *L, SuperMatrix *U, mem_usage_t *mem_usage,
	   doublecomplex rpg, doublecomplex rcond, doublecomplex *ferr,
	   doublecomplex *berr, char *equed)
{
    SCformat *Lstore;
    NCformat *Ustore;
    extern SuperLUStat_t SuperLUStat;
    double   *utime;
    flops_t  *ops;
    
    utime = SuperLUStat.utime;
    ops   = SuperLUStat.ops;
    
    if ( utime[FACT] != 0. )
	printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT],
	       ops[FACT]*1e-6/utime[FACT]);
    printf("Identify relaxed snodes	= %8.2f\n", utime[RELAX]);
    if ( utime[SOLVE] != 0. )
	printf("Solve flops = %.0f, Mflops = %8.2f\n", ops[SOLVE],
	       ops[SOLVE]*1e-6/utime[SOLVE]);
    
    Lstore = (SCformat *) L->Store;
    Ustore = (NCformat *) U->Store;
    printf("\tNo of nonzeros in factor L = %d\n", Lstore->nnz);
    printf("\tNo of nonzeros in factor U = %d\n", Ustore->nnz);
    printf("\tNo of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
	
    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	   mem_usage->for_lu/1e6, mem_usage->total/1e6,
	   mem_usage->expansions);
	
    printf("\tFactor\tMflops\tSolve\tMflops\tEtree\tEquil\tRcond\tRefine\n");
    printf("PERF:%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
	   utime[FACT], ops[FACT]*1e-6/utime[FACT],
	   utime[SOLVE], ops[SOLVE]*1e-6/utime[SOLVE],
	   utime[ETREE], utime[EQUIL], utime[RCOND], utime[REFINE]);
    
    printf("\tRpg\t\tRcond\t\tFerr\t\tBerr\t\tEquil?\n");
    printf("NUM:\t%e\t%e\t%e\t%e\t%s\n",
	   rpg, rcond, ferr[0], berr[0], equed);
    
}


void PrintDoublecomplex(char *name, int_t len, doublecomplex *x)
{
    register int_t i;
    
    printf("%10s:\tReal\tImag\n", name);
    for (i = 0; i < len; ++i)
	printf("\t%d\t%.4f\t%.4f\n", i, x[i].r, x[i].i);
}

/* 
 * Print the blocks in the factored matrix L.
 */
void zPrintLblocks(int_t iam, int_t nsupers, gridinfo_t *grid,
		   Glu_persist_t *Glu_persist, LocalLU_t *Llu)
{
    register int_t c, extra, gb, j, lb, nsupc, nsupr, len, nb, ncb;
    register int_t k, mycol, r;
    int_t *xsup = Glu_persist->xsup;
    int_t *index;
    doublecomplex *nzval;

    printf("\n(%d) L BLOCKS IN COLUMN-MAJOR ORDER -->\n", iam);
    ncb = nsupers / grid->npcol;
    extra = nsupers % grid->npcol;
    mycol = MYCOL( iam, grid );
    if ( mycol < extra ) ++ncb;
    for (lb = 0; lb < ncb; ++lb) {
	index = Llu->Lrowind_bc_ptr[lb];
	if ( index ) { /* Not an empty column */
	    nzval = Llu->Lnzval_bc_ptr[lb];
	    nb = index[0];
	    nsupr = index[1];
	    gb = lb * grid->npcol + mycol;
	    nsupc = SuperSize( gb );
	    printf("(%d) block column (local) %d, # row blocks %d\n",
		   iam, lb, nb);
	    for (c = 0, k = BC_HEADER, r = 0; c < nb; ++c) {
		len = index[k+1];
		printf("(%d) row-block %d: block # %d\tlength %d\n", 
		       iam, c, index[k], len);
		PrintInt10("lsub", len, &index[k+LB_DESCRIPTOR]);
		for (j = 0; j < nsupc; ++j) {
		    PrintDoublecomplex("nzval", len, &nzval[r + j*nsupr]);
		}
		k += LB_DESCRIPTOR + len;
		r += len;
	    }
	}
	printf("(%d)", iam);
 	PrintInt10("ToSendR[]", grid->npcol, Llu->ToSendR[lb]);
	PrintInt10("fsendx_plist[]", grid->nprow, Llu->fsendx_plist[lb]);
    }
    printf("nfrecvx %4d\n", Llu->nfrecvx);
    k = CEILING( nsupers, grid->nprow );
    PrintInt10("fmod", k, Llu->fmod);
    
} /* PRINTLBLOCKS */


/* 
 * Print the blocks in the factored matrix U.
 */
void zPrintUblocks(int_t iam, int_t nsupers, gridinfo_t *grid, 
		   Glu_persist_t *Glu_persist, LocalLU_t *Llu)
{
    register int_t c, extra, jb, k, lb, len, nb, nrb, nsupc;
    register int_t myrow, r;
    int_t *xsup = Glu_persist->xsup;
    int_t *index;
    doublecomplex *nzval;

    printf("\n(%d) U BLOCKS IN ROW-MAJOR ORDER -->\n", iam);
    nrb = nsupers / grid->nprow;
    extra = nsupers % grid->nprow;
    myrow = MYROW( iam, grid );
    if ( myrow < extra ) ++nrb;
    for (lb = 0; lb < nrb; ++lb) {
	index = Llu->Ufstnz_br_ptr[lb];
	if ( index ) { /* Not an empty row */
	    nzval = Llu->Unzval_br_ptr[lb];
	    nb = index[0];
	    printf("(%d) block row (local) %d, # column blocks %d\n",
		   iam, lb, nb);
	    r  = 0;
	    for (c = 0, k = BR_HEADER; c < nb; ++c) {
		jb = index[k];
		len = index[k+1];
		printf("(%d) col-block %d: block # %d\tlength %d\n", 
		       iam, c, jb, index[k+1]);
		nsupc = SuperSize( jb );
		PrintInt10("fstnz", nsupc, &index[k+UB_DESCRIPTOR]);
		PrintDoublecomplex("nzval", len, &nzval[r]);
		k += UB_DESCRIPTOR + nsupc;
		r += len;
	    }

	    printf("(%d) ToSendD[] %d\n", iam, Llu->ToSendD[lb]);
	}
    }
} /* PRINTUBLOCKS */

/* cg5.cua
            b = A*x           y = L\b
   0        1 + 4.0000i       1.0000 + 4.0000i
   1        0 + 5.0000i	      1.3529 + 5.4118i
   2        1 + 4.0000i	      1.0000 + 4.0000i
   3        2 + 3.0000i	      2.0000 + 3.0000i
   4        1 + 4.0000i	      3.5882 + 4.3529i
   5        1 + 4.0000i	      4.1250 + 3.3202i
   6          + 5.0000i	      4.4640 + 3.8632i
   7        2 + 3.0000i	      2.0000 + 3.0000i
   8        2 + 3.0000i	      2.0000 + 3.0000i
   9        1 + 4.0000i	      1.0000 + 4.0000i
  10        1 + 4.0000i	      3.5882 + 4.3529i
  11          + 5.0000i	           0 + 5.0000i
  12        1 + 4.0000i	      5.1793 + 4.6604i
  13        2 + 3.0000i	      2.0000 + 3.0000i
  14        1 + 4.0000i	      1.0000 + 4.0000i
  15          + 5.0000i	      1.3529 + 5.4118i
  16        1 + 4.0000i	      4.0045 + 3.8950i
  17          + 5.0000i	      3.0338 + 4.6248i
  18        1 + 4.0000i	      5.4495 + 2.2703i
  19          + 5.0000i	      4.0980 + 3.7290i
  20          + 5.0000i	      4.2680 + 3.7739i
  21          + 5.0000i	      5.3514 + 2.9480i
  22        1 + 4.0000i	      4.4178 + 2.0476i
  23        1 + 4.0000i	      3.5615 + 2.8322i
  24          + 5.0000i	      4.7526 + 2.2605i
*/
