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
#include "superlu_ddefs.h"

void
dCreate_CompCol_Matrix(SuperMatrix *A, int_t m, int_t n, int_t nnz, 
			double *nzval, int_t *rowind, int_t *colptr,
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
dCopy_CompCol_Matrix(SuperMatrix *A, SuperMatrix *B)
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
	((double *)Bstore->nzval)[i] = ((double *)Astore->nzval)[i];
    for (i = 0; i < nnz; ++i) Bstore->rowind[i] = Astore->rowind[i];
    for (i = 0; i <= ncol; ++i) Bstore->colptr[i] = Astore->colptr[i];
}


void dPrint_CompCol_Matrix(SuperMatrix *A)
{
    NCformat     *Astore;
    register int i;
    double       *dp;
    
    printf("\nCompCol matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (NCformat *) A->Store;
    dp = (double *) Astore->nzval;
    printf("nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
    printf("\nnzval: ");
    for (i = 0; i < Astore->nnz; ++i) printf("%f  ", dp[i]);
    printf("\nrowind: ");
    for (i = 0; i < Astore->nnz; ++i) printf("%d  ", Astore->rowind[i]);
    printf("\ncolptr: ");
    for (i = 0; i <= A->ncol; ++i) printf("%d  ", Astore->colptr[i]);
    printf("\nend CompCol matrix.\n");
}

void dPrint_Dense_Matrix(SuperMatrix *A)
{
    DNformat     *Astore;
    register int i;
    double       *dp;
    
    printf("\nDense matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (DNformat *) A->Store;
    dp = (double *) Astore->nzval;
    printf("nrow %d, ncol %d, lda %d\n", A->nrow,A->ncol,Astore->lda);
    printf("\nnzval: ");
    for (i = 0; i < A->nrow; ++i) printf("%f  ", dp[i]);
    printf("\nend Dense matrix.\n");
}

void
dCreate_Dense_Matrix(SuperMatrix *X, int_t m, int_t n, double *x, int_t ldx,
		     Stype_t stype, Dtype_t dtype, Mtype_t mtype)
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
    Xstore->nzval = (double *) x;
}

void
dCopy_Dense_Matrix(int_t M, int_t N, double *X, int_t ldx,
		   double *Y, int_t ldy)
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
dCreate_SuperNode_Matrix(SuperMatrix *L, int_t m, int_t n, int_t nnz, 
			 double *nzval, int_t *nzval_colptr,
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
void dcheck_tempv(int n, double *tempv)
{
    int i;
	
    for (i = 0; i < n; i++) {
	if (tempv[i] != 0.0) 
	{
	    fprintf(stderr,"tempv[%d] = %f\n", i,tempv[i]);
	    ABORT("dcheck_tempv");
	}
    }
}


void
dGenXtrue(int_t n, int_t nrhs, double *x, int_t ldx)
{
    int  i, j;
    for (j = 0; j < nrhs; ++j)
	for (i = 0; i < n; ++i) {
	    if ( i % 2 ) x[i + j*ldx] = 1.0;/* + (double)(i+1.)/n;*/
	    else x[i + j*ldx] = 1.0;
	}
}

/*
 * Let rhs[i] = sum of i-th row of A, so the solution vector is all 1's
 */
void
dFillRHS(char *trans, int_t nrhs, double *x, int_t ldx,
	 SuperMatrix *A, double *rhs, int_t ldb)
{
    double one = 1.0;
    double zero = 0.0;

    sp_dgemm(trans, "N", A->nrow, nrhs, A->ncol, one, A,
	     x, ldx, zero, rhs, ldb);

}

/* 
 * Fills a double precision array with a given value.
 */
void 
dfill(double *a, int_t alen, double dval)
{
    register int_t i;
    for (i = 0; i < alen; i++) a[i] = dval;
}



/* 
 * Check the inf-norm of the error vector 
 */
void dinf_norm_error(int_t n, int_t nrhs, double *x, int_t ldx,
		     double *xtrue, int_t ldxtrue)
{
    double err, xnorm;
    double *x_work, *xtrue_work;
    int i, j;

    for (j = 0; j < nrhs; j++) {
      x_work = &x[j*ldx];
      xtrue_work = &xtrue[j*ldxtrue];
      err = xnorm = 0.0;
      for (i = 0; i < n; i++) {
	err = MAX(err, fabs(x_work[i] - xtrue_work[i]));
	xnorm = MAX(xnorm, fabs(x_work[i]));
      }
      err = err / xnorm;
      printf(".. ||X - Xtrue||/||X|| = %e\n", err);
    }
}



/* Print performance of the code. */
void
dPrintPerf(SuperMatrix *L, SuperMatrix *U, mem_usage_t *mem_usage,
	   double rpg, double rcond, double *ferr,
	   double *berr, char *equed)
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


void PrintDouble5(char *name, int_t len, double *x)
{
    register int_t i;
    
    printf("%10s:", name);
    for (i = 0; i < len; ++i) {
	if ( i % 5 == 0 ) printf("\n[%2d-%2d] ", i, i+4);
	printf("%14e", x[i]);
    }
    printf("\n");
}

/* 
 * Print the blocks in the factored matrix L.
 */
void PrintLblocks(int_t iam, int_t nsupers, gridinfo_t *grid,
		  Glu_persist_t *Glu_persist, LocalLU_t *Llu)
{
    register int_t c, extra, gb, j, lb, nsupc, nsupr, len, nb, ncb;
    register int_t k, mycol, r;
    int_t *xsup = Glu_persist->xsup;
    int_t *index;
    double *nzval;

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
		    PrintDouble5("nzval", len, &nzval[r + j*nsupr]);
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
void PrintUblocks(int_t iam, int_t nsupers, gridinfo_t *grid, 
		  Glu_persist_t *Glu_persist, LocalLU_t *Llu)
{
    register int_t c, extra, jb, k, lb, len, nb, nrb, nsupc;
    register int_t myrow, r;
    int_t *xsup = Glu_persist->xsup;
    int_t *index;
    double *nzval;

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
		PrintDouble5("nzval", len, &nzval[r]);
		k += UB_DESCRIPTOR + nsupc;
		r += len;
	    }

	    printf("(%d) ToSendD[] %d\n", iam, Llu->ToSendD[lb]);
	}
    }
} /* PRINTUBLOCKS */


void
GenXtrueRHS(int nrhs, SuperMatrix *A, Glu_persist_t *Glu_persist,
	    gridinfo_t *grid, double **xact, int *ldx, double **b, int *ldb)
{
    int_t gb, gbrow, i, iam, irow, j, lb, lsup, myrow, n, nlrows,
          nsupr, nsupers, rel;
    int_t *supno, *xsup, *lxsup;
    double *x, *bb;
    NCformat *Astore;
    double   *Aval;

    n = A->ncol;
    *ldb = 0;
    supno = Glu_persist->supno;
    xsup = Glu_persist->xsup;
    nsupers = supno[n-1] + 1;
    iam = grid->iam;
    myrow = MYROW( iam, grid );
    Astore = A->Store;
    Aval = Astore->nzval;
    lb = CEILING( nsupers, grid->nprow ) + 1;
    if ( !(lxsup = intMalloc(lb)) )
	ABORT("Malloc fails for lxsup[].");

    lsup = 0;
    nlrows = 0;
    for (j = 0; j < nsupers; ++j) {
	i = PROW( j, grid );
	if ( myrow == i ) {
	    nsupr = SuperSize( j );
	    *ldb += nsupr;
	    lxsup[lsup++] = nlrows;
	    nlrows += nsupr;
	}
    }
    *ldx = n;
    if ( !(x = doubleMalloc(*ldx * nrhs)) )
	ABORT("Malloc fails for x[].");
    if ( !(bb = doubleCalloc(*ldb * nrhs)) )
	ABORT("Calloc fails for bb[].");
    for (j = 0; j < nrhs; ++j)
	for (i = 0; i < n; ++i) x[i + j*(*ldx)] = 1.0;

    /* Form b = A*x. */
    for (j = 0; j < n; ++j)
	for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; ++i) {
	    irow = Astore->rowind[i];
	    gb = supno[irow];
	    gbrow = PROW( gb, grid );
	    if ( myrow == gbrow ) {
		rel = irow - xsup[gb];
		lb = LBi( gb, grid );
		bb[lxsup[lb] + rel] += Aval[i] * x[j];
	    }
	}

    /* Memory allocated but not freed: xact, b */
    *xact = x;
    *b = bb;

    SUPERLU_FREE(lxsup);

#if ( PRNTlevel>=2 )
    for (i = 0; i < grid->nprow*grid->npcol; ++i) {
	if ( iam == i ) {
	    printf("\n(%d)\n", iam);
	    PrintDouble5("rhs", *ldb, *b);
	}
	MPI_Barrier( grid->comm );
    }
#endif

} /* GENXTRUERHS */

/* g5.rua
          b = A*x    y = L\b
   0      1	     1.0000
   1      0	     0.2500
   2      1	     1.0000
   3      2	     2.0000
   4      1	     1.7500
   5      1	     1.8917
   6      0	     1.1879
   7      2	     2.0000
   8      2	     2.0000
   9      1	     1.0000
   10     1	     1.7500
   11     0	          0
   12     1	     1.8750
   13     2	     2.0000
   14     1	     1.0000
   15     0	     0.2500
   16     1	     1.7667
   17     0	     0.6419
   18     1	     2.2504
   19     0	     1.1563
   20     0	     0.9069
   21     0	     1.4269
   22     1	     2.7510
   23     1	     2.2289
   24     0	     2.4332

   g6.rua
       b=A*x  y=L\b
    0    0         0
    1    1    1.0000
    2    1    1.0000
    3    2    2.5000
    4    0         0
    5    2    2.0000
    6    1    1.0000
    7    1    1.7500
    8    1    1.0000
    9    0    0.2500
   10    0    0.5667
   11    1    2.0787
   12    0    0.8011
   13    1    1.9838
   14    1    1.0000
   15    1    1.0000
   16    2    2.5000
   17    0    0.8571
   18    0         0
   19    1    1.0000
   20    0    0.2500
   21    1    1.0000
   22    2    2.0000
   23    1    1.7500
   24    1    1.8917
   25    0    1.1879
   26    0    0.8011
   27    1    1.9861
   28    1    2.0199
   29    0    1.3620
   30    0    0.6136
   31    1    2.3677
   32    0    1.1011
   33    0    1.5258
   34    0    1.7628
   35    0    2.1658
*/

