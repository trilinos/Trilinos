#include <math.h>
#include "zsuperlu_ddefs.h"


void pzgssv(SuperMatrix *A, doublecomplex *B, int ldb, int nrhs,
	    gridinfo_t *grid, GlobalLU_t *Glu, LocalLU_t *Llu, int *info)
{
/* 
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 *
 * Purpose
 * =======
 *
 * PZGSSV solves a system of distributed linear equations A*X=B,
 * using the LU factorization from PDGSTRF.
 * It performs the following steps:
 *
 *   1. If fact = 'E', scaling factors are computed to equilibrate the system:
 *        trans = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
 *        trans = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
 *        trans = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
 *      Whether or not the system will be equilibrated depends on the
 *      scaling of the matrix A, but if equilibration is used, A is
 *      overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if trans='N')
 *      or diag(C)*B (if trans = 'T' or 'C').
 *
 *   2. Permute rows and columns of A, forming Pc*A*Pc^T, where Pc is a
 *      permutation matrix that usually preserves sparisity.
 *      For more details of this step, see sp_colorder.c.
 *
 *   3, Permute rows of A, forming Pr*A, where Pr is a permutation
 *      matrix that tries to put large entries on the diagonal.
 *
 *   4. The LU decomposition is used to factor the matrix A (after
 *      equilibration if fact = 'E') as Pr*Pc*A*Pc^T = L*U,
 *      with Pc and Pr determined in Steps 2 and 3, respectively.
 *
 *   5. The system of equations is solved for X using the factored form of A.
 *      The rows of B (after scaling if equilibration was used) are
 *      permuted by Pr*Pc.
 *
 *   6. Iterative refinement is applied to improve the computed solution
 *      matrix and calculate error bounds and backward error estimates
 *      for it.
 *
 *   7. If equilibration was used, the matrix X is premultiplied by
 *      diag(C) (if trans = 'N') or diag(R) (if trans = 'T' or 'C') so
 *      that it solves the original system before equilibration.
 *
 *
 * Arguments
 * =========
 *
 * A      (input) SuperMatrix*
 *        Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *        of linear equations is A->nrow. The type of A can be:
 *        Stype = NC; Dtype = _Z; Mtype = GE.
 *
 *        NOTE: Currently, A must reside on all processes when calling
 *              this routine.
 *
 * B      (input/output) doublecomplex*
 *        On entry, the right-hand side matrix of dimension (A->nrow, nrhs).
 *        On exit, the solution matrix if info = 0;
 *
 *        NOTE: Currently, B must reside on all processes when calling
 *              this routine.
 *
 * ldb    (input) int (global)
 *        The leading dimension of matrix B.
 *
 * nrhs   (input) int (global)
 *        The number of right-hand sides.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Glu    (output) GlobalLU_t*
 *        Global data structure (Glu->xsup, Glu->supno) replicated on 
 *        all processes, describing the supernode partition in the
 *        factored matrices L and U.
 *	      xsup[s] is the leading column of the s-th supernode.
 *            supno[i] is the supernode number to which column i belongs;
 *        See zsp_defs.h for the definition of 'GlobalLU_t' structure.
 *
 * Llu    (output) LocalLU_t*
 *        Local data structures to store distributed L and U matrices.
 *        See zsuperlu_ddefs.h for the definition of 'LocalLU_t' structure.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        > 0: if info = i, and i is
 *            <= A->ncol: U(i,i) is exactly zero. The factorization has
 *               been completed, but the factor U is exactly singular,
 *               so the solution could not be computed.
 *            > A->ncol: number of bytes allocated when memory allocation
 *               failure occurred, plus A->ncol.
 *
 */
    SuperMatrix AC;
    NCformat *Astore;
    NCPformat *ACstore;
    doublecomplex *a;
    int_t    *perm_c; /* column permutation vector */
    int_t    *perm_r; /* row permutations from partial pivoting */
    int_t    *etree;  /* column elimination tree */
    int_t    *colcnt, *colptr, *rowind; /* Adjacency structure for A_bar. */
    int_t    colequ, equil, notran, rowequ, static_pivot;
    int_t    i, ii, iinfo, j, irow, m, n, nnz, permc_spec;
    int      iam;
    int      l_ldb, l_ldx;  /* LDA for matrices B and X (Local). */
    char     fact[1], equed[1], norm[1], refact[1], trans[1];
    char     **cpp, c;
    void     *work;
    int_t    lwork;
    double   *C, *R, amax, anorm, colcnd, rowcnd;
    double   droptol, eta1, eta2;
    doublecomplex *b, *x, *b_col, *b_work, *x_col;
    double   t;
    SuperLUStat_t stat;
    mem_usage_t num_mem_usage, symb_mem_usage;

    *trans = 'N';
    notran = lsame_(trans, "N");

    /* Test input parameters. */
    *info = 0;
    if ( A->nrow != A->ncol || A->nrow < 0 ||
         A->Stype != NC || A->Dtype != _Z || A->Mtype != GE )
	*info = -1;
    else if ( ldb < A->nrow ) 
	*info = -3;
    else if ( nrhs <= 0 )
	*info = -4;
    if ( *info ) {
	i = -(*info);
	pxerbla("PZGSSV", grid, -*info);
	return;
    }

    /* Initialization */
    iam = grid->iam;
    *fact = 'E';
    *refact = 'N';
    *trans = 'N';
    work = NULL;
    lwork = 0;
    eta1 = 1e-26;
    Astore = A->Store;
    m = A->nrow;
    n = A->ncol;
    nnz = Astore->nnz;
    a = Astore->nzval;

    /* Initialize the statistics variables. */
    PStatInit(&stat);

    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");

    /* Equilibration. */
    equil = (*fact == 'E' || *fact == 'e');
    if ( equil ) {
#if ( DEBUGlevel>=1 )
	CHECK_MALLOC(iam, "Enter equil");
#endif
	t = SuperLU_timer_();

	if ( !(C = (double *) SUPERLU_MALLOC(A->ncol*sizeof(double))) )
	    ABORT("Malloc fails for C[].");
	if ( !(R = (double *) SUPERLU_MALLOC(A->nrow*sizeof(double))) )
	    ABORT("Malloc fails for R[].");

	/* Compute row and column scalings to equilibrate matrix A. */
	zgsequ(A, R, C, &rowcnd, &colcnd, &amax, &iinfo);
	
	if ( iinfo == 0 ) {
	    /* Equilibrate matrix A. */
	    zlaqgs(A, R, C, rowcnd, colcnd, amax, equed);
	    rowequ = lsame_(equed, "R") || lsame_(equed, "B");
	    colequ = lsame_(equed, "C") || lsame_(equed, "B");
	}
#if ( PRNTlevel>=1 )
	if ( !iam ) {
	    printf(".. equilibrated? *equed = %c\n", *equed);
	    /*fflush(stdout);*/
	}
#endif
	*fact = 'N';
	stat.utime[EQUIL] = SuperLU_timer_() - t;
#if ( DEBUGlevel>=1 )
	CHECK_MALLOC(iam, "Exit equil");
#endif
    }
    
    /* Compute norm(A), which will be used to drop small entries. */
    if ( notran ) *(unsigned char *)norm = '1';
    else *(unsigned char *)norm = 'I';
    anorm = zlangs(norm, A);

    t = SuperLU_timer_();
    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering 
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     */    	
    permc_spec = 2;
    get_perm_c(iam, permc_spec, A, perm_c);

    /* Permute columns of A to form A*Pc', where Pc is a permutation matrix. */
    sp_colorder(refact, A, perm_c, etree, &AC);

    /* Form Pc*A*Pc' to preserve the diagonal of the original matrix. */
    ACstore = AC.Store;
    for (j = 0; j < n; ++j) 
	for (i = ACstore->colbeg[j]; i < ACstore->colend[j]; ++i) {
	    irow = ACstore->rowind[i];
	    ACstore->rowind[i] = perm_c[irow];
	}
    stat.utime[COLPERM] = SuperLU_timer_() - t;


    /* Check whether diagonal is zero-free. */
    if ( !(colcnt = intMalloc(n)) ) ABORT("Malloc fails for colcnt[].");
    for (j = 0; j < n; ++j) 
	colcnt[j] = ACstore->colend[j] - ACstore->colbeg[j];
    static_pivot = CheckZeroDiagonal(n, ACstore->rowind,
				     ACstore->colbeg, colcnt);
#if ( PRNTlevel>=1 )
    if ( !iam )
	printf(".. static pivoting? %d\n", static_pivot);
#endif
    if ( static_pivot ) { /* Do this only when diagonal has zeros. */

	droptol = eta1 * anorm;
#if ( PRNTlevel>=1 )
	if ( !iam ) {
	    printf(".. # zeros on diagonal = %d\n", static_pivot);
	    printf(".. anorm = %e, droptol = %4.2e * anorm = %e\n", 
		   anorm, eta1, droptol);
	}
#endif
        /* ------------------------------------------------------------
           Permute rows of A so that it does not have zeros on its
           diagonal.
           ------------------------------------------------------------*/
        t = SuperLU_timer_();

        /* Remove the tiny elements before performing maximum matching. */
        if ( !(colptr = intMalloc(n+1)) ) ABORT("Malloc fails for colptr[].");
        if ( !(rowind = intMalloc(nnz)) ) ABORT("Malloc fails for rowind[].");;
        for (colptr[0] = 0, ii = 0, j = 0; j < n; ++j) {
            colcnt[j] = 0;
    	    for (i = ACstore->colbeg[j]; i < ACstore->colend[j]; ++i) {
                irow = ACstore->rowind[i];
	        if ( z_abs( &a[i] ) > droptol ) {
		    rowind[ii++] = irow;
		    ++colcnt[j];
		}
	    }
	    colptr[j+1] = ii;
	}
#if ( PRNTlevel>=1 )
	if ( !iam ) printf(".. Removed small entries -- nnz(A) = %d\n", ii);
#endif
	zfdperm(m, nnz, ACstore->rowind, ACstore->colbeg, colcnt, perm_r);

	/* Apply perm_r[] to matrix AC. */
	for (j = 0; j < n; ++j) {
	    for (i = ACstore->colbeg[j]; i < ACstore->colend[j]; ++i) {
		irow = ACstore->rowind[i];
		ACstore->rowind[i] = perm_r[irow];
	    }
	}
	
	SUPERLU_FREE(colptr);
	SUPERLU_FREE(rowind);
	stat.utime[ZFDPERM] = SuperLU_timer_() - t;

    } /* if static_pivot ... */

#if ( DEBUGlevel>=1 )
    /* Check whether diagonal is zero-free. */
    CheckZeroDiagonal(n, ACstore->rowind, ACstore->colbeg, colcnt);
#endif
#if ( DEBUGlevel>=2 )
    {
	/* Write the permuted matrix to a file. */
	doublecomplex *nzval = ACstore->nzval;
	FILE *fp, *fopen();

	fp = fopen("AC.triple", "w");
	for (j = 0; j < n; ++j) 
	    for (i = ACstore->colbeg[j]; i < ACstore->colend[j]; ++i) {
		irow = ACstore->rowind[i];
		fprintf(fp, "%8d%8d%20.8e\n", irow+1, j+1, nzval[i]);
	    }
	fclose(fp);
    }
#endif

    SUPERLU_FREE(colcnt);


    /* Allocate storage common to the symbolic factor routines */
    iinfo = LUSubInit(refact, work, lwork, m, n, nnz, Glu);

    /* Perform a symbolic factorization on matrix A and sets up the nonzero
       data structures which are suitable for supernodal GENP. */
#if ( PRNTlevel>=1 ) 
    if ( !iam )	printf(".. symbfact(): relax %4d, maxsuper %4d, fill %4d\n",
		       sp_ienv(2), sp_ienv(3), sp_ienv(6));
#endif
    t = SuperLU_timer_();
    iinfo = symbfact(iam, &AC, perm_c, etree, perm_r, Glu);
    stat.utime[SYMBFAC] = SuperLU_timer_() - t;
    SUPERLU_FREE(etree);

    if ( !iam ) {
	if ( iinfo == 0 ) {
	    printf("\tNo of supers        %8d\n", Glu->supno[n-1]+1);
	    printf("\tSize of super G(L)  %8d\n", Glu->xlsub[n]);
	    printf("\tSize of G(U)        %8d\n", Glu->xusub[n]);
	    QuerySpace(n, Glu, &symb_mem_usage);
	    printf("\tSYMBfact:\tL\\U MB %.2f\ttotal MB %.2f\texpansions %d\n",
		   symb_mem_usage.for_lu*1e-6, 
		   symb_mem_usage.total*1e-6,
		   symb_mem_usage.expansions);
	} else {
	    printf("symbfact() error returns %d\n", iinfo);
	}
    }

    /* Distribute the L and U factors onto the process grid. */
    t = SuperLU_timer_();
    zdistribute(n, &AC, Glu, grid, Llu);
    stat.utime[DIST] = SuperLU_timer_() - t;

    /* Deallocate storage used in the symbolic factor routines. */
    iinfo = LUSubFree(work, lwork, Glu);

    /* Perform numerical factorization. */
    t = SuperLU_timer_();
    pzgstrf(&AC, anorm, Glu, grid, Llu, &stat, info);
    stat.utime[FACT] = SuperLU_timer_() - t;

    {
	float for_lu, total;
	zQuerySpace(n, Glu, Llu, grid, &num_mem_usage);
	MPI_Reduce( &num_mem_usage.for_lu, &for_lu,
		   1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
	MPI_Reduce( &num_mem_usage.total, &total,
		   1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
	if ( !iam )
	    printf("\tNUMfact:\tL\\U MB %.2f\ttotal MB %.2f\n",
		   for_lu*1e-6, total*1e-6);
    }
    
#if ( PRNTlevel>=1 )
    if ( !iam ) printf(".. pzgstrf INFO = %d\n", *info);
#endif

    /* Scale the right-hand side. */
    if ( notran ) {
	if ( rowequ ) {
	    b_col = B;
	    for (j = 0; j < nrhs; ++j) {
		for (i = 0; i < m; ++i) zd_mult(&b_col[i], &b_col[i], R[i]);
		b_col += ldb;
	    }
	}
    } else if ( colequ ) {
	b_col = B;
	for (j = 0; j < nrhs; ++j) {
	    for (i = 0; i < m; ++i) zd_mult(&b_col[i], &b_col[i], C[i]);
	    b_col += ldb;
	}
    }

    /* Permute the right-hand side to form Pr*Pc*B. */
    if ( !(b_work = doublecomplexMalloc(n)) )
	ABORT("Malloc fails for b_work[]");
    if ( notran ) {
	b_col = B;
	for (j = 0; j < nrhs; ++j) {
	    for (i = 0; i < n; ++i) b_work[perm_c[i]] = b_col[i];
	    if ( static_pivot )
		for (i = 0; i < n; ++i) b_col[perm_r[i]] = b_work[i];
	    else
		for (i = 0; i < n; ++i) b_col[i] = b_work[i];
	    b_col += ldb;
	}
    }

#if 0
    /* Distribute the right-hand side according to supernode partition. */
    {
	int_t jj, k, knsupc, myrow, nsupers, *ilsum, *xsup;
	doublecomplex *b_col, *B_tmp;
	
	myrow = MYROW( iam, grid );
	xsup = Glu->xsup;
	nsupers = Glu->supno[n-1]+1;
	ilsum = Llu->ilsum;
	l_ldb = Llu->ldalsum;  /* Computed from ZDISTRIBUTE. */
	if ( !(b = doublecomplexMalloc(l_ldb * nrhs)) )
	    ABORT("Malloc fails for b[]");
	jj = 0;
	for (k = 0; k < nsupers; ++k) {
	    knsupc = SuperSize( k );
	    if ( myrow == PROW( k, grid ) ) {
		ii = ilsum[LBi(k,grid)];
		b_col = b + ii;
		B_tmp = B + jj;
		for (j = 0; j < nrhs; ++j) {
		    for (i = 0; i < knsupc; ++i) b_col[i] = B_tmp[i];
		    b_col += l_ldb;  B_tmp += ldb;
		}
	    }
	    jj += knsupc;
	}
#if ( DEBUGlevel>=1 )
	PrintDoublecomplex("Distributed b", l_ldb, b);
#endif
    }
#endif


#if 0
    /* Save a copy of the right-hand side. */
    l_ldx = l_ldb;
    if ( !(x = doublecomplexMalloc(l_ldx * nrhs)) )
	ABORT("Malloc fails for x[]");
    x_col = x;  b_col = b;
    for (j = 0; j < nrhs; ++j) {
	for (i = 0; i < l_ldb; ++i) x_col[i] = b_col[i];
	x_col += l_ldx;  b_col += l_ldb;
    }
#endif

#if ( DEBUGlevel>=1 )
    PrintDoublecomplex("B", n, B);
#endif

    /* Solve lower and upper triangular systems. */
    t = SuperLU_timer_();
    pzgstrs(n, perm_r, perm_c, Glu, grid, Llu, B, ldb, nrhs,
	    &stat, info);
    stat.utime[SOLVE] = SuperLU_timer_() - t;

#if 0
    /* Improve the solution by iterative refinement. */
    t = SuperLU_timer_();
    pzgsrfs(n, &AC, anorm, perm_r, perm_c, Glu, Llu, grid, b, l_ldb,
	    x, l_ldx, nrhs, &stat, info);
    stat.utime[REFINE] = SuperLU_timer_() - t;
#endif


    /* Compute the final solution X <= Pc'*X. */
    for (i = 0; i < nrhs; i++) {
	b_col = &B[i*ldb];
	for (j = 0; j < n; ++j) b_work[j] = b_col[perm_c[j]];
	for (j = 0; j < n; ++j) b_col[j] = b_work[j];
    }
	
    /* Transform the solution matrix X to a solution of the original system
       before the equilibration. */
    if ( notran ) {
	if ( colequ ) {
	    b_col = B;
	    for (j = 0; j < nrhs; ++j) {
		for (i = 0; i < n; ++i)
                    zd_mult(&b_col[i], &b_col[i], C[i]);
		b_col += ldb;
	    }
	}
    } else if ( rowequ ) {
	b_col = B;
	for (j = 0; j < nrhs; ++j) {
	    for (i = 0; i < n; ++i)
		zd_mult(&b_col[i], &b_col[i], R[i]);
            b_col += ldb;
	}
    }


    PStatPrint(&stat, grid);

    PStatFree(&stat);
    SUPERLU_FREE(perm_c);
    SUPERLU_FREE(perm_r);
    if ( equil ) {
	SUPERLU_FREE(C);
	SUPERLU_FREE(R);
    }
#if 0
    SUPERLU_FREE(b);
    SUPERLU_FREE(x);
#endif
    SUPERLU_FREE(b_work);
    Destroy_CompCol_Permuted(&AC);

}
