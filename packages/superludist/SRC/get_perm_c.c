/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include "superlu_ddefs.h"
#include "colamd.h"

void
get_colamd_dist(
		const int_t m,  /* number of rows in matrix A. */
		const int_t n,  /* number of columns in matrix A. */
		const int_t nnz,/* number of nonzeros in matrix A. */
		int_t *colptr,  /* column pointer of size n+1 for matrix A. */
		int_t *rowind,  /* row indices of size nz for matrix A. */
		int_t *perm_c   /* out - the column permutation vector. */
		)
{
    int_t Alen, *A, i, info, *p;
    double *knobs;

    Alen = colamd_recommended(nnz, m, n);

    if ( !(knobs = (double *) SUPERLU_MALLOC(COLAMD_KNOBS * sizeof(double))) )
        ABORT("Malloc fails for knobs");
    colamd_set_defaults(knobs);

    if (!(A = (int_t *) SUPERLU_MALLOC(Alen * sizeof(int_t))) )
        ABORT("Malloc fails for A[]");
    if (!(p = (int_t *) SUPERLU_MALLOC((n+1) * sizeof(int_t))) )
        ABORT("Malloc fails for p[]");
    for (i = 0; i <= n; ++i) p[i] = colptr[i];
    for (i = 0; i < nnz; ++i) A[i] = rowind[i];
    info = colamd(m, n, Alen, A, p, knobs);
    if ( info == FALSE ) 
        ABORT("COLAMD failed");

    for (i = 0; i < n; ++i) perm_c[p[i]] = i;

    SUPERLU_FREE(knobs);
    SUPERLU_FREE(A);
    SUPERLU_FREE(p);
}


void
getata_dist(
	    const int_t m,    /* number of rows in matrix A. */
	    const int_t n,    /* number of columns in matrix A. */
	    const int_t nz,   /* number of nonzeros in matrix A */
	    int_t *colptr,    /* column pointer of size n+1 for matrix A. */
	    int_t *rowind,    /* row indices of size nz for matrix A. */
	    int_t *atanz,     /* out - on exit, returns the actual number of
				 nonzeros in matrix A'*A. */
	    int_t **ata_colptr, /* out - size n+1 */
	    int_t **ata_rowind  /* out - size *atanz */
	    )
{
/*
 * Purpose
 * =======
 *
 * Form the structure of A'*A. A is an m-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'*A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (ata_colptr, ata_rowind).
 *
 * This routine is modified from GETATA routine by Tim Davis.
 * The complexity of this algorithm is: SUM_{i=1,m} r(i)^2,
 * i.e., the sum of the square of the row counts.
 *
 * Questions
 * =========
 *     o  Do I need to withhold the *dense* rows?
 *     o  How do I know the number of nonzeros in A'*A?
 * 
 */
    register int_t i, j, k, col, num_nz, ti, trow;
    int_t *marker, *b_colptr, *b_rowind;
    int_t *t_colptr, *t_rowind; /* a column oriented form of T = A' */

    if ( !(marker = (int_t*) SUPERLU_MALLOC( (SUPERLU_MAX(m,n)+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (int_t*) SUPERLU_MALLOC( (m+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC t_colptr[]");
    if ( !(t_rowind = (int_t*) SUPERLU_MALLOC( nz * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for t_rowind[]");

    
    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < m; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i)
	    ++marker[rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < m; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose the matrix from A to T */
    for (j = 0; j < n; ++j)
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    col = rowind[i];
	    t_rowind[marker[col]] = j;
	    ++marker[col];
	}

    
    /* ----------------------------------------------------------------
       compute B = T * A, where column j of B is:

       Struct (B_*j) =    UNION   ( Struct (T_*k) )
                        A_kj != 0

       do not include the diagonal entry
   
       ( Partition A as: A = (A_*1, ..., A_*n)
         Then B = T * A = (T * A_*1, ..., T * A_*n), where
         T * A_*j = (T_*1, ..., T_*m) * A_*j.  )
       ---------------------------------------------------------------- */

    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;

    /* First pass determines number of nonzeros in B */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
	    k = rowind[i];
	    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
		trow = t_rowind[ti];
		if ( marker[trow] != j ) {
		    marker[trow] = j;
		    num_nz++;
		}
	    }
	}
    }
    *atanz = num_nz;
    
    /* Allocate storage for A'*A */
    if ( !(*ata_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for ata_colptr[]");
    if ( *atanz ) {
	if ( !(*ata_rowind = (int_t*)SUPERLU_MALLOC(*atanz*sizeof(int_t)) ) ) {
	    fprintf(stderr, ".. atanz = %ld\n", *atanz);
	    ABORT("SUPERLU_MALLOC fails for ata_rowind[]");
	}
    }
    b_colptr = *ata_colptr; /* aliasing */
    b_rowind = *ata_rowind;
    
    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;
    
    /* Compute each column of B, one at a time */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	b_colptr[j] = num_nz;
	
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    /* A_kj is nonzero, add pattern of column T_*k to B_*j */
	    k = rowind[i];
	    for (ti = t_colptr[k]; ti < t_colptr[k+1]; ++ti) {
		trow = t_rowind[ti];
		if ( marker[trow] != j ) {
		    marker[trow] = j;
		    b_rowind[num_nz++] = trow;
		}
	    }
	}
    }
    b_colptr[n] = num_nz;
       
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_colptr);
    SUPERLU_FREE(t_rowind);
}


void
a_plus_at_dist(
	       const int_t n,    /* number of columns in matrix A. */
	       const int_t nz,   /* number of nonzeros in matrix A */
	       int_t *colptr,    /* column pointer of size n+1 for matrix A. */
	       int_t *rowind,    /* row indices of size nz for matrix A. */
	       int_t *bnz,       /* out - on exit, returns the actual number of
				    nonzeros in matrix A'+A. */
	       int_t **b_colptr, /* out - size n+1 */
	       int_t **b_rowind  /* out - size *bnz */
	       )
{
/*
 * Purpose
 * =======
 *
 * Form the structure of A'+A. A is an n-by-n matrix in column oriented
 * format represented by (colptr, rowind). The output A'+A is in column
 * oriented format (symmetrically, also row oriented), represented by
 * (b_colptr, b_rowind).
 *
 */
    register int_t i, j, k, col, num_nz;
    int_t *t_colptr, *t_rowind; /* a column oriented form of T = A' */
    int_t *marker;

    if ( !(marker = (int_t*) SUPERLU_MALLOC( n * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for marker[]");
    if ( !(t_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for t_colptr[]");
    if ( !(t_rowind = (int_t*) SUPERLU_MALLOC( nz * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails t_rowind[]");

    
    /* Get counts of each column of T, and set up column pointers */
    for (i = 0; i < n; ++i) marker[i] = 0;
    for (j = 0; j < n; ++j) {
	for (i = colptr[j]; i < colptr[j+1]; ++i)
	    ++marker[rowind[i]];
    }
    t_colptr[0] = 0;
    for (i = 0; i < n; ++i) {
	t_colptr[i+1] = t_colptr[i] + marker[i];
	marker[i] = t_colptr[i];
    }

    /* Transpose the matrix from A to T */
    for (j = 0; j < n; ++j)
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    col = rowind[i];
	    t_rowind[marker[col]] = j;
	    ++marker[col];
	}


    /* ----------------------------------------------------------------
       compute B = A + T, where column j of B is:

       Struct (B_*j) = Struct (A_*k) UNION Struct (T_*k)

       do not include the diagonal entry
       ---------------------------------------------------------------- */

    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;

    /* First pass determines number of nonzeros in B */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	/* Add pattern of column A_*k to B_*j */
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    k = rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		++num_nz;
	    }
	}

	/* Add pattern of column T_*k to B_*j */
	for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	    k = t_rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		++num_nz;
	    }
	}
    }
    *bnz = num_nz;
    
    /* Allocate storage for A+A' */
    if ( !(*b_colptr = (int_t*) SUPERLU_MALLOC( (n+1) * sizeof(int_t)) ) )
	ABORT("SUPERLU_MALLOC fails for b_colptr[]");
    if ( *bnz ) {
	if ( !(*b_rowind = (int_t*) SUPERLU_MALLOC( *bnz * sizeof(int_t)) ) )
	    ABORT("SUPERLU_MALLOC fails for b_rowind[]");
    }
    
    /* Zero the diagonal flag */
    for (i = 0; i < n; ++i) marker[i] = -1;
    
    /* Compute each column of B, one at a time */
    num_nz = 0;
    for (j = 0; j < n; ++j) {
	(*b_colptr)[j] = num_nz;
	
	/* Flag the diagonal so it's not included in the B matrix */
	marker[j] = j;

	/* Add pattern of column A_*k to B_*j */
	for (i = colptr[j]; i < colptr[j+1]; ++i) {
	    k = rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		(*b_rowind)[num_nz++] = k;
	    }
	}

	/* Add pattern of column T_*k to B_*j */
	for (i = t_colptr[j]; i < t_colptr[j+1]; ++i) {
	    k = t_rowind[i];
	    if ( marker[k] != j ) {
		marker[k] = j;
		(*b_rowind)[num_nz++] = k;
	    }
	}
    }
    (*b_colptr)[n] = num_nz;
       
    SUPERLU_FREE(marker);
    SUPERLU_FREE(t_colptr);
    SUPERLU_FREE(t_rowind);
} /* a_plus_at_dist */


void
get_perm_c_dist(int_t pnum, int_t ispec, SuperMatrix *A, int_t *perm_c)
/*
 * Purpose
 * =======
 *
 * GET_PERM_C_DIST obtains a permutation matrix Pc, by applying the multiple
 * minimum degree ordering code by Joseph Liu to matrix A'*A or A+A',
 * or using approximate minimum degree column ordering by Davis et. al.
 * The LU factorization of A*Pc tends to have less fill than the LU 
 * factorization of A.
 *
 * Arguments
 * =========
 *
 * ispec   (input) colperm_t
 *         Specifies what type of column permutation to use to reduce fill.
 *         = NATURAL: natural ordering (i.e., Pc = I)
 *         = MMD_AT_PLUS_A: minimum degree ordering on structure of A'+A
 *         = MMD_ATA: minimum degree ordering on structure of A'*A
 *         = COLAMD: approximate minimum degree for unsymmetric matrices
 * 
 * A       (input) SuperMatrix*
 *         Matrix A in A*X=B, of dimension (A->nrow, A->ncol). The number
 *         of the linear equations is A->nrow. Currently, the type of A 
 *         can be: Stype = NC; Dtype = _D; Mtype = GE. In the future,
 *         more general A can be handled.
 *
 * perm_c  (output) int*
 *	   Column permutation vector of size A->ncol, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 *
 */
{
    NCformat *Astore = A->Store;
    int_t m, n, bnz, *b_colptr, i;
    int_t delta, maxint, nofsub, *invp;
    int_t *b_rowind, *dhead, *qsize, *llist, *marker;
    double t, SuperLU_timer_();

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(pnum, "Enter get_perm_c_dist()");
#endif

    m = A->nrow;
    n = A->ncol;

    t = SuperLU_timer_();
    switch ( ispec ) {
        case NATURAL: /* Natural ordering */
	      for (i = 0; i < n; ++i) perm_c[i] = i;
#if ( PRNTlevel>=1 )
	      if ( !pnum ) printf(".. Use natural column ordering.\n");
#endif
	      return;
        case MMD_AT_PLUS_A: /* Minimum degree ordering on A'+A */
	      if ( m != n ) ABORT("Matrix is not square");
	      a_plus_at_dist(n, Astore->nnz, Astore->colptr, Astore->rowind,
			     &bnz, &b_colptr, &b_rowind);
	      t = SuperLU_timer_() - t;
	      /*printf("Form A'+A time = %8.3f\n", t);*/
#if ( PRNTlevel>=1 )
	      if ( !pnum ) printf(".. Use minimum degree ordering on A'+A.\n");
#endif
	      break;
        case MMD_ATA: /* Minimum degree ordering on A'*A */
	      getata_dist(m, n, Astore->nnz, Astore->colptr, Astore->rowind,
			  &bnz, &b_colptr, &b_rowind);
	      t = SuperLU_timer_() - t;
	      /*printf("Form A'*A time = %8.3f\n", t);*/
#if ( PRNTlevel>=1 )
	      if ( !pnum ) printf(".. Use minimum degree ordering on A'*A.\n");
#endif
	      break;
        case COLAMD: /* Approximate minimum degree column ordering. */
	      get_colamd_dist(m, n, Astore->nnz, Astore->colptr,
			      Astore->rowind, perm_c);
#if ( PRNTlevel>=1 )
	      if ( !pnum )
		  printf(".. Use approximate minimum degree column ordering.\n");
#endif
	      return; 
        default:
	      ABORT("Invalid ISPEC");
    }

    if ( bnz ) {
	t = SuperLU_timer_();

	/* Initialize and allocate storage for GENMMD. */
	delta = 1; /* DELTA is a parameter to allow the choice of nodes
		      whose degree <= min-degree + DELTA. */
	maxint = 2147483647; /* 2**31 - 1 */
	invp = (int_t *) SUPERLU_MALLOC((n+delta)*sizeof(int_t));
	if ( !invp ) ABORT("SUPERLU_MALLOC fails for invp.");
	dhead = (int_t *) SUPERLU_MALLOC((n+delta)*sizeof(int_t));
	if ( !dhead ) ABORT("SUPERLU_MALLOC fails for dhead.");
	qsize = (int_t *) SUPERLU_MALLOC((n+delta)*sizeof(int_t));
	if ( !qsize ) ABORT("SUPERLU_MALLOC fails for qsize.");
	llist = (int_t *) SUPERLU_MALLOC(n*sizeof(int_t));
	if ( !llist ) ABORT("SUPERLU_MALLOC fails for llist.");
	marker = (int_t *) SUPERLU_MALLOC(n*sizeof(int_t));
	if ( !marker ) ABORT("SUPERLU_MALLOC fails for marker.");
	
	/* Transform adjacency list into 1-based indexing required by GENMMD.*/
	for (i = 0; i <= n; ++i) ++b_colptr[i];
	for (i = 0; i < bnz; ++i) ++b_rowind[i];
	
	genmmd_dist_(&n, b_colptr, b_rowind, perm_c, invp, &delta, dhead, 
		     qsize, llist, marker, &maxint, &nofsub);

	/* Transform perm_c into 0-based indexing. */
	for (i = 0; i < n; ++i) --perm_c[i];

	SUPERLU_FREE(b_colptr);
	SUPERLU_FREE(b_rowind);
	SUPERLU_FREE(invp);
	SUPERLU_FREE(dhead);
	SUPERLU_FREE(qsize);
	SUPERLU_FREE(llist);
	SUPERLU_FREE(marker);

	t = SuperLU_timer_() - t;
	/*    printf("call GENMMD time = %8.3f\n", t);*/

    } else { /* Empty adjacency structure */
	for (i = 0; i < n; ++i) perm_c[i] = i;
    }
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(pnum, "Exit get_perm_c_dist()");
#endif
} /* get_perm_c_dist */
