#include <math.h>
#include "superlu_ddefs.h"


void
pdgssvx_ABglobal(superlu_options_t *options, SuperMatrix *A, 
		 ScalePermstruct_t *ScalePermstruct,
		 double B[], int ldb, int nrhs, gridinfo_t *grid,
		 LUstruct_t *LUstruct, double *berr,
		 SuperLUStat_t *stat, int *info)
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
 * pdgssvx_ABglobal solves a system of linear equations A*X=B,
 * by using Gaussian elimination with "static pivoting" to
 * compute the LU factorization of A.
 *
 * Static pivoting is a technique that combines the numerical stability
 * of partial pivoting with the scalability of Cholesky (no pivoting),
 * to run accurately and efficiently on large numbers of processors.
 *
 * See our paper at http://www.nersc.gov/~xiaoye/SuperLU/ for a detailed
 * description of the parallel algorithms.
 *
 * Here are the options for using this code:
 *
 *   1. Independent of all the other options specified below, the
 *      user must supply
 *
 *      -  B, the matrix of right hand sides, and its dimensions ldb and nrhs
 *      -  grid, a structure describing the 2D processor mesh
 *      -  options->IterRefine, which determines whether or not to
 *            improve the accuracy of the computed solution using 
 *            iterative refinement
 *
 *      On output, B is overwritten with the solution X.
 *
 *   2. Depending on options->Fact, the user has several options
 *      for solving A*X=B. The standard option is for factoring
 *      A "from scratch". (The other options, described below,
 *      are used when A is sufficiently similar to a previously 
 *      solved problem to save time by reusing part or all of 
 *      the previous factorization.)
 *
 *      -  options->Fact = DOFACT: A is factored "from scratch"
 *
 *      In this case the user must also supply
 *
 *      -  A, the input matrix
 *
 *      as well as the following options, which are described in more 
 *      detail below:
 *
 *      -  options->Equil,   to specify how to scale the rows and columns
 *                           of A to "equilibrate" it (to try to reduce its
 *                           condition number and so improve the
 *                           accuracy of the computed solution)
 *
 *      -  options->RowPerm, to specify how to permute the rows of A
 *                           (typically to control numerical stability)
 *
 *      -  options->ColPerm, to specify how to permute the columns of A
 *                           (typically to control fill-in and enhance
 *                           parallelism during factorization)
 *
 *      -  options->ReplaceTinyPivot, to specify how to deal with tiny
 *                           pivots encountered during factorization
 *                           (to control numerical stability)
 *
 *      The outputs returned include
 *         
 *      -  ScalePermstruct,  modified to describe how the input matrix A
 *                           was equilibrated and permuted:
 *         -  ScalePermstruct->DiagScale, indicates whether the rows and/or
 *                                        columns of A were scaled
 *         -  ScalePermstruct->R, array of row scale factors
 *         -  ScalePermstruct->C, array of column scale factors
 *         -  ScalePermstruct->perm_r, row permutation vector
 *         -  ScalePermstruct->perm_c, column permutation vector
 *
 *            (part of ScalePermstruct may also need to be supplied on input,
 *             depending on options->RowPerm and options->ColPerm as described 
 *             later).
 *
 *      -  A, the input matrix A overwritten by the scaled and permuted matrix
 *                Pc*Pr*diag(R)*A*diag(C)
 *             where 
 *                Pr and Pc are row and columns permutation matrices determined
 *                  by ScalePermstruct->perm_r and ScalePermstruct->perm_c, 
 *                  respectively, and 
 *                diag(R) and diag(C) are diagonal scaling matrices determined
 *                  by ScalePermstruct->DiagScale, ScalePermstruct->R and 
 *                  ScalePermstruct->C
 *
 *      -  LUstruct, which contains the L and U factorization of A1 where
 *
 *                A1 = Pc*Pr*diag(R)*A*diag(C)*Pc^T = L*U
 *
 *              (Note that A1 = Aout * Pc^T, where Aout is the matrix stored
 *               in A on output.)
 *
 *   3. The second value of options->Fact assumes that a matrix with the same
 *      sparsity pattern as A has already been factored:
 *     
 *      -  options->Fact = SamePattern: A is factored, assuming that it has
 *            the same nonzero pattern as a previously factored matrix. In this
 *            case the algorithm saves time by reusing the previously computed
 *            column permutation vector stored in ScalePermstruct->perm_c
 *            and the "elimination tree" of A stored in LUstruct->etree
 *
 *      In this case the user must still specify the following options
 *      as before:
 *
 *      -  options->Equil
 *      -  options->RowPerm
 *      -  options->ReplaceTinyPivot
 *
 *      but not options->ColPerm, whose value is ignored. This is because the
 *      previous column permutation from ScalePermstruct->perm_c is used as
 *      input. The user must also supply 
 *
 *      -  A, the input matrix
 *      -  ScalePermstruct->perm_c, the column permutation
 *      -  LUstruct->etree, the elimination tree
 *
 *      The outputs returned include
 *         
 *      -  A, the input matrix A overwritten by the scaled and permuted matrix
 *            as described above
 *      -  ScalePermstruct,  modified to describe how the input matrix A was
 *                           equilibrated and row permuted
 *      -  LUstruct, modified to contain the new L and U factors
 *
 *   4. The third value of options->Fact assumes that a matrix B with the same
 *      sparsity pattern as A has already been factored, and where the
 *      row permutation of B can be reused for A. This is useful when A and B
 *      have similar numerical values, so that the same row permutation
 *      will make both factorizations numerically stable. This lets us reuse
 *      all of the previously computed structure of L and U.
 *
 *      -  options->Fact = SamePattern_SameRowPerm: A is factored,
 *            assuming not only the same nonzero pattern as the previously
 *            factored matrix B, but reusing B's row permutation.
 *
 *      In this case the user must still specify the following options
 *      as before:
 *
 *      -  options->Equil
 *      -  options->ReplaceTinyPivot
 *
 *      but not options->RowPerm or options->ColPerm, whose values are ignored.
 *      This is because the permutations from ScalePermstruct->perm_r and
 *      ScalePermstruct->perm_c are used as input.
 *
 *      The user must also supply 
 *
 *      -  A, the input matrix
 *      -  ScalePermstruct->DiagScale, how the previous matrix was row and/or
 *                                     column scaled
 *      -  ScalePermstruct->R, the row scalings of the previous matrix, if any
 *      -  ScalePermstruct->C, the columns scalings of the previous matrix, 
 *                             if any
 *      -  ScalePermstruct->perm_r, the row permutation of the previous matrix
 *      -  ScalePermstruct->perm_c, the column permutation of the previous 
 *                                  matrix
 *      -  all of LUstruct, the previously computed information about L and U
 *                (the actual numerical values of L and U stored in
 *                 LUstruct->Llu are ignored)
 *
 *      The outputs returned include
 *         
 *      -  A, the input matrix A overwritten by the scaled and permuted matrix
 *            as described above
 *      -  ScalePermstruct,  modified to describe how the input matrix A was
 *                           equilibrated 
 *                  (thus ScalePermstruct->DiagScale, R and C may be modified)
 *      -  LUstruct, modified to contain the new L and U factors
 *
 *   5. The fourth and last value of options->Fact assumes that A is
 *      identical to a matrix that has already been factored on a previous 
 *      call, and reuses its entire LU factorization
 *
 *      -  options->Fact = Factored: A is identical to a previously
 *            factorized matrix, so the entire previous factorization
 *            can be reused.
 *
 *      In this case all the other options mentioned above are ignored
 *      (options->Equil, options->RowPerm, options->ColPerm, 
 *       options->ReplaceTinyPivot)
 *       options->Iterefine is not ignored.
 *
 *      The user must also supply 
 *
 *      -  A, the unfactored matrix, only in the case that iterative refinment
 *            is to be done (specifically A must be the output A from 
 *            the previous call, so that it has been scaled and permuted)
 *      -  all of ScalePermstruct
 *      -  all of LUstruct, including the actual numerical values of L and U
 *
 *      all of which are unmodified on output.
 *         
 * Arguments
 * =========
 *
 * options (input) superlu_options_t*
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed.
 *         The following fields should be defined for this structure:
 *         
 *         o Fact (fact_t)
 *           Specifies whether or not the factored form of the matrix
 *           A is supplied on entry, and if not, how the matrix A should
 *           be factorized based on the previous history.
 *
 *           = DOFACT: The matrix A will be factorized from scratch.
 *                 Inputs:  A
 *                          options->Equil, RowPerm, ColPerm, ReplaceTinyPivot
 *                 Outputs: modified A
 *                             (possibly row and/or column scaled and/or 
 *                              permuted)
 *                          all of ScalePermstruct
 *                          all of LUstruct
 *
 *           = SamePattern: the matrix A will be factorized assuming
 *             that a factorization of a matrix with the same sparsity
 *             pattern was performed prior to this one. Therefore, this
 *             factorization will reuse column permutation vector 
 *             ScalePermstruct->perm_c and the elimination tree
 *             LUstruct->etree
 *                 Inputs:  A
 *                          options->Equil, RowPerm, ReplaceTinyPivot
 *                          ScalePermstruct->perm_c
 *                          LUstruct->etree
 *                 Outputs: modified A
 *                             (possibly row and/or column scaled and/or 
 *                              permuted)
 *                          rest of ScalePermstruct (DiagScale, R, C, perm_r)
 *                          rest of LUstruct (GLU_persist, Llu)
 *
 *           = SamePattern_SameRowPerm: the matrix A will be factorized
 *             assuming that a factorization of a matrix with the same
 *             sparsity	pattern and similar numerical values was performed
 *             prior to this one. Therefore, this factorization will reuse
 *             both row and column scaling factors R and C, and the
 *             both row and column permutation vectors perm_r and perm_c,
 *             distributed data structure set up from the previous symbolic
 *             factorization.
 *                 Inputs:  A
 *                          options->Equil, ReplaceTinyPivot
 *                          all of ScalePermstruct
 *                          all of LUstruct
 *                 Outputs: modified A
 *                             (possibly row and/or column scaled and/or 
 *                              permuted)
 *                          modified LUstruct->Llu
 *           = FACTORED: the matrix A is already factored.
 *                 Inputs:  all of ScalePermstruct
 *                          all of LUstruct
 *
 *         o Equil (equi_t)
 *           Specifies whether to equilibrate the system.
 *           = NEQU:  no equilibration.
 *           = EQUI: scaling factors are computed to equilibrate the system:
 *                      diag(R)*A*diag(C)*inv(diag(C))*X = diag(R)*B.
 *                  Whether or not the system will be equilibrated depends
 *                  on the scaling of the matrix A, but if equilibration is
 *                  used, A is overwritten by diag(R)*A*diag(C) and B by
 *                  diag(R)*B.
 *
 *         o RowPerm (rowperm_t)
 *           Specifies how to permute rows of the matrix A.
 *           = NATURAL:   use the natural ordering.
 *           = LargeDiag: use the Duff/Koster algorithm to permute rows of
 *                        the original matrix to make the diagonal large
 *                        relative to the off-diagonal.
 *           = MY_PERMR:  use the ordering given in ScalePermstruct->perm_r
 *                        input by the user.
 *           
 *         o ColPerm (colperm_t)
 *           Specifies what type of column permutation to use to reduce fill.
 *           = NATURAL:       natural ordering.
 *           = MMD_ATA:       minimum degree ordering on structure of A'*A.
 *           = MMD_AT_PLUS_A: minimum degree ordering on structure of A'+A.
 *           = COLAMD:        approximate minimum degree column ordering.
 *           = MY_PERMC:      the ordering given in ScalePermstruct->perm_c.
 *         
 *         o ReplaceTinyPivot (place_t)
 *           = NO:  do not modify pivots
 *           = YES: replace tiny pivots by sqrt(epsilon)*norm(A) during 
 *                  LU factorization.
 *
 *         o IterRefine (IterRefine_t)
 *           Specifies how to perform iterative refinement.
 *           = NO:     no iterative refinement.
 *           = DOUBLE: accumulate residual in double precision.
 *           = EXTRA:  accumulate residual in extra precision.
 *
 *         NOTE: all options must be indentical on all processes when
 *               calling this routine.
 *
 * A (input/output) SuperMatrix*
 *         On entry, matrix A in A*X=B, of dimension (A->nrow, A->ncol).
 *         The number of linear equations is A->nrow. The type of A must be:
 *         Stype = NC; Dtype = _D; Mtype = GE. That is, A is stored in
 *         compressed column format (also known as Harwell-Boeing format).
 *         See supermatrix.h for the definition of 'SuperMatrix'.
 *         This routine only handles square A, however, the LU factorization
 *         routine pdgstrf can factorize rectangular matrices.
 *         On exit, A may be overwtirren by Pc*Pr*diag(R)*A*diag(C),
 *         depending on ScalePermstruct->DiagScale, options->RowPerm and
 *         options->colpem:
 *             if ScalePermstruct->DiagScale != NOEQUIL, A is overwritten by
 *                diag(R)*A*diag(C).
 *             if options->RowPerm != NATURAL, A is further overwritten by
 *                Pr*diag(R)*A*diag(C).
 *             if options->ColPerm != NATURAL, A is further overwritten by
 *                Pc*Pr*diag(R)*A*diag(C).
 *         If all the above condition are true, the LU decomposition is
 *         performed on the matrix Pc*Pr*diag(R)*A*diag(C)*Pc^T.
 *
 *         NOTE: Currently, A must reside in all processes when calling
 *               this routine.
 *
 * ScalePermstruct (input/output) ScalePermstruct_t*
 *         The data structure to store the scaling and permutation vectors
 *         describing the transformations performed to the matrix A.
 *         It contains the following fields:
 *
 *         o DiagScale (DiagScale_t)
 *           Specifies the form of equilibration that was done.
 *           = NOEQUIL: no equilibration.
 *           = ROW:     row equilibration, i.e., A was premultiplied by
 *                      diag(R).
 *           = COL:     Column equilibration, i.e., A was postmultiplied
 *                      by diag(C).
 *           = BOTH:    both row and column equilibration, i.e., A was 
 *                      replaced by diag(R)*A*diag(C).
 *           If options->Fact = FACTORED or SamePattern_SameRowPerm,
 *           DiagScale is an input argument; otherwise it is an output
 *           argument.
 *
 *         o perm_r (int*)
 *           Row permutation vector, which defines the permutation matrix Pr;
 *           perm_r[i] = j means row i of A is in position j in Pr*A.
 *           If options->RowPerm = MY_PERMR, or
 *           options->Fact = SamePattern_SameRowPerm, perm_r is an
 *           input argument; otherwise it is an output argument.
 *
 *         o perm_c (int*)
 *           Column permutation vector, which defines the 
 *           permutation matrix Pc; perm_c[i] = j means column i of A is 
 *           in position j in A*Pc.
 *           If options->ColPerm = MY_PERMC or options->Fact = SamePattern
 *           or options->Fact = SamePattern_SameRowPerm, perm_c is an
 *           input argument; otherwise, it is an output argument.
 *           On exit, perm_c may be overwritten by the product of the input
 *           perm_c and a permutation that postorders the elimination tree
 *           of Pc*A'*A*Pc'; perm_c is not changed if the elimination tree
 *           is already in postorder.
 *
 *         o R (double*) dimension (A->nrow)
 *           The row scale factors for A.
 *           If DiagScale = ROW or BOTH, A is multiplied on the left by 
 *                          diag(R).
 *           If DiagScale = NOEQUIL or COL, R is not defined.
 *           If options->Fact = FACTORED or SamePattern_SameRowPerm, R is
 *           an input argument; otherwise, R is an output argument.
 *
 *         o C (double*) dimension (A->ncol)
 *           The column scale factors for A.
 *           If DiagScale = COL or BOTH, A is multiplied on the right by 
 *                          diag(C).
 *           If DiagScale = NOEQUIL or ROW, C is not defined.
 *           If options->Fact = FACTORED or SamePattern_SameRowPerm, C is
 *           an input argument; otherwise, C is an output argument.
 *         
 * B       (input/output) double*
 *         On entry, the right-hand side matrix of dimension (A->nrow, nrhs).
 *         On exit, the solution matrix if info = 0;
 *
 *         NOTE: Currently, B must reside in all processes when calling
 *               this routine.
 *
 * ldb     (input) int (global)
 *         The leading dimension of matrix B.
 *
 * nrhs    (input) int (global)
 *         The number of right-hand sides.
 *         If nrhs = 0, only LU decomposition is performed, the forward
 *         and back substitution are skipped.
 *
 * grid    (input) gridinfo_t*
 *         The 2D process mesh. It contains the MPI communicator, the number
 *         of process rows (NPROW), the number of process columns (NPCOL),
 *         and my process rank. It is an input argument to all the
 *         parallel routines.
 *         Grid can be initialized by subroutine SUPERLU_GRIDINIT.
 *         See superlu_ddefs.h for the definition of 'gridinfo_t'.
 *
 * LUstruct (input/output) LUstruct_t*
 *         The data structures to store the distributed L and U factors.
 *         It contains the following fields:
 *
 *         o etree (int*) dimension (A->ncol)
 *           Elimination tree of Pc*(A'+A)*Pc' or Pc*A'*A*Pc', dimension A->ncol.
 *           It is computed in sp_colorder() during the first factorization,
 *           and is reused in the subsequent factorizations of the matrices
 *           with the same nonzero pattern.
 *           On exit of sp_colorder(), the columns of A are permuted so that
 *           the etree is in a certain postorder. This postorder is reflected
 *           in ScalePermstruct->perm_c.
 *           NOTE:
 *           Etree is a vector of parent pointers for a forest whose vertices
 *           are the integers 0 to A->ncol-1; etree[root]==A->ncol.
 *
 *         o Glu_persist (Glu_persist_t*)
 *           Global data structure (xsup, supno) replicated on all processes,
 *           describing the supernode partition in the factored matrices
 *           L and U:
 *	       xsup[s] is the leading column of the s-th supernode,
 *             supno[i] is the supernode number to which column i belongs.
 *
 *         o Llu (LocalLU_t*)
 *           The distributed data structures to store L and U factors.
 *           See superlu_ddefs.h for the definition of 'LocalLU_t'.
 *
 * berr    (output) double*, dimension (nrhs)
 *         The componentwise relative backward error of each solution   
 *         vector X(j) (i.e., the smallest relative change in   
 *         any element of A or B that makes X(j) an exact solution).
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics on runtime and floating-point operation count.
 *        See util.h for the definition of 'SuperLUStat_t'.
 *
 * info    (output) int*
 *         = 0: successful exit
 *         > 0: if info = i, and i is
 *             <= A->ncol: U(i,i) is exactly zero. The factorization has
 *                been completed, but the factor U is exactly singular,
 *                so the solution could not be computed.
 *             > A->ncol: number of bytes allocated when memory allocation
 *                failure occurred, plus A->ncol.
 *
 *
 * See superlu_ddefs.h for the definitions of varioous data types.
 *
 */
    SuperMatrix AC;
    NCformat *Astore;
    NCPformat *ACstore;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    Glu_freeable_t *Glu_freeable;
            /* The nonzero structures of L and U factors, which are
	       replicated on all processrs.
	           (lsub, xlsub) contains the compressed subscript of
		                 supernodes in L.
          	   (usub, xusub) contains the compressed subscript of
		                 nonzero segments in U.
	      If options->Fact != SamePattern_SameRowPerm, they are 
	      computed by SYMBFACT routine, and then used by DDISTRIBUTE
	      routine. They will be freed after DDISTRIBUTE routine.
	      If options->Fact == SamePattern_SameRowPerm, these
	      structures are not used.                                  */
    fact_t   Fact;
    double   *a;
    int_t    *perm_r; /* row permutations from partial pivoting */
    int_t    *perm_c; /* column permutation vector */
    int_t    *etree;  /* elimination tree */
    int_t    *colptr, *rowind;
    int_t    colequ, Equil, factored, job, notran, rowequ;
    int_t    i, iinfo, j, irow, m, n, nnz, permc_spec;
    int      iam;
    int      ldx;  /* LDA for matrix X (global). */
    char     equed[1], norm[1];
    double   *C, *R, *C1, *R1, amax, anorm, colcnd, rowcnd;
    double   *X, *b_col, *b_work, *x_col;
    double   t;
    mem_usage_t             symb_mem_usage;
#if ( PRNTlevel>=1 )
    mem_usage_t             num_mem_usage;
#endif
#if ( PRNTlevel>= 2 )
    double                  dmin, dsum, dprod;
#endif

    perm_r = ScalePermstruct->perm_r;
    perm_c = ScalePermstruct->perm_c;

    /* Test input parameters. */
    *info = 0;
    Fact = options->Fact;
    if ( Fact < 0 || Fact > FACTORED )
	*info = -10;
    else if ( options->RowPerm < 0 || options->RowPerm > MY_PERMR )
	*info = -11;
    else if ( options->ColPerm < 0 || options->ColPerm > MY_PERMC )
	*info = -12;
    else if ( options->IterRefine < 0 || options->IterRefine > EXTRA )
	*info = -13;
    else if ( options->IterRefine == EXTRA ) {
	*info = -14;
	fprintf(stderr, "Extra precise iterative refinement: yet to support.");
    } else if ( A->nrow != A->ncol || A->nrow < 0 ||
         A->Stype != NC || A->Dtype != _D || A->Mtype != GE )
	*info = -2;
    else if ( ldb < A->nrow )
	*info = -4;
    else if ( nrhs <= 0 )
	*info = -5;
    if ( *info ) {
	i = -(*info);
	pxerbla("pdgssvx_ABglobal", grid, -*info);
	return;
    }

    /* Initialization */
    factored = (Fact == FACTORED);
    Equil = (!factored && options->Equil == EQUI);
    notran = (options->Trans == NOTRANS);
    iam = grid->iam;
    job = 5;
    m = A->nrow;
    n = A->ncol;
    Astore = (NCformat *) A->Store;
    nnz = Astore->nnz;
    a      = (double *) Astore->nzval;
    colptr = Astore->colptr;
    rowind = Astore->rowind;
    if ( factored || (Fact == SamePattern_SameRowPerm && Equil) ) {
	rowequ = (ScalePermstruct->DiagScale == ROW) ||
	         (ScalePermstruct->DiagScale == BOTH);
	colequ = (ScalePermstruct->DiagScale == COL) ||
	         (ScalePermstruct->DiagScale == BOTH);
    } else rowequ = colequ = FALSE;
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter pdgssvx_ABglobal()");
#endif
    perm_r = ScalePermstruct->perm_r;
    perm_c = ScalePermstruct->perm_c;
    etree = LUstruct->etree;
    R = ScalePermstruct->R;
    C = ScalePermstruct->C;
    if ( Equil ) {
	/* Allocate storage if not done so before. */
	switch ( ScalePermstruct->DiagScale ) {
	    case NOEQUIL:
		if ( !(R = (double *) doubleMalloc(m)) )
		    ABORT("Malloc fails for R[].");
	        if ( !(C = (double *) doubleMalloc(n)) )
		    ABORT("Malloc fails for C[].");
		ScalePermstruct->R = R;
		ScalePermstruct->C = C;
		break;
	    case ROW: 
	        if ( !(C = (double *) doubleMalloc(n)) )
		    ABORT("Malloc fails for C[].");
		ScalePermstruct->C = C;
		break;
	    case COL: 
		if ( !(R = (double *) doubleMalloc(m)) )
		    ABORT("Malloc fails for R[].");
		ScalePermstruct->R = R;
		break;
	}
    }

    /* ------------------------------------------------------------
       Diagonal scaling to equilibrate the matrix.
       ------------------------------------------------------------*/
    if ( Equil ) {
#if ( DEBUGlevel>=1 )
	CHECK_MALLOC(iam, "Enter equil");
#endif
	t = SuperLU_timer_();

	if ( Fact == SamePattern_SameRowPerm ) {
	    /* Reuse R and C. */
	    switch ( ScalePermstruct->DiagScale ) {
	      case NOEQUIL:
		break;
	      case ROW:
		for (j = 0; j < n; ++j) {
		    for (i = colptr[j]; i < colptr[j+1]; ++i) {
			irow = rowind[i];
			a[i] *= R[irow];       /* Scale rows. */
		    }
		}
		break;
	      case COL:
		for (j = 0; j < n; ++j)
		    for (i = colptr[j]; i < colptr[j+1]; ++i)
			a[i] *= C[j];          /* Scale columns. */
		break;
	      case BOTH: 
		for (j = 0; j < n; ++j) {
		    for (i = colptr[j]; i < colptr[j+1]; ++i) {
			irow = rowind[i];
			a[i] *= R[irow] * C[j]; /* Scale rows and columns. */
		    }
		}
	        break;
	    }
	} else {
	    if ( !iam ) {
		/* Compute row and column scalings to equilibrate matrix A. */
		dgsequ(A, R, C, &rowcnd, &colcnd, &amax, &iinfo);
	    
		MPI_Bcast( &iinfo, 1, mpi_int_t, 0, grid->comm );
		if ( iinfo == 0 ) {
		    MPI_Bcast( R,       m, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( C,       n, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( &rowcnd, 1, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( &colcnd, 1, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( &amax,   1, MPI_DOUBLE, 0, grid->comm );
		} else {
		    if ( iinfo > 0 ) {
			if ( iinfo <= m )
			    fprintf(stderr, "The %d-th row of A is exactly zero\n", 
				    iinfo);
			else fprintf(stderr, "The %d-th column of A is exactly zero\n", 
				     iinfo-n);
			exit(-1);
		    }
		}
	    } else {
		MPI_Bcast( &iinfo, 1, mpi_int_t, 0, grid->comm );
		if ( iinfo == 0 ) {
		    MPI_Bcast( R,       m, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( C,       n, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( &rowcnd, 1, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( &colcnd, 1, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( &amax,   1, MPI_DOUBLE, 0, grid->comm );
		} else {
		    ABORT("DGSEQU failed\n");
		}
	    }
	
	    /* Equilibrate matrix A. */
	    dlaqgs(A, R, C, rowcnd, colcnd, amax, equed);
	    if ( lsame_(equed, "R") ) {
		ScalePermstruct->DiagScale = ROW;
		rowequ = ROW;
	    } else if ( lsame_(equed, "C") ) {
		ScalePermstruct->DiagScale = COL;
		colequ = COL;
	    } else if ( lsame_(equed, "B") ) {
		ScalePermstruct->DiagScale = BOTH;
		rowequ = ROW;
		colequ = COL;
	    } else ScalePermstruct->DiagScale = NOEQUIL;

#if ( PRNTlevel>=1 )
	    if ( !iam ) {
		printf(".. equilibrated? *equed = %c\n", *equed);
		/*fflush(stdout);*/
	    }
#endif
	} /* if Fact ... */

	stat->utime[EQUIL] = SuperLU_timer_() - t;
#if ( DEBUGlevel>=1 )
	CHECK_MALLOC(iam, "Exit equil");
#endif
    } /* if Equil ... */
    
    /* ------------------------------------------------------------
       Permute rows of A. 
       ------------------------------------------------------------*/
    if ( options->RowPerm != NOROWPERM ) {
	t = SuperLU_timer_();

	if ( Fact == SamePattern_SameRowPerm /* Reuse perm_r. */
	    || options->RowPerm == MY_PERMR ) { /* Use my perm_r. */
            if ( !factored ) {
	       for (i = 0; i < colptr[n]; ++i) {
		    irow = rowind[i]; 
		    rowind[i] = perm_r[irow];
	       }
	    }
	} else if ( !factored ) {
	    if ( job == 5 ) {
		/* Allocate storage for scaling factors. */
		if ( !(R1 = (double *) SUPERLU_MALLOC(m * sizeof(double))) ) 
		    ABORT("SUPERLU_MALLOC fails for R1[]");
		if ( !(C1 = (double *) SUPERLU_MALLOC(n * sizeof(double))) )
		    ABORT("SUPERLU_MALLOC fails for C1[]");
	    }

	    if ( !iam ) {
		/* Process 0 finds a row permutation for large diagonal. */
		dldperm(job, m, nnz, colptr, rowind, a, perm_r, R1, C1);
		
		MPI_Bcast( perm_r, m, mpi_int_t, 0, grid->comm );
		if ( job == 5 && Equil ) {
		    MPI_Bcast( R1, m, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( C1, n, MPI_DOUBLE, 0, grid->comm );
		}
	    } else {
		MPI_Bcast( perm_r, m, mpi_int_t, 0, grid->comm );
		if ( job == 5 && Equil ) {
		    MPI_Bcast( R1, m, MPI_DOUBLE, 0, grid->comm );
		    MPI_Bcast( C1, n, MPI_DOUBLE, 0, grid->comm );
		}
	    }

#if ( PRNTlevel>=2 )
	    dmin = dlamch_("Overflow");
	    dsum = 0.0;
	    dprod = 1.0;
#endif
	    if ( job == 5 ) {
		if ( Equil ) {
		    for (i = 0; i < n; ++i) {
			R1[i] = exp(R1[i]);
			C1[i] = exp(C1[i]);
		    }
		    for (j = 0; j < n; ++j) {
			for (i = colptr[j]; i < colptr[j+1]; ++i) {
			    irow = rowind[i];
			    a[i] *= R1[irow] * C1[j]; /* Scale the matrix. */
			    rowind[i] = perm_r[irow];
#if ( PRNTlevel>=2 )
			    if ( rowind[i] == j ) /* New diagonal */
				dprod *= fabs(a[i]);
#endif
			}
		    }

		    /* Multiply together the scaling factors. */
		    if ( rowequ ) for (i = 0; i < m; ++i) R[i] *= R1[i];
		    else for (i = 0; i < m; ++i) R[i] = R1[i];
		    if ( colequ ) for (i = 0; i < n; ++i) C[i] *= C1[i];
		    else for (i = 0; i < n; ++i) C[i] = C1[i];
		    
		    ScalePermstruct->DiagScale = BOTH;
		    rowequ = colequ = 1;
		} else { /* No equilibration. */
		    for (i = colptr[0]; i < colptr[n]; ++i) {
			    irow = rowind[i];
			    rowind[i] = perm_r[irow];
		    }
		}
		SUPERLU_FREE (R1);
		SUPERLU_FREE (C1);
	    } else { /* job = 2,3,4 */
		for (j = 0; j < n; ++j) {
		    for (i = colptr[j]; i < colptr[j+1]; ++i) {
			irow = rowind[i];
			rowind[i] = perm_r[irow];
#if ( PRNTlevel>=2 )
			if ( rowind[i] == j ) { /* New diagonal */
			    if ( job == 2 || job == 3 )
				dmin = MIN(dmin, fabs(a[i]));
			    else if ( job == 4 )
				dsum += fabs(a[i]);
			    else if ( job == 5 )
				dprod *= fabs(a[i]);
			}
#endif
		    }
		}
	    }

#if ( PRNTlevel>=2 )
	    if ( job == 2 || job == 3 ) {
		if ( !iam ) printf("\tsmallest diagonal %e\n", dmin);
	    } else if ( job == 4 ) {
		if ( !iam ) printf("\tsum of diagonal %e\n", dsum);
	    } else if ( job == 5 ) {
		if ( !iam ) printf("\t product of diagonal %e\n", dprod);
	    }
#endif
	    
        } /* else !factored */

	t = SuperLU_timer_() - t;
	stat->utime[ROWPERM] = t;
#if ( PRNTlevel>=1 )
	if ( !iam ) printf(".. LDPERM job %d\t time: %.2f\n", job, t);
#endif
    
    } /* if options->RowPerm ... */

    if ( !factored || options->IterRefine ) {
	/* Compute norm(A), which will be used to adjust small diagonal. */
	if ( notran ) *(unsigned char *)norm = '1';
	else *(unsigned char *)norm = 'I';
	anorm = dlangs(norm, A);
    }

    /* ------------------------------------------------------------
       Perform the LU factorization.
       ------------------------------------------------------------*/
    if ( !factored ) {
	t = SuperLU_timer_();
	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = NATURAL:  natural ordering 
	 *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
	 *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
	 *   permc_spec = COLAMD:   approximate minimum degree column ordering
	 *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
	 */
	permc_spec = options->ColPerm;
	if ( permc_spec != MY_PERMC && Fact == DOFACT )
	    /* Use an ordering provided by SuperLU */
	    get_perm_c(iam, permc_spec, A, perm_c);

        /* Compute the elimination tree of Pc*(A'+A)*Pc' or Pc*A'*A*Pc'
         * (a.k.a. column etree), depending on the choice of ColPerm.
         * Adjust perm_c[] to be consistent with a postorder of etree.
         * Permute columns of A to form A*Pc'.
         */
        sp_colorder(options, A, perm_c, etree, &AC);

	/* Form Pc*A*Pc' to preserve the diagonal of the matrix Pr*A. */
        ACstore = (NCPformat *) AC.Store;
	for (j = 0; j < n; ++j) 
	    for (i = ACstore->colbeg[j]; i < ACstore->colend[j]; ++i) {
		irow = ACstore->rowind[i];
		ACstore->rowind[i] = perm_c[irow];
	    }
	stat->utime[COLPERM] = SuperLU_timer_() - t;

	/* Perform a symbolic factorization on matrix A and set up the
	   nonzero data structures which are suitable for supernodal GENP. */
	if ( Fact != SamePattern_SameRowPerm ) {
#if ( PRNTlevel>=1 ) 
	    if ( !iam ) 
		printf(".. symbfact(): relax %4d, maxsuper %4d, fill %4d\n",
		       sp_ienv(2), sp_ienv(3), sp_ienv(6));
#endif
	    t = SuperLU_timer_();
	    if ( !(Glu_freeable = (Glu_freeable_t *)
		   SUPERLU_MALLOC(sizeof(Glu_freeable_t))) )
		ABORT("Malloc fails for Glu_freeable.");

	    iinfo = symbfact(iam, &AC, perm_c, etree, 
			     Glu_persist, Glu_freeable);

	    stat->utime[SYMBFAC] = SuperLU_timer_() - t;

	    if ( !iam ) {
		if ( iinfo < 0 ) {
#if ( PRNTlevel>=1 )
		    printf("\tNo of supers %8d\n", Glu_persist->supno[n-1]+1);
		    printf("\tSize of G(L) %8d\n", Glu_freeable->xlsub[n]);
		    printf("\tSize of G(U) %8d\n", Glu_freeable->xusub[n]);
#endif
                    QuerySpace(n, -iinfo, Glu_freeable, &symb_mem_usage);
		    printf("\tSYMBfact:\tL\\U MB %.2f\ttotal MB %.2f\texpansions %d\n",
			   symb_mem_usage.for_lu*1e-6, 
			   symb_mem_usage.total*1e-6,
			   symb_mem_usage.expansions);
		} else {
		    fprintf(stderr, "symbfact() error returns %d\n", iinfo);
		    exit(-1);
		}
	    }
	}

	/* Distribute the L and U factors onto the process grid. */
	t = SuperLU_timer_();
	ddistribute(Fact, n, &AC, Glu_freeable, LUstruct, grid);
	stat->utime[DIST] = SuperLU_timer_() - t;

	/* Deallocate storage used in symbolic factor. */
	if ( Fact != SamePattern_SameRowPerm ) {
	    iinfo = symbfact_SubFree(Glu_freeable);
	    SUPERLU_FREE(Glu_freeable);
	}

	/* Perform numerical factorization in parallel. */
	t = SuperLU_timer_();
	pdgstrf(options, m, n, anorm, LUstruct, grid, stat, info);
	stat->utime[FACT] = SuperLU_timer_() - t;

#if ( PRNTlevel>=1 )
	{
          int_t TinyPivots;
          float for_lu, total;
          dQuerySpace(n, LUstruct, grid, &num_mem_usage);
          MPI_Reduce( &num_mem_usage.for_lu, &for_lu,
                     1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
          MPI_Reduce( &num_mem_usage.total, &total,
                     1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
          MPI_Allreduce( &stat->TinyPivots, &TinyPivots, 1, mpi_int_t,
                        MPI_SUM, grid->comm );
          stat->TinyPivots = TinyPivots;
          if ( !iam ) {
              printf("\tNUMfact space for all PEs:"
                     "\tL\\U MB %.2f\ttotal MB %.2f\n"
                     "\tNumber of tiny pivots: %10d\n",
                     for_lu*1e-6, total*1e-6, stat->TinyPivots);
          }
        }
#endif
    
#if ( PRNTlevel>=2 )
	if ( !iam ) printf(".. pdgstrf INFO = %d\n", *info);
#endif
 
    } else if ( options->IterRefine ) { /* options->Fact==FACTORED */
        /* Permute columns of A to form A*Pc' using the existing perm_c.
         * NOTE: rows of A were previously permuted to Pc*A.
         */
        sp_colorder(options, A, perm_c, NULL, &AC);
    } /* if !factored ... */
	
    /* ------------------------------------------------------------
       Compute the solution matrix X.
       ------------------------------------------------------------*/
    if ( nrhs ) {

	if ( !(b_work = doubleMalloc(n)) )
	    ABORT("Malloc fails for b_work[]");

	/* ------------------------------------------------------------
	   Scale the right-hand side if equilibration was performed. 
	   ------------------------------------------------------------*/
	if ( notran ) {
	    if ( rowequ ) {
		b_col = B;
		for (j = 0; j < nrhs; ++j) {
		    for (i = 0; i < m; ++i) b_col[i] *= R[i];
		    b_col += ldb;
		}
	    }
	} else if ( colequ ) {
	    b_col = B;
	    for (j = 0; j < nrhs; ++j) {
		for (i = 0; i < m; ++i) b_col[i] *= C[i];
		b_col += ldb;
	    }
	}

	/* ------------------------------------------------------------
	   Permute the right-hand side to form Pr*B.
	   ------------------------------------------------------------*/
	if ( options->RowPerm != NOROWPERM ) {
	    if ( notran ) {
		b_col = B;
		for (j = 0; j < nrhs; ++j) {
		    for (i = 0; i < m; ++i) b_work[perm_r[i]] = b_col[i];
		    for (i = 0; i < m; ++i) b_col[i] = b_work[i];
		    b_col += ldb;
		}
	    }
	}


	/* ------------------------------------------------------------
	   Permute the right-hand side to form Pc*B.
	   ------------------------------------------------------------*/
	if ( notran ) {
	    b_col = B;
	    for (j = 0; j < nrhs; ++j) {
		for (i = 0; i < m; ++i) b_work[perm_c[i]] = b_col[i];
		for (i = 0; i < m; ++i) b_col[i] = b_work[i];
		b_col += ldb;
	    }
	}


	/* Save a copy of the right-hand side. */
	ldx = ldb;
	if ( !(X = doubleMalloc(ldx * nrhs)) )
	    ABORT("Malloc fails for X[]");
	x_col = X;  b_col = B;
	for (j = 0; j < nrhs; ++j) {
	    for (i = 0; i < ldb; ++i) x_col[i] = b_col[i];
	    x_col += ldx;  b_col += ldb;
	}

	/* ------------------------------------------------------------
	   Solve the linear system.
	   ------------------------------------------------------------*/
	pdgstrs_Bglobal(n, LUstruct, grid, X, ldb, nrhs, stat, info);

	/* ------------------------------------------------------------
	   Use iterative refinement to improve the computed solution and
	   compute error bounds and backward error estimates for it.
	   ------------------------------------------------------------*/
	if ( options->IterRefine ) {
	    /* Improve the solution by iterative refinement. */
	    t = SuperLU_timer_();
	    pdgsrfs_ABXglobal(n, &AC, anorm, LUstruct, grid, B, ldb,
			      X, ldx, nrhs, berr, stat, info);
	    stat->utime[REFINE] = SuperLU_timer_() - t;
	}

	/* Permute the solution matrix X <= Pc'*X. */
	for (j = 0; j < nrhs; j++) {
	    b_col = &B[j*ldb];
	    x_col = &X[j*ldx];
	    for (i = 0; i < n; ++i) b_col[i] = x_col[perm_c[i]];
	}
	
	/* Transform the solution matrix X to a solution of the original system
	   before the equilibration. */
	if ( notran ) {
	    if ( colequ ) {
		b_col = B;
		for (j = 0; j < nrhs; ++j) {
		    for (i = 0; i < n; ++i) b_col[i] *= C[i];
		    b_col += ldb;
		}
	    }
	} else if ( rowequ ) {
	    b_col = B;
	    for (j = 0; j < nrhs; ++j) {
		for (i = 0; i < n; ++i) b_col[i] *= R[i];
		b_col += ldb;
	    }
	}

	SUPERLU_FREE(b_work);
	SUPERLU_FREE(X);

    } /* if nrhs != 0 */

#if ( PRNTlevel>=1 )
    if ( !iam ) printf(".. DiagScale = %d\n", ScalePermstruct->DiagScale);
#endif

    /* Deallocate storage. */
    if ( Equil && Fact != SamePattern_SameRowPerm ) {
	switch ( ScalePermstruct->DiagScale ) {
	    case NOEQUIL:
	        SUPERLU_FREE(R);
		SUPERLU_FREE(C);
		break;
	    case ROW: 
		SUPERLU_FREE(C);
		break;
	    case COL: 
		SUPERLU_FREE(R);
		break;
	}
    }
    if ( !factored || (factored && options->IterRefine) )
	Destroy_CompCol_Permuted(&AC);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit pdgssvx_ABglobal()");
#endif
}
