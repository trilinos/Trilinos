

#include <math.h>
#include "superlu_ddefs.h"

void
pdgssvx(superlu_options_t *options, SuperMatrix *A, 
	ScalePermstruct_t *ScalePermstruct,
	double B[], int ldb, int nrhs, gridinfo_t *grid,
	LUstruct_t *LUstruct, SOLVEstruct_t *SOLVEstruct, double *berr,
	SuperLUStat_t *stat, int *info)
{
/* 
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 *
 * Purpose
 * =======
 *
 * PDGSSVX solves a system of linear equations A*X=B,
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
 *      -  B, the matrix of right-hand sides, distributed by block rows,
 *            and its dimensions ldb (local) and nrhs (global)
 *      -  grid, a structure describing the 2D processor mesh
 *      -  options->IterRefine, which determines whether or not to
 *            improve the accuracy of the computed solution using 
 *            iterative refinement
 *
 *      On output, B is overwritten with the solution X.
 *
 *   2. Depending on options->Fact, the user has four options
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
 *        o  A, the input matrix
 *
 *        as well as the following options to determine what matrix to
 *        factorize.
 *
 *        o  options->Equil,   to specify how to scale the rows and columns
 *                             of A to "equilibrate" it (to try to reduce its
 *                             condition number and so improve the
 *                             accuracy of the computed solution)
 *
 *        o  options->RowPerm, to specify how to permute the rows of A
 *                             (typically to control numerical stability)
 *
 *        o  options->ColPerm, to specify how to permute the columns of A
 *                             (typically to control fill-in and enhance
 *                             parallelism during factorization)
 *
 *        o  options->ReplaceTinyPivot, to specify how to deal with tiny
 *                             pivots encountered during factorization
 *                             (to control numerical stability)
 *
 *      The outputs returned include
 *         
 *        o  ScalePermstruct,  modified to describe how the input matrix A
 *                             was equilibrated and permuted:
 *          .  ScalePermstruct->DiagScale, indicates whether the rows and/or
 *                                         columns of A were scaled
 *          .  ScalePermstruct->R, array of row scale factors
 *          .  ScalePermstruct->C, array of column scale factors
 *          .  ScalePermstruct->perm_r, row permutation vector
 *          .  ScalePermstruct->perm_c, column permutation vector
 *
 *          (part of ScalePermstruct may also need to be supplied on input,
 *           depending on options->RowPerm and options->ColPerm as described 
 *           later).
 *
 *        o  A, the input matrix A overwritten by the scaled and permuted
 *              matrix diag(R)*A*diag(C)*Pc^T, where 
 *              Pc is the row permutation matrix determined by
 *                  ScalePermstruct->perm_c
 *              diag(R) and diag(C) are diagonal scaling matrices determined
 *                  by ScalePermstruct->DiagScale, ScalePermstruct->R and 
 *                  ScalePermstruct->C
 *
 *        o  LUstruct, which contains the L and U factorization of A1 where
 *
 *                A1 = Pc*Pr*diag(R)*A*diag(C)*Pc^T = L*U
 *
 *               (Note that A1 = Pc*Pr*Aout, where Aout is the matrix stored
 *                in A on output.)
 *
 *   3. The second value of options->Fact assumes that a matrix with the same
 *      sparsity pattern as A has already been factored:
 *     
 *      -  options->Fact = SamePattern: A is factored, assuming that it has
 *            the same nonzero pattern as a previously factored matrix. In
 *            this case the algorithm saves time by reusing the previously
 *            computed column permutation vector stored in
 *            ScalePermstruct->perm_c and the "elimination tree" of A
 *            stored in LUstruct->etree
 *
 *      In this case the user must still specify the following options
 *      as before:
 *
 *        o  options->Equil
 *        o  options->RowPerm
 *        o  options->ReplaceTinyPivot
 *
 *      but not options->ColPerm, whose value is ignored. This is because the
 *      previous column permutation from ScalePermstruct->perm_c is used as
 *      input. The user must also supply 
 *
 *        o  A, the input matrix
 *        o  ScalePermstruct->perm_c, the column permutation
 *        o  LUstruct->etree, the elimination tree
 *
 *      The outputs returned include
 *         
 *        o  A, the input matrix A overwritten by the scaled and permuted
 *              matrix as described above
 *        o  ScalePermstruct, modified to describe how the input matrix A was
 *                            equilibrated and row permuted
 *        o  LUstruct, modified to contain the new L and U factors
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
 *        o  options->Equil
 *        o  options->ReplaceTinyPivot
 *
 *      but not options->RowPerm or options->ColPerm, whose values are
 *      ignored. This is because the permutations from ScalePermstruct->perm_r
 *      and ScalePermstruct->perm_c are used as input.
 *
 *      The user must also supply 
 *
 *        o  A, the input matrix
 *        o  ScalePermstruct->DiagScale, how the previous matrix was row
 *                                       and/or column scaled
 *        o  ScalePermstruct->R, the row scalings of the previous matrix,
 *                               if any
 *        o  ScalePermstruct->C, the columns scalings of the previous matrix, 
 *                               if any
 *        o  ScalePermstruct->perm_r, the row permutation of the previous
 *                                    matrix
 *        o  ScalePermstruct->perm_c, the column permutation of the previous 
 *                                    matrix
 *        o  all of LUstruct, the previously computed information about
 *                            L and U (the actual numerical values of L and U
 *                            stored in LUstruct->Llu are ignored)
 *
 *      The outputs returned include
 *         
 *        o  A, the input matrix A overwritten by the scaled and permuted
 *              matrix as described above
 *        o  ScalePermstruct,  modified to describe how the input matrix A was
 *                             equilibrated (thus ScalePermstruct->DiagScale,
 *                             R and C may be modified)
 *        o  LUstruct, modified to contain the new L and U factors
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
 *
 *      The user must also supply 
 *
 *        o  A, the unfactored matrix, only in the case that iterative
 *              refinment is to be done (specifically A must be the output
 *              A from the previous call, so that it has been scaled and permuted)
 *        o  all of ScalePermstruct
 *        o  all of LUstruct, including the actual numerical values of
 *           L and U
 *
 *      all of which are unmodified on output.
 *         
 * Arguments
 * =========
 *
 * options (input) superlu_options_t* (global)
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
 *         o Equil (yes_no_t)
 *           Specifies whether to equilibrate the system.
 *           = NO:  no equilibration.
 *           = YES: scaling factors are computed to equilibrate the system:
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
 *           = MMD_AT_PLUS_A: minimum degree ordering on structure of A'+A.
 *           = MMD_ATA:       minimum degree ordering on structure of A'*A.
 *           = COLAMD:        approximate minimum degree column ordering.
 *           = MY_PERMC:      the ordering given in ScalePermstruct->perm_c.
 *         
 *         o ReplaceTinyPivot (yes_no_t)
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
 * A (input/output) SuperMatrix* (local)
 *         On entry, matrix A in A*X=B, of dimension (A->nrow, A->ncol).
 *           The number of linear equations is A->nrow. The type of A must be:
 *           Stype = SLU_NR_loc; Dtype = SLU_D; Mtype = SLU_GE.
 *           That is, A is stored in distributed compressed row format.
 *           See supermatrix.h for the definition of 'SuperMatrix'.
 *           This routine only handles square A, however, the LU factorization
 *           routine PDGSTRF can factorize rectangular matrices.
 *         On exit, A may be overwtirren by diag(R)*A*diag(C)*Pc^T,
 *           depending on ScalePermstruct->DiagScale and options->ColPerm:
 *             if ScalePermstruct->DiagScale != NOEQUIL, A is overwritten by
 *                diag(R)*A*diag(C).
 *             if options->ColPerm != NATURAL, A is further overwritten by
 *                diag(R)*A*diag(C)*Pc^T.
 *           If all the above condition are true, the LU decomposition is
 *           performed on the matrix Pc*Pr*diag(R)*A*diag(C)*Pc^T.
 *
 * ScalePermstruct (input/output) ScalePermstruct_t* (global)
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
 * B       (input/output) double* (local)
 *         On entry, the right-hand side matrix of dimension (m_loc, nrhs),
 *           where, m_loc is the number of rows stored locally on my
 *           process and is defined in the data structure of matrix A.
 *         On exit, the solution matrix if info = 0;
 *
 * ldb     (input) int (local)
 *         The leading dimension of matrix B.
 *
 * nrhs    (input) int (global)
 *         The number of right-hand sides.
 *         If nrhs = 0, only LU decomposition is performed, the forward
 *         and back substitutions are skipped.
 *
 * grid    (input) gridinfo_t* (global)
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
 *         o etree (int*) dimension (A->ncol) (global)
 *           Elimination tree of Pc*(A'+A)*Pc' or Pc*A'*A*Pc'.
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
 *         o Glu_persist (Glu_persist_t*) (global)
 *           Global data structure (xsup, supno) replicated on all processes,
 *           describing the supernode partition in the factored matrices
 *           L and U:
 *	       xsup[s] is the leading column of the s-th supernode,
 *             supno[i] is the supernode number to which column i belongs.
 *
 *         o Llu (LocalLU_t*) (local)
 *           The distributed data structures to store L and U factors.
 *           See superlu_ddefs.h for the definition of 'LocalLU_t'.
 *
 * SOLVEstruct (input/output) SOLVEstruct_t*
 *         The data structure to hold the communication pattern used
 *         in the phases of triangular solution and iterative refinement.
 *         This pattern should be intialized only once for repeated solutions.
 *         If options->SolveInitialized = YES, it is an input argument.
 *         If options->SolveInitialized = NO and nrhs != 0, it is an output
 *         argument. See superlu_ddefs.h for the definition of 'SOLVEstruct_t'.
 *
 * berr    (output) double*, dimension (nrhs) (global)
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
 * See superlu_ddefs.h for the definitions of varioous data types.
 *
 */
    NRformat_loc *Astore;
    SuperMatrix GA;      /* Global A in NC format */
    NCformat *GAstore;
    double   *a_GA;
    SuperMatrix GAC;      /* Global A in NCP format (add n end pointers) */
    NCPformat *GACstore;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    Glu_freeable_t *Glu_freeable;
            /* The nonzero structures of L and U factors, which are
	       replicated on all processrs.
	           (lsub, xlsub) contains the compressed subscript of
		                 supernodes in L.
          	   (usub, xusub) contains the compressed subscript of
		                 nonzero segments in U.
	      If options->Fact != SamePattern_SameRowPerm, they are 
	      computed by SYMBFACT routine, and then used by PDDISTRIBUTE
	      routine. They will be freed after PDDISTRIBUTE routine.
	      If options->Fact == SamePattern_SameRowPerm, these
	      structures are not used.                                  */
    fact_t   Fact;
    double   *a;
    int_t    *colptr, *rowind;
    int_t    *perm_r; /* row permutations from partial pivoting */
    int_t    *perm_c; /* column permutation vector */
    int_t    *etree;  /* elimination tree */
    int_t    *rowptr, *colind;  /* Local A in NR*/
    int_t    *rowind_loc, *colptr_loc;
    int_t    colequ, Equil, factored, job, notran, rowequ, need_value;
    int_t    i, iinfo, j, irow, m, n, nnz, permc_spec, dist_mem_use;
    int_t    nnz_loc, m_loc, fst_row, icol;
    int      iam;
    int      ldx;  /* LDA for matrix X (local). */
    char     equed[1], norm[1];
    double   *C, *R, *C1, *R1, amax, anorm, colcnd, rowcnd;
    double   *X, *b_col, *b_work, *x_col;
    double   t;
    static mem_usage_t num_mem_usage, symb_mem_usage;
#if ( PRNTlevel>= 2 )
    double   dmin, dsum, dprod;
#endif
    int_t procs;

    /* Initialization. */
    m = A->nrow;
    n = A->ncol;
    Astore = (NRformat_loc *) A->Store;
    nnz_loc = Astore->nnz_loc;
    m_loc = Astore->m_loc;
    fst_row = Astore->fst_row;
    a = Astore->nzval;
    rowptr = Astore->rowptr;
    colind = Astore->colind;

    /* Test the input parameters. */
    *info = 0;
    Fact = options->Fact;
    if ( Fact < 0 || Fact > FACTORED )
	*info = -1;
    else if ( options->RowPerm < 0 || options->RowPerm > MY_PERMR )
	*info = -1;
    else if ( options->ColPerm < 0 || options->ColPerm > MY_PERMC )
	*info = -1;
    else if ( options->IterRefine < 0 || options->IterRefine > EXTRA )
	*info = -1;
    else if ( options->IterRefine == EXTRA ) {
	*info = -1;
	fprintf(stderr, "Extra precise iterative refinement yet to support.");
    } else if ( A->nrow != A->ncol || A->nrow < 0 || A->Stype != SLU_NR_loc
		|| A->Dtype != SLU_D || A->Mtype != SLU_GE )
	*info = -2;
    else if ( ldb < m_loc )
	*info = -5;
    else if ( nrhs < 0 )
	*info = -6;
    if ( *info ) {
	i = -(*info);
	pxerbla("pdgssvx", grid, -*info);
	return;
    }

    factored = (Fact == FACTORED);
    Equil = (!factored && options->Equil == YES);
    notran = (options->Trans == NOTRANS);
    iam = grid->iam;
    job = 5;
    if ( factored || (Fact == SamePattern_SameRowPerm && Equil) ) {
	rowequ = (ScalePermstruct->DiagScale == ROW) ||
	         (ScalePermstruct->DiagScale == BOTH);
	colequ = (ScalePermstruct->DiagScale == COL) ||
	         (ScalePermstruct->DiagScale == BOTH);
    } else rowequ = colequ = FALSE;

    /* The following arrays are replicated on all processes. */
    perm_r = ScalePermstruct->perm_r;
    perm_c = ScalePermstruct->perm_c;
    etree = LUstruct->etree;
    R = ScalePermstruct->R;
    C = ScalePermstruct->C;
    /********/

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter pdgssvx()");
#endif

    if ( Equil ) { /* Not factored & ask for equilibration */
	/* Allocate storage if not done so before. */
	switch ( ScalePermstruct->DiagScale ) {
	    case NOEQUIL:
		if ( !(R = (double *) doubleMalloc_dist(m)) )
		    ABORT("Malloc fails for R[].");
	        if ( !(C = (double *) doubleMalloc_dist(n)) )
		    ABORT("Malloc fails for C[].");
		ScalePermstruct->R = R;
		ScalePermstruct->C = C;
		break;
	    case ROW: 
	        if ( !(C = (double *) doubleMalloc_dist(n)) )
		    ABORT("Malloc fails for C[].");
		ScalePermstruct->C = C;
		break;
	    case COL: 
		if ( !(R = (double *) doubleMalloc_dist(m)) )
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
		irow = fst_row;
		for (j = 0; j < m_loc; ++j) {
		    for (i = rowptr[j]; i < rowptr[j+1]; ++i) {
			a[i] *= R[irow];       /* Scale rows. */
		    }
		    ++irow;
		}
		break;
	      case COL:
		for (j = 0; j < m_loc; ++j)
		    for (i = rowptr[j]; i < rowptr[j+1]; ++i){
		        icol = colind[i];
			a[i] *= C[icol];          /* Scale columns. */
		    }
		break;
	      case BOTH:
		irow = fst_row;
		for (j = 0; j < m_loc; ++j) {
		    for (i = rowptr[j]; i < rowptr[j+1]; ++i) {
			icol = colind[i];
			a[i] *= R[irow] * C[icol]; /* Scale rows and cols. */
		    }
		    ++irow;
		}
	        break;
	    }
	} else { /* Compute R & C from scratch */
            /* Compute the row and column scalings. */
	    pdgsequ(A, R, C, &rowcnd, &colcnd, &amax, &iinfo, grid);

	    /* Equilibrate matrix A if it is badly-scaled. */
	    pdlaqgs(A, R, C, rowcnd, colcnd, amax, equed);

	    if ( lsame_(equed, "R") ) {
		ScalePermstruct->DiagScale = rowequ = ROW;
	    } else if ( lsame_(equed, "C") ) {
		ScalePermstruct->DiagScale = colequ = COL;
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

    if ( !factored ) { /* Skip this if already factored. */
        /*
         * Gather A from the distributed compressed row format to
         * global A in compressed column format.
         * Numerical values are gathered only when a row permutation
         * for large diagonal is sought after.
         */
	if ( Fact != SamePattern_SameRowPerm ) {
            need_value = (options->RowPerm == LargeDiag);
            pdCompRow_loc_to_CompCol_global(need_value, A, grid, &GA);
            GAstore = (NCformat *) GA.Store;
            colptr = GAstore->colptr;
            rowind = GAstore->rowind;
            nnz = GAstore->nnz;
            if ( need_value ) a_GA = GAstore->nzval;
            else assert(GAstore->nzval == NULL);
	}

        /* ------------------------------------------------------------
           Find the row permutation for A.
           ------------------------------------------------------------*/
        if ( options->RowPerm != NO ) {
	    t = SuperLU_timer_();
	    if ( Fact != SamePattern_SameRowPerm ) {
	        if ( options->RowPerm == MY_PERMR ) { /* Use user's perm_r. */
	            /* Permute the global matrix GA for symbfact() */
	            for (i = 0; i < colptr[n]; ++i) {
	            	irow = rowind[i]; 
		    	rowind[i] = perm_r[irow];
	            }
	        } else { /* options->RowPerm == LargeDiag */
	            /* Get a new perm_r[] */
	            if ( job == 5 ) {
		        /* Allocate storage for scaling factors. */
		        if ( !(R1 = doubleMalloc_dist(m)) )
		            ABORT("SUPERLU_MALLOC fails for R1[]");
		    	if ( !(C1 = doubleMalloc_dist(n)) )
		            ABORT("SUPERLU_MALLOC fails for C1[]");
	            }

	            if ( !iam ) {
		        /* Process 0 finds a row permutation */
		        dldperm(job, m, nnz, colptr, rowind, a_GA,
		                perm_r, R1, C1);
		
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

		            /* Scale the distributed matrix */
		            irow = fst_row;
		            for (j = 0; j < m_loc; ++j) {
			        for (i = rowptr[j]; i < rowptr[j+1]; ++i) {
			            icol = colind[i];
			            a[i] *= R1[irow] * C1[icol];
			        }
			        ++irow;
		            }

		            /* Multiply together the scaling factors. */
		            if ( rowequ ) for (i = 0; i < m; ++i) R[i] *= R1[i];
		            else for (i = 0; i < m; ++i) R[i] = R1[i];
		            if ( colequ ) for (i = 0; i < n; ++i) C[i] *= C1[i];
		            else for (i = 0; i < n; ++i) C[i] = C1[i];
		    
		            ScalePermstruct->DiagScale = BOTH;
		            rowequ = colequ = 1;

		        } /* end Equil */

                        /* Now permute global A to prepare for symbfact() */
                        for (j = 0; j < n; ++j) {
		            for (i = colptr[j]; i < colptr[j+1]; ++i) {
	                        irow = rowind[i];
		                rowind[i] = perm_r[irow];
#if ( PRNTlevel>=2 )
			        if ( rowind[i] == j ) /* New diagonal */
			            dprod *= fabs(a[i]);
#endif
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
				        dmin = SUPERLU_MIN(dmin, fabs(a[i]));
			            else if ( job == 4 )
				        dsum += fabs(a[i]);
			            else if ( job == 5 )
				        dprod *= fabs(a[i]);
			        }
#endif
		            } /* end for i ... */
		        } /* end for j ... */
	            } /* end else job ... */

#if ( PRNTlevel>=2 )
	            if ( job == 2 || job == 3 ) {
		        if ( !iam ) printf("\tsmallest diagonal %e\n", dmin);
	            } else if ( job == 4 ) {
		        if ( !iam ) printf("\tsum of diagonal %e\n", dsum);
	            } else if ( job == 5 ) {
		        if ( !iam ) printf("\t product of diagonal %e\n", dprod);
	            }
#endif
	    
                } /* end if options->RowPerm ... */

	        t = SuperLU_timer_() - t;
	        stat->utime[ROWPERM] = t;
#if ( PRNTlevel>=1 )
	        if ( !iam ) printf(".. LDPERM job %d\t time: %.2f\n", job, t);
#endif
            } /* end if Fact ... */
        } else { /* options->RowPerm == NOROWPERM */
            for (i = 0; i <m; ++i) perm_r[i] = i;
        }

#if ( DEBUGlevel>=1 )
        if ( !iam ) PrintInt10("perm_r",  m, perm_r);
#endif
    } /* end if (!factored) */

    if ( !factored || options->IterRefine ) {
	/* Compute norm(A), which will be used to adjust small diagonal. */
	if ( notran ) *(unsigned char *)norm = '1';
	else *(unsigned char *)norm = 'I';
	anorm = pdlangs(norm, A, grid);
#if ( PRNTlevel>=1 )
	if ( !iam ) printf(".. anorm %e\n", anorm);
#endif
    }

    /* ------------------------------------------------------------
       Perform the LU factorization.
       ------------------------------------------------------------*/
    if ( !factored ) {
	t = SuperLU_timer_();
	/*
	 * Get column permutation vector perm_c[], according to permc_spec:
	 *   permc_spec = NATURAL:  natural ordering 
	 *   permc_spec = MMD_AT_PLUS_A: minimum degree on structure of A'+A
	 *   permc_spec = MMD_ATA:  minimum degree on structure of A'*A
	 *   permc_spec = COLAMD:   approximate minimum degree column ordering
	 *   permc_spec = MY_PERMC: the ordering already supplied in perm_c[]
	 */
	permc_spec = options->ColPerm;
	if ( permc_spec != MY_PERMC && Fact == DOFACT )
	    get_perm_c_dist(iam, permc_spec, &GA, perm_c);

	stat->utime[COLPERM] = SuperLU_timer_() - t;

	/* Compute the elimination tree of Pc*(A'+A)*Pc' or Pc*A'*A*Pc'
	   (a.k.a. column etree), depending on the choice of ColPerm.
	   Adjust perm_c[] to be consistent with a postorder of etree.
	   Permute columns of A to form A*Pc'. */
	if ( Fact != SamePattern_SameRowPerm ) {
	    int_t *GACcolbeg, *GACcolend, *GACrowind;

	    sp_colorder(options, &GA, perm_c, etree, &GAC); 

	    /* Form Pc*A*Pc' to preserve the diagonal of the matrix GAC. */
	    GACstore = GAC.Store;
	    GACcolbeg = GACstore->colbeg;
	    GACcolend = GACstore->colend;
	    GACrowind = GACstore->rowind;
	    for (j = 0; j < n; ++j) {
	        for (i = GACcolbeg[j]; i < GACcolend[j]; ++i) {
		    irow = GACrowind[i];
		    GACrowind[i] = perm_c[irow];
	        }
	    }

	    /* Perform a symbolic factorization on Pc*Pr*A*Pc' and set up the
	       nonzero data structures for L & U. */
#if ( PRNTlevel>=1 ) 
            if ( !iam ) 
		printf(".. symbfact(): relax %4d, maxsuper %4d, fill %4d\n",
		       sp_ienv_dist(2), sp_ienv_dist(3), sp_ienv_dist(6));
#endif
	    t = SuperLU_timer_();
	    if ( !(Glu_freeable = (Glu_freeable_t *)
		   SUPERLU_MALLOC(sizeof(Glu_freeable_t))) )
		ABORT("Malloc fails for Glu_freeable.");

	    /* Every process does this. */
	    iinfo = symbfact(iam, &GAC, perm_c, etree, 
			     Glu_persist, Glu_freeable);

	    stat->utime[SYMBFAC] = SuperLU_timer_() - t;
	    if ( iinfo < 0 ) { /* Successful return */
		QuerySpace_dist(n, -iinfo, Glu_freeable, &symb_mem_usage);
#if ( PRNTlevel>=1 )
		if ( !iam ) {
		    printf("\tNo of supers %ld\n", Glu_persist->supno[n-1]+1);
		    printf("\tSize of G(L) %ld\n", Glu_freeable->xlsub[n]);
		    printf("\tSize of G(U) %ld\n", Glu_freeable->xusub[n]);
		    printf("\tint %d, short %d, float %d, double %d\n", 
			   sizeof(int_t), sizeof(short), sizeof(float),
			   sizeof(double));
		    printf("\tSYMBfact (MB):\tL\\U %.2f\ttotal %.2f\texpansions %d\n",
			   symb_mem_usage.for_lu*1e-6, 
			   symb_mem_usage.total*1e-6,
			   symb_mem_usage.expansions);
		}
#endif
	    } else {
		if ( !iam ) {
		    fprintf(stderr, "symbfact() error returns %d\n", iinfo);
		    exit(-1);
		}
	    }
	} /* end if Fact ... */

	/* Apply column permutation to the original distributed A */
	for (j = 0; j < nnz_loc; ++j) colind[j] = perm_c[colind[j]];

	/* Distribute Pc*Pr*diag(R)*A*diag(C)*Pc' into L and U storage. 
	   NOTE: the row permutation Pc*Pr is applied internally in the
	   distribution routine. */
	t = SuperLU_timer_();
	dist_mem_use = pddistribute(Fact, n, A, ScalePermstruct,
                                  Glu_freeable, LUstruct, grid);
	stat->utime[DIST] = SuperLU_timer_() - t;

	/* Deallocate storage used in symbolic factorization. */
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
	    float for_lu, total, max, avg, temp;
	    dQuerySpace_dist(n, LUstruct, grid, &num_mem_usage);
	    MPI_Reduce( &num_mem_usage.for_lu, &for_lu,
		       1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
	    MPI_Reduce( &num_mem_usage.total, &total,
		       1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
	    temp = SUPERLU_MAX(symb_mem_usage.total,
			       symb_mem_usage.for_lu +
			       (float)dist_mem_use + num_mem_usage.for_lu);
	    temp = SUPERLU_MAX(temp, num_mem_usage.total);
	    MPI_Reduce( &temp, &max,
		       1, MPI_FLOAT, MPI_MAX, 0, grid->comm );
	    MPI_Reduce( &temp, &avg,
		       1, MPI_FLOAT, MPI_SUM, 0, grid->comm );
	    MPI_Allreduce( &stat->TinyPivots, &TinyPivots, 1, mpi_int_t,
			  MPI_SUM, grid->comm );
	    stat->TinyPivots = TinyPivots;
	    if ( !iam ) {
		printf("\tNUMfact (MB) all PEs:\tL\\U\t%.2f\tall\t%.2f\n",
		       for_lu*1e-6, total*1e-6);
		printf("\tAll space (MB):"
		       "\t\ttotal\t%.2f\tAvg\t%.2f\tMax\t%.2f\n",
		       avg*1e-6, avg/grid->nprow/grid->npcol*1e-6, max*1e-6);
		printf("\tNumber of tiny pivots: %10d\n", stat->TinyPivots);
	    }
	}
#endif
    
        /* Destroy GA */
        if ( Fact != SamePattern_SameRowPerm )
            Destroy_CompCol_Matrix_dist(&GA);
    } /* end if (!factored) */
	
    /* ------------------------------------------------------------
       Compute the solution matrix X.
       ------------------------------------------------------------*/
    if ( nrhs ) {

	if ( !(b_work = doubleMalloc_dist(n)) )
	    ABORT("Malloc fails for b_work[]");

	/* ------------------------------------------------------------
	   Scale the right-hand side if equilibration was performed. 
	   ------------------------------------------------------------*/
	if ( notran ) {
	    if ( rowequ ) {
		b_col = B;
		for (j = 0; j < nrhs; ++j) {
		    irow = fst_row;
		    for (i = 0; i < m_loc; ++i) {
		        b_col[i] *= R[irow];
		        ++irow;
		    }
		    b_col += ldb;
		}
	    }
	} else if ( colequ ) {
	    b_col = B;
	    for (j = 0; j < nrhs; ++j) {
	        irow = fst_row;
		for (i = 0; i < m_loc; ++i) {
		    b_col[i] *= C[irow];
		    ++irow;
		}
		b_col += ldb;
	    }
	}

	/* Save a copy of the right-hand side. */
	ldx = ldb;
	if ( !(X = doubleMalloc_dist(((size_t)ldx) * nrhs)) )
	    ABORT("Malloc fails for X[]");
	x_col = X;  b_col = B;
	for (j = 0; j < nrhs; ++j) {
	    for (i = 0; i < m_loc; ++i) x_col[i] = b_col[i];
	    x_col += ldx;  b_col += ldb;
	}

	/* ------------------------------------------------------------
	   Solve the linear system.
	   ------------------------------------------------------------*/
	if ( options->SolveInitialized == NO ) {
	    dSolveInit(options, A, perm_r, perm_c, nrhs, LUstruct, grid,
		       SOLVEstruct);
	}

	pdgstrs(n, LUstruct, ScalePermstruct, grid, X, m_loc, 
		fst_row, ldb, nrhs, SOLVEstruct, stat, info);

	/* ------------------------------------------------------------
	   Use iterative refinement to improve the computed solution and
	   compute error bounds and backward error estimates for it.
	   ------------------------------------------------------------*/
	if ( options->IterRefine ) {
	    /* Improve the solution by iterative refinement. */
	    int_t *it, *colind_gsmv = SOLVEstruct->A_colind_gsmv;

	    t = SuperLU_timer_();
	    if ( options->RefineInitialized == NO || Fact == DOFACT ) {
	        /* All these cases need to re-initialize gsmv structure */
	        if ( options->RefineInitialized )
		    pdgsmv_finalize(SOLVEstruct->gsmv_comm);
	        pdgsmv_init(A, SOLVEstruct->row_to_proc, grid,
			    SOLVEstruct->gsmv_comm);
	       
                /* Save a copy of the transformed local col indices. */
	        if ( colind_gsmv ) SUPERLU_FREE(colind_gsmv);
	        if ( !(it = intMalloc_dist(nnz_loc)) )
		    ABORT("Malloc fails for colind_gsmv[]");
	        colind_gsmv = SOLVEstruct->A_colind_gsmv = it;
	        for (i = 0; i < nnz_loc; ++i) colind_gsmv[i] = colind[i];
	        options->RefineInitialized = YES;
	    } else if ( Fact == SamePattern ||
			Fact == SamePattern_SameRowPerm ) {
	        double at;
	        int_t k, jcol, p;
	        /* Swap to beginning the part of A corresponding to the
		   local part of X, as was done in pdgsmv_init() */
	        for (i = 0; i < m_loc; ++i) { /* Loop through each row */
		    k = rowptr[i];
		    for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
		        jcol = colind[j];
		        p = SOLVEstruct->row_to_proc[jcol];
		        if ( p == iam ) { /* Local */
		            at = a[k]; a[k] = a[j]; a[j] = at;
		            ++k;
		        }
		    }
	        }
	      
	        /* Re-use the local col indices of A obtained from the
		   previous call to pdgsmv_init() */
	        for (i = 0; i < nnz_loc; ++i) colind[i] = colind_gsmv[i];
	    }

	    pdgsrfs(n, A, anorm, LUstruct, ScalePermstruct, grid,
		    B, ldb, X, ldx, nrhs, SOLVEstruct, berr, stat, info);
	    stat->utime[REFINE] = SuperLU_timer_() - t;
	}

	/* Permute the solution matrix B <= Pc'*X. */
	pdPermute_Dense_Matrix(fst_row, m_loc, SOLVEstruct->row_to_proc,
			       SOLVEstruct->inv_perm_c,
			       X, ldx, B, ldb, nrhs, grid);
#if ( DEBUGlevel>=2 )
	printf("\n (%d) .. After pdPermute_Dense_Matrix(): b =\n", iam);
	for (i = 0; i < m_loc; ++i)
	  printf("\t(%d)\t%4d\t%.10f\n", iam, i+fst_row, B[i]);
#endif
	
	/* Transform the solution matrix X to a solution of the original
	   system before the equilibration. */
	if ( notran ) {
	    if ( colequ ) {
		b_col = B;
		for (j = 0; j < nrhs; ++j) {
		    irow = fst_row;
		    for (i = 0; i < m_loc; ++i) {
		        b_col[i] *= C[irow];
		        ++irow;
		    }
		    b_col += ldb;
		}
	    }
	} else if ( rowequ ) {
	    b_col = B;
	    for (j = 0; j < nrhs; ++j) {
	        irow = fst_row;
		for (i = 0; i < m_loc; ++i) {
		    b_col[i] *= R[irow];
		    ++irow;
		}
		b_col += ldb;
	    }
	}

	SUPERLU_FREE(b_work);
	SUPERLU_FREE(X);

    } /* end if nrhs != 0 */

#if ( PRNTlevel>=1 )
    if ( !iam ) printf(".. DiagScale = %d\n", ScalePermstruct->DiagScale);
#endif

    /* Deallocate R and/or C if it was not used. */
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
    if ( !factored && Fact != SamePattern_SameRowPerm )
 	Destroy_CompCol_Permuted_dist(&GAC);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit pdgssvx()");
#endif

}
