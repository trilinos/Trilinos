/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include <math.h>
#include "superlu_ddefs.h"
#if ( VAMPIR>=1 )
#include <VT.h>
#endif

/*
 * Internal prototypes
 */
static void pdgstrf2(superlu_options_t *, int_t, double, Glu_persist_t *,
		     gridinfo_t *, LocalLU_t *, SuperLUStat_t *, int *);
#ifdef _CRAY
static void pdgstrs2(int_t, int_t, Glu_persist_t *, gridinfo_t *,
		     LocalLU_t *, SuperLUStat_t *, _fcd, _fcd, _fcd);
#else
static void pdgstrs2(int_t, int_t, Glu_persist_t *, gridinfo_t *,
		     LocalLU_t *, SuperLUStat_t *);
#endif

/* 
 * Sketch of the algorithm
 * =======================
 *
 * The following relations hold:
 *     * A_kk = L_kk * U_kk
 *     * L_ik = Aik * U_kk^(-1)
 *     * U_kj = L_kk^(-1) * A_kj
 *
 *              ----------------------------------
 *              |   |                            |
 *              ----|-----------------------------
 *              |   | \ U_kk|                    |
 *              |   |   \   |        U_kj        |
 *              |   |L_kk \ |         ||         |
 *              ----|-------|---------||----------
 *              |   |       |         \/         |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   | L_ik ==>       A_ij        |
 *              |   |       |                    |
 *              |   |       |                    |
 *              |   |       |                    |
 *              ----------------------------------
 *
 * Handle the first block of columns separately.
 *     * Factor diagonal and subdiagonal blocks and test for exact
 *       singularity. ( pdgstrf2(0), one column at a time )
 *     * Compute block row of U
 *     * Update trailing matrix
 * 
 * Loop over the remaining blocks of columns.
 *   mycol = MYCOL( iam, grid );
 *   myrow = MYROW( iam, grid );
 *   N = nsupers;
 *   For (k = 1; k < N; ++k) {
 *       krow = PROW( k, grid );
 *       kcol = PCOL( k, grid );
 *       Pkk = PNUM( krow, kcol, grid );
 *
 *     * Factor diagonal and subdiagonal blocks and test for exact
 *       singularity.
 *       if ( mycol == kcol ) {
 *           pdgstrf2(k), one column at a time 
 *       }
 *
 *     * Parallel triangular solve
 *       if ( iam == Pkk ) multicast L_k,k to this process row;
 *       if ( myrow == krow && mycol != kcol ) {
 *          Recv L_k,k from process Pkk;
 *          for (j = k+1; j < N; ++j) 
 *              if ( PCOL( j, grid ) == mycol && A_k,j != 0 )
 *                 U_k,j = L_k,k \ A_k,j;
 *       }
 *
 *     * Parallel rank-k update
 *       if ( myrow == krow ) multicast U_k,k+1:N to this process column;
 *       if ( mycol == kcol ) multicast L_k+1:N,k to this process row;
 *       if ( myrow != krow ) {
 *          Pkj = PNUM( krow, mycol, grid );
 *          Recv U_k,k+1:N from process Pkj;
 *       }
 *       if ( mycol != kcol ) {
 *          Pik = PNUM( myrow, kcol, grid );
 *          Recv L_k+1:N,k from process Pik;
 *       }
 *       for (j = k+1; k < N; ++k) {
 *          for (i = k+1; i < N; ++i) 
 *              if ( myrow == PROW( i, grid ) && mycol == PCOL( j, grid )
 *                   && L_i,k != 0 && U_k,j != 0 )
 *                 A_i,j = A_i,j - L_i,k * U_k,j;
 *       }
 *  }
 *
 *
 * Remaining issues
 *   (1) Use local indices for L subscripts and SPA.  [DONE]
 *
 */
/************************************************************************/
void pdgstrf
/************************************************************************/
(
 superlu_options_t *options, int m, int n, double anorm,
 LUstruct_t *LUstruct, gridinfo_t *grid, SuperLUStat_t *stat, int *info
 )
/* 
 * Purpose
 * =======
 *
 *  PDGSTRF performs the LU factorization in parallel.
 *
 * Arguments
 * =========
 * 
 * options (input) superlu_options_t*
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed.
 *         The following field should be defined:
 *         o ReplaceTinyPivot (yes_no_t)
 *           Specifies whether to replace the tiny diagonals by
 *           sqrt(epsilon)*norm(A) during LU factorization.
 *
 * m      (input) int
 *        Number of rows in the matrix.
 *
 * n      (input) int
 *        Number of columns in the matrix.
 *
 * anorm  (input) double
 *        The norm of the original matrix A, or the scaled A if
 *        equilibration was done.
 *
 * LUstruct (input/output) LUstruct_t*
 *         The data structures to store the distributed L and U factors.
 *         The following fields should be defined:
 *
 *         o Glu_persist (input) Glu_persist_t*
 *           Global data structure (xsup, supno) replicated on all processes,
 *           describing the supernode partition in the factored matrices
 *           L and U:
 *	       xsup[s] is the leading column of the s-th supernode,
 *             supno[i] is the supernode number to which column i belongs.
 *
 *         o Llu (input/output) LocalLU_t*
 *           The distributed data structures to store L and U factors.
 *           See superlu_ddefs.h for the definition of 'LocalLU_t'.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh. It contains the MPI communicator, the number
 *        of process rows (NPROW), the number of process columns (NPCOL),
 *        and my process rank. It is an input argument to all the
 *        parallel routines.
 *        Grid can be initialized by subroutine SUPERLU_GRIDINIT.
 *        See superlu_ddefs.h for the definition of 'gridinfo_t'.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics on runtime and floating-point operation count.
 *        See util.h for the definition of 'SuperLUStat_t'.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, U(i,i) is exactly zero. The factorization has
 *             been completed, but the factor U is exactly singular,
 *             and division by zero will occur if it is used to solve a
 *             system of equations.
 *
 */
{
#ifdef _CRAY
    _fcd ftcs = _cptofcd("N", strlen("N"));
    _fcd ftcs1 = _cptofcd("L", strlen("L"));
    _fcd ftcs2 = _cptofcd("N", strlen("N"));
    _fcd ftcs3 = _cptofcd("U", strlen("U"));
#endif
    double alpha = 1.0, beta = 0.0;
    int_t *xsup;
    int_t *lsub, *lsub1, *usub, *Lsub_buf, *Usub_buf;
    int_t fnz, i, ib, ijb, ilst, it, iukp, jb, jj, klst, knsupc,
          lb, lib, ldv, ljb, lptr, lptr0, lptrj, luptr, luptr0, luptrj,
          nlb, nub, nsupc, rel, rukp;
    int_t Pc, Pr;
    int   iam, kcol, krow, mycol, myrow, pi, pj;
    int   j, k, lk, nsupers;
    int   nsupr, nbrow, segsize;
    int   msgcnt[4]; /* Count the size of the message xfer'd in each buffer:
		      *     0 : transferred in Lsub_buf[]
		      *     1 : transferred in Lval_buf[]
		      *     2 : transferred in Usub_buf[] 
		      *     3 : transferred in Uval_buf[]
		      */
    int_t  msg0, msg2;
    double *lusup, *lusup1, *uval, *Lval_buf, *Uval_buf;
    int_t  **Ufstnz_br_ptr, **Lrowind_bc_ptr;
    double **Unzval_br_ptr, **Lnzval_bc_ptr;
    int_t  *index;
    double *nzval;
    int_t  *iuip, *ruip;/* Pointers to U index/nzval; size ceil(NSUPERS/Pr). */
    double *ucol;
    int_t  *indirect;
    double *tempv, *tempv2d;
    int_t iinfo;
    int_t *ToRecv, *ToSendD, **ToSendR;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    LocalLU_t *Llu = LUstruct->Llu;
    superlu_scope_t *scp;
    MPI_Status status;
    double s_eps, thresh;
    double *tempU2d, *tempu;
    int    full, ldt, ldu, lead_zero, ncols;
#if ( DEBUGlevel>=1 ) 
    int_t num_copy=0, num_update=0;
#endif
    extern double slamch_(char *);
#if ( PRNTlevel==3 )
    int_t  zero_msg = 0, total_msg = 0;
#endif
#if ( PROFlevel>=1 )
    double t1, t2;
#endif

    /* Test the input parameters. */
    *info = 0;
    if ( m < 0 ) *info = -2;
    else if ( n < 0 ) *info = -3;
    if ( *info ) {
	pxerbla("pdgstrf", grid, -*info);
	return;
    }

    /* Quick return if possible. */
    if ( m == 0 || n == 0 ) return;

    /*
     * Initialization.
     */
    iam = grid->iam;
    Pc = grid->npcol;
    Pr = grid->nprow;
    myrow = MYROW( iam, grid );
    mycol = MYCOL( iam, grid );
    nsupers = Glu_persist->supno[n-1] + 1;
    xsup = Glu_persist->xsup;
    s_eps = slamch_("Epsilon");
    thresh = s_eps * anorm;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter pdgstrf()");
#endif

    if ( !(Llu->Lsub_buf = intMalloc(Llu->bufmax[0])) )
	ABORT("Malloc fails for Lsub_buf[].");
    if ( !(Llu->Lval_buf = doubleMalloc(Llu->bufmax[1])) )
	ABORT("Malloc fails for Lval_buf[].");
    if ( !(Llu->Usub_buf = intMalloc(Llu->bufmax[2])) )
	ABORT("Malloc fails for Usub_buf[].");
    if ( !(Llu->Uval_buf = doubleMalloc(Llu->bufmax[3])) )
	ABORT("Malloc fails for Uval_buf[].");
    if ( !(Llu->ujrow = doubleMalloc(sp_ienv(3))) )
	ABORT("Malloc fails for ujrow[].");

#if ( PRNTlevel>=1 )
    if ( !iam ) {
	printf(".. thresh = s_eps %e * anorm %e = %e\n", s_eps, anorm, thresh);
	printf(".. Buffer size: Lsub %d\tLval %d\tUsub %d\tUval %d\tLDA %d\n",
	       Llu->bufmax[0], Llu->bufmax[1], 
	       Llu->bufmax[2], Llu->bufmax[3], Llu->bufmax[4]);
    }
#endif

    Lval_buf = Llu->Lval_buf;
    Lsub_buf = Llu->Lsub_buf;
    Uval_buf = Llu->Uval_buf;
    Usub_buf = Llu->Usub_buf;
    Lrowind_bc_ptr = Llu->Lrowind_bc_ptr;
    Lnzval_bc_ptr = Llu->Lnzval_bc_ptr;
    Ufstnz_br_ptr = Llu->Ufstnz_br_ptr;
    Unzval_br_ptr = Llu->Unzval_br_ptr;
    ToRecv = Llu->ToRecv;
    ToSendD = Llu->ToSendD;
    ToSendR = Llu->ToSendR;

    ldt = sp_ienv(3); /* Size of maximum supernode */
    if ( !(tempv2d = doubleCalloc(2*ldt*ldt)) )
	ABORT("Calloc fails for tempv2d[].");
    tempU2d = tempv2d + ldt*ldt;
    if ( !(indirect = intMalloc(ldt)) )
	ABORT("Malloc fails for indirect[].");
    k = CEILING( nsupers, Pr ); /* Number of local block rows */
    if ( !(iuip = intMalloc(k)) )
	ABORT("Malloc fails for iuip[].");
    if ( !(ruip = intMalloc(k)) )
	ABORT("Malloc fails for ruip[].");

#if ( VAMPIR>=1 )
    VT_symdef(1, "Send-L", "Comm");
    VT_symdef(2, "Recv-L", "Comm");
    VT_symdef(3, "Send-U", "Comm");
    VT_symdef(4, "Recv-U", "Comm");
    VT_symdef(5, "TRF2", "Factor");
    VT_symdef(100, "Factor", "Factor");
    VT_traceon();
    VT_begin(100);
#endif

    /* Handle the first block column separately. */
    if ( mycol == 0 ) {
#if ( VAMPIR>=1 )
	VT_begin(5);
#endif
	pdgstrf2(options, 0, thresh, Glu_persist, grid, Llu, stat, info);
#if ( VAMPIR>=1 )
	VT_end(5);
#endif

	scp = &grid->rscp; /* The scope of process row. */

	/* Process column *kcol* multicasts numeric values of L(:,k) 
	   to process rows. */
	lsub = Lrowind_bc_ptr[0];
	lusup = Lnzval_bc_ptr[0];
	if ( lsub ) {
	    msgcnt[0] = lsub[1] + BC_HEADER + lsub[0]*LB_DESCRIPTOR;
	    msgcnt[1] = lsub[1] * SuperSize( 0 );
	} else {
	    msgcnt[0] = msgcnt[1] = 0;
	}
	
	for (pj = 0; pj < Pc; ++pj) {
	    if ( ToSendR[0][pj] != EMPTY ) {
#if ( PROFlevel>=1 )
		TIC(t1);
#endif
#if ( VAMPIR>=1 )
		VT_begin(1);
#endif
		MPI_Send( lsub, msgcnt[0], mpi_int_t, pj, 0, scp->comm );
		MPI_Send( lusup, msgcnt[1], MPI_DOUBLE, pj, 1, scp->comm );
#if ( DEBUGlevel>=2 )
		printf("(%d) Send L(:,%4d): lsub %4d, lusup %4d to Pc %2d\n",
		       iam, 0, msgcnt[0], msgcnt[1], pj);
#endif
#if ( VAMPIR>=1 )
		VT_end(1);
#endif
#if ( PROFlevel>=1 )
		TOC(t2, t1);
		stat->utime[COMM] += t2;
#endif
	    }
	} /* for pj ... */
    } /* if mycol == 0 */

    /*
     * Loop through all block columns.
     */
    for (k = 0; k < nsupers; ++k) {

	knsupc = SuperSize( k );
	krow = PROW( k, grid );
	kcol = PCOL( k, grid );

	if ( mycol == kcol ) {
	    lk = LBj( k, grid ); /* Local block number. */
	    lsub = Lrowind_bc_ptr[lk];
	    lusup = Lnzval_bc_ptr[lk];
	} else {
	    if ( ToRecv[k] >= 1 ) { /* Recv block column L(:,k). */
		scp = &grid->rscp; /* The scope of process row. */
#if ( PROFlevel>=1 )
		TIC(t1);
#endif
#if ( VAMPIR>=1 )
		VT_begin(2);
#endif
		/*probe_recv(iam, kcol, (4*k)%NTAGS, mpi_int_t, scp->comm, 
		  Llu->bufmax[0]);*/
		MPI_Recv( Lsub_buf, Llu->bufmax[0], mpi_int_t, kcol, 
			 (4*k)%NTAGS, scp->comm, &status );
		MPI_Get_count( &status, mpi_int_t, &msgcnt[0] );
		/*probe_recv(iam, kcol, (4*k+1)%NTAGS, MPI_DOUBLE, scp->comm, 
		  Llu->bufmax[1]);*/
		MPI_Recv( Lval_buf, Llu->bufmax[1], MPI_DOUBLE, kcol, 
			 (4*k+1)%NTAGS, scp->comm, &status );
		MPI_Get_count( &status, MPI_DOUBLE, &msgcnt[1] );
#if ( VAMPIR>=1 )
		VT_end(2);
#endif
#if ( PROFlevel>=1 )
		TOC(t2, t1);
		stat->utime[COMM] += t2;
#endif
#if ( DEBUGlevel>=2 )
		printf("(%d) Recv L(:,%4d): lsub %4d, lusup %4d from Pc %2d\n",
		       iam, k, msgcnt[0], msgcnt[1], kcol);
#endif
		lsub = Lsub_buf;
		lusup = Lval_buf;
#if ( PRNTlevel==3 )
		++total_msg;
		if ( !msgcnt[0] ) ++zero_msg;
#endif
	    } else msgcnt[0] = 0;
	} /* if mycol = Pc(k) */

	scp = &grid->cscp; /* The scope of process column. */

	if ( myrow == krow ) {
	    /* Parallel triangular solve across process row *krow* --
	       U(k,j) = L(k,k) \ A(k,j).  */
#ifdef _CRAY
	    pdgstrs2(n, k, Glu_persist, grid, Llu, stat, ftcs1, ftcs2, ftcs3);
#else
	    pdgstrs2(n, k, Glu_persist, grid, Llu, stat);
#endif

	    /* Multicasts U(k,:) to process columns. */
	    lk = LBi( k, grid );
	    usub = Ufstnz_br_ptr[lk];
	    uval = Unzval_br_ptr[lk];
	    if ( usub )	{
		msgcnt[2] = usub[2];
		msgcnt[3] = usub[1];
	    } else {
		msgcnt[2] = msgcnt[3] = 0;
	    }

	    if ( ToSendD[lk] == YES ) {
		for (pi = 0; pi < Pr; ++pi) {
		    if ( pi != myrow ) {
#if ( PROFlevel>=1 )
			TIC(t1);
#endif
#if ( VAMPIR>=1 )
			VT_begin(3);
#endif
			MPI_Send( usub, msgcnt[2], mpi_int_t, pi,
				 (4*k+2)%NTAGS, scp->comm);
			MPI_Send( uval, msgcnt[3], MPI_DOUBLE, pi,
				 (4*k+3)%NTAGS, scp->comm);
#if ( VAMPIR>=1 )
			VT_end(3);
#endif
#if ( PROFlevel>=1 )
			TOC(t2, t1);
			stat->utime[COMM] += t2;
#endif
#if ( DEBUGlevel>=2 )
			printf("(%d) Send U(%4d,:) to Pr %2d\n", iam, k, pi);
#endif
		    } /* if pi ... */
		} /* for pi ... */
	    } /* if ToSendD ... */
	} else { /* myrow != krow */
	    if ( ToRecv[k] == 2 ) { /* Recv block row U(k,:). */
#if ( PROFlevel>=1 )
		TIC(t1);
#endif
#if ( VAMPIR>=1 )
		VT_begin(4);
#endif
		/*probe_recv(iam, krow, (4*k+2)%NTAGS, mpi_int_t, scp->comm, 
		  Llu->bufmax[2]);*/
		MPI_Recv( Usub_buf, Llu->bufmax[2], mpi_int_t, krow,
			 (4*k+2)%NTAGS, scp->comm, &status );
		MPI_Get_count( &status, mpi_int_t, &msgcnt[2] );
		/*probe_recv(iam, krow, (4*k+3)%NTAGS, MPI_DOUBLE, scp->comm, 
		  Llu->bufmax[3]);*/
		MPI_Recv( Uval_buf, Llu->bufmax[3], MPI_DOUBLE, krow, 
			 (4*k+3)%NTAGS, scp->comm, &status );
		MPI_Get_count( &status, MPI_DOUBLE, &msgcnt[3] );
#if ( VAMPIR>=1 )
		VT_end(4);
#endif
#if ( PROFlevel>=1 )
		TOC(t2, t1);
		stat->utime[COMM] += t2;
#endif
		usub = Usub_buf;
		uval = Uval_buf;
#if ( DEBUGlevel>=2 )
		printf("(%d) Recv U(%4d,:) from Pr %2d\n", iam, k, krow);
#endif
#if ( PRNTlevel==3 )
		++total_msg;
		if ( !msgcnt[2] ) ++zero_msg;
#endif
	    } else msgcnt[2] = 0;
	} /* if myrow == Pr(k) */
	  
	/* 
	 * Parallel rank-k update; pair up blocks L(i,k) and U(k,j).
	 *  for (j = k+1; k < N; ++k) {
	 *     for (i = k+1; i < N; ++i) 
	 *         if ( myrow == PROW( i, grid ) && mycol == PCOL( j, grid )
	 *              && L(i,k) != 0 && U(k,j) != 0 )
	 *             A(i,j) = A(i,j) - L(i,k) * U(k,j);
	 */
	msg0 = msgcnt[0];
	msg2 = msgcnt[2];
	if ( msg0 && msg2 ) { /* L(:,k) and U(k,:) are not empty. */
	    nsupr = lsub[1]; /* LDA of lusup. */
	    if ( myrow == krow ) { /* Skip diagonal block L(k,k). */
		lptr0 = BC_HEADER + LB_DESCRIPTOR + lsub[BC_HEADER+1];
		luptr0 = knsupc;
		nlb = lsub[0] - 1;
	    } else {
		lptr0 = BC_HEADER;
		luptr0 = 0;
		nlb = lsub[0];
	    }
	    lptr = lptr0;
	    for (lb = 0; lb < nlb; ++lb) { /* Initialize block row pointers. */
		ib = lsub[lptr];
		lib = LBi( ib, grid );
		iuip[lib] = BR_HEADER;
		ruip[lib] = 0;
		lptr += LB_DESCRIPTOR + lsub[lptr+1];
	    }
	    nub = usub[0];    /* Number of blocks in the block row U(k,:) */
	    iukp = BR_HEADER; /* Skip header; Pointer to index[] of U(k,:) */
	    rukp = 0;         /* Pointer to nzval[] of U(k,:) */
	    klst = FstBlockC( k+1 );
	    
	    /* 
	     * Update the first block column A(:,k+1).
	     */
	    jb = usub[iukp];   /* Global block number of block U(k,j). */
	    if ( jb == k+1 ) { /* First update (k+1)-th block. */
		--nub;
		lptr = lptr0;
		luptr = luptr0;
		ljb = LBj( jb, grid ); /* Local block number of U(k,j). */
		nsupc = SuperSize( jb );
		iukp += UB_DESCRIPTOR; /* Start fstnz of block U(k,j). */

		/* Prepare to call DGEMM. */
		jj = iukp;
		while ( usub[jj] == klst ) ++jj;
		ldu = klst - usub[jj++];
		ncols = 1;
		full = 1;
		for (; jj < iukp+nsupc; ++jj) {
		    segsize = klst - usub[jj];
		    if ( segsize ) {
		        ++ncols;
			if ( segsize != ldu ) full = 0;
		        if ( segsize > ldu ) ldu = segsize;
		    }
		}
#if ( DEBUGlevel>=1 )
		++num_update;
#endif
		if ( full ) {
		    tempu = &uval[rukp];
		} else { /* Copy block U(k,j) into tempU2d. */
#if ( DEBUGlevel>=1 )
		  printf("(%d) full=%d,k=%d,jb=%d,ldu=%d,ncols=%d,nsupc=%d\n",
			 iam, full, k, jb, ldu, ncols, nsupc);
		  ++num_copy;
#endif
		    tempu = tempU2d;
		    for (jj = iukp; jj < iukp+nsupc; ++jj) {
		        segsize = klst - usub[jj];
			if ( segsize ) {
			    lead_zero = ldu - segsize;
			    for (i = 0; i < lead_zero; ++i) tempu[i] = 0.0;
			    tempu += lead_zero;
			    for (i = 0; i < segsize; ++i)
				tempu[i] = uval[rukp+i];
			    rukp += segsize;
			    tempu += segsize;
			}
		    }
		    tempu = tempU2d;
		    rukp -= usub[iukp - 1]; /* Return to start of U(k,j). */
		} /* if full ... */

		for (lb = 0; lb < nlb; ++lb) { 
		    ib = lsub[lptr]; /* Row block L(i,k). */
		    nbrow = lsub[lptr+1];  /* Number of full rows. */
		    lptr += LB_DESCRIPTOR; /* Skip descriptor. */
		    tempv = tempv2d;
#ifdef _CRAY
		    SGEMM(ftcs, ftcs, &nbrow, &ncols, &ldu, &alpha, 
			  &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			  tempu, &ldu, &beta, tempv, &ldt);
#else
		    dgemm_("N", "N", &nbrow, &ncols, &ldu, &alpha, 
			   &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			   tempu, &ldu, &beta, tempv, &ldt);
#endif
		    stat->ops[FACT] += 2 * nbrow * ldu * ncols;

		    /* Now gather the result into the destination block. */
		    if ( ib < jb ) { /* A(i,j) is in U. */
			ilst = FstBlockC( ib+1 );
			lib = LBi( ib, grid );
			index = Ufstnz_br_ptr[lib];
			ijb = index[iuip[lib]];
			while ( ijb < jb ) { /* Search for dest block. */
			    ruip[lib] += index[iuip[lib]+1];
			    iuip[lib] += UB_DESCRIPTOR + SuperSize( ijb );
			    ijb = index[iuip[lib]];
			}
			iuip[lib] += UB_DESCRIPTOR; /* Skip descriptor. */

			tempv = tempv2d;
			for (jj = 0; jj < nsupc; ++jj) {
			    segsize = klst - usub[iukp + jj];
			    fnz = index[iuip[lib]++];
			    if ( segsize ) { /* Nonzero segment in U(k.j). */
				ucol = &Unzval_br_ptr[lib][ruip[lib]];
				for (i = 0, it = 0; i < nbrow; ++i) {
				    rel = lsub[lptr + i] - fnz;
				    ucol[rel] -= tempv[it++];
				}
				tempv += ldt;
			    }
			    ruip[lib] += ilst - fnz;
			}
		    } else { /* A(i,j) is in L. */
			index = Lrowind_bc_ptr[ljb];
			ldv = index[1];   /* LDA of the dest lusup. */
			lptrj = BC_HEADER;
			luptrj = 0;
			ijb = index[lptrj];
			while ( ijb != ib ) { /* Search for dest block -- 
						 blocks are not ordered! */
			    luptrj += index[lptrj+1];
			    lptrj += LB_DESCRIPTOR + index[lptrj+1];
			    ijb = index[lptrj];
			}
			/*
			 * Build indirect table. This is needed because the
			 * indices are not sorted.
			 */
			fnz = FstBlockC( ib );
			lptrj += LB_DESCRIPTOR;
			for (i = 0; i < index[lptrj-1]; ++i) {
			    rel = index[lptrj + i] - fnz;
			    indirect[rel] = i;
			}
			nzval = Lnzval_bc_ptr[ljb] + luptrj;
			tempv = tempv2d;
			for (jj = 0; jj < nsupc; ++jj) {
			    segsize = klst - usub[iukp + jj];
			    if ( segsize ) {
/*#pragma _CRI cache_bypass nzval,tempv*/
				for (it = 0, i = 0; i < nbrow; ++i) {
				    rel = lsub[lptr + i] - fnz;
				    nzval[indirect[rel]] -= tempv[it++];
				}
				tempv += ldt;
			    }
			    nzval += ldv;
			}
		    } /* if ib < jb ... */
		    lptr += nbrow;
		    luptr += nbrow;
		} /* for lb ... */
		rukp += usub[iukp - 1]; /* Move to block U(k,j+1) */
		iukp += nsupc;
	    }  /* if jb == k+1 */
	} /* if L(:,k) and U(k,:) not empty */


	if ( k+1 < nsupers && mycol == PCOL( k+1, grid ) ) {
#if ( VAMPIR>=1 )
	VT_begin(5);
#endif
	    /* Factor diagonal and subdiagonal blocks and test for exact
	       singularity.  */
	    pdgstrf2(options, k+1, thresh, Glu_persist, grid, Llu, stat, info);
#if ( VAMPIR>=1 )
	VT_end(5);
#endif

	    /* Process column *kcol* multicasts numeric values of L(:,k) 
	       to process rows. */
	    lk = LBj( k+1, grid ); /* Local block number. */
	    lsub1 = Lrowind_bc_ptr[lk];
 	    if ( lsub1 ) {
		msgcnt[0] = lsub1[1] + BC_HEADER + lsub1[0]*LB_DESCRIPTOR;
		msgcnt[1] = lsub1[1] * SuperSize( k+1 );
	    } else {
		msgcnt[0] = 0;
		msgcnt[1] = 0;
	    }
	    scp = &grid->rscp; /* The scope of process row. */
	    for (pj = 0; pj < Pc; ++pj) {
		if ( ToSendR[lk][pj] != EMPTY ) {
		    lusup1 = Lnzval_bc_ptr[lk];
#if ( PROFlevel>=1 )
		    TIC(t1);
#endif
#if ( VAMPIR>=1 )
		    VT_begin(1);
#endif
#if 1
		    MPI_Send( lsub1, msgcnt[0], mpi_int_t, pj,
			     (4*(k+1))%NTAGS, scp->comm );
		    MPI_Send( lusup1, msgcnt[1], MPI_DOUBLE, pj,
			     (4*(k+1)+1)%NTAGS, scp->comm );
#else
/* This is worse than blocking send. */
		    MPI_Isend( lsub1, msgcnt[0], mpi_int_t, pj,
			      (4*(k+1))%NTAGS, scp->comm, &request[0] );
		    MPI_Isend( lusup1, msgcnt[1], MPI_DOUBLE, pj,
			     (4*(k+1)+1)%NTAGS, scp->comm, &request[1] );
#endif
#if ( VAMPIR>=1 )
		    VT_end(1);
#endif
#if ( PROFlevel>=1 )
		    TOC(t2, t1);
		    stat->utime[COMM] += t2;
#endif
#if ( DEBUGlevel>=2 )
		    printf("(%d) Send L(:,%4d): lsub %4d, lusup %4d to Pc %2d\n",
			   iam, k+1, msgcnt[0], msgcnt[1], pj);
#endif
		}
	    } /* for pj ... */
	} /* if mycol == Pc(k+1) */

	if ( msg0 && msg2 ) { /* L(:,k) and U(k,:) are not empty. */
	    /* 
	     * Update all other blocks using block row U(k,:)
	     */
	    for (j = 0; j < nub; ++j) { 
		lptr = lptr0;
		luptr = luptr0;
		jb = usub[iukp];  /* Global block number of block U(k,j). */
		ljb = LBj( jb, grid ); /* Local block number of U(k,j). */
		nsupc = SuperSize( jb );
		iukp += UB_DESCRIPTOR; /* Start fstnz of block U(k,j). */

		/* Prepare to call DGEMM. */
		jj = iukp;
		while ( usub[jj] == klst ) ++jj;
		ldu = klst - usub[jj++];
		ncols = 1;
		full = 1;
		for (; jj < iukp+nsupc; ++jj) {
		    segsize = klst - usub[jj];
		    if ( segsize ) {
		        ++ncols;
			if ( segsize != ldu ) full = 0;
		        if ( segsize > ldu ) ldu = segsize;
		    }
		}
#if ( DEBUGlevel>=2 )
		printf("(%d) full=%d,k=%d,jb=%d,ldu=%d,ncols=%d,nsupc=%d\n",
		       iam, full, k, jb, ldu, ncols, nsupc);
		++num_update;
#endif
		if ( full ) {
		    tempu = &uval[rukp];
		} else { /* Copy block U(k,j) into tempU2d. */
#if ( DEBUGlevel>=1 ) 		  
		    ++num_copy;
#endif
		    tempu = tempU2d;
		    for (jj = iukp; jj < iukp+nsupc; ++jj) {
		        segsize = klst - usub[jj];
			if ( segsize ) {
			    lead_zero = ldu - segsize;
			    for (i = 0; i < lead_zero; ++i) tempu[i] = 0.0;
			    tempu += lead_zero;
			    for (i = 0; i < segsize; ++i)
			        tempu[i] = uval[rukp+i];
			    rukp += segsize;
			    tempu += segsize;
			}
		    }
		    tempu = tempU2d;
		    rukp -= usub[iukp - 1]; /* Return to start of U(k,j). */
		} /* if full ... */

		for (lb = 0; lb < nlb; ++lb) { 
		    ib = lsub[lptr];       /* Row block L(i,k). */
		    nbrow = lsub[lptr+1];  /* Number of full rows. */
		    lptr += LB_DESCRIPTOR; /* Skip descriptor. */
		    tempv = tempv2d;
#ifdef _CRAY
		    SGEMM(ftcs, ftcs, &nbrow, &ncols, &ldu, &alpha, 
			  &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			  tempu, &ldu, &beta, tempv, &ldt);
#else
		    dgemm_("N", "N", &nbrow, &ncols, &ldu, &alpha, 
			   &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr, 
			   tempu, &ldu, &beta, tempv, &ldt);
#endif
		    stat->ops[FACT] += 2 * nbrow * ldu * ncols;

		    /* Now gather the result into the destination block. */
		    if ( ib < jb ) { /* A(i,j) is in U. */
			ilst = FstBlockC( ib+1 );
			lib = LBi( ib, grid );
			index = Ufstnz_br_ptr[lib];
			ijb = index[iuip[lib]];
			while ( ijb < jb ) { /* Search for dest block. */
			    ruip[lib] += index[iuip[lib]+1];
			    iuip[lib] += UB_DESCRIPTOR + SuperSize( ijb );
			    ijb = index[iuip[lib]];
			}
			iuip[lib] += UB_DESCRIPTOR; /* Skip descriptor. */

			tempv = tempv2d;
			for (jj = 0; jj < nsupc; ++jj) {
			    segsize = klst - usub[iukp + jj];
			    fnz = index[iuip[lib]++];
			    if ( segsize ) { /* Nonzero segment in U(k.j). */
				ucol = &Unzval_br_ptr[lib][ruip[lib]];
				for (i = 0, it = 0; i < nbrow; ++i) {
				    rel = lsub[lptr + i] - fnz;
				    ucol[rel] -= tempv[it++];
				}
				tempv += ldt;
			    }
			    ruip[lib] += ilst - fnz;
			}
		    } else { /* A(i,j) is in L. */
			index = Lrowind_bc_ptr[ljb];
			ldv = index[1];   /* LDA of the dest lusup. */
			lptrj = BC_HEADER;
			luptrj = 0;
			ijb = index[lptrj];
			while ( ijb != ib ) { /* Search for dest block -- 
						 blocks are not ordered! */
			    luptrj += index[lptrj+1];
			    lptrj += LB_DESCRIPTOR + index[lptrj+1];
			    ijb = index[lptrj];
			}
			/*
			 * Build indirect table. This is needed because the
			 * indices are not sorted.
			 */
			fnz = FstBlockC( ib );
			lptrj += LB_DESCRIPTOR;
			for (i = 0; i < index[lptrj-1]; ++i) {
			    rel = index[lptrj + i] - fnz;
			    indirect[rel] = i;
			}
			nzval = Lnzval_bc_ptr[ljb] + luptrj;
			tempv = tempv2d;
			for (jj = 0; jj < nsupc; ++jj) {
			    segsize = klst - usub[iukp + jj];
			    if ( segsize ) {
/*#pragma _CRI cache_bypass nzval,tempv*/
				for (it = 0, i = 0; i < nbrow; ++i) {
				    rel = lsub[lptr + i] - fnz;
				    nzval[indirect[rel]] -= tempv[it++];
				}
				tempv += ldt;
			    }
			    nzval += ldv;
			}
		    } /* if ib < jb ... */
		    lptr += nbrow;
		    luptr += nbrow;
		} /* for lb ... */
		rukp += usub[iukp - 1]; /* Move to block U(k,j+1) */
		iukp += nsupc;
	    } /* for j ... */
	} /* if  k L(:,k) and U(k,:) are not empty */

    } /* for k ... */

#if ( VAMPIR>=1 )
    VT_end(100);
    VT_traceoff();
#endif

    SUPERLU_FREE(Llu->Lsub_buf);
    SUPERLU_FREE(Llu->Lval_buf);
    SUPERLU_FREE(Llu->Usub_buf);
    SUPERLU_FREE(Llu->Uval_buf);
    SUPERLU_FREE(Llu->ujrow);

    SUPERLU_FREE(tempv2d);
    SUPERLU_FREE(indirect);
    SUPERLU_FREE(iuip);
    SUPERLU_FREE(ruip);

    /* Prepare error message. */
    if ( *info == 0 ) *info = n + 1;
#if ( PROFlevel>=1 )
    TIC(t1);
#endif
    MPI_Allreduce( info, &iinfo, 1, mpi_int_t, MPI_MIN, grid->comm );
#if ( PROFlevel>=1 )
    TOC(t2, t1);
    stat->utime[COMM] += t2;
#endif
    if ( iinfo == n + 1 ) *info = 0;
    else *info = iinfo;


#if ( PRNTlevel==3 )
    MPI_Allreduce( &zero_msg, &iinfo, 1, mpi_int_t, MPI_SUM, grid->comm );
    if ( !iam ) printf(".. # msg of zero size\t%d\n", iinfo);
    MPI_Allreduce( &total_msg, &iinfo, 1, mpi_int_t, MPI_SUM, grid->comm );
    if ( !iam ) printf(".. # total msg\t%d\n", iinfo);
#endif

#if ( PRNTlevel==2 )
    for (i = 0; i < Pr * Pc; ++i) {
	if ( iam == i ) {
	    PrintLblocks(iam, nsupers, grid, Glu_persist, Llu);
	    PrintUblocks(iam, nsupers, grid, Glu_persist, Llu);
	    printf("(%d)\n", iam);
	    PrintInt10("Recv", nsupers, Llu->ToRecv);
	}
	MPI_Barrier( grid->comm );
    }
#endif

#if ( DEBUGlevel>=1 )
    printf("(%d) num_copy=%d, num_update=%d\n", iam, num_copy, num_update);
    CHECK_MALLOC(iam, "Exit pdgstrf()");
#endif

} /* PDGSTRF */


/************************************************************************/
static void pdgstrf2
/************************************************************************/
(
 superlu_options_t *options,
 int_t k, double thresh, Glu_persist_t *Glu_persist, gridinfo_t *grid,
 LocalLU_t *Llu, SuperLUStat_t *stat, int* info
 )
/* 
 * Purpose
 * =======
 *   Factor diagonal and subdiagonal blocks and test for exact singularity.
 *   Only the process column that owns block column *k* participates
 *   in the work.
 * 
 * Arguments
 * =========
 *
 * k      (input) int (global)
 *        The column number of the block column to be factorized.
 *
 * thresh (input) double (global)
 *        The threshold value = s_eps * anorm.
 *
 * Glu_persist (input) Glu_persist_t*
 *        Global data structures (xsup, supno) replicated on all processes.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Llu    (input/output) LocalLU_t*
 *        Local data structures to store distributed L and U matrices.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics about the factorization.
 *        See SuperLUStat_t structure defined in util.h.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, U(i,i) is exactly zero. The factorization has
 *             been completed, but the factor U is exactly singular,
 *             and division by zero will occur if it is used to solve a
 *             system of equations.
 *
 */
{
    int    c, iam, l, pkk;
    int    incx = 1, incy = 1;
    int    nsupr; /* number of rows in the block (LDA) */
    int    luptr;
    int_t  i, krow, j, jfst, jlst;
    int_t  nsupc; /* number of columns in the block */
    int_t  *xsup = Glu_persist->xsup;
    double *lusup, temp, alpha = -1;
    double *ujrow;

    *info = 0;

    /* Quick return. */

    /* Initialization. */
    iam   = grid->iam;
    krow  = PROW( k, grid );
    pkk   = PNUM( PROW(k, grid), PCOL(k, grid), grid );
    j     = LBj( k, grid ); /* Local block number */
    jfst  = FstBlockC( k );
    jlst  = FstBlockC( k+1 );
    lusup = Llu->Lnzval_bc_ptr[j];
    nsupc = SuperSize( k );
    if ( Llu->Lrowind_bc_ptr[j] ) nsupr = Llu->Lrowind_bc_ptr[j][1];
    ujrow = Llu->ujrow;

    luptr = 0; /* Point to the diagonal entries. */
    c = nsupc;
    for (j = 0; j < jlst - jfst; ++j) {
	/* Broadcast the j-th row (nsupc - j) elements to
	   the process column. */
	if ( iam == pkk ) { /* Diagonal process. */
	    i = luptr;
	    if ( options->ReplaceTinyPivot == YES || lusup[i] == 0.0 ) {
		if ( fabs(lusup[i]) < thresh ) { /* Diagonal */
#if ( PRNTlevel>=2 )
		    printf("(%d) .. col %d, tiny pivot %e  ",
			   iam, jfst+j, lusup[i]);
#endif
		    /* Keep the replaced diagonal with the same sign. */
		    if ( lusup[i] < 0 ) lusup[i] = -thresh;
		    else lusup[i] = thresh;
#if ( PRNTlevel>=2 )
		    printf("replaced by %e\n", lusup[i]);
#endif
		    ++(stat->TinyPivots);
		}
	    }
	    for (l = 0; l < c; ++l, i += nsupr)	ujrow[l] = lusup[i];
	}
#if 0
	dbcast_col(ujrow, c, pkk, UjROW, grid, &c);
#else
	/*bcast_tree(ujrow, c, MPI_DOUBLE, krow, (24*k+j)%NTAGS,
		   grid, COMM_COLUMN, &c);*/
	MPI_Bcast(ujrow, c, MPI_DOUBLE, krow, (grid->cscp).comm);
#endif

#if ( DEBUGlevel>=2 )
if ( k == 3329 && j == 2 ) {
	if ( iam == pkk ) {
	    printf("..(%d) k %d, j %d: Send ujrow[0] %e\n",iam,k,j,ujrow[0]);
	} else {
	    printf("..(%d) k %d, j %d: Recv ujrow[0] %e\n",iam,k,j,ujrow[0]);
	}
}
#endif

	if ( !lusup ) { /* Empty block column. */
	    --c;
	    if ( ujrow[0] == 0.0 ) *info = j+jfst+1;
	    continue;
	}

	/* Test for singularity. */
	if ( ujrow[0] == 0.0 ) {
	    *info = j+jfst+1;
	} else {
	    /* Scale the j-th column of the matrix. */
	    temp = 1.0 / ujrow[0];
	    if ( iam == pkk ) {
		for (i = luptr+1; i < luptr-j+nsupr; ++i) lusup[i] *= temp;
		stat->ops[FACT] += nsupr-j-1;
	    } else {
		for (i = luptr; i < luptr+nsupr; ++i) lusup[i] *= temp;
		stat->ops[FACT] += nsupr;
	    }
	}
	    
	/* Rank-1 update of the trailing submatrix. */
	if ( --c ) {
	    if ( iam == pkk ) {
		l = nsupr - j - 1;
#ifdef _CRAY
		SGER(&l, &c, &alpha, &lusup[luptr+1], &incx,
		     &ujrow[1], &incy, &lusup[luptr+nsupr+1], &nsupr);
#else
		dger_(&l, &c, &alpha, &lusup[luptr+1], &incx,
		      &ujrow[1], &incy, &lusup[luptr+nsupr+1], &nsupr);
#endif
		stat->ops[FACT] += 2 * l * c;
	    } else {
#ifdef _CRAY
		SGER(&nsupr, &c, &alpha, &lusup[luptr], &incx, 
		     &ujrow[1], &incy, &lusup[luptr+nsupr], &nsupr);
#else
		dger_(&nsupr, &c, &alpha, &lusup[luptr], &incx, 
		      &ujrow[1], &incy, &lusup[luptr+nsupr], &nsupr);
#endif
		stat->ops[FACT] += 2 * nsupr * c;
	    }
	}
	
	/* Move to the next column. */
	if ( iam == pkk ) luptr += nsupr + 1;
	else luptr += nsupr;

    } /* for j ... */

} /* PDGSTRF2 */


/************************************************************************/
static void pdgstrs2
/************************************************************************/
#ifdef _CRAY
(
 int_t m, int_t k, Glu_persist_t *Glu_persist, gridinfo_t *grid,
 LocalLU_t *Llu, SuperLUStat_t *stat, _fcd ftcs1, _fcd ftcs2, _fcd ftcs3
 )
#else
(
 int_t m, int_t k, Glu_persist_t *Glu_persist, gridinfo_t *grid,
 LocalLU_t *Llu, SuperLUStat_t *stat
 )
#endif
/* 
 * Purpose
 * =======
 *   Perform parallel triangular solves
 *           U(k,:) := A(k,:) \ L(k,k). 
 *   Only the process column that owns block column *k* participates
 *   in the work.
 * 
 * Arguments
 * =========
 *
 * m      (input) int (global)
 *        Number of rows in the matrix.
 *
 * k      (input) int (global)
 *        The row number of the block row to be factorized.
 *
 * Glu_persist (input) Glu_persist_t*
 *        Global data structures (xsup, supno) replicated on all processes.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Llu    (input/output) LocalLU_t*
 *        Local data structures to store distributed L and U matrices.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics about the factorization; 
 *        See SuperLUStat_t structure defined in util.h.
 *
 */
{
    int    iam, pkk;
    int    incx = 1;
    int    nsupr; /* number of rows in the block L(:,k) (LDA) */
    int    segsize;
    int_t  nsupc; /* number of columns in the block */
    int_t  luptr, iukp, rukp;
    int_t  b, gb, j, klst, knsupc, lk, nb;
    int_t  *xsup = Glu_persist->xsup;
    int_t  *usub;
    double *lusup, *uval;

    /* Quick return. */
    lk = LBi( k, grid ); /* Local block number */
    if ( !Llu->Unzval_br_ptr[lk] ) return;

    /* Initialization. */
    iam  = grid->iam;
    pkk  = PNUM( PROW(k, grid), PCOL(k, grid), grid );
    klst = FstBlockC( k+1 );
    knsupc = SuperSize( k );
    usub = Llu->Ufstnz_br_ptr[lk]; /* index[] of block row U(k,:) */
    uval = Llu->Unzval_br_ptr[lk];
    nb = usub[0];
    iukp = BR_HEADER;
    rukp = 0;
    if ( iam == pkk ) {
	lk = LBj( k, grid );
	nsupr = Llu->Lrowind_bc_ptr[lk][1]; /* LDA of lusup[] */
	lusup = Llu->Lnzval_bc_ptr[lk];
    } else {
	nsupr = Llu->Lsub_buf[1]; /* LDA of lusup[] */
	lusup = Llu->Lval_buf;
    }

    /* Loop through all the row blocks. */
    for (b = 0; b < nb; ++b) {
	gb = usub[iukp];
	nsupc = SuperSize( gb );
	iukp += UB_DESCRIPTOR;

	/* Loop through all the segments in the block. */
	for (j = 0; j < nsupc; ++j) {
	    segsize = klst - usub[iukp++]; 
	    if ( segsize ) { /* Nonzero segment. */
		luptr = (knsupc - segsize) * (nsupr + 1);
#ifdef _CRAY
		STRSV(ftcs1, ftcs2, ftcs3, &segsize, &lusup[luptr], &nsupr, 
		      &uval[rukp], &incx);
#else
		dtrsv_("L", "N", "U", &segsize, &lusup[luptr], &nsupr, 
		       &uval[rukp], &incx);
#endif
		stat->ops[FACT] += segsize * (segsize + 1);
		rukp += segsize;
	    }
	}
    } /* for b ... */

} /* PDGSTRS2 */

int
probe_recv(int iam, int source, int tag, MPI_Datatype datatype, MPI_Comm comm,
	   int buf_size)
{
    MPI_Status status;
    int count; 

    MPI_Probe( source, tag, comm, &status );
    MPI_Get_count( &status, datatype, &count );
    if ( count > buf_size ) {
        printf("(%d) Recv'ed count %d > buffer size $d\n",
	       iam, count, buf_size);
	exit(-1);
    }
    return 0;
}
