

/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 */

#include "superlu_ddefs.h"

#define ISEND_IRECV

/*
 * Function prototypes
 */
#ifdef _CRAY
fortran void STRSM(_fcd, _fcd, _fcd, _fcd, int*, int*, double*,
		   double*, int*, double*, int*);
fortran void SGEMM(_fcd, _fcd, int*, int*, int*, double*, double*, 
		   int*, double*, int*, double*, double*, int*);
_fcd ftcs1;
_fcd ftcs2;
_fcd ftcs3;
#endif


int_t
pdReDistribute_B_to_X(double *B, int_t m_loc, int nrhs, int_t ldb,
                      int_t fst_row, int_t *ilsum, double *x,
		      ScalePermstruct_t *ScalePermstruct,
		      Glu_persist_t *Glu_persist,
		      gridinfo_t *grid, SOLVEstruct_t *SOLVEstruct)
{
/*
 * Purpose
 * =======
 *   Re-distribute B on the diagonal processes of the 2D process mesh.
 * 
 * Note
 * ====
 *   This routine can only be called after the routine pxgstrs_init(),
 *   in which the structures of the send and receive buffers are set up.
 *
 * Arguments
 * =========
 * 
 * B      (input) double*
 *        The distributed right-hand side matrix of the possibly
 *        equilibrated system.
 *
 * m_loc  (input) int (local)
 *        The local row dimension of matrix B.
 *
 * nrhs   (input) int (global)
 *        Number of right-hand sides.
 *
 * ldb    (input) int (local)
 *        Leading dimension of matrix B.
 *
 * fst_row (input) int (global)
 *        The row number of B's first row in the global matrix.
 *
 * ilsum  (input) int* (global)
 *        Starting position of each supernode in a full array.
 *
 * x      (output) double*
 *        The solution vector. It is valid only on the diagonal processes.
 *
 * ScalePermstruct (input) ScalePermstruct_t*
 *        The data structure to store the scaling and permutation vectors
 *        describing the transformations performed to the original matrix A.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * SOLVEstruct (input) SOLVEstruct_t*
 *        Contains the information for the communication during the
 *        solution phase.
 *
 * Return value
 * ============
 *
 */
    int  *SendCnt, *SendCnt_nrhs, *RecvCnt, *RecvCnt_nrhs;
    int  *sdispls, *sdispls_nrhs, *rdispls, *rdispls_nrhs;
    int  *ptr_to_ibuf, *ptr_to_dbuf;
    int_t  *perm_r, *perm_c; /* row and column permutation vectors */
    int_t  *send_ibuf, *recv_ibuf;
    double *send_dbuf, *recv_dbuf;
    int_t  *xsup, *supno;
    int_t  i, ii, irow, gbi, j, jj, k, knsupc, l, lk;
    int    p, procs;
    pxgstrs_comm_t *gstrs_comm = SOLVEstruct->gstrs_comm;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Enter pdReDistribute_B_to_X()");
#endif

    /* ------------------------------------------------------------
       INITIALIZATION.
       ------------------------------------------------------------*/
    perm_r = ScalePermstruct->perm_r;
    perm_c = ScalePermstruct->perm_c;
    procs = grid->nprow * grid->npcol;
    xsup = Glu_persist->xsup;
    supno = Glu_persist->supno;
    SendCnt      = gstrs_comm->B_to_X_SendCnt;
    SendCnt_nrhs = gstrs_comm->B_to_X_SendCnt +   procs;
    RecvCnt      = gstrs_comm->B_to_X_SendCnt + 2*procs;
    RecvCnt_nrhs = gstrs_comm->B_to_X_SendCnt + 3*procs;
    sdispls      = gstrs_comm->B_to_X_SendCnt + 4*procs;
    sdispls_nrhs = gstrs_comm->B_to_X_SendCnt + 5*procs;
    rdispls      = gstrs_comm->B_to_X_SendCnt + 6*procs;
    rdispls_nrhs = gstrs_comm->B_to_X_SendCnt + 7*procs;
    ptr_to_ibuf  = gstrs_comm->ptr_to_ibuf;
    ptr_to_dbuf  = gstrs_comm->ptr_to_dbuf;

    /* ------------------------------------------------------------
       NOW COMMUNICATE THE ACTUAL DATA.
       ------------------------------------------------------------*/
    k = sdispls[procs-1] + SendCnt[procs-1]; /* Total number of sends */
    l = rdispls[procs-1] + RecvCnt[procs-1]; /* Total number of receives */
    if ( !(send_ibuf = intMalloc_dist(k + l)) )
        ABORT("Malloc fails for send_ibuf[].");
    recv_ibuf = send_ibuf + k;
    if ( !(send_dbuf = doubleMalloc_dist((k + l)* (size_t)nrhs)) )
        ABORT("Malloc fails for send_dbuf[].");
    recv_dbuf = send_dbuf + k * nrhs;
    
    for (p = 0; p < procs; ++p) {
        ptr_to_ibuf[p] = sdispls[p];
        ptr_to_dbuf[p] = sdispls[p] * nrhs;
    }

    /* Copy the row indices and values to the send buffer. */
    for (i = 0, l = fst_row; i < m_loc; ++i, ++l) {
        irow = perm_c[perm_r[l]]; /* Row number in Pc*Pr*A */
	gbi = BlockNum( irow );
	p = PNUM( PROW(gbi,grid), PCOL(gbi,grid), grid ); /* Diagonal process */
	k = ptr_to_ibuf[p];
	send_ibuf[k] = irow;
	k = ptr_to_dbuf[p];
	RHS_ITERATE(j) { /* RHS is stored in row major in the buffer. */
	    send_dbuf[k++] = B[i + j*ldb];
	}
	++ptr_to_ibuf[p];
	ptr_to_dbuf[p] += nrhs;
    }

    /* Communicate the (permuted) row indices. */
    MPI_Alltoallv(send_ibuf, SendCnt, sdispls, mpi_int_t,
		  recv_ibuf, RecvCnt, rdispls, mpi_int_t, grid->comm);

    /* Communicate the numerical values. */
    MPI_Alltoallv(send_dbuf, SendCnt_nrhs, sdispls_nrhs, MPI_DOUBLE,
		  recv_dbuf, RecvCnt_nrhs, rdispls_nrhs, MPI_DOUBLE,
		  grid->comm);
    
    /* ------------------------------------------------------------
       Copy buffer into X on the diagonal processes.
       ------------------------------------------------------------*/
    ii = 0;
    for (p = 0; p < procs; ++p) {
        jj = rdispls_nrhs[p];
        for (i = 0; i < RecvCnt[p]; ++i) {
	    /* Only the diagonal processes do this; the off-diagonal processes
	       have 0 RecvCnt. */
	    irow = recv_ibuf[ii]; /* The permuted row index. */
	    k = BlockNum( irow );
	    knsupc = SuperSize( k );
	    lk = LBi( k, grid );  /* Local block number. */
	    l = X_BLK( lk );
	    x[l - XK_H] = k;      /* Block number prepended in the header. */
	    irow = irow - FstBlockC(k); /* Relative row number in X-block */
	    RHS_ITERATE(j) {
	        x[l + irow + j*knsupc] = recv_dbuf[jj++];
	    }
	    ++ii;
	}
    }

    SUPERLU_FREE(send_ibuf);
    SUPERLU_FREE(send_dbuf);
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Exit pdReDistribute_B_to_X()");
#endif
    return 0;
} /* pdReDistribute_B_to_X */

int_t
pdReDistribute_X_to_B(int_t n, double *B, int_t m_loc, int_t ldb, int_t fst_row,
		      int_t nrhs, double *x, int_t *ilsum,
		      ScalePermstruct_t *ScalePermstruct,
		      Glu_persist_t *Glu_persist, gridinfo_t *grid,
		      SOLVEstruct_t *SOLVEstruct)
{
/*
 * Purpose
 * =======
 *   Re-distribute X on the diagonal processes to B distributed on all
 *   the processes.
 *
 * Note
 * ====
 *   This routine can only be called after the routine pxgstrs_init(),
 *   in which the structures of the send and receive buffers are set up.
 *
 */
    int_t  i, ii, irow, j, jj, k, knsupc, nsupers, l, lk;
    int_t  *xsup, *supno;
    int  *SendCnt, *SendCnt_nrhs, *RecvCnt, *RecvCnt_nrhs;
    int  *sdispls, *rdispls, *sdispls_nrhs, *rdispls_nrhs;
    int  *ptr_to_ibuf, *ptr_to_dbuf;
    int_t  *send_ibuf, *recv_ibuf;
    double *send_dbuf, *recv_dbuf;
    int_t  *row_to_proc = SOLVEstruct->row_to_proc; /* row-process mapping */
    pxgstrs_comm_t *gstrs_comm = SOLVEstruct->gstrs_comm;
    int  iam, p, q, pkk, procs;
    int_t  num_diag_procs, *diag_procs;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Enter pdReDistribute_X_to_B()");
#endif

    /* ------------------------------------------------------------
       INITIALIZATION.
       ------------------------------------------------------------*/
    xsup = Glu_persist->xsup;
    supno = Glu_persist->supno;
    nsupers = Glu_persist->supno[n-1] + 1;
    iam = grid->iam;
    procs = grid->nprow * grid->npcol;
 
    SendCnt      = gstrs_comm->X_to_B_SendCnt;
    SendCnt_nrhs = gstrs_comm->X_to_B_SendCnt +   procs;
    RecvCnt      = gstrs_comm->X_to_B_SendCnt + 2*procs;
    RecvCnt_nrhs = gstrs_comm->X_to_B_SendCnt + 3*procs;
    sdispls      = gstrs_comm->X_to_B_SendCnt + 4*procs;
    sdispls_nrhs = gstrs_comm->X_to_B_SendCnt + 5*procs;
    rdispls      = gstrs_comm->X_to_B_SendCnt + 6*procs;
    rdispls_nrhs = gstrs_comm->X_to_B_SendCnt + 7*procs;
    ptr_to_ibuf  = gstrs_comm->ptr_to_ibuf;
    ptr_to_dbuf  = gstrs_comm->ptr_to_dbuf;

    k = sdispls[procs-1] + SendCnt[procs-1]; /* Total number of sends */
    l = rdispls[procs-1] + RecvCnt[procs-1]; /* Total number of receives */
    if ( !(send_ibuf = intMalloc_dist(k + l)) )
        ABORT("Malloc fails for send_ibuf[].");
    recv_ibuf = send_ibuf + k;
    if ( !(send_dbuf = doubleMalloc_dist((k + l)*nrhs)) )
        ABORT("Malloc fails for send_dbuf[].");
    recv_dbuf = send_dbuf + k * nrhs;
    for (p = 0; p < procs; ++p) {
        ptr_to_ibuf[p] = sdispls[p];
        ptr_to_dbuf[p] = sdispls_nrhs[p];
    }
    num_diag_procs = SOLVEstruct->num_diag_procs;
    diag_procs = SOLVEstruct->diag_procs;

    for (p = 0; p < num_diag_procs; ++p) {  /* For all diagonal processes. */
	pkk = diag_procs[p];
	if ( iam == pkk ) {
	    for (k = p; k < nsupers; k += num_diag_procs) {
		knsupc = SuperSize( k );
		lk = LBi( k, grid ); /* Local block number */
		irow = FstBlockC( k );
		l = X_BLK( lk );
		for (i = 0; i < knsupc; ++i) {
#if 0
		    ii = inv_perm_c[irow]; /* Apply X <== Pc'*Y */
#else
		    ii = irow;
#endif
		    q = row_to_proc[ii];
		    jj = ptr_to_ibuf[q];
		    send_ibuf[jj] = ii;
		    jj = ptr_to_dbuf[q];
		    RHS_ITERATE(j) { /* RHS stored in row major in buffer. */
		        send_dbuf[jj++] = x[l + i + j*knsupc];
		    }
		    ++ptr_to_ibuf[q];
		    ptr_to_dbuf[q] += nrhs;
		    ++irow;
		}
	    }
	}
    }
    
    /* ------------------------------------------------------------
        COMMUNICATE THE (PERMUTED) ROW INDICES AND NUMERICAL VALUES.
       ------------------------------------------------------------*/
    MPI_Alltoallv(send_ibuf, SendCnt, sdispls, mpi_int_t,
		  recv_ibuf, RecvCnt, rdispls, mpi_int_t, grid->comm);
    MPI_Alltoallv(send_dbuf, SendCnt_nrhs, sdispls_nrhs, MPI_DOUBLE, 
		  recv_dbuf, RecvCnt_nrhs, rdispls_nrhs, MPI_DOUBLE,
		  grid->comm);

    /* ------------------------------------------------------------
       COPY THE BUFFER INTO B.
       ------------------------------------------------------------*/
    for (i = 0, k = 0; i < m_loc; ++i) {
	irow = recv_ibuf[i];
	irow -= fst_row; /* Relative row number */
	RHS_ITERATE(j) { /* RHS is stored in row major in the buffer. */
	    B[irow + j*ldb] = recv_dbuf[k++];
	}
    }

    SUPERLU_FREE(send_ibuf);
    SUPERLU_FREE(send_dbuf);
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Exit pdReDistribute_X_to_B()");
#endif
    return 0;

} /* pdReDistribute_X_to_B */


void
pdgstrs(int_t n, LUstruct_t *LUstruct, 
	ScalePermstruct_t *ScalePermstruct,
	gridinfo_t *grid, double *B,
	int_t m_loc, int_t fst_row, int_t ldb, int nrhs,
	SOLVEstruct_t *SOLVEstruct,
	SuperLUStat_t *stat, int *info)
{
/*
 * Purpose
 * =======
 *
 * PDGSTRS solves a system of distributed linear equations
 * A*X = B with a general N-by-N matrix A using the LU factorization
 * computed by PDGSTRF.
 * If the equilibration, and row and column permutations were performed,
 * the LU factorization was performed for A1 where
 *     A1 = Pc*Pr*diag(R)*A*diag(C)*Pc^T = L*U
 * and the linear system solved is
 *     A1 * Y = Pc*Pr*B1, where B was overwritten by B1 = diag(R)*B, and
 * the permutation to B1 by Pc*Pr is applied internally in this routine.
 * 
 * Arguments
 * =========
 *
 * n      (input) int (global)
 *        The order of the system of linear equations.
 *
 * LUstruct (input) LUstruct_t*
 *        The distributed data structures storing L and U factors.
 *        The L and U factors are obtained from PDGSTRF for
 *        the possibly scaled and permuted matrix A.
 *        See superlu_ddefs.h for the definition of 'LUstruct_t'.
 *        A may be scaled and permuted into A1, so that
 *        A1 = Pc*Pr*diag(R)*A*diag(C)*Pc^T = L*U
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh. It contains the MPI communicator, the number
 *        of process rows (NPROW), the number of process columns (NPCOL),
 *        and my process rank. It is an input argument to all the
 *        parallel routines.
 *        Grid can be initialized by subroutine SUPERLU_GRIDINIT.
 *        See superlu_defs.h for the definition of 'gridinfo_t'.
 *
 * B      (input/output) double*
 *        On entry, the distributed right-hand side matrix of the possibly
 *        equilibrated system. That is, B may be overwritten by diag(R)*B.
 *        On exit, the distributed solution matrix Y of the possibly
 *        equilibrated system if info = 0, where Y = Pc*diag(C)^(-1)*X,
 *        and X is the solution of the original system.
 *
 * m_loc  (input) int (local)
 *        The local row dimension of matrix B.
 *
 * fst_row (input) int (global)
 *        The row number of B's first row in the global matrix.
 *
 * ldb    (input) int (local)
 *        The leading dimension of matrix B.
 *
 * nrhs   (input) int (global)
 *        Number of right-hand sides.
 * 
 * SOLVEstruct (output) SOLVEstruct_t* (global)
 *        Contains the information for the communication during the
 *        solution phase.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics about the triangular solves.
 *        See util.h for the definition of 'SuperLUStat_t'.
 *
 * info   (output) int*
 * 	   = 0: successful exit
 *	   < 0: if info = -i, the i-th argument had an illegal value
 *        
 */
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    LocalLU_t *Llu = LUstruct->Llu;
    double alpha = 1.0;
    double zero = 0.0;
    double *lsum;  /* Local running sum of the updates to B-components */
    double *x;     /* X component at step k. */
		    /* NOTE: x and lsum are of same size. */
    double *lusup, *dest;
    double *recvbuf, *tempv;
    double *rtemp; /* Result of full matrix-vector multiply. */
    int_t  **Ufstnz_br_ptr = Llu->Ufstnz_br_ptr;
    int_t  *Urbs, *Urbs1; /* Number of row blocks in each block column of U. */
    Ucb_indptr_t **Ucb_indptr;/* Vertical linked list pointing to Uindex[] */
    int_t  **Ucb_valptr;      /* Vertical linked list pointing to Unzval[] */
    int_t  iam, kcol, krow, mycol, myrow;
    int_t  i, ii, il, j, jj, k, lb, ljb, lk, lptr, luptr;
    int_t  nb, nlb, nub, nsupers;
    int_t  *xsup, *supno, *lsub, *usub;
    int_t  *ilsum;    /* Starting position of each supernode in lsum (LOCAL)*/
    int_t  Pc, Pr;
    int    knsupc, nsupr;
    int    ldalsum;   /* Number of lsum entries locally owned. */
    int    maxrecvsz, p, pi;
    int_t  **Lrowind_bc_ptr;
    double **Lnzval_bc_ptr;
    MPI_Status status;
#ifdef ISEND_IRECV
    MPI_Request *send_req, recv_req;
    int test_flag;
#endif
    pxgstrs_comm_t *gstrs_comm = SOLVEstruct->gstrs_comm;

    /*-- Counts used for L-solve --*/
    int_t  *fmod;         /* Modification count for L-solve --
                             Count the number of local block products to
                             be summed into lsum[lk]. */
    int_t  **fsendx_plist = Llu->fsendx_plist;
    int_t  nfrecvx = Llu->nfrecvx; /* Number of X components to be recv'd. */
    int_t  *frecv;        /* Count of lsum[lk] contributions to be received
                             from processes in this row. 
                             It is only valid on the diagonal processes. */
    int_t  nfrecvmod = 0; /* Count of total modifications to be recv'd. */
    int_t  nleaf = 0, nroot = 0;

    /*-- Counts used for U-solve --*/
    int_t  *bmod;         /* Modification count for L-solve. */
    int_t  **bsendx_plist = Llu->bsendx_plist;
    int_t  nbrecvx = Llu->nbrecvx; /* Number of X components to be recv'd. */
    int_t  *brecv;        /* Count of modifications to be recv'd from
			     processes in this row. */
    int_t  nbrecvmod = 0; /* Count of total modifications to be recv'd. */
    double t;
#if ( DEBUGlevel>=2 )
    int_t Ublocks = 0;
#endif

    t = SuperLU_timer_();

    /* Test input parameters. */
    *info = 0;
    if ( n < 0 ) *info = -1;
    else if ( nrhs < 0 ) *info = -9;
    if ( *info ) {
	pxerbla("PDGSTRS", grid, -*info);
	return;
    }
	
    /*
     * Initialization.
     */
    iam = grid->iam;
    Pc = grid->npcol;
    Pr = grid->nprow;
    myrow = MYROW( iam, grid );
    mycol = MYCOL( iam, grid );
    xsup = Glu_persist->xsup;
    supno = Glu_persist->supno;
    nsupers = supno[n-1] + 1;
    Lrowind_bc_ptr = Llu->Lrowind_bc_ptr;
    Lnzval_bc_ptr = Llu->Lnzval_bc_ptr;
    nlb = CEILING( nsupers, Pr ); /* Number of local block rows. */

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter pdgstrs()");
#endif

    stat->ops[SOLVE] = 0.0;

    /* Save the count to be altered so it can be used by
       subsequent call to PDGSTRS. */
    if ( !(fmod = intMalloc_dist(nlb)) )
	ABORT("Calloc fails for fmod[].");
    for (i = 0; i < nlb; ++i) fmod[i] = Llu->fmod[i];
    if ( !(frecv = intMalloc_dist(nlb)) )
	ABORT("Malloc fails for frecv[].");
    Llu->frecv = frecv;

#ifdef ISEND_IRECV
    if ( !(send_req = (MPI_Request*) SUPERLU_MALLOC(Pr*sizeof(MPI_Request))) )
	ABORT("Malloc fails for send_req[].");
    for (i = 0; i < Pr; ++i) send_req[i] = MPI_REQUEST_NULL;
#endif

#ifdef _CRAY
    ftcs1 = _cptofcd("L", strlen("L"));
    ftcs2 = _cptofcd("N", strlen("N"));
    ftcs3 = _cptofcd("U", strlen("U"));
#endif


    /* Obtain ilsum[] and ldalsum for process column 0. */
    ilsum = Llu->ilsum;
    ldalsum = Llu->ldalsum;

    /* Allocate working storage. */
    knsupc = sp_ienv_dist(3);
    maxrecvsz = knsupc * nrhs + SUPERLU_MAX( XK_H, LSUM_H );
    if ( !(lsum = doubleCalloc_dist(((size_t)ldalsum)*nrhs + nlb*LSUM_H)) )
	ABORT("Calloc fails for lsum[].");
    if ( !(x = doubleMalloc_dist(ldalsum * nrhs + nlb * XK_H)) )
	ABORT("Malloc fails for x[].");
    if ( !(recvbuf = doubleMalloc_dist(maxrecvsz)) )
	ABORT("Malloc fails for recvbuf[].");
    if ( !(rtemp = doubleMalloc_dist(maxrecvsz)) )
	ABORT("Malloc fails for rtemp[].");

    
    /*---------------------------------------------------
     * Forward solve Ly = b.
     *---------------------------------------------------*/
    /* Redistribute B into X on the diagonal processes. */
    pdReDistribute_B_to_X(B, m_loc, nrhs, ldb, fst_row, ilsum, x, 
			  ScalePermstruct, Glu_persist, grid, SOLVEstruct);

    /* Set up the headers in lsum[]. */
    ii = 0;
    for (k = 0; k < nsupers; ++k) {
	knsupc = SuperSize( k );
	krow = PROW( k, grid );
	if ( myrow == krow ) {
	    lk = LBi( k, grid );   /* Local block number. */
	    il = LSUM_BLK( lk );
	    lsum[il - LSUM_H] = k; /* Block number prepended in the header. */
	}
	ii += knsupc;
    }

    /*
     * Compute frecv[] and nfrecvmod counts on the diagonal processes.
     */
    {
	superlu_scope_t *scp = &grid->rscp;

	for (k = 0; k < nsupers; ++k) {
	    krow = PROW( k, grid );
	    if ( myrow == krow ) {
		lk = LBi( k, grid );    /* Local block number. */
		kcol = PCOL( k, grid ); /* Root process in this row scope. */
		if ( mycol != kcol && fmod[lk] )
		    i = 1;  /* Contribution from non-diagonal process. */
		else i = 0;
		MPI_Reduce( &i, &frecv[lk], 1, mpi_int_t,
			   MPI_SUM, kcol, scp->comm );
		if ( mycol == kcol ) { /* Diagonal process. */
		    nfrecvmod += frecv[lk];
		    if ( !frecv[lk] && !fmod[lk] ) ++nleaf;
#if ( DEBUGlevel>=2 )
		    printf("(%2d) frecv[%4d]  %2d\n", iam, k, frecv[lk]);
		    assert( frecv[lk] < Pc );
#endif
		}
	    }
	}
    }

    /* ---------------------------------------------------------
       Solve the leaf nodes first by all the diagonal processes.
       --------------------------------------------------------- */
#if ( DEBUGlevel>=1 )
    printf("(%2d) nleaf %4d\n", iam, nleaf);
#endif
    for (k = 0; k < nsupers && nleaf; ++k) {
	krow = PROW( k, grid );
	kcol = PCOL( k, grid );
	if ( myrow == krow && mycol == kcol ) { /* Diagonal process */
	    knsupc = SuperSize( k );
	    lk = LBi( k, grid );
	    if ( frecv[lk]==0 && fmod[lk]==0 ) {
		fmod[lk] = -1;  /* Do not solve X[k] in the future. */
		ii = X_BLK( lk );
		lk = LBj( k, grid ); /* Local block number, column-wise. */
		lsub = Lrowind_bc_ptr[lk];
		lusup = Lnzval_bc_ptr[lk];
		nsupr = lsub[1];
#ifdef _CRAY
		STRSM(ftcs1, ftcs1, ftcs2, ftcs3, &knsupc, &nrhs, &alpha,
		      lusup, &nsupr, &x[ii], &knsupc);
#else
		dtrsm_("L", "L", "N", "U", &knsupc, &nrhs, &alpha, 
		       lusup, &nsupr, &x[ii], &knsupc);
#endif
		stat->ops[SOLVE] += knsupc * (knsupc - 1) * nrhs;
		--nleaf;
#if ( DEBUGlevel>=2 )
		printf("(%2d) Solve X[%2d]\n", iam, k);
#endif
		
		/*
		 * Send Xk to process column Pc[k].
		 */
		for (p = 0; p < Pr; ++p)
		    if ( fsendx_plist[lk][p] != EMPTY ) {
			pi = PNUM( p, kcol, grid );
#ifdef ISEND_IRECV
#if 1
			MPI_Test( &send_req[p], &test_flag, &status );
#else
			if ( send_req[p] != MPI_REQUEST_NULL ) 
			    MPI_Wait( &send_req[p], &status );
#endif
			MPI_Isend( &x[ii - XK_H], knsupc * nrhs + XK_H,
				 MPI_DOUBLE, pi, Xk,
                                 grid->comm, &send_req[p]);
#else
			MPI_Send( &x[ii - XK_H], knsupc * nrhs + XK_H,
				 MPI_DOUBLE, pi, Xk, grid->comm );
#endif
#if ( DEBUGlevel>=2 )
			printf("(%2d) Sent X[%2.0f] to P %2d\n",
			       iam, x[ii-XK_H], pi);
#endif
		    }
		
		/*
		 * Perform local block modifications: lsum[i] -= L_i,k * X[k]
		 */
		nb = lsub[0] - 1;
		lptr = BC_HEADER + LB_DESCRIPTOR + knsupc;
		luptr = knsupc; /* Skip diagonal block L(k,k). */
		
		dlsum_fmod(lsum, x, &x[ii], rtemp, nrhs, knsupc, k,
			   fmod, nb, lptr, luptr, xsup, grid, Llu, 
			   send_req, stat);
#ifdef ISEND_IRECV
		/* Wait for previous Isends to complete. */
		for (p = 0; p < Pr; ++p) {
		    if ( fsendx_plist[lk][p] != EMPTY )
			/*MPI_Wait( &send_req[p], &status );*/
			MPI_Test( &send_req[p], &test_flag, &status );
		}
#endif
	    }
	} /* if diagonal process ... */
    } /* for k ... */

    /* -----------------------------------------------------------
       Compute the internal nodes asynchronously by all processes.
       ----------------------------------------------------------- */
#if ( DEBUGlevel>=1 )
    printf("(%2d) nfrecvx %4d,  nfrecvmod %4d,  nleaf %4d\n",
	   iam, nfrecvx, nfrecvmod, nleaf);
#endif

    while ( nfrecvx || nfrecvmod ) { /* While not finished. */

	/* Receive a message. */
#ifdef ISEND_IRECV
	/* -MPI- FATAL: Remote protocol queue full */
	MPI_Irecv( recvbuf, maxrecvsz, MPI_DOUBLE,
                 MPI_ANY_SOURCE, MPI_ANY_TAG, grid->comm, &recv_req );
	MPI_Wait( &recv_req, &status );
#else
	MPI_Recv( recvbuf, maxrecvsz, MPI_DOUBLE,
                  MPI_ANY_SOURCE, MPI_ANY_TAG, grid->comm, &status );
#endif

	k = *recvbuf;

#if ( DEBUGlevel>=2 )
	printf("(%2d) Recv'd block %d, tag %2d\n", iam, k, status.MPI_TAG);
#endif
	
	switch ( status.MPI_TAG ) {
	  case Xk:
	      --nfrecvx;
	      lk = LBj( k, grid ); /* Local block number, column-wise. */
	      lsub = Lrowind_bc_ptr[lk];
	      lusup = Lnzval_bc_ptr[lk];
	      if ( lsub ) {
		  nb   = lsub[0];
		  lptr = BC_HEADER;
		  luptr = 0;
		  knsupc = SuperSize( k );

		  /*
		   * Perform local block modifications: lsum[i] -= L_i,k * X[k]
		   */
		  dlsum_fmod(lsum, x, &recvbuf[XK_H], rtemp, nrhs, knsupc, k,
			     fmod, nb, lptr, luptr, xsup, grid, Llu, 
			     send_req, stat);
	      } /* if lsub */

	      break;

	  case LSUM:
	      --nfrecvmod;
	      lk = LBi( k, grid ); /* Local block number, row-wise. */
	      ii = X_BLK( lk );
	      knsupc = SuperSize( k );
	      tempv = &recvbuf[LSUM_H];
	      RHS_ITERATE(j) {
		  for (i = 0; i < knsupc; ++i)
		      x[i + ii + j*knsupc] += tempv[i + j*knsupc];
	      }

	      if ( (--frecv[lk])==0 && fmod[lk]==0 ) {
		  fmod[lk] = -1; /* Do not solve X[k] in the future. */
		  lk = LBj( k, grid ); /* Local block number, column-wise. */
		  lsub = Lrowind_bc_ptr[lk];
		  lusup = Lnzval_bc_ptr[lk];
		  nsupr = lsub[1];
#ifdef _CRAY
		  STRSM(ftcs1, ftcs1, ftcs2, ftcs3, &knsupc, &nrhs, &alpha,
			lusup, &nsupr, &x[ii], &knsupc);
#else
		  dtrsm_("L", "L", "N", "U", &knsupc, &nrhs, &alpha, 
			 lusup, &nsupr, &x[ii], &knsupc);
#endif
		  stat->ops[SOLVE] += knsupc * (knsupc - 1) * nrhs;
#if ( DEBUGlevel>=2 )
		  printf("(%2d) Solve X[%2d]\n", iam, k);
#endif
		
		  /*
		   * Send Xk to process column Pc[k].
		   */
		  kcol = PCOL( k, grid );
		  for (p = 0; p < Pr; ++p)
		      if ( fsendx_plist[lk][p] != EMPTY ) {
			  pi = PNUM( p, kcol, grid );
#ifdef ISEND_IRECV
#if 1
			  MPI_Test( &send_req[p], &test_flag, &status );
#else
			  if ( send_req[p] != MPI_REQUEST_NULL )
			    MPI_Wait( &send_req[p], &status );
#endif
			  MPI_Isend( &x[ii-XK_H], knsupc * nrhs + XK_H,
                                     MPI_DOUBLE, pi, Xk,
                                     grid->comm, &send_req[p]);
#else
			  MPI_Send( &x[ii - XK_H], knsupc * nrhs + XK_H,
				    MPI_DOUBLE, pi, Xk, grid->comm );
#endif
#if ( DEBUGlevel>=2 )
			  printf("(%2d) Sent X[%2.0f] to P %2d\n",
				 iam, x[ii-XK_H], pi);
#endif
		      }

		  /*
		   * Perform local block modifications.
		   */
		  nb = lsub[0] - 1;
		  lptr = BC_HEADER + LB_DESCRIPTOR + knsupc;
		  luptr = knsupc; /* Skip diagonal block L(k,k). */

		  dlsum_fmod(lsum, x, &x[ii], rtemp, nrhs, knsupc, k,
			     fmod, nb, lptr, luptr, xsup, grid, Llu,
			     send_req, stat);
#ifdef ISEND_IRECV
		  /* Wait for the previous Isends to complete. */
		  for (p = 0; p < Pr; ++p) {
		      if ( fsendx_plist[lk][p] != EMPTY )
			  MPI_Test( &send_req[p], &test_flag, &status );
		  }
#endif
	      } /* if */

	      break;

#if ( DEBUGlevel>=1 )	      
	    default:
	      printf("(%2d) Recv'd wrong message tag %4d\n", status.MPI_TAG);
	      break;
#endif
	  } /* switch */

    } /* while not finished ... */


#if ( PRNTlevel>=2 )
    t = SuperLU_timer_() - t;
    if ( !iam ) printf(".. L-solve time\t%8.2f\n", t);
    t = SuperLU_timer_();
#endif

#if ( DEBUGlevel==2 )
    {
      printf("(%d) .. After L-solve: y =\n", iam);
      for (i = 0, k = 0; k < nsupers; ++k) {
	  krow = PROW( k, grid );
	  kcol = PCOL( k, grid );
	  if ( myrow == krow && mycol == kcol ) { /* Diagonal process */
	      knsupc = SuperSize( k );
	      lk = LBi( k, grid );
	      ii = X_BLK( lk );
	      for (j = 0; j < knsupc; ++j)
		printf("\t(%d)\t%4d\t%.10f\n", iam, xsup[k]+j, x[ii+j]);
	      fflush(stdout);
	  }
	  MPI_Barrier( grid->comm );
      }
    }
#endif

    SUPERLU_FREE(fmod);
    SUPERLU_FREE(frecv);
    SUPERLU_FREE(rtemp);

    /* MPI_Barrier( grid->comm );  Drain messages in the forward solve. */


    /*---------------------------------------------------
     * Back solve Ux = y.
     *
     * The Y components from the forward solve is already
     * on the diagonal processes.
     *---------------------------------------------------*/

    /* Save the count to be altered so it can be used by
       subsequent call to PDGSTRS. */
    if ( !(bmod = intMalloc_dist(nlb)) )
	ABORT("Calloc fails for bmod[].");
    for (i = 0; i < nlb; ++i) bmod[i] = Llu->bmod[i];
    if ( !(brecv = intMalloc_dist(nlb)) )
	ABORT("Malloc fails for brecv[].");
    Llu->brecv = brecv;

    /*
     * Compute brecv[] and nbrecvmod counts on the diagonal processes.
     */
    {
	superlu_scope_t *scp = &grid->rscp;

	for (k = 0; k < nsupers; ++k) {
	    krow = PROW( k, grid );
	    if ( myrow == krow ) {
		lk = LBi( k, grid );    /* Local block number. */
		kcol = PCOL( k, grid ); /* Root process in this row scope. */
		if ( mycol != kcol && bmod[lk] )
		    i = 1;  /* Contribution from non-diagonal process. */
		else i = 0;
		MPI_Reduce( &i, &brecv[lk], 1, mpi_int_t,
			   MPI_SUM, kcol, scp->comm );
		if ( mycol == kcol ) { /* Diagonal process. */
		    nbrecvmod += brecv[lk];
		    if ( !brecv[lk] && !bmod[lk] ) ++nroot;
#if ( DEBUGlevel>=2 )
		    printf("(%2d) brecv[%4d]  %2d\n", iam, k, brecv[lk]);
		    assert( brecv[lk] < Pc );
#endif
		}
	    }
	}
    }

    /* Re-initialize lsum to zero. Each block header is already in place. */
    for (k = 0; k < nsupers; ++k) {
	krow = PROW( k, grid );
	if ( myrow == krow ) {
	    knsupc = SuperSize( k );
	    lk = LBi( k, grid );
	    il = LSUM_BLK( lk );
	    dest = &lsum[il];
	    RHS_ITERATE(j) {
		for (i = 0; i < knsupc; ++i) dest[i + j*knsupc] = zero;
	    }
	}
    }

    /* Set up additional pointers for the index and value arrays of U.
       nlb is the number of local block rows. */
    nub = CEILING( nsupers, Pc ); /* Number of local block columns. */
    if ( !(Urbs = (int_t *) intCalloc_dist(2*nub)) )
	ABORT("Malloc fails for Urbs[]"); /* Record number of nonzero
					     blocks in a block column. */
    Urbs1 = Urbs + nub;
    if ( !(Ucb_indptr = SUPERLU_MALLOC(nub * sizeof(Ucb_indptr_t *))) )
        ABORT("Malloc fails for Ucb_indptr[]");
    if ( !(Ucb_valptr = SUPERLU_MALLOC(nub * sizeof(int_t *))) )
        ABORT("Malloc fails for Ucb_valptr[]");

    /* Count number of row blocks in a block column. 
       One pass of the skeleton graph of U. */
    for (lk = 0; lk < nlb; ++lk) {
	usub = Ufstnz_br_ptr[lk];
	if ( usub ) { /* Not an empty block row. */
	    /* usub[0] -- number of column blocks in this block row. */
#if ( DEBUGlevel>=2 )
	    Ublocks += usub[0];
#endif
	    i = BR_HEADER; /* Pointer in index array. */
	    for (lb = 0; lb < usub[0]; ++lb) { /* For all column blocks. */
		k = usub[i];            /* Global block number */
		++Urbs[LBj(k,grid)];
		i += UB_DESCRIPTOR + SuperSize( k );
	    }
	}
    }

    /* Set up the vertical linked lists for the row blocks.
       One pass of the skeleton graph of U. */
    for (lb = 0; lb < nub; ++lb)
	if ( Urbs[lb] ) { /* Not an empty block column. */
	    if ( !(Ucb_indptr[lb]
		   = SUPERLU_MALLOC(Urbs[lb] * sizeof(Ucb_indptr_t))) )
		ABORT("Malloc fails for Ucb_indptr[lb][]");
	    if ( !(Ucb_valptr[lb] = (int_t *) intMalloc_dist(Urbs[lb])) )
		ABORT("Malloc fails for Ucb_valptr[lb][]");
	}
    for (lk = 0; lk < nlb; ++lk) { /* For each block row. */
	usub = Ufstnz_br_ptr[lk];
	if ( usub ) { /* Not an empty block row. */
	    i = BR_HEADER; /* Pointer in index array. */
	    j = 0;         /* Pointer in nzval array. */
	    for (lb = 0; lb < usub[0]; ++lb) { /* For all column blocks. */
		k = usub[i];          /* Global block number, column-wise. */
		ljb = LBj( k, grid ); /* Local block number, column-wise. */
		Ucb_indptr[ljb][Urbs1[ljb]].lbnum = lk;
		Ucb_indptr[ljb][Urbs1[ljb]].indpos = i;
		Ucb_valptr[ljb][Urbs1[ljb]] = j;
		++Urbs1[ljb];
		j += usub[i+1];
		i += UB_DESCRIPTOR + SuperSize( k );
	    }
	}
    }

#if ( DEBUGlevel>=2 )
    for (p = 0; p < Pr*Pc; ++p) {
	if (iam == p) {
	    printf("(%2d) .. Ublocks %d\n", iam, Ublocks);
	    for (lb = 0; lb < nub; ++lb) {
		printf("(%2d) Local col %2d: # row blocks %2d\n",
		       iam, lb, Urbs[lb]);
		if ( Urbs[lb] ) {
		    for (i = 0; i < Urbs[lb]; ++i)
			printf("(%2d) .. row blk %2d:\
                               lbnum %d, indpos %d, valpos %d\n",
			       iam, i, 
			       Ucb_indptr[lb][i].lbnum,
			       Ucb_indptr[lb][i].indpos,
			       Ucb_valptr[lb][i]);
		}
	    }
	}
	MPI_Barrier( grid->comm );
    }
    for (p = 0; p < Pr*Pc; ++p) {
	if ( iam == p ) {
	    printf("\n(%d) bsendx_plist[][]", iam);
	    for (lb = 0; lb < nub; ++lb) {
		printf("\n(%d) .. local col %2d: ", iam, lb);
		for (i = 0; i < Pr; ++i)
		    printf("%4d", bsendx_plist[lb][i]);
	    }
	    printf("\n");
	}
	MPI_Barrier( grid->comm );
    }
#endif /* DEBUGlevel */


#if ( PRNTlevel>=3 )
    t = SuperLU_timer_() - t;
    if ( !iam) printf(".. Setup U-solve time\t%8.2f\n", t);
    t = SuperLU_timer_();
#endif

    /*
     * Solve the roots first by all the diagonal processes.
     */
#if ( DEBUGlevel>=1 )
    printf("(%2d) nroot %4d\n", iam, nroot);
#endif
    for (k = nsupers-1; k >= 0 && nroot; --k) {
	krow = PROW( k, grid );
	kcol = PCOL( k, grid );
	if ( myrow == krow && mycol == kcol ) { /* Diagonal process. */
	    knsupc = SuperSize( k );
	    lk = LBi( k, grid ); /* Local block number, row-wise. */
	    if ( brecv[lk]==0 && bmod[lk]==0 ) {
		bmod[lk] = -1;       /* Do not solve X[k] in the future. */
		ii = X_BLK( lk );
		lk = LBj( k, grid ); /* Local block number, column-wise */
		lsub = Lrowind_bc_ptr[lk];
		lusup = Lnzval_bc_ptr[lk];
		nsupr = lsub[1];
#ifdef _CRAY
		STRSM(ftcs1, ftcs3, ftcs2, ftcs2, &knsupc, &nrhs, &alpha,
		      lusup, &nsupr, &x[ii], &knsupc);
#else
		dtrsm_("L", "U", "N", "N", &knsupc, &nrhs, &alpha, 
		       lusup, &nsupr, &x[ii], &knsupc);
#endif
		stat->ops[SOLVE] += knsupc * (knsupc + 1) * nrhs;
		--nroot;
#if ( DEBUGlevel>=2 )
		printf("(%2d) Solve X[%2d]\n", iam, k);
#endif
		/*
		 * Send Xk to process column Pc[k].
		 */
		for (p = 0; p < Pr; ++p)
		    if ( bsendx_plist[lk][p] != EMPTY ) {
			pi = PNUM( p, kcol, grid );
#ifdef ISEND_IRECV
#if 1
			MPI_Test( &send_req[p], &test_flag, &status );
#else
			if ( send_req[p] != MPI_REQUEST_NULL )
			  MPI_Wait( &send_req[p], &status );
#endif
			MPI_Isend( &x[ii - XK_H], knsupc * nrhs + XK_H,
                                   MPI_DOUBLE, pi, Xk,
                                   grid->comm, &send_req[p]);
#else
			MPI_Send( &x[ii - XK_H], knsupc * nrhs + XK_H,
                                  MPI_DOUBLE, pi, Xk,
                                  grid->comm );
#endif
#if ( DEBUGlevel>=2 )
			printf("(%2d) Sent X[%2.0f] to P %2d\n",
			       iam, x[ii-XK_H], pi);
#endif
		    }
		
		/*
		 * Perform local block modifications: lsum[i] -= U_i,k * X[k]
		 */
		if ( Urbs[lk] ) 
		    dlsum_bmod(lsum, x, &x[ii], nrhs, k, bmod, Urbs,
			       Ucb_indptr, Ucb_valptr, xsup, grid, Llu,
			       send_req, stat);
#ifdef ISEND_IRECV
		/* Wait for the previous Isends to complete. */
		for (p = 0; p < Pr; ++p) {
		    if ( bsendx_plist[lk][p] != EMPTY )
			MPI_Test( &send_req[p], &test_flag, &status );
		}
#endif
	    } /* if root ... */
	} /* if diagonal process ... */
    } /* for k ... */


    /*
     * Compute the internal nodes asychronously by all processes.
     */
    while ( nbrecvx || nbrecvmod ) { /* While not finished. */

	/* Receive a message. */
	MPI_Recv( recvbuf, maxrecvsz, MPI_DOUBLE,
                  MPI_ANY_SOURCE, MPI_ANY_TAG, grid->comm, &status );
	k = *recvbuf;

#if ( DEBUGlevel>=2 )
	printf("(%2d) Recv'd block %d, tag %2d\n", iam, k, status.MPI_TAG);
#endif

	switch ( status.MPI_TAG ) {
	    case Xk:
	        --nbrecvx;
		lk = LBj( k, grid ); /* Local block number, column-wise. */
		/*
		 * Perform local block modifications:
		 *         lsum[i] -= U_i,k * X[k]
		 */
		dlsum_bmod(lsum, x, &recvbuf[XK_H], nrhs, k, bmod, Urbs,
			   Ucb_indptr, Ucb_valptr, xsup, grid, Llu, 
			   send_req, stat);

	        break;

	    case LSUM:
		--nbrecvmod;
		lk = LBi( k, grid ); /* Local block number, row-wise. */
		ii = X_BLK( lk );
		knsupc = SuperSize( k );
		tempv = &recvbuf[LSUM_H];
		RHS_ITERATE(j) {
		    for (i = 0; i < knsupc; ++i)
			x[i + ii + j*knsupc] += tempv[i + j*knsupc];
		}

		if ( (--brecv[lk])==0 && bmod[lk]==0 ) {
		    bmod[lk] = -1; /* Do not solve X[k] in the future. */
		    lk = LBj( k, grid ); /* Local block number, column-wise. */
		    lsub = Lrowind_bc_ptr[lk];
		    lusup = Lnzval_bc_ptr[lk];
		    nsupr = lsub[1];
#ifdef _CRAY
		    STRSM(ftcs1, ftcs3, ftcs2, ftcs2, &knsupc, &nrhs, &alpha,
			  lusup, &nsupr, &x[ii], &knsupc);
#else
		    dtrsm_("L", "U", "N", "N", &knsupc, &nrhs, &alpha, 
			   lusup, &nsupr, &x[ii], &knsupc);
#endif
		    stat->ops[SOLVE] += knsupc * (knsupc + 1) * nrhs;
#if ( DEBUGlevel>=2 )
		    printf("(%2d) Solve X[%2d]\n", iam, k);
#endif
		    /*
		     * Send Xk to process column Pc[k].
		     */
		    kcol = PCOL( k, grid );
		    for (p = 0; p < Pr; ++p)
			if ( bsendx_plist[lk][p] != EMPTY ) {
			    pi = PNUM( p, kcol, grid );
#ifdef ISEND_IRECV
#if 1
			    MPI_Test( &send_req[p], &test_flag, &status );
#else
			    if ( send_req[p] != MPI_REQUEST_NULL )
			        MPI_Wait( &send_req[p], &status );
#endif
			    MPI_Isend( &x[ii - XK_H], knsupc * nrhs + XK_H,
                                       MPI_DOUBLE, pi, Xk,
                                       grid->comm, &send_req[p] );
#else
			    MPI_Send( &x[ii - XK_H], knsupc * nrhs + XK_H,
                                      MPI_DOUBLE, pi, Xk,
                                      grid->comm );
#endif
#if ( DEBUGlevel>=2 )
			    printf("(%2d) Sent X[%2.0f] to P %2d\n",
				   iam, x[ii - XK_H], pi);
#endif
			}
		
		    /*
		     * Perform local block modifications: 
		     *         lsum[i] -= U_i,k * X[k]
		     */
		    if ( Urbs[lk] )
			dlsum_bmod(lsum, x, &x[ii], nrhs, k, bmod, Urbs,
				   Ucb_indptr, Ucb_valptr, xsup, grid, Llu,
				   send_req, stat);
#ifdef ISEND_IRECV
		    /* Wait for the previous Isends to complete. */
		    for (p = 0; p < Pr; ++p) {
			if ( bsendx_plist[lk][p] != EMPTY )
			    /*MPI_Wait( &send_req[p], &status );*/
			    MPI_Test( &send_req[p], &test_flag, &status );
		    }
#endif
		} /* if becomes solvable */
		
		break;

#if ( DEBUGlevel>=1 )
	      default:
		printf("(%2d) Recv'd wrong message tag %4d\n", status.MPI_TAG);
		break;
#endif		

	} /* switch */

    } /* while not finished ... */

#if ( PRNTlevel>=3 )
    t = SuperLU_timer_() - t;
    if ( !iam ) printf(".. U-solve time\t%8.2f\n", t);
#endif

#if ( DEBUGlevel>=2 )
    {
	double *x_col;
	int diag;
	printf("\n(%d) .. After U-solve: x (ON DIAG PROCS) = \n", iam);
	ii = 0;
	for (k = 0; k < nsupers; ++k) {
	    knsupc = SuperSize( k );
	    krow = PROW( k, grid );
	    kcol = PCOL( k, grid );
	    diag = PNUM( krow, kcol, grid);
	    if ( iam == diag ) { /* Diagonal process. */
		lk = LBi( k, grid );
		jj = X_BLK( lk );
		x_col = &x[jj];
		RHS_ITERATE(j) {
		    for (i = 0; i < knsupc; ++i) { /* X stored in blocks */
			printf("\t(%d)\t%4d\t%.10f\n",
			       iam, xsup[k]+i, x_col[i]);
		    }
		    x_col += knsupc;
		}
	    }
	    ii += knsupc;
	} /* for k ... */
    }
#endif

    pdReDistribute_X_to_B(n, B, m_loc, ldb, fst_row, nrhs, x, ilsum,
			  ScalePermstruct, Glu_persist, grid, SOLVEstruct);


    /* Deallocate storage. */
    SUPERLU_FREE(lsum);
    SUPERLU_FREE(x);
    SUPERLU_FREE(recvbuf);
    for (i = 0; i < nub; ++i) {
	if ( Urbs[i] ) {
	    SUPERLU_FREE(Ucb_indptr[i]);
	    SUPERLU_FREE(Ucb_valptr[i]);
	}
    }
    SUPERLU_FREE(Ucb_indptr);
    SUPERLU_FREE(Ucb_valptr);
    SUPERLU_FREE(Urbs);
    SUPERLU_FREE(bmod);
    SUPERLU_FREE(brecv);
#ifdef ISEND_IRECV
    for (p = 0; p < Pr; ++p) {
        if ( send_req[p] != MPI_REQUEST_NULL )
	    MPI_Wait( &send_req[p], &status );
    }
    SUPERLU_FREE(send_req);
#endif

    stat->utime[SOLVE] = SuperLU_timer_() - t;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit pdgstrs()");
#endif

} /* PDGSTRS */

