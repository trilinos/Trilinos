

/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 */

#include <math.h>
#include "superlu_ddefs.h"

/*
 * Gather A from the distributed compressed row format to
 * global A in compressed column format.
 */
int pdCompRow_loc_to_CompCol_global
(
 int_t need_value, /* Input. Whether need to gather numerical values */
 SuperMatrix *A,   /* Input. Distributed matrix in NRformat_loc format. */
 gridinfo_t *grid, /* Input */
 SuperMatrix *GA   /* Output */
)
{
    NRformat_loc *Astore;
    NCformat *GAstore;
    double *a, *a_loc;
    int_t *colind, *rowptr;
    int_t *colptr_loc, *rowind_loc;
    int_t m_loc, n, i, j, k, l;
    int_t colnnz, fst_row, m_loc_max, nnz_loc, nnz_max, nnz;
    double *a_recv;  /* Buffer to receive the blocks of values. */
    double *a_buf;   /* Buffer to merge blocks into block columns. */
    int_t *colcnt, *itemp;
    int_t *colptr_send; /* Buffer to redistribute the column pointers of the 
			   local block rows.
			   Use n_loc+1 pointers for each block. */
    int_t *colptr_blk;  /* The column pointers for each block, after
			   redistribution to the local block columns. 
			   Use n_loc+1 pointers for each block. */
    int_t *rowind_recv; /* Buffer to receive the blocks of row indices. */
    int_t *rowind_buf;  /* Buffer to merge blocks into block columns. */
    int_t *fst_rows, *n_locs;
    int   *sendcnts, *sdispls, *recvcnts, *rdispls, *itemp_32;
    int   it, n_loc, procs;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Enter pdCompRow_loc_to_CompCol_global");
#endif

    /* Initialization. */
    n = A->ncol;
    Astore = (NRformat_loc *) A->Store;
    nnz_loc = Astore->nnz_loc;
    m_loc = Astore->m_loc;
    fst_row = Astore->fst_row;
    a = Astore->nzval;
    rowptr = Astore->rowptr;
    colind = Astore->colind;
    n_loc = m_loc; /* NOTE: CURRENTLY ONLY WORK FOR SQUARE MATRIX */

    /* ------------------------------------------------------------
       FIRST PHASE: TRANSFORM A INTO DISTRIBUTED COMPRESSED COLUMN.
       ------------------------------------------------------------*/
    dCompRow_to_CompCol_dist(m_loc, n, nnz_loc, a, colind, rowptr, &a_loc,
                             &rowind_loc, &colptr_loc);
    /* Change local row index numbers to global numbers. */
    for (i = 0; i < nnz_loc; ++i) rowind_loc[i] += fst_row;

#if ( DEBUGlevel>=2 )
    printf("Proc %d\n", grid->iam);
    PrintInt10("rowind_loc", nnz_loc, rowind_loc);
    PrintInt10("colptr_loc", n+1, colptr_loc);
#endif

    procs = grid->nprow * grid->npcol;
    if ( !(fst_rows = (int_t *) intMalloc_dist(2*procs)) )
	  ABORT("Malloc fails for fst_rows[]");
    n_locs = fst_rows + procs;
    MPI_Allgather(&fst_row, 1, mpi_int_t, fst_rows, 1, mpi_int_t,
		  grid->comm);
    for (i = 0; i < procs-1; ++i) n_locs[i] = fst_rows[i+1] - fst_rows[i];
    n_locs[procs-1] = n - fst_rows[procs-1];
    if ( !(recvcnts = SUPERLU_MALLOC(5*procs * sizeof(int))) )
	  ABORT("Malloc fails for recvcnts[]");
    sendcnts = recvcnts + procs;
    rdispls = sendcnts + procs;
    sdispls = rdispls + procs;
    itemp_32 = sdispls + procs;

    /* All-to-all transfer column pointers of each block.
       Now the matrix view is P-by-P block-partition. */
    /* n column starts for each column, and procs column ends for each block */
    if ( !(colptr_send = intMalloc_dist(n + procs)) )
	   ABORT("Malloc fails for colptr_send[]");
    if ( !(colptr_blk = intMalloc_dist( (((size_t) n_loc)+1)*procs)) )
	   ABORT("Malloc fails for colptr_blk[]");
    for (i = 0, j = 0; i < procs; ++i) {
        for (k = j; k < j + n_locs[i]; ++k) colptr_send[i+k] = colptr_loc[k];
	colptr_send[i+k] = colptr_loc[k]; /* Add an END marker */
	sendcnts[i] = n_locs[i] + 1;
	assert(j == fst_rows[i]);
	sdispls[i] = j + i;
	recvcnts[i] = n_loc + 1;
	rdispls[i] = i * (n_loc + 1);
	j += n_locs[i]; /* First column of next block in colptr_loc[] */
    }
    MPI_Alltoallv(colptr_send, sendcnts, sdispls, mpi_int_t,
		  colptr_blk, recvcnts, rdispls, mpi_int_t, grid->comm);

    /* Adjust colptr_blk[] so that they contain the local indices of the
       column pointers in the receive buffer. */
    nnz = 0; /* The running sum of the nonzeros counted by far */
    k = 0;
    for (i = 0; i < procs; ++i) {
	for (j = rdispls[i]; j < rdispls[i] + n_loc; ++j) {
	    colnnz = colptr_blk[j+1] - colptr_blk[j];
	    /*assert(k<=j);*/
	    colptr_blk[k] = nnz;
	    nnz += colnnz; /* Start of the next column */
	    ++k;
	}
	colptr_blk[k++] = nnz; /* Add an END marker for each block */
    }
    /*assert(k == (n_loc+1)*procs);*/

    /* Now prepare to transfer row indices and values. */
    sdispls[0] = 0;
    for (i = 0; i < procs-1; ++i) {
        sendcnts[i] = colptr_loc[fst_rows[i+1]] - colptr_loc[fst_rows[i]];
	sdispls[i+1] = sdispls[i] + sendcnts[i];
    }
    sendcnts[procs-1] = colptr_loc[n] - colptr_loc[fst_rows[procs-1]];
    for (i = 0; i < procs; ++i) {
        j = rdispls[i]; /* Point to this block in colptr_blk[]. */
	recvcnts[i] = colptr_blk[j+n_loc] - colptr_blk[j];
    }
    rdispls[0] = 0; /* Recompute rdispls[] for row indices. */
    for (i = 0; i < procs-1; ++i) rdispls[i+1] = rdispls[i] + recvcnts[i];

    k = rdispls[procs-1] + recvcnts[procs-1]; /* Total received */
    if ( !(rowind_recv = (int_t *) intMalloc_dist(2*k)) )
        ABORT("Malloc fails for rowind_recv[]");
    rowind_buf = rowind_recv + k;
    MPI_Alltoallv(rowind_loc, sendcnts, sdispls, mpi_int_t,
		  rowind_recv, recvcnts, rdispls, mpi_int_t, grid->comm);
    if ( need_value ) {
        if ( !(a_recv = (double *) doubleMalloc_dist(2*k)) )
	    ABORT("Malloc fails for rowind_recv[]");
	a_buf = a_recv + k;
	MPI_Alltoallv(a_loc, sendcnts, sdispls, MPI_DOUBLE,
                      a_recv, recvcnts, rdispls, MPI_DOUBLE,
                      grid->comm);
    }
      
    /* Reset colptr_loc[] to point to the n_loc global columns. */
    colptr_loc[0] = 0;
    itemp = colptr_send;
    for (j = 0; j < n_loc; ++j) {
        colnnz = 0;
	for (i = 0; i < procs; ++i) {
	    k = i * (n_loc + 1) + j; /* j-th column in i-th block */
	    colnnz += colptr_blk[k+1] - colptr_blk[k];
	}
	colptr_loc[j+1] = colptr_loc[j] + colnnz;
	itemp[j] = colptr_loc[j]; /* Save a copy of the column starts */
    }
    itemp[n_loc] = colptr_loc[n_loc];
      
    /* Merge blocks of row indices into columns of row indices. */
    for (i = 0; i < procs; ++i) {
        k = i * (n_loc + 1);
	for (j = 0; j < n_loc; ++j) { /* i-th block */
	    for (l = colptr_blk[k+j]; l < colptr_blk[k+j+1]; ++l) {
	        rowind_buf[itemp[j]] = rowind_recv[l];
		++itemp[j];
	    }
	}
    }

    if ( need_value ) {
        for (j = 0; j < n_loc+1; ++j) itemp[j] = colptr_loc[j];
        for (i = 0; i < procs; ++i) {
	    k = i * (n_loc + 1);
	    for (j = 0; j < n_loc; ++j) { /* i-th block */
	        for (l = colptr_blk[k+j]; l < colptr_blk[k+j+1]; ++l) {
		    a_buf[itemp[j]] = a_recv[l];
		    ++itemp[j];
		}
	    }
	}
    }

    /* ------------------------------------------------------------
       SECOND PHASE: GATHER TO GLOBAL A IN COMPRESSED COLUMN FORMAT.
       ------------------------------------------------------------*/
    GA->nrow  = A->nrow;
    GA->ncol  = A->ncol;
    GA->Stype = SLU_NC;
    GA->Dtype = A->Dtype;
    GA->Mtype = A->Mtype;
    GAstore = GA->Store = (NCformat *) SUPERLU_MALLOC ( sizeof(NCformat) );
    if ( !GAstore ) ABORT ("SUPERLU_MALLOC fails for GAstore");

    /* First gather the size of each piece. */
    nnz_loc = colptr_loc[n_loc];
    MPI_Allgather(&nnz_loc, 1, mpi_int_t, itemp, 1, mpi_int_t, grid->comm);
    for (i = 0, nnz = 0; i < procs; ++i) nnz += itemp[i];
    GAstore->nnz = nnz;
    
    if ( !(GAstore->rowind = (int_t *) intMalloc_dist (nnz)) )
        ABORT ("SUPERLU_MALLOC fails for GAstore->rowind[]");
    if ( !(GAstore->colptr = (int_t *) intMalloc_dist (n+1)) )
        ABORT ("SUPERLU_MALLOC fails for GAstore->colptr[]");
      
    /* Allgatherv for row indices. */
    rdispls[0] = 0;
    for (i = 0; i < procs-1; ++i) {
        rdispls[i+1] = rdispls[i] + itemp[i];
        itemp_32[i] = itemp[i];
    }
    itemp_32[procs-1] = itemp[procs-1];
    it = nnz_loc;
    MPI_Allgatherv(rowind_buf, it, mpi_int_t, GAstore->rowind, 
		   itemp_32, rdispls, mpi_int_t, grid->comm);
    if ( need_value ) {
      if ( !(GAstore->nzval = (double *) doubleMalloc_dist (nnz)) )
          ABORT ("SUPERLU_MALLOC fails for GAstore->rnzval[]");
      MPI_Allgatherv(a_buf, it, MPI_DOUBLE, GAstore->nzval, 
		     itemp_32, rdispls, MPI_DOUBLE, grid->comm);
    } else GAstore->nzval = NULL;

    /* Now gather the column pointers. */
    rdispls[0] = 0;
    for (i = 0; i < procs-1; ++i) {
        rdispls[i+1] = rdispls[i] + n_locs[i];
        itemp_32[i] = n_locs[i];
    }
    itemp_32[procs-1] = n_locs[procs-1];
    MPI_Allgatherv(colptr_loc, n_loc, mpi_int_t, GAstore->colptr, 
		   itemp_32, rdispls, mpi_int_t, grid->comm);

    /* Recompute column pointers. */
    for (i = 1; i < procs; ++i) {
        k = rdispls[i];
	for (j = 0; j < n_locs[i]; ++j) GAstore->colptr[k++] += itemp[i-1];
	itemp[i] += itemp[i-1]; /* prefix sum */
    }
    GAstore->colptr[n] = nnz;

#if ( DEBUGlevel>=2 )
    if ( !grid->iam ) {
        printf("After pdCompRow_loc_to_CompCol_global()\n");
	dPrint_CompCol_Matrix_dist(GA);
    }
#endif

    SUPERLU_FREE(a_loc);
    SUPERLU_FREE(rowind_loc);
    SUPERLU_FREE(colptr_loc);
    SUPERLU_FREE(fst_rows);
    SUPERLU_FREE(recvcnts);
    SUPERLU_FREE(colptr_send);
    SUPERLU_FREE(colptr_blk);
    SUPERLU_FREE(rowind_recv);
    if ( need_value) SUPERLU_FREE(a_recv);
#if ( DEBUGlevel>=1 )
    if ( !grid->iam ) printf("sizeof(NCformat) %d\n", sizeof(NCformat));
    CHECK_MALLOC(grid->iam, "Exit pdCompRow_loc_to_CompCol_global");
#endif
    return 0;
} /* pdCompRow_loc_to_CompCol_global */


/*
 * Permute the distributed dense matrix: B <= perm(X).
 * perm[i] = j means the i-th row of X is in the j-th row of B.
 */
int pdPermute_Dense_Matrix
(
 int_t fst_row,
 int_t m_loc,
 int_t row_to_proc[],
 int_t perm[],
 double X[], int ldx,
 double B[], int ldb,
 int nrhs,
 gridinfo_t *grid
)
{
    int_t i, j, k, l;
    int p, procs;
    int *sendcnts, *sendcnts_nrhs, *recvcnts, *recvcnts_nrhs;
    int *sdispls, *sdispls_nrhs, *rdispls, *rdispls_nrhs;
    int *ptr_to_ibuf, *ptr_to_dbuf;
    int_t *send_ibuf, *recv_ibuf;
    double *send_dbuf, *recv_dbuf;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Enter pdPermute_Dense_Matrix()");
#endif

    procs = grid->nprow * grid->npcol;
    if ( !(sendcnts = SUPERLU_MALLOC(10*procs * sizeof(int))) )
        ABORT("Malloc fails for sendcnts[].");
    sendcnts_nrhs = sendcnts + procs;
    recvcnts = sendcnts_nrhs + procs;
    recvcnts_nrhs = recvcnts + procs;
    sdispls = recvcnts_nrhs + procs;
    sdispls_nrhs = sdispls + procs;
    rdispls = sdispls_nrhs + procs;
    rdispls_nrhs = rdispls + procs;
    ptr_to_ibuf = rdispls_nrhs + procs;
    ptr_to_dbuf = ptr_to_ibuf + procs;

    for (i = 0; i < procs; ++i) sendcnts[i] = 0;

    /* Count the number of X entries to be sent to each process.*/
    for (i = fst_row; i < fst_row + m_loc; ++i) {
        p = row_to_proc[perm[i]];
	++sendcnts[p];
    }
    MPI_Alltoall(sendcnts, 1, MPI_INT, recvcnts, 1, MPI_INT, grid->comm);
    sdispls[0] = rdispls[0] = 0;
    sdispls_nrhs[0] = rdispls_nrhs[0] = 0;
    sendcnts_nrhs[0] = sendcnts[0] * nrhs;
    recvcnts_nrhs[0] = recvcnts[0] * nrhs;
    for (i = 1; i < procs; ++i) {
        sdispls[i] = sdispls[i-1] + sendcnts[i-1];
	sdispls_nrhs[i] = sdispls[i] * nrhs;
	rdispls[i] = rdispls[i-1] + recvcnts[i-1];
	rdispls_nrhs[i] = rdispls[i] * nrhs;
	sendcnts_nrhs[i] = sendcnts[i] * nrhs;
	recvcnts_nrhs[i] = recvcnts[i] * nrhs;
    }
    k = sdispls[procs-1] + sendcnts[procs-1];/* Total number of sends */
    l = rdispls[procs-1] + recvcnts[procs-1];/* Total number of recvs */
    /*assert(k == m_loc);*/
    /*assert(l == m_loc);*/
    if ( !(send_ibuf = intMalloc_dist(k + l)) )
        ABORT("Malloc fails for send_ibuf[].");
    recv_ibuf = send_ibuf + k;
    if ( !(send_dbuf = doubleMalloc_dist((k + l)*nrhs)) )
        ABORT("Malloc fails for send_dbuf[].");
    recv_dbuf = send_dbuf + k * nrhs;

    for (i = 0; i < procs; ++i) {
        ptr_to_ibuf[i] = sdispls[i];
	ptr_to_dbuf[i] = sdispls_nrhs[i];
    }

    /* Fill in the send buffers: send_ibuf[] and send_dbuf[]. */
    for (i = fst_row; i < fst_row + m_loc; ++i) {
        j = perm[i];
	p = row_to_proc[j];
	send_ibuf[ptr_to_ibuf[p]] = j;
	j = ptr_to_dbuf[p];
	RHS_ITERATE(k) { /* RHS stored in row major in the buffer */
	    send_dbuf[j++] = X[i-fst_row + k*ldx];
	}
	++ptr_to_ibuf[p];
	ptr_to_dbuf[p] += nrhs;
    }
	  
    /* Transfer the (permuted) row indices and numerical values. */
    MPI_Alltoallv(send_ibuf, sendcnts, sdispls, mpi_int_t,
		  recv_ibuf, recvcnts, rdispls, mpi_int_t, grid->comm);
    MPI_Alltoallv(send_dbuf, sendcnts_nrhs, sdispls_nrhs, MPI_DOUBLE,
		  recv_dbuf, recvcnts_nrhs, rdispls_nrhs, MPI_DOUBLE,
		  grid->comm);

    /* Copy the buffer into b. */
    for (i = 0, l = 0; i < m_loc; ++i) {
        j = recv_ibuf[i] - fst_row; /* Relative row number */
	RHS_ITERATE(k) { /* RHS stored in row major in the buffer */
	    B[j + k*ldb] = recv_dbuf[l++];
	}
    }

    SUPERLU_FREE(sendcnts);
    SUPERLU_FREE(send_ibuf);
    SUPERLU_FREE(send_dbuf);
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Exit pdPermute_Dense_Matrix()");
#endif
    return 0;
} /* pdPermute_Dense_Matrix */


/*
 * Initialize the data structure for the solution phase.
 */
int dSolveInit(superlu_options_t *options, SuperMatrix *A, 
	       int_t perm_r[], int_t perm_c[], int_t nrhs,
	       LUstruct_t *LUstruct, gridinfo_t *grid,
	       SOLVEstruct_t *SOLVEstruct)
{
    int_t *row_to_proc, *inv_perm_c, *itemp;
    NRformat_loc *Astore;
    int_t        i, fst_row, m_loc, p;
    SuperMatrix  A_loc;
    int          procs;

    Astore = (NRformat_loc *) A->Store;
    fst_row = Astore->fst_row;
    m_loc = Astore->m_loc;
    procs = grid->nprow * grid->npcol;
    
    if ( !(row_to_proc = intMalloc_dist(A->nrow)) )
	ABORT("Malloc fails for row_to_proc[]");
    SOLVEstruct->row_to_proc = row_to_proc;
    if ( !(inv_perm_c = intMalloc_dist(A->ncol)) )
        ABORT("Malloc fails for inv_perm_c[].");
    for (i = 0; i < A->ncol; ++i) inv_perm_c[perm_c[i]] = i;
    SOLVEstruct->inv_perm_c = inv_perm_c;

    /* ------------------------------------------------------------
       EVERY PROCESS NEEDS TO KNOW GLOBAL PARTITION.
       SET UP THE MAPPING BETWEEN ROWS AND PROCESSES.
       
       NOTE: For those processes that do not own any row, it must
             must be set so that fst_row == A->nrow. 
       ------------------------------------------------------------*/
    if ( !(itemp = intMalloc_dist(procs+1)) )
        ABORT("Malloc fails for itemp[]");
    MPI_Allgather(&fst_row, 1, mpi_int_t, itemp, 1, mpi_int_t,
		  grid->comm);
    itemp[procs] = A->nrow;
    for (p = 0; p < procs; ++p) {
        for (i = itemp[p] ; i < itemp[p+1]; ++i) row_to_proc[i] = p;
    }
#if ( DEBUGlevel>=1 )
    if ( !grid->iam ) {
      printf("fst_row = %d\n", fst_row);
      PrintInt10("row_to_proc", A->nrow, row_to_proc);
      PrintInt10("inv_perm_c", A->ncol, inv_perm_c);
    }
#endif
    SUPERLU_FREE(itemp);

#if 0
    /* Compute the mapping between rows and processes. */
    /* XSL NOTE: What happens if # of mapped processes is smaller
       than total Procs?  For the processes without any row, let
       fst_row be EMPTY (-1). Make sure this case works! */
    MPI_Allgather(&fst_row, 1, mpi_int_t, itemp, 1, mpi_int_t,
		  grid->comm);
    itemp[procs] = n;
    for (p = 0; p < procs; ++p) {
        j = itemp[p];
	if ( j != EMPTY ) {
	    k = itemp[p+1];
	    if ( k == EMPTY ) k = n;
	    for (i = j ; i < k; ++i) row_to_proc[i] = p;
	}
    }
#endif    

    get_diag_procs(A->ncol, LUstruct->Glu_persist, grid,
		   &SOLVEstruct->num_diag_procs,
		   &SOLVEstruct->diag_procs,
		   &SOLVEstruct->diag_len);

    if ( !(SOLVEstruct->gstrs_comm = (pxgstrs_comm_t *)
	   SUPERLU_MALLOC(sizeof(pxgstrs_comm_t))) )
        ABORT("Malloc fails for gstrs_comm[]");
    pxgstrs_init(A->ncol, m_loc, nrhs, fst_row, perm_r, perm_c, grid, 
		 LUstruct->Glu_persist, SOLVEstruct);

    if ( !(SOLVEstruct->gsmv_comm = (pdgsmv_comm_t *)
           SUPERLU_MALLOC(sizeof(pdgsmv_comm_t))) )
        ABORT("Malloc fails for gsmv_comm[]");
    SOLVEstruct->A_colind_gsmv = NULL;
    
    options->SolveInitialized = YES;
    return 0;
} /* dSolveInit */

/*
 * Release the resources used for the solution phase.
 */
void dSolveFinalize(superlu_options_t *options, SOLVEstruct_t *SOLVEstruct)
{
    int_t *it;
    pxgstrs_finalize(SOLVEstruct->gstrs_comm);
    SUPERLU_FREE(SOLVEstruct->gstrs_comm);
    if ( options->RefineInitialized ) {
        pdgsmv_finalize(SOLVEstruct->gsmv_comm);
	options->RefineInitialized = NO;
    }
    SUPERLU_FREE(SOLVEstruct->gsmv_comm);
    SUPERLU_FREE(SOLVEstruct->row_to_proc);
    SUPERLU_FREE(SOLVEstruct->inv_perm_c);
    SUPERLU_FREE(SOLVEstruct->diag_procs);
    SUPERLU_FREE(SOLVEstruct->diag_len);
    if ( it = SOLVEstruct->A_colind_gsmv ) SUPERLU_FREE(it);
    options->SolveInitialized = NO;
} /* dSolveFinalize */

/* 
 * Check the inf-norm of the error vector 
 */
void pdinf_norm_error(int iam, int_t n, int_t nrhs, double x[], int_t ldx,
		      double xtrue[], int_t ldxtrue, gridinfo_t *grid) 
{
    double err, xnorm, temperr, tempxnorm;
    double *x_work, *xtrue_work;
    int i, j;

    for (j = 0; j < nrhs; j++) {
      x_work = &x[j*ldx];
      xtrue_work = &xtrue[j*ldxtrue];
      err = xnorm = 0.0;
      for (i = 0; i < n; i++) {
	err = SUPERLU_MAX(err, fabs(x_work[i] - xtrue_work[i]));
	xnorm = SUPERLU_MAX(xnorm, fabs(x_work[i]));
      }

      /* get the golbal max err & xnrom */
      temperr = err;
      tempxnorm = xnorm;
      MPI_Allreduce( &temperr, &err, 1, MPI_DOUBLE, MPI_MAX, grid->comm);
      MPI_Allreduce( &tempxnorm, &xnorm, 1, MPI_DOUBLE, MPI_MAX, grid->comm);

      err = err / xnorm;
      if ( !iam ) printf("(%d) .. ||X-Xtrue||/||X|| = %e\n", grid->iam, err);
    }
}

