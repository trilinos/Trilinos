

/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 */

#include <math.h>
#include "superlu_ddefs.h"

void pdgsmv_init
(
 SuperMatrix *A,       /* Matrix A permuted by columns (input/output).
			  The type of A can be:
			  Stype = NR_loc; Dtype = D; Mtype = GE. */
 int_t *row_to_proc,   /* Input. Mapping between rows and processes. */
 gridinfo_t *grid,     /* Input */
 pdgsmv_comm_t *gsmv_comm /* Output. The data structure for communication. */
 )
{
    NRformat_loc *Astore;
    int iam, p, procs;
    int *SendCounts, *RecvCounts;
    int_t i, j, k, l, m, m_loc, n, fst_row, jcol;
    int_t TotalIndSend, TotalValSend;
    int_t *colind, *rowptr;
    int_t *ind_tosend = NULL, *ind_torecv = NULL;
    int_t *ptr_ind_tosend, *ptr_ind_torecv;
    int_t *extern_start, *spa, *itemp;
    double *nzval, *val_tosend = NULL, *val_torecv = NULL, t;
    MPI_Request *send_req, *recv_req;
    MPI_Status status;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Enter pdgsmv_init()");
#endif

    /* ------------------------------------------------------------
       INITIALIZATION.
       ------------------------------------------------------------*/
    iam = grid->iam;
    procs = grid->nprow * grid->npcol;
    Astore = (NRformat_loc *) A->Store;
    m = A->nrow;
    n = A->ncol;
    m_loc = Astore->m_loc;
    fst_row = Astore->fst_row;
    colind = Astore->colind;
    rowptr = Astore->rowptr;
    nzval = Astore->nzval;
    if ( !(SendCounts = SUPERLU_MALLOC(2*procs * sizeof(int))) )
        ABORT("Malloc fails for SendCounts[]");
    /*for (i = 0; i < 2*procs; ++i) SendCounts[i] = 0;*/
    RecvCounts = SendCounts + procs;
    if ( !(ptr_ind_tosend = intMalloc_dist(2*(procs+1))) )
        ABORT("Malloc fails for ptr_ind_tosend[]");
    ptr_ind_torecv = ptr_ind_tosend + procs + 1;
    if ( !(extern_start = intMalloc_dist(m_loc)) )
        ABORT("Malloc fails for extern_start[]");
    for (i = 0; i < m_loc; ++i) extern_start[i] = rowptr[i];

    /* ------------------------------------------------------------
       COUNT THE NUMBER OF X ENTRIES TO BE SENT TO EACH PROCESS.
       THIS IS THE UNION OF THE COLUMN INDICES OF MY ROWS.
       SWAP TO THE BEGINNING THE PART OF A CORRESPONDING TO THE
       LOCAL PART OF X.
       THIS ACCOUNTS FOR THE FIRST PASS OF ACCESSING MATRIX A.
       ------------------------------------------------------------*/
    if ( !(spa = intCalloc_dist(n)) ) /* Aid in global to local translation */
        ABORT("Malloc fails for spa[]");
    for (p = 0; p < procs; ++p) SendCounts[p] = 0;
    for (i = 0; i < m_loc; ++i) { /* Loop through each row */
        k = extern_start[i];
        for (j = rowptr[i]; j < rowptr[i+1]; ++j) {/* Each nonzero in row i */
	    jcol = colind[j];
            p = row_to_proc[jcol];
	    if ( p != iam ) { /* External */
	        if ( spa[jcol] == 0 ) { /* First time see this index */
		    ++SendCounts[p];
		    spa[jcol] = 1;
                }
	    } else { /* Swap to beginning the part of A corresponding
			to the local part of X */
		l = colind[k];
		t = nzval[k];
		colind[k] = jcol;
		nzval[k] = nzval[j];
		colind[j] = l;
		nzval[j] = t;
		++k;
	    }
	}
	extern_start[i] = k;
    }

    /* ------------------------------------------------------------
       LOAD THE X-INDICES TO BE SENT TO THE OTHER PROCESSES.
       THIS ACCOUNTS FOR THE SECOND PASS OF ACCESSING MATRIX A.
       ------------------------------------------------------------*/
    /* Build pointers to ind_tosend[]. */
    ptr_ind_tosend[0] = 0;
    for (p = 0, TotalIndSend = 0; p < procs; ++p) {
        TotalIndSend += SendCounts[p]; /* Total to send. */
	ptr_ind_tosend[p+1] = ptr_ind_tosend[p] + SendCounts[p];
    }
#if 0
    ptr_ind_tosend[iam] = 0; /* Local part of X */
#endif
    if ( TotalIndSend ) {
        if ( !(ind_tosend = intMalloc_dist(TotalIndSend)) )
	    ABORT("Malloc fails for ind_tosend[]"); /* Exclude local part of X */
    }

    /* Build SPA to aid global to local translation. */
    for (i = 0; i < n; ++i) spa[i] = EMPTY;
    for (i = 0; i < m_loc; ++i) { /* Loop through each row of A */
        for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
	    jcol = colind[j];
	    if ( spa[jcol] == EMPTY ) { /* First time see this index */
	        p = row_to_proc[jcol];
		if ( p == iam ) { /* Local */
		  /*assert(jcol>=fst_row);*/
		  spa[jcol] = jcol - fst_row; /* Relative position in local X */
		} else {          /* External */
		  ind_tosend[ptr_ind_tosend[p]] = jcol; /* Still global */
		  spa[jcol] = ptr_ind_tosend[p]; /* Position in ind_tosend[] */
		  ++ptr_ind_tosend[p];
		}
	    }
	}
    }
    
    /* ------------------------------------------------------------
       TRANSFORM THE COLUMN INDICES OF MATRIX A INTO LOCAL INDICES.
       THIS ACCOUNTS FOR THE THIRD PASS OF ACCESSING MATRIX A.
       ------------------------------------------------------------*/
    for (i = 0; i < m_loc; ++i) {
        for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
	    jcol = colind[j];
	    colind[j] = spa[jcol];
	}
    }

    /* ------------------------------------------------------------
       COMMUNICATE THE EXTERNAL INDICES OF X.
       ------------------------------------------------------------*/
    MPI_Alltoall(SendCounts, 1, MPI_INT, RecvCounts, 1, MPI_INT,
		 grid->comm);

    /* Build pointers to ind_torecv[]. */
    ptr_ind_torecv[0] = 0;
    for (p = 0, TotalValSend = 0; p < procs; ++p) {
        TotalValSend += RecvCounts[p]; /* Total to receive. */
	ptr_ind_torecv[p+1] = ptr_ind_torecv[p] + RecvCounts[p];
    }
    if ( TotalValSend ) {
        if ( !(ind_torecv = intMalloc_dist(TotalValSend)) )
	    ABORT("Malloc fails for ind_torecv[]");
    }

    if ( !(send_req = (MPI_Request *)
	   SUPERLU_MALLOC(2*procs *sizeof(MPI_Request))))
        ABORT("Malloc fails for recv_req[].");
    recv_req = send_req + procs;
    for (p = 0; p < procs; ++p) {
        ptr_ind_tosend[p] -= SendCounts[p]; /* Reset pointer to beginning */
        if ( SendCounts[p] ) {
	    MPI_Isend(&ind_tosend[ptr_ind_tosend[p]], SendCounts[p],
		      mpi_int_t, p, iam, grid->comm, &send_req[p]);
	}
	if ( RecvCounts[p] ) {
	    MPI_Irecv(&ind_torecv[ptr_ind_torecv[p]], RecvCounts[p],
		      mpi_int_t, p, p, grid->comm, &recv_req[p]);
	}
    }
    for (p = 0; p < procs; ++p) {
        if ( SendCounts[p] ) MPI_Wait(&send_req[p], &status);
	if ( RecvCounts[p] ) MPI_Wait(&recv_req[p], &status);
    }

    /* Allocate storage for the X values to to transferred. */
    if ( TotalIndSend &&
         !(val_torecv = doubleMalloc_dist(TotalIndSend)) )
        ABORT("Malloc fails for val_torecv[].");
    if ( TotalValSend &&
         !(val_tosend = doubleMalloc_dist(TotalValSend)) )
        ABORT("Malloc fails for val_tosend[].");

    gsmv_comm->extern_start = extern_start;
    gsmv_comm->ind_tosend = ind_tosend;
    gsmv_comm->ind_torecv = ind_torecv;
    gsmv_comm->ptr_ind_tosend = ptr_ind_tosend;
    gsmv_comm->ptr_ind_torecv = ptr_ind_torecv;
    gsmv_comm->SendCounts = SendCounts;
    gsmv_comm->RecvCounts = RecvCounts;
    gsmv_comm->val_tosend = val_tosend;
    gsmv_comm->val_torecv = val_torecv;
    gsmv_comm->TotalIndSend = TotalIndSend;
    gsmv_comm->TotalValSend = TotalValSend;
    
    SUPERLU_FREE(spa);
    SUPERLU_FREE(send_req);

#if ( DEBUGlevel>=1 )
    PrintInt10("pdgsmv_init::rowptr", m_loc+1, rowptr);
    PrintInt10("pdgsmv_init::extern_start", m_loc, extern_start);
    CHECK_MALLOC(iam, "Exit pdgsmv_init()");
#endif

} /* PDGSMV_INIT */


/*
 * Performs sparse matrix-vector multiplication.
 */
void
pdgsmv
(
 int_t  abs,               /* Input. Do abs(A)*abs(x). */
 SuperMatrix *A_internal,  /* Input. Matrix A permuted by columns.
			      The column indices are translated into
			      the relative positions in the gathered x-vector.
			      The type of A can be:
			      Stype = NR_loc; Dtype = SLU_D; Mtype = GE. */
 gridinfo_t *grid,         /* Input */
 pdgsmv_comm_t *gsmv_comm, /* Input. The data structure for communication. */
 double x[],       /* Input. The distributed source vector */
 double ax[]       /* Output. The distributed destination vector */
)
{
    NRformat_loc *Astore;
    int iam, procs;
    int_t i, j, p, m, m_loc, n, fst_row, jcol;
    int_t *colind, *rowptr;
    int   *SendCounts, *RecvCounts;
    int_t *ind_tosend, *ind_torecv, *ptr_ind_tosend, *ptr_ind_torecv;
    int_t *extern_start, TotalValSend;
    double *nzval, *val_tosend, *val_torecv;
    double zero = 0.0;
    MPI_Request *send_req, *recv_req;
    MPI_Status status;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(grid->iam, "Enter pdgsmv()");
#endif

    /* ------------------------------------------------------------
       INITIALIZATION.
       ------------------------------------------------------------*/
    iam = grid->iam;
    procs = grid->nprow * grid->npcol;
    Astore = (NRformat_loc *) A_internal->Store;
    m = A_internal->nrow;
    n = A_internal->ncol;
    m_loc = Astore->m_loc;
    fst_row = Astore->fst_row;
    colind = Astore->colind;
    rowptr = Astore->rowptr;
    nzval = (double *) Astore->nzval;
    extern_start = gsmv_comm->extern_start;
    ind_torecv = gsmv_comm->ind_torecv;
    ptr_ind_tosend = gsmv_comm->ptr_ind_tosend;
    ptr_ind_torecv = gsmv_comm->ptr_ind_torecv;
    SendCounts = gsmv_comm->SendCounts;
    RecvCounts = gsmv_comm->RecvCounts;
    val_tosend = (double *) gsmv_comm->val_tosend;
    val_torecv = (double *) gsmv_comm->val_torecv;
    TotalValSend = gsmv_comm->TotalValSend;

    /* ------------------------------------------------------------
       COPY THE X VALUES INTO THE SEND BUFFER.
       ------------------------------------------------------------*/
    for (i = 0; i < TotalValSend; ++i) {
        j = ind_torecv[i] - fst_row; /* Relative index in x[] */
	val_tosend[i] = x[j];
    }

    /* ------------------------------------------------------------
       COMMUNICATE THE X VALUES.
       ------------------------------------------------------------*/
    if ( !(send_req = (MPI_Request *)
	   SUPERLU_MALLOC(2*procs *sizeof(MPI_Request))))
        ABORT("Malloc fails for recv_req[].");
    recv_req = send_req + procs;
    for (p = 0; p < procs; ++p) {
        if ( RecvCounts[p] ) {
	    MPI_Isend(&val_tosend[ptr_ind_torecv[p]], RecvCounts[p],
                      MPI_DOUBLE, p, iam,
                      grid->comm, &send_req[p]);
	}
	if ( SendCounts[p] ) {
	    MPI_Irecv(&val_torecv[ptr_ind_tosend[p]], SendCounts[p],
                      MPI_DOUBLE, p, p,
                      grid->comm, &recv_req[p]);
	}
    }
    for (p = 0; p < procs; ++p) {
        if ( RecvCounts[p] ) MPI_Wait(&send_req[p], &status);
	if ( SendCounts[p] ) MPI_Wait(&recv_req[p], &status);
    }
    
    /* ------------------------------------------------------------
       PERFORM THE ACTUAL MULTIPLICATION.
       ------------------------------------------------------------*/
    if ( abs ) { /* Perform abs(A)*abs(x) */
        for (i = 0; i < m_loc; ++i) { /* Loop through each row */
	    ax[i] = 0.0;

	    /* Multiply the local part. */
	    for (j = rowptr[i]; j < extern_start[i]; ++j) {
	        jcol = colind[j];
		ax[i] += fabs(nzval[j]) * fabs(x[jcol]);
	    }

	    /* Multiply the external part. */
	    for (; j < rowptr[i+1]; ++j) {
	        jcol = colind[j];
	        ax[i] += fabs(nzval[j]) * fabs(val_torecv[jcol]);
	    }
	}
    } else {
        for (i = 0; i < m_loc; ++i) { /* Loop through each row */
	    ax[i] = zero;

	    /* Multiply the local part. */
	    for (j = rowptr[i]; j < extern_start[i]; ++j) {
	        jcol = colind[j];
		ax[i] += nzval[j] * x[jcol];
	    }

	    /* Multiply the external part. */
	    for (; j < rowptr[i+1]; ++j) {
	        jcol = colind[j];
	        ax[i] += nzval[j] * val_torecv[jcol];
	    }
	}
    }

    SUPERLU_FREE(send_req);
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit pdgsmv()");
#endif

} /* PDGSMV */

void pdgsmv_finalize(pdgsmv_comm_t *gsmv_comm)
{
    int_t *it;
    double *dt;
    SUPERLU_FREE(gsmv_comm->extern_start);
    if ( it = gsmv_comm->ind_tosend ) SUPERLU_FREE(it);
    if ( it = gsmv_comm->ind_torecv ) SUPERLU_FREE(it);
    SUPERLU_FREE(gsmv_comm->ptr_ind_tosend);
    SUPERLU_FREE(gsmv_comm->SendCounts);
    if ( dt = gsmv_comm->val_tosend ) SUPERLU_FREE(dt);
    if ( dt = gsmv_comm->val_torecv ) SUPERLU_FREE(dt);
}

