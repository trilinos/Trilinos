
#include "superlu_zdefs.h"

int_t
zReDistribute_A(SuperMatrix *A, ScalePermstruct_t *ScalePermstruct,
                Glu_freeable_t *Glu_freeable, int_t *xsup, int_t *supno,
                gridinfo_t *grid, int_t *colptr[], int_t *rowind[],
                doublecomplex *a[])
{
/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 * Purpose
 * =======
 *   Re-distribute A on the 2D process mesh.
 * 
 * Arguments
 * =========
 * 
 * A      (input) SuperMatrix*
 *	  The distributed input matrix A of dimension (A->nrow, A->ncol).
 *        A may be overwritten by diag(R)*A*diag(C)*Pc^T.
 *        The type of A can be: Stype = SLU_NR_loc; Dtype = SLU_Z; Mtype = SLU_GE.
 *
 * ScalePermstruct (input) ScalePermstruct_t*
 *        The data structure to store the scaling and permutation vectors
 *        describing the transformations performed to the original matrix A.
 *
 * Glu_freeable (input) *Glu_freeable_t
 *        The global structure describing the graph of L and U.
 * 
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * colptr (output) int*
 *
 * rowind (output) int*
 *
 * a      (output) doublecomplex*
 *
 * Return value
 * ============
 *
 */
    NRformat_loc *Astore;
    int_t  *perm_r; /* row permutation vector */
    int_t  *perm_c; /* column permutation vector */
    int_t  i, irow, fst_row, j, jcol, k, gbi, gbj, n, m_loc, jsize;
    int_t  nnz_loc;    /* number of local nonzeros */
    int_t  nnz_remote; /* number of remote nonzeros to be sent */
    int_t  SendCnt; /* number of remote nonzeros to be sent */
    int_t  RecvCnt; /* number of remote nonzeros to be sent */
    int_t  *nnzToSend, *nnzToRecv, maxnnzToRecv;
    int_t  *ia, *ja, **ia_send, *index, *itemp;
    int_t  *ptr_to_send;
    doublecomplex *aij, **aij_send, *nzval, *dtemp;
    doublecomplex *nzval_a;
    int    iam, it, p, procs;
    MPI_Request *send_req;
    MPI_Status  status;
    

    /* ------------------------------------------------------------
       INITIALIZATION.
       ------------------------------------------------------------*/
    iam = grid->iam;
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter zReDistribute_A()");
#endif
    perm_r = ScalePermstruct->perm_r;
    perm_c = ScalePermstruct->perm_c;
    procs = grid->nprow * grid->npcol;
    Astore = (NRformat_loc *) A->Store;
    n = A->ncol;
    m_loc = Astore->m_loc;
    fst_row = Astore->fst_row;
    nnzToRecv = intCalloc_dist(2*procs);
    nnzToSend = nnzToRecv + procs;


    /* ------------------------------------------------------------
       COUNT THE NUMBER OF NONZEROS TO BE SENT TO EACH PROCESS,
       THEN ALLOCATE SPACE.
       THIS ACCOUNTS FOR THE FIRST PASS OF A.
       ------------------------------------------------------------*/
    for (i = 0; i < m_loc; ++i) {
        for (j = Astore->rowptr[i]; j < Astore->rowptr[i+1]; ++j) {
  	    irow = perm_c[perm_r[i+fst_row]];  /* Row number in Pc*Pr*A */
	    jcol = Astore->colind[j];
	    gbi = BlockNum( irow );
	    gbj = BlockNum( jcol );
	    p = PNUM( PROW(gbi,grid), PCOL(gbj,grid), grid );
	    ++nnzToSend[p]; 
	}
    }

    /* All-to-all communication */
    MPI_Alltoall( nnzToSend, 1, mpi_int_t, nnzToRecv, 1, mpi_int_t,
		  grid->comm);

    maxnnzToRecv = 0;
    nnz_loc = SendCnt = RecvCnt = 0;

    for (p = 0; p < procs; ++p) {
	if ( p != iam ) {
	    SendCnt += nnzToSend[p];
	    RecvCnt += nnzToRecv[p];
	    maxnnzToRecv = SUPERLU_MAX( nnzToRecv[p], maxnnzToRecv );
	} else {
	    nnz_loc += nnzToRecv[p];
	    /*assert(nnzToSend[p] == nnzToRecv[p]);*/
	}
    }
    k = nnz_loc + RecvCnt; /* Total nonzeros ended up in my process. */

    /* Allocate space for storing the triplets after redistribution. */
    if ( !(ia = intMalloc_dist(2*k)) )
        ABORT("Malloc fails for ia[].");
    ja = ia + k;
    if ( !(aij = doublecomplexMalloc_dist(k)) )
        ABORT("Malloc fails for aij[].");

    /* Allocate temporary storage for sending/receiving the A triplets. */
    if ( procs > 1 ) {
      if ( !(send_req = (MPI_Request *)
	     SUPERLU_MALLOC(2*procs *sizeof(MPI_Request))) )
	ABORT("Malloc fails for send_req[].");
      if ( !(ia_send = (int_t **) SUPERLU_MALLOC(procs*sizeof(int_t*))) )
        ABORT("Malloc fails for ia_send[].");
      if ( !(aij_send = (doublecomplex **)SUPERLU_MALLOC(procs*sizeof(doublecomplex*))) )
        ABORT("Malloc fails for aij_send[].");
      if ( !(index = intMalloc_dist(2*SendCnt)) )
        ABORT("Malloc fails for index[].");
      if ( !(nzval = doublecomplexMalloc_dist(SendCnt)) )
        ABORT("Malloc fails for nzval[].");
      if ( !(ptr_to_send = intCalloc_dist(procs)) )
        ABORT("Malloc fails for ptr_to_send[].");
      if ( !(itemp = intMalloc_dist(2*maxnnzToRecv)) )
        ABORT("Malloc fails for itemp[].");
      if ( !(dtemp = doublecomplexMalloc_dist(maxnnzToRecv)) )
        ABORT("Malloc fails for dtemp[].");

      for (i = 0, j = 0, p = 0; p < procs; ++p) {
          if ( p != iam ) {
	      ia_send[p] = &index[i];
	      i += 2 * nnzToSend[p]; /* ia/ja indices alternate */
	      aij_send[p] = &nzval[j];
	      j += nnzToSend[p];
	  }
      }
    } /* if procs > 1 */
      
    if ( !(*colptr = intCalloc_dist(n+1)) )
        ABORT("Malloc fails for *colptr[].");

    /* ------------------------------------------------------------
       LOAD THE ENTRIES OF A INTO THE (IA,JA,AIJ) STRUCTURES TO SEND.
       THIS ACCOUNTS FOR THE SECOND PASS OF A.
       ------------------------------------------------------------*/
    nnz_loc = 0; /* Reset the local nonzero count. */
    nzval_a = Astore->nzval;
    for (i = 0; i < m_loc; ++i) {
        for (j = Astore->rowptr[i]; j < Astore->rowptr[i+1]; ++j) {
  	    irow = perm_c[perm_r[i+fst_row]];  /* Row number in Pc*Pr*A */
	    jcol = Astore->colind[j];
	    gbi = BlockNum( irow );
	    gbj = BlockNum( jcol );
	    p = PNUM( PROW(gbi,grid), PCOL(gbj,grid), grid );

	    if ( p != iam ) { /* remote */
	        k = ptr_to_send[p];
	        ia_send[p][k] = irow;
	        ia_send[p][k + nnzToSend[p]] = jcol;
		aij_send[p][k] = nzval_a[j];
		++ptr_to_send[p]; 
	    } else {          /* local */
	        ia[nnz_loc] = irow;
	        ja[nnz_loc] = jcol;
		aij[nnz_loc] = nzval_a[j];
		++nnz_loc;
		++(*colptr)[jcol]; /* Count nonzeros in each column */
	    }
	}
    }

    /* ------------------------------------------------------------
       PERFORM REDISTRIBUTION. THIS INVOLVES ALL-TO-ALL COMMUNICATION.
       NOTE: Can possibly use MPI_Alltoallv.
       ------------------------------------------------------------*/
    for (p = 0; p < procs; ++p) {
        if ( p != iam ) {
	    it = 2*nnzToSend[p];
	    MPI_Isend( ia_send[p], it, mpi_int_t,
		       p, iam, grid->comm, &send_req[p] );
	    it = nnzToSend[p];
	    MPI_Isend( aij_send[p], it, SuperLU_MPI_DOUBLE_COMPLEX,
	               p, iam+procs, grid->comm, &send_req[procs+p] ); 
	}
    }

    for (p = 0; p < procs; ++p) {
        if ( p != iam ) {
	    it = 2*nnzToRecv[p];
	    MPI_Recv( itemp, it, mpi_int_t, p, p, grid->comm, &status ); 
	    it = nnzToRecv[p];
            MPI_Recv( dtemp, it, SuperLU_MPI_DOUBLE_COMPLEX, p, p+procs,
		      grid->comm, &status );
	    for (i = 0; i < nnzToRecv[p]; ++i) {
	        ia[nnz_loc] = itemp[i];
		jcol = itemp[i + nnzToRecv[p]];
		/*assert(jcol<n);*/
	        ja[nnz_loc] = jcol;
		aij[nnz_loc] = dtemp[i];
		++nnz_loc;
		++(*colptr)[jcol]; /* Count nonzeros in each column */ 
	    }
	}
    }

    for (p = 0; p < procs; ++p) {
        if ( p != iam ) {
	    MPI_Wait( &send_req[p], &status);
	    MPI_Wait( &send_req[procs+p], &status);
	}
    }

    /* ------------------------------------------------------------
       DEALLOCATE TEMPORARY STORAGE
       ------------------------------------------------------------*/

    SUPERLU_FREE(nnzToRecv);

    if ( procs > 1 ) {
	SUPERLU_FREE(send_req);
	SUPERLU_FREE(ia_send);
	SUPERLU_FREE(aij_send);
	SUPERLU_FREE(index);
	SUPERLU_FREE(nzval);
	SUPERLU_FREE(ptr_to_send);
	SUPERLU_FREE(itemp);
	SUPERLU_FREE(dtemp);
    }

    /* ------------------------------------------------------------
       CONVERT THE TRIPLET FORMAT INTO THE CCS FORMAT.
       ------------------------------------------------------------*/
    if ( !(*rowind = intMalloc_dist(nnz_loc)) )
        ABORT("Malloc fails for *rowind[].");
    if ( !(*a = doublecomplexMalloc_dist(nnz_loc)) )
        ABORT("Malloc fails for *a[].");

    /* Initialize the array of column pointers */
    k = 0;
    jsize = (*colptr)[0];
    (*colptr)[0] = 0;
    for (j = 1; j < n; ++j) {
	k += jsize;
	jsize = (*colptr)[j];
	(*colptr)[j] = k;
    }
    
    /* Copy the triplets into the column oriented storage */
    for (i = 0; i < nnz_loc; ++i) {
	j = ja[i];
	k = (*colptr)[j];
	(*rowind)[k] = ia[i];
	(*a)[k] = aij[i];
	++(*colptr)[j];
    }

    /* Reset the column pointers to the beginning of each column */
    for (j = n; j > 0; --j) (*colptr)[j] = (*colptr)[j-1];
    (*colptr)[0] = 0;

    SUPERLU_FREE(ia);
    SUPERLU_FREE(aij);

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit zReDistribute_A()");
#endif
 
} /* zReDistribute_A */

int_t
pzdistribute(fact_t fact, int_t n, SuperMatrix *A,
	     ScalePermstruct_t *ScalePermstruct,
	     Glu_freeable_t *Glu_freeable, LUstruct_t *LUstruct,
	     gridinfo_t *grid)
/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 *
 * Purpose
 * =======
 *   Distribute the matrix onto the 2D process mesh.
 * 
 * Arguments
 * =========
 * 
 * fact (input) fact_t
 *        Specifies whether or not the L and U structures will be re-used.
 *        = SamePattern_SameRowPerm: L and U structures are input, and
 *                                   unchanged on exit.
 *        = DOFACT or SamePattern: L and U structures are computed and output.
 *
 * n      (input) int
 *        Dimension of the matrix.
 *
 * A      (input) SuperMatrix*
 *	  The distributed input matrix A of dimension (A->nrow, A->ncol).
 *        A may be overwritten by diag(R)*A*diag(C)*Pc^T.
 *        The type of A can be: Stype = NR; Dtype = SLU_Z; Mtype = GE.
 *
 * ScalePermstruct (input) ScalePermstruct_t*
 *        The data structure to store the scaling and permutation vectors
 *        describing the transformations performed to the original matrix A.
 *
 * Glu_freeable (input) *Glu_freeable_t
 *        The global structure describing the graph of L and U.
 * 
 * LUstruct (input) LUstruct_t*
 *        Data structures for L and U factors.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Return value
 * ============
 *   > 0, working storage required (in bytes).
 *
 */
{
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    LocalLU_t *Llu = LUstruct->Llu;
    int_t bnnz, fsupc, i, irow, istart, j, jb, jj, k, len, len1, nsupc;
    int_t ljb;  /* local block column number */
    int_t nrbl; /* number of L blocks in current block column */
    int_t nrbu; /* number of U blocks in current block column */
    int_t gb;   /* global block number; 0 < gb <= nsuper */
    int_t lb;   /* local block number; 0 < lb <= ceil(NSUPERS/Pr) */
    int iam, jbrow, kcol, mycol, myrow, pc, pr;
    int_t mybufmax[NBUFFERS];
#if 0
    NCPformat *Astore;
#else /* XSL ==> */
    NRformat_loc *Astore;
#endif
    doublecomplex *a;
    int_t *asub, *xa;
#if 0
    int_t *xa_begin, *xa_end;
#endif
    int_t *xsup = Glu_persist->xsup;    /* supernode and column mapping */
    int_t *supno = Glu_persist->supno;   
    int_t *lsub, *xlsub, *usub, *xusub;
    int_t nsupers;
    int_t next_lind;      /* next available position in index[*] */
    int_t next_lval;      /* next available position in nzval[*] */
    int_t *index;         /* indices consist of headers and row subscripts */
    doublecomplex *lusup, *uval; /* nonzero values in L and U */
    doublecomplex **Lnzval_bc_ptr;  /* size ceil(NSUPERS/Pc) */
    int_t  **Lrowind_bc_ptr; /* size ceil(NSUPERS/Pc) */
    doublecomplex **Unzval_br_ptr;  /* size ceil(NSUPERS/Pr) */
    int_t  **Ufstnz_br_ptr;  /* size ceil(NSUPERS/Pr) */

    /*-- Counts to be used in factorization. --*/
    int_t  *ToRecv, *ToSendD, **ToSendR;

    /*-- Counts to be used in lower triangular solve. --*/
    int_t  *fmod;          /* Modification count for L-solve.        */
    int_t  **fsendx_plist; /* Column process list to send down Xk.   */
    int_t  nfrecvx = 0;    /* Number of Xk I will receive.           */
    int_t  kseen;

    /*-- Counts to be used in upper triangular solve. --*/
    int_t  *bmod;          /* Modification count for U-solve.        */
    int_t  **bsendx_plist; /* Column process list to send down Xk.   */
    int_t  nbrecvx = 0;    /* Number of Xk I will receive.           */
    int_t  *ilsum;         /* starting position of each supernode in 
			      the full array (local)                 */

    /*-- Auxiliary arrays; freed on return --*/
    int_t *rb_marker;  /* block hit marker; size ceil(NSUPERS/Pr)           */
    int_t *Urb_length; /* U block length; size ceil(NSUPERS/Pr)             */
    int_t *Urb_indptr; /* pointers to U index[]; size ceil(NSUPERS/Pr)      */
    int_t *Urb_fstnz;  /* # of fstnz in a block row; size ceil(NSUPERS/Pr)  */
    int_t *Ucbs;       /* number of column blocks in a block row            */
    int_t *Lrb_length; /* L block length; size ceil(NSUPERS/Pr)             */
    int_t *Lrb_number; /* global block number; size ceil(NSUPERS/Pr)        */
    int_t *Lrb_indptr; /* pointers to L index[]; size ceil(NSUPERS/Pr)      */
    int_t *Lrb_valptr; /* pointers to L nzval[]; size ceil(NSUPERS/Pr)      */
    doublecomplex *dense, *dense_col; /* SPA */
    doublecomplex zero = {0.0, 0.0};
    int_t ldaspa;     /* LDA of SPA */
    int_t mem_use = 0, iword, dword;

#if ( PRNTlevel>=1 )
    int_t nLblocks = 0, nUblocks = 0;
#endif
#if ( PROFlevel>=1 ) 
    double t, t_u, t_l;
    int_t u_blks;
#endif

    /* Initialization. */
    iam = grid->iam;
    myrow = MYROW( iam, grid );
    mycol = MYCOL( iam, grid );
    for (i = 0; i < NBUFFERS; ++i) mybufmax[i] = 0;
    nsupers  = supno[n-1] + 1;
    Astore   = (NRformat_loc *) A->Store;

#if ( PRNTlevel>=1 )
    iword = sizeof(int_t);
    dword = sizeof(doublecomplex);
#endif

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter pzdistribute()");
#endif

    zReDistribute_A(A, ScalePermstruct, Glu_freeable, xsup, supno,
		      grid, &xa, &asub, &a);

    if ( fact == SamePattern_SameRowPerm ) {
#if ( PROFlevel>=1 )
	t_l = t_u = 0; u_blks = 0;
#endif
	/* We can propagate the new values of A into the existing
	   L and U data structures.            */
	ilsum = Llu->ilsum;
	ldaspa = Llu->ldalsum;
	if ( !(dense = doublecomplexCalloc_dist(ldaspa * sp_ienv_dist(3))) )
	    ABORT("Calloc fails for SPA dense[].");
	nrbu = CEILING( nsupers, grid->nprow ); /* Number of local block rows */
	if ( !(Urb_length = intCalloc_dist(nrbu)) )
	    ABORT("Calloc fails for Urb_length[].");
	if ( !(Urb_indptr = intMalloc_dist(nrbu)) )
	    ABORT("Malloc fails for Urb_indptr[].");
	for (lb = 0; lb < nrbu; ++lb) 
	    Urb_indptr[lb] = BR_HEADER; /* Skip header in U index[]. */
	Lrowind_bc_ptr = Llu->Lrowind_bc_ptr;
	Lnzval_bc_ptr = Llu->Lnzval_bc_ptr;
	Ufstnz_br_ptr = Llu->Ufstnz_br_ptr;
	Unzval_br_ptr = Llu->Unzval_br_ptr;
#if ( PRNTlevel>=1 )
	mem_use += 2*nrbu*iword + ldaspa*sp_ienv_dist(3)*dword;
#endif
	for (jb = 0; jb < nsupers; ++jb) { /* Loop through each block column */
	    pc = PCOL( jb, grid );
	    if ( mycol == pc ) { /* Block column jb in my process column */
		fsupc = FstBlockC( jb );
		nsupc = SuperSize( jb );

		/* Scatter A into SPA. */
		for (j = fsupc, dense_col = dense; j < FstBlockC(jb+1); ++j) {
		    for (i = xa[j]; i < xa[j+1]; ++i) {
			irow = asub[i];
			gb = BlockNum( irow );
			if ( myrow == PROW( gb, grid ) ) {
			    lb = LBi( gb, grid );
			    irow = ilsum[lb] + irow - FstBlockC( gb );
			    dense_col[irow] = a[i];
			}
		    }
		    dense_col += ldaspa;
		}

#if ( PROFlevel>=1 )
		t = SuperLU_timer_();
#endif
		/* Gather the values of A from SPA into Unzval[]. */
		for (lb = 0; lb < nrbu; ++lb) {
		    index = Ufstnz_br_ptr[lb];
		    if ( index && index[Urb_indptr[lb]] == jb ) {
			uval = Unzval_br_ptr[lb];
			len = Urb_indptr[lb] + UB_DESCRIPTOR;
			gb = lb * grid->nprow + myrow;/* Global block number */
			k = FstBlockC( gb+1 );
			irow = ilsum[lb] - FstBlockC( gb );
			for (jj = 0, dense_col = dense; jj < nsupc; ++jj) {
			    j = index[len+jj];
			    for (i = j; i < k; ++i) {
				uval[Urb_length[lb]++] = dense_col[irow+i];
				dense_col[irow+i] = zero;
			    }
			    dense_col += ldaspa;
			}
			Urb_indptr[lb] += UB_DESCRIPTOR + nsupc;
		    } /* if index != NULL */
		} /* for lb ... */
#if ( PROFlevel>=1 )
		t_u += SuperLU_timer_() - t;
		t = SuperLU_timer_();
#endif
		/* Gather the values of A from SPA into Lnzval[]. */
		ljb = LBj( jb, grid ); /* Local block number */
		index = Lrowind_bc_ptr[ljb];
		if ( index ) {
		    nrbl = index[0];   /* Number of row blocks. */
		    len = index[1];    /* LDA of lusup[]. */
		    lusup = Lnzval_bc_ptr[ljb];
		    next_lind = BC_HEADER;
		    next_lval = 0;
		    for (jj = 0; jj < nrbl; ++jj) {
			gb = index[next_lind++];
			len1 = index[next_lind++]; /* Rows in the block. */
			lb = LBi( gb, grid );
			for (bnnz = 0; bnnz < len1; ++bnnz) {
			    irow = index[next_lind++]; /* Global index. */
			    irow = ilsum[lb] + irow - FstBlockC( gb );
			    k = next_lval++;
			    for (j = 0, dense_col = dense; j < nsupc; ++j) {
				lusup[k] = dense_col[irow];
				dense_col[irow] = zero;
				k += len;
				dense_col += ldaspa;
			    }
			} /* for bnnz ... */
		    } /* for jj ... */
		} /* if index ... */
#if ( PROFlevel>=1 )
		t_l += SuperLU_timer_() - t;
#endif
	    } /* if mycol == pc */
	} /* for jb ... */

	SUPERLU_FREE(dense);
	SUPERLU_FREE(Urb_length);
	SUPERLU_FREE(Urb_indptr);
#if ( PROFlevel>=1 )
	if ( !iam ) printf(".. 2nd distribute time: L %.2f\tU %.2f\tu_blks %d\tnrbu %d\n",
			   t_l, t_u, u_blks, nrbu);
#endif

    } else {
        /* ------------------------------------------------------------
	   FIRST TIME CREATING THE L AND U DATA STRUCTURES.
	   ------------------------------------------------------------*/

#if ( PROFlevel>=1 )
	t_l = t_u = 0; u_blks = 0;
#endif
	/* We first need to set up the L and U data structures and then
	 * propagate the values of A into them.
	 */
	lsub = Glu_freeable->lsub;    /* compressed L subscripts */
	xlsub = Glu_freeable->xlsub;
	usub = Glu_freeable->usub;    /* compressed U subscripts */
	xusub = Glu_freeable->xusub;
    
	if ( !(ToRecv = intCalloc_dist(nsupers)) )
	    ABORT("Calloc fails for ToRecv[].");

	k = CEILING( nsupers, grid->npcol );/* Number of local column blocks */
	if ( !(ToSendR = (int_t **) SUPERLU_MALLOC(k*sizeof(int_t*))) )
	    ABORT("Malloc fails for ToSendR[].");
	j = k * grid->npcol;
	if ( !(index = intMalloc_dist(j)) )
	    ABORT("Malloc fails for index[].");
#if ( PRNTlevel>=1 )
	mem_use += k*sizeof(int_t*) + (j + nsupers)*iword;
#endif
	for (i = 0; i < j; ++i) index[i] = EMPTY;
	for (i = 0,j = 0; i < k; ++i, j += grid->npcol) ToSendR[i] = &index[j];
	k = CEILING( nsupers, grid->nprow ); /* Number of local block rows */

	/* Pointers to the beginning of each block row of U. */
	if ( !(Unzval_br_ptr = 
              (doublecomplex**)SUPERLU_MALLOC(k * sizeof(doublecomplex*))) )
	    ABORT("Malloc fails for Unzval_br_ptr[].");
	if ( !(Ufstnz_br_ptr = (int_t**)SUPERLU_MALLOC(k * sizeof(int_t*))) )
	    ABORT("Malloc fails for Ufstnz_br_ptr[].");
	
	if ( !(ToSendD = intCalloc_dist(k)) )
	    ABORT("Malloc fails for ToSendD[].");
	if ( !(ilsum = intMalloc_dist(k+1)) )
	    ABORT("Malloc fails for ilsum[].");

	/* Auxiliary arrays used to set up U block data structures.
	   They are freed on return. */
	if ( !(rb_marker = intCalloc_dist(k)) )
	    ABORT("Calloc fails for rb_marker[].");
	if ( !(Urb_length = intCalloc_dist(k)) )
	    ABORT("Calloc fails for Urb_length[].");
	if ( !(Urb_indptr = intMalloc_dist(k)) )
	    ABORT("Malloc fails for Urb_indptr[].");
	if ( !(Urb_fstnz = intCalloc_dist(k)) )
	    ABORT("Calloc fails for Urb_fstnz[].");
	if ( !(Ucbs = intCalloc_dist(k)) )
	    ABORT("Calloc fails for Ucbs[].");
#if ( PRNTlevel>=1 )	
	mem_use += 2*k*sizeof(int_t*) + (7*k+1)*iword;
#endif
	/* Compute ldaspa and ilsum[]. */
	ldaspa = 0;
	ilsum[0] = 0;
	for (gb = 0; gb < nsupers; ++gb) {
	    if ( myrow == PROW( gb, grid ) ) {
		i = SuperSize( gb );
		ldaspa += i;
		lb = LBi( gb, grid );
		ilsum[lb + 1] = ilsum[lb] + i;
	    }
	}
	
            
	/* ------------------------------------------------------------
	   COUNT NUMBER OF ROW BLOCKS AND THE LENGTH OF EACH BLOCK IN U.
	   THIS ACCOUNTS FOR ONE-PASS PROCESSING OF G(U).
	   ------------------------------------------------------------*/
	
	/* Loop through each supernode column. */
	for (jb = 0; jb < nsupers; ++jb) {
	    pc = PCOL( jb, grid );
	    fsupc = FstBlockC( jb );
	    nsupc = SuperSize( jb );
	    /* Loop through each column in the block. */
	    for (j = fsupc; j < fsupc + nsupc; ++j) {
		/* usub[*] contains only "first nonzero" in each segment. */
		for (i = xusub[j]; i < xusub[j+1]; ++i) {
		    irow = usub[i]; /* First nonzero of the segment. */
		    gb = BlockNum( irow );
		    kcol = PCOL( gb, grid );
		    ljb = LBj( gb, grid );
		    if ( mycol == kcol && mycol != pc ) ToSendR[ljb][pc] = YES;
		    pr = PROW( gb, grid );
		    lb = LBi( gb, grid );
		    if ( mycol == pc ) {
			if  ( myrow == pr ) {
			    ToSendD[lb] = YES;
			    /* Count nonzeros in entire block row. */
			    Urb_length[lb] += FstBlockC( gb+1 ) - irow;
			    if (rb_marker[lb] <= jb) {/* First see the block */
				rb_marker[lb] = jb + 1;
				Urb_fstnz[lb] += nsupc;
				++Ucbs[lb]; /* Number of column blocks
					       in block row lb. */
#if ( PRNTlevel>=1 )
				++nUblocks;
#endif
			    }
			    ToRecv[gb] = 1;
			} else ToRecv[gb] = 2; /* Do I need 0, 1, 2 ? */
		    }
		} /* for i ... */
	    } /* for j ... */
	} /* for jb ... */
	
	/* Set up the initial pointers for each block row in U. */
	nrbu = CEILING( nsupers, grid->nprow );/* Number of local block rows */
	for (lb = 0; lb < nrbu; ++lb) {
	    len = Urb_length[lb];
	    rb_marker[lb] = 0; /* Reset block marker. */
	    if ( len ) {
		/* Add room for descriptors */
		len1 = Urb_fstnz[lb] + BR_HEADER + Ucbs[lb] * UB_DESCRIPTOR;
		if ( !(index = intMalloc_dist(len1+1)) )
		    ABORT("Malloc fails for Uindex[].");
		Ufstnz_br_ptr[lb] = index;
		if ( !(Unzval_br_ptr[lb] = doublecomplexMalloc_dist(len)) )
		    ABORT("Malloc fails for Unzval_br_ptr[*][].");
		mybufmax[2] = SUPERLU_MAX( mybufmax[2], len1 );
		mybufmax[3] = SUPERLU_MAX( mybufmax[3], len );
		index[0] = Ucbs[lb]; /* Number of column blocks */
		index[1] = len;      /* Total length of nzval[] */
		index[2] = len1;     /* Total length of index[] */
		index[len1] = -1;    /* End marker */
	    } else {
		Ufstnz_br_ptr[lb] = NULL;
		Unzval_br_ptr[lb] = NULL;
	    }
	    Urb_length[lb] = 0; /* Reset block length. */
	    Urb_indptr[lb] = BR_HEADER; /* Skip header in U index[]. */
	} /* for lb ... */

	SUPERLU_FREE(Urb_fstnz);
	SUPERLU_FREE(Ucbs);
#if ( PRNTlevel>=1 )
        mem_use -= 2*k * iword;
#endif
	/* Auxiliary arrays used to set up L block data structures.
	   They are freed on return.
	   k is the number of local row blocks.   */
	if ( !(Lrb_length = intCalloc_dist(k)) )
	    ABORT("Calloc fails for Lrb_length[].");
	if ( !(Lrb_number = intMalloc_dist(k)) )
	    ABORT("Malloc fails for Lrb_number[].");
	if ( !(Lrb_indptr = intMalloc_dist(k)) )
	    ABORT("Malloc fails for Lrb_indptr[].");
	if ( !(Lrb_valptr = intMalloc_dist(k)) )
	    ABORT("Malloc fails for Lrb_valptr[].");
	if ( !(dense = doublecomplexCalloc_dist(ldaspa * sp_ienv_dist(3))) )
	    ABORT("Calloc fails for SPA dense[].");

	/* These counts will be used for triangular solves. */
	if ( !(fmod = intCalloc_dist(k)) )
	    ABORT("Calloc fails for fmod[].");
	if ( !(bmod = intCalloc_dist(k)) )
	    ABORT("Calloc fails for bmod[].");
	/* ------------------------------------------------ */
#if ( PRNTlevel>=1 )	
	mem_use += 6*k*iword + ldaspa*sp_ienv_dist(3)*dword;
#endif
	k = CEILING( nsupers, grid->npcol );/* Number of local block columns */

	/* Pointers to the beginning of each block column of L. */
	if ( !(Lnzval_bc_ptr = 
              (doublecomplex**)SUPERLU_MALLOC(k * sizeof(doublecomplex*))) )
	    ABORT("Malloc fails for Lnzval_bc_ptr[].");
	if ( !(Lrowind_bc_ptr = (int_t**)SUPERLU_MALLOC(k * sizeof(int_t*))) )
	    ABORT("Malloc fails for Lrowind_bc_ptr[].");
	Lrowind_bc_ptr[k-1] = NULL;

	/* These lists of processes will be used for triangular solves. */
	if ( !(fsendx_plist = (int_t **) SUPERLU_MALLOC(k*sizeof(int_t*))) )
	    ABORT("Malloc fails for fsendx_plist[].");
	len = k * grid->nprow;
	if ( !(index = intMalloc_dist(len)) )
	    ABORT("Malloc fails for fsendx_plist[0]");
	for (i = 0; i < len; ++i) index[i] = EMPTY;
	for (i = 0, j = 0; i < k; ++i, j += grid->nprow)
	    fsendx_plist[i] = &index[j];
	if ( !(bsendx_plist = (int_t **) SUPERLU_MALLOC(k*sizeof(int_t*))) )
	    ABORT("Malloc fails for bsendx_plist[].");
	if ( !(index = intMalloc_dist(len)) )
	    ABORT("Malloc fails for bsendx_plist[0]");
	for (i = 0; i < len; ++i) index[i] = EMPTY;
	for (i = 0, j = 0; i < k; ++i, j += grid->nprow)
	    bsendx_plist[i] = &index[j];
	/* -------------------------------------------------------------- */
#if ( PRNTlevel>=1 )
	mem_use += 4*k*sizeof(int_t*) + 2*len*iword;
#endif

	/*------------------------------------------------------------
	  PROPAGATE ROW SUBSCRIPTS AND VALUES OF A INTO L AND U BLOCKS.
	  THIS ACCOUNTS FOR ONE-PASS PROCESSING OF A, L AND U.
	  ------------------------------------------------------------*/

	for (jb = 0; jb < nsupers; ++jb) {
	    pc = PCOL( jb, grid );
	    if ( mycol == pc ) { /* Block column jb in my process column */
		fsupc = FstBlockC( jb );
		nsupc = SuperSize( jb );
		ljb = LBj( jb, grid ); /* Local block number */
		
		/* Scatter A into SPA. */
		for (j = fsupc, dense_col = dense; j < FstBlockC(jb+1); ++j) {
		    for (i = xa[j]; i < xa[j+1]; ++i) {
			irow = asub[i];
			gb = BlockNum( irow );
			if ( myrow == PROW( gb, grid ) ) {
			    lb = LBi( gb, grid );
			    irow = ilsum[lb] + irow - FstBlockC( gb );
			    dense_col[irow] = a[i];
			}
		    }
		    dense_col += ldaspa;
		}

		jbrow = PROW( jb, grid );

#if ( PROFlevel>=1 )
		t = SuperLU_timer_();
#endif
		/*------------------------------------------------
		 * SET UP U BLOCKS.
		 *------------------------------------------------*/
		kseen = 0;
		/* Loop through each column in the block column. */
		for (j = fsupc; j < FstBlockC( jb+1 ); ++j) {
		    istart = xusub[j];
		    for (i = istart; i < xusub[j+1]; ++i) {
			irow = usub[i]; /* First nonzero in the segment. */
			gb = BlockNum( irow );
			pr = PROW( gb, grid );
			if ( pr != jbrow ) 
			    bsendx_plist[ljb][pr] = YES;
			if ( myrow == pr ) {
			    lb = LBi( gb, grid ); /* Local block number */
			    index = Ufstnz_br_ptr[lb];
			    if (rb_marker[lb] <= jb) {/* First see the block */
				rb_marker[lb] = jb + 1;
				index[Urb_indptr[lb]] = jb; /* Descriptor */
				Urb_indptr[lb] += UB_DESCRIPTOR;
				len = Urb_indptr[lb];
				for (k = 0; k < nsupc; ++k)
				    index[len+k] = FstBlockC( gb+1 );
				if ( gb != jb )/* Exclude diagonal block. */
				    ++bmod[lb];/* Mod. count for back solve */
				if ( kseen == 0 && myrow != jbrow ) {
				    ++nbrecvx;
				    kseen = 1;
				}
			    } else {
				len = Urb_indptr[lb];/* Start fstnz in index */
			    }
			    jj = j - fsupc;
			    index[len+jj] = irow;
			} /* if myrow == pr ... */
		    } /* for i ... */
		} /* for j ... */

		/* Figure out how many nonzeros in each block, and gather
		   the initial values of A from SPA into Uval. */
		for (lb = 0; lb < nrbu; ++lb) {
		    if ( rb_marker[lb] == jb + 1 ) { /* Not an empty block. */
			index = Ufstnz_br_ptr[lb];
			uval = Unzval_br_ptr[lb];
			len = Urb_indptr[lb];
			gb = lb * grid->nprow + myrow;/* Global block number */
			k = FstBlockC( gb+1 );
			irow = ilsum[lb] - FstBlockC( gb );
			for (jj=0, bnnz=0, dense_col=dense; jj < nsupc; ++jj) {
			    j = index[len+jj];  /* First nonzero in segment. */
			    bnnz += k - j;
			    for (i = j; i < k; ++i) {
				uval[Urb_length[lb]++] = dense_col[irow + i];
				dense_col[irow + i] = zero;
			    }
			    dense_col += ldaspa;
			}
			index[len-1] = bnnz; /* Set block length in Descriptor */
			Urb_indptr[lb] += nsupc;
		    }
		} /* for lb ... */
#if ( PROFlevel>=1 )
		t_u += SuperLU_timer_() - t;
		t = SuperLU_timer_();
#endif		
		/*------------------------------------------------
		 * SET UP L BLOCKS.
		 *------------------------------------------------*/

		/* Count number of blocks and length of each block. */
		nrbl = 0;
		len = 0; /* Number of row subscripts I own. */
		kseen = 0;
		istart = xlsub[fsupc];
		for (i = istart; i < xlsub[fsupc+1]; ++i) {
		    irow = lsub[i];
		    gb = BlockNum( irow ); /* Global block number */
		    pr = PROW( gb, grid ); /* Process row owning this block */
		    if ( pr != jbrow )
			fsendx_plist[ljb][pr] = YES;
		    if ( myrow == pr ) {
			lb = LBi( gb, grid );  /* Local block number */
			if (rb_marker[lb] <= jb) { /* First see this block */
			    rb_marker[lb] = jb + 1;
			    Lrb_length[lb] = 1;
			    Lrb_number[nrbl++] = gb;
			    if ( gb != jb ) /* Exclude diagonal block. */
				++fmod[lb]; /* Mod. count for forward solve */
			    if ( kseen == 0 && myrow != jbrow ) {
				++nfrecvx;
				kseen = 1;
			    }
#if ( PRNTlevel>=1 )
			    ++nLblocks;
#endif
			} else {
			    ++Lrb_length[lb];
			}
			++len;
		    }
		} /* for i ... */

		if ( nrbl ) { /* Do not ensure the blocks are sorted! */
		    /* Set up the initial pointers for each block in 
		       index[] and nzval[]. */
		    /* Add room for descriptors */
		    len1 = len + BC_HEADER + nrbl * LB_DESCRIPTOR;
		    if ( !(index = intMalloc_dist(len1)) ) 
			ABORT("Malloc fails for index[]");
		    Lrowind_bc_ptr[ljb] = index;
		    if (!(Lnzval_bc_ptr[ljb] = 
                         doublecomplexMalloc_dist(len*nsupc))) {
			fprintf(stderr, "col block %d ", jb);
			ABORT("Malloc fails for Lnzval_bc_ptr[*][]");
		    }
		    mybufmax[0] = SUPERLU_MAX( mybufmax[0], len1 );
		    mybufmax[1] = SUPERLU_MAX( mybufmax[1], len*nsupc );
		    mybufmax[4] = SUPERLU_MAX( mybufmax[4], len );
		    index[0] = nrbl;  /* Number of row blocks */
		    index[1] = len;   /* LDA of the nzval[] */
		    next_lind = BC_HEADER;
		    next_lval = 0;
		    for (k = 0; k < nrbl; ++k) {
			gb = Lrb_number[k];
			lb = LBi( gb, grid );
			len = Lrb_length[lb];
			Lrb_length[lb] = 0;  /* Reset vector of block length */
			index[next_lind++] = gb; /* Descriptor */
			index[next_lind++] = len; 
			Lrb_indptr[lb] = next_lind;
			Lrb_valptr[lb] = next_lval;
			next_lind += len;
			next_lval += len;
		    }
		    /* Propagate the compressed row subscripts to Lindex[],
                       and the initial values of A from SPA into Lnzval[]. */
		    lusup = Lnzval_bc_ptr[ljb];
		    len = index[1];  /* LDA of lusup[] */
		    for (i = istart; i < xlsub[fsupc+1]; ++i) {
			irow = lsub[i];
			gb = BlockNum( irow );
			if ( myrow == PROW( gb, grid ) ) {
			    lb = LBi( gb, grid );
			    k = Lrb_indptr[lb]++; /* Random access a block */
			    index[k] = irow;
			    k = Lrb_valptr[lb]++;
			    irow = ilsum[lb] + irow - FstBlockC( gb );
			    for (j = 0, dense_col = dense; j < nsupc; ++j) {
				lusup[k] = dense_col[irow];
				dense_col[irow] = zero;
				k += len;
				dense_col += ldaspa;
			    }
			}
		    } /* for i ... */
		} else {
		    Lrowind_bc_ptr[ljb] = NULL;
		    Lnzval_bc_ptr[ljb] = NULL;
		} /* if nrbl ... */
#if ( PROFlevel>=1 )
		t_l += SuperLU_timer_() - t;
#endif
	    } /* if mycol == pc */

	} /* for jb ... */

	Llu->Lrowind_bc_ptr = Lrowind_bc_ptr;
	Llu->Lnzval_bc_ptr = Lnzval_bc_ptr;
	Llu->Ufstnz_br_ptr = Ufstnz_br_ptr;
	Llu->Unzval_br_ptr = Unzval_br_ptr;
	Llu->ToRecv = ToRecv;
	Llu->ToSendD = ToSendD;
	Llu->ToSendR = ToSendR;
	Llu->fmod = fmod;
	Llu->fsendx_plist = fsendx_plist;
	Llu->nfrecvx = nfrecvx;
	Llu->bmod = bmod;
	Llu->bsendx_plist = bsendx_plist;
	Llu->nbrecvx = nbrecvx;
	Llu->ilsum = ilsum;
	Llu->ldalsum = ldaspa;
	
#if ( PRNTlevel>=1 )
	if ( !iam ) printf(".. # L blocks %d\t# U blocks %d\n",
			   nLblocks, nUblocks);
#endif

	SUPERLU_FREE(rb_marker);
	SUPERLU_FREE(Urb_length);
	SUPERLU_FREE(Urb_indptr);
	SUPERLU_FREE(Lrb_length);
	SUPERLU_FREE(Lrb_number);
	SUPERLU_FREE(Lrb_indptr);
	SUPERLU_FREE(Lrb_valptr);
	SUPERLU_FREE(dense);

	/* Find the maximum buffer size. */
	MPI_Allreduce(mybufmax, Llu->bufmax, NBUFFERS, mpi_int_t, 
		      MPI_MAX, grid->comm);
#if ( PROFlevel>=1 )
	if ( !iam ) printf(".. 1st distribute time: L %.2f\tU %.2f\tu_blks %d\tnrbu %d\n",
			   t_l, t_u, u_blks, nrbu);
#endif

    } /* else fact != SamePattern_SameRowPerm */

    SUPERLU_FREE(xa);
    SUPERLU_FREE(asub);
    SUPERLU_FREE(a);

#if ( DEBUGlevel>=1 )
    /* Memory allocated but not freed:
       ilsum, fmod, fsendx_plist, bmod, bsendx_plist  */
    CHECK_MALLOC(iam, "Exit pzdistribute()");
#endif
    
    return (mem_use);
} /* PZDISTRIBUTE */
