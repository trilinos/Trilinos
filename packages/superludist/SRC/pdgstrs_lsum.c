

/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 * Modified:
 *     Feburary 7, 2001    use MPI_Isend/MPI_Irecv
 *     October 2, 2001     use MPI_Isend/MPI_Irecv with MPI_Test
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

/************************************************************************/
void dlsum_fmod
/************************************************************************/
(
 double *lsum,    /* Sum of local modifications.                        */
 double *x,       /* X array (local)                                    */
 double *xk,      /* X[k].                                              */
 double *rtemp,   /* Result of full matrix-vector multiply.             */
 int   nrhs,      /* Number of right-hand sides.                        */
 int   knsupc,    /* Size of supernode k.                               */
 int_t k,         /* The k-th component of X.                           */
 int_t *fmod,     /* Modification count for L-solve.                    */
 int_t nlb,       /* Number of L blocks.                                */
 int_t lptr,      /* Starting position in lsub[*].                      */
 int_t luptr,     /* Starting position in lusup[*].                     */
 int_t *xsup,
 gridinfo_t *grid,
 LocalLU_t *Llu,
 MPI_Request send_req[],
 SuperLUStat_t *stat
)
{
/*
 * Purpose
 * =======
 *   Perform local block modifications: lsum[i] -= L_i,k * X[k].
 */
    double alpha = 1.0, beta = 0.0;
    double *lusup, *lusup1;
    double *dest;
    int    iam, iknsupc, myrow, nbrow, nsupr, nsupr1, p, pi;
    int_t  i, ii, ik, il, ikcol, irow, j, lb, lk, rel;
    int_t  *lsub, *lsub1, nlb1, lptr1, luptr1;
    int_t  *ilsum = Llu->ilsum; /* Starting position of each supernode in lsum.   */
    int_t  *frecv = Llu->frecv;
    int_t  **fsendx_plist = Llu->fsendx_plist;
    MPI_Status status;
    int test_flag;

    iam = grid->iam;
    myrow = MYROW( iam, grid );
    lk = LBj( k, grid ); /* Local block number, column-wise. */
    lsub = Llu->Lrowind_bc_ptr[lk];
    lusup = Llu->Lnzval_bc_ptr[lk];
    nsupr = lsub[1];

    for (lb = 0; lb < nlb; ++lb) {
	ik = lsub[lptr]; /* Global block number, row-wise. */
	nbrow = lsub[lptr+1];
#ifdef _CRAY
	SGEMM( ftcs2, ftcs2, &nbrow, &nrhs, &knsupc,
	      &alpha, &lusup[luptr], &nsupr, xk,
	      &knsupc, &beta, rtemp, &nbrow );
#else
	dgemm_( "N", "N", &nbrow, &nrhs, &knsupc,
	       &alpha, &lusup[luptr], &nsupr, xk,
	       &knsupc, &beta, rtemp, &nbrow );
#endif
	stat->ops[SOLVE] += 2 * nbrow * nrhs * knsupc + nbrow * nrhs;
   
	lk = LBi( ik, grid ); /* Local block number, row-wise. */
	iknsupc = SuperSize( ik );
	il = LSUM_BLK( lk );
	dest = &lsum[il];
	lptr += LB_DESCRIPTOR;
	rel = xsup[ik]; /* Global row index of block ik. */
	for (i = 0; i < nbrow; ++i) {
	    irow = lsub[lptr++] - rel; /* Relative row. */
	    RHS_ITERATE(j)
		dest[irow + j*iknsupc] -= rtemp[i + j*nbrow];
	}
	luptr += nbrow;
		    
	if ( (--fmod[lk])==0 ) { /* Local accumulation done. */
	    ikcol = PCOL( ik, grid );
	    p = PNUM( myrow, ikcol, grid );
	    if ( iam != p ) {
#ifdef ISEND_IRECV
#if 1
	        MPI_Test( &send_req[myrow], &test_flag, &status );
#else
	        if ( send_req[myrow] != MPI_REQUEST_NULL ) 
		    MPI_Wait( &send_req[myrow], &status );
#endif
		MPI_Isend( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			 MPI_DOUBLE, p, LSUM, grid->comm, &send_req[myrow] );
#else
		MPI_Send( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			 MPI_DOUBLE, p, LSUM, grid->comm );
#endif
#if ( DEBUGlevel>=2 )
		printf("(%2d) Sent LSUM[%2.0f], size %2d, to P %2d\n",
		       iam, lsum[il-LSUM_H], iknsupc*nrhs+LSUM_H, p);
#endif
	    } else { /* Diagonal process: X[i] += lsum[i]. */
		ii = X_BLK( lk );
		RHS_ITERATE(j)
		    for (i = 0; i < iknsupc; ++i)
			x[i + ii + j*iknsupc] += lsum[i + il + j*iknsupc];
		if ( frecv[lk]==0 ) { /* Becomes a leaf node. */
		    fmod[lk] = -1; /* Do not solve X[k] in the future. */
		    lk = LBj( ik, grid );/* Local block number, column-wise. */
		    lsub1 = Llu->Lrowind_bc_ptr[lk];
		    lusup1 = Llu->Lnzval_bc_ptr[lk];
		    nsupr1 = lsub1[1];
#ifdef _CRAY
		    STRSM(ftcs1, ftcs1, ftcs2, ftcs3, &iknsupc, &nrhs, &alpha,
			  lusup1, &nsupr1, &x[ii], &iknsupc);
#else
		    dtrsm_("L", "L", "N", "U", &iknsupc, &nrhs, &alpha, 
			   lusup1, &nsupr1, &x[ii], &iknsupc);
#endif
		    stat->ops[SOLVE] += iknsupc * (iknsupc - 1) * nrhs;
#if ( DEBUGlevel>=2 )
		    printf("(%2d) Solve X[%2d]\n", iam, ik);
#endif
		
		    /*
		     * Send Xk to process column Pc[k].
		     */
		    for (p = 0; p < grid->nprow; ++p)
			if ( fsendx_plist[lk][p] != EMPTY ) {
			    pi = PNUM( p, ikcol, grid );
#ifdef ISEND_IRECV
#if 1	      
			    MPI_Test( &send_req[p], &test_flag, &status );
#else
			    if ( send_req[p] != MPI_REQUEST_NULL ) 
			        MPI_Wait( &send_req[p], &status );
#endif
			    MPI_Isend( &x[ii - XK_H], iknsupc * nrhs + XK_H,
				      MPI_DOUBLE, pi, Xk, grid->comm,
				      &send_req[p] );
#else
			    MPI_Send( &x[ii - XK_H], iknsupc * nrhs + XK_H,
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
		    nlb1 = lsub1[0] - 1;
		    lptr1 = BC_HEADER + LB_DESCRIPTOR + iknsupc;
		    luptr1 = iknsupc; /* Skip diagonal block L(I,I). */

		    dlsum_fmod(lsum, x, &x[ii], rtemp, nrhs, iknsupc, ik,
			       fmod, nlb1, lptr1, luptr1, xsup,
			       grid, Llu, send_req, stat);
#ifdef ISEND_IRECV
		    /* Wait for previous Isends to complete. */
		    for (p = 0; p < grid->nprow; ++p) {
			if ( fsendx_plist[lk][p] != EMPTY )
			    /*MPI_Wait( &send_req[p], &status );*/
			    MPI_Test( &send_req[p], &test_flag, &status );
		    }
#endif
		} /* if frecv[lk] == 0 */
	    } /* if iam == p */
	} /* if fmod[lk] == 0 */

    } /* for lb ... */

} /* dLSUM_FMOD */


/************************************************************************/
void dlsum_bmod
/************************************************************************/
(
 double *lsum,        /* Sum of local modifications.                    */
 double *x,           /* X array (local).                               */
 double *xk,          /* X[k].                                          */
 int    nrhs,	      /* Number of right-hand sides.                    */
 int_t  k,            /* The k-th component of X.                       */
 int_t  *bmod,        /* Modification count for L-solve.                */
 int_t  *Urbs,        /* Number of row blocks in each block column of U.*/
 Ucb_indptr_t **Ucb_indptr,/* Vertical linked list pointing to Uindex[].*/
 int_t  **Ucb_valptr, /* Vertical linked list pointing to Unzval[].     */
 int_t  *xsup,
 gridinfo_t *grid,
 LocalLU_t *Llu,
 MPI_Request send_req[],
 SuperLUStat_t *stat
 )
{
/*
 * Purpose
 * =======
 *   Perform local block modifications: lsum[i] -= U_i,k * X[k].
 */
    double alpha = 1.0;
    int    iam, iknsupc, knsupc, myrow, nsupr, p, pi;
    int_t  fnz, gik, gikcol, i, ii, ik, ikfrow, iklrow, il, irow,
           j, jj, lk, lk1, nub, ub, uptr;
    int_t  *usub;
    double *uval, *dest, *y;
    int_t  *lsub;
    double *lusup;
    int_t  *ilsum = Llu->ilsum; /* Starting position of each supernode in lsum.   */
    int_t  *brecv = Llu->brecv;
    int_t  **bsendx_plist = Llu->bsendx_plist;
    MPI_Status status;
    int test_flag;

    iam = grid->iam;
    myrow = MYROW( iam, grid );
    knsupc = SuperSize( k );
    lk = LBj( k, grid ); /* Local block number, column-wise. */
    nub = Urbs[lk];

    for (ub = 0; ub < nub; ++ub) {
	ik = Ucb_indptr[lk][ub].lbnum; /* Local block number, row-wise. */
	usub = Llu->Ufstnz_br_ptr[ik];
	uval = Llu->Unzval_br_ptr[ik];
	i = Ucb_indptr[lk][ub].indpos; /* Start of the block in usub[]. */
	i += UB_DESCRIPTOR;
	il = LSUM_BLK( ik );
	gik = ik * grid->nprow + myrow;   /* Global block number, row-wise. */
	iknsupc = SuperSize( gik );
	ikfrow = FstBlockC( gik );
	iklrow = FstBlockC( gik+1 );

	RHS_ITERATE(j) {
	    dest = &lsum[il + j*iknsupc];
	    y = &xk[j*knsupc];
	    uptr = Ucb_valptr[lk][ub]; /* Start of the block in uval[]. */
	    for (jj = 0; jj < knsupc; ++jj) {
		fnz = usub[i + jj];
		if ( fnz < iklrow ) { /* Nonzero segment. */
		    /* AXPY */
		    for (irow = fnz; irow < iklrow; ++irow)
			dest[irow - ikfrow] -= uval[uptr++] * y[jj];
		    stat->ops[SOLVE] += 2 * (iklrow - fnz);
		}
	    } /* for jj ... */
	}

	--bmod[ik];
	if ( !(bmod[ik]) ) { /* Local accumulation done. */
	    gikcol = PCOL( gik, grid );
	    p = PNUM( myrow, gikcol, grid );
	    if ( iam != p ) {
#ifdef ISEND_IRECV
#if 1
	        MPI_Test( &send_req[myrow], &test_flag, &status );
#else
	        if ( send_req[myrow] != MPI_REQUEST_NULL ) 
		    MPI_Wait( &send_req[myrow], &status );
#endif
		MPI_Isend( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			  MPI_DOUBLE, p, LSUM, grid->comm, &send_req[myrow] );
#else
		MPI_Send( &lsum[il - LSUM_H], iknsupc * nrhs + LSUM_H,
			 MPI_DOUBLE, p, LSUM, grid->comm );
#endif
#if ( DEBUGlevel>=2 )
		printf("(%2d) Sent LSUM[%2.0f], size %2d, to P %2d\n",
		       iam, lsum[il-LSUM_H], iknsupc*nrhs+LSUM_H, p);
#endif
	    } else { /* Diagonal process: X[i] += lsum[i]. */
		ii = X_BLK( ik );
		dest = &x[ii];
		RHS_ITERATE(j)
		    for (i = 0; i < iknsupc; ++i)
			dest[i + j*iknsupc] += lsum[i + il + j*iknsupc];
		if ( !brecv[ik] ) { /* Becomes a leaf node. */
		    bmod[ik] = -1; /* Do not solve X[k] in the future. */
		    lk1 = LBj( gik, grid ); /* Local block number. */
		    lsub = Llu->Lrowind_bc_ptr[lk1];
		    lusup = Llu->Lnzval_bc_ptr[lk1];
		    nsupr = lsub[1];
#ifdef _CRAY
		    STRSM(ftcs1, ftcs3, ftcs2, ftcs2, &iknsupc, &nrhs, &alpha,
			  lusup, &nsupr, &x[ii], &iknsupc);
#else
		    dtrsm_("L", "U", "N", "N", &iknsupc, &nrhs, &alpha, 
			   lusup, &nsupr, &x[ii], &iknsupc);
#endif
		    stat->ops[SOLVE] += iknsupc * (iknsupc + 1) * nrhs;
#if ( DEBUGlevel>=2 )
		    printf("(%2d) Solve X[%2d]\n", iam, gik);
#endif

		    /*
		     * Send Xk to process column Pc[k].
		     */
		    for (p = 0; p < grid->nprow; ++p)
			if ( bsendx_plist[lk1][p] != EMPTY ) {
			    pi = PNUM( p, gikcol, grid );
#ifdef ISEND_IRECV
#if 1
			    MPI_Test( &send_req[p], &test_flag, &status );
#else
			    if ( send_req[p] != MPI_REQUEST_NULL ) 
				MPI_Wait( &send_req[p], &status );
#endif
			    MPI_Isend( &x[ii - XK_H], iknsupc * nrhs + XK_H,
				     MPI_DOUBLE, pi, Xk, grid->comm,
				     &send_req[p] );
#else
			    MPI_Send( &x[ii - XK_H], iknsupc * nrhs + XK_H,
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
		    if ( Urbs[lk1] )
			dlsum_bmod(lsum, x, &x[ii], nrhs, gik, bmod, Urbs,
				   Ucb_indptr, Ucb_valptr, xsup, grid, Llu,
				   send_req, stat);
#ifdef ISEND_IRECV
		    /* Wait for the previous Isends to complete. */
		    for (p = 0; p < grid->nprow; ++p) {
			if ( bsendx_plist[lk1][p] != EMPTY )
			    /*MPI_Wait( &send_req[p], &status );*/
			    MPI_Test( &send_req[p], &test_flag, &status );
		    }
#endif
		} /* if brecv[ik] == 0 */
	    }
	} /* if bmod[ik] == 0 */

    } /* for ub ... */

} /* dlSUM_BMOD */

