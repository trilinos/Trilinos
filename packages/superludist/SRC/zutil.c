
/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 */

/*
  Copyright (c) 1994 by Xerox Corporation.  All rights reserved.
 
  THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
  EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 
  Permission is hereby granted to use or copy this program for any
  purpose, provided the above notices are retained on all copies.
  Permission to modify the code and to distribute modified code is
  granted, provided the above notices are retained, and a notice that
  the code was modified is included with the above copyright notice.
*/

#include <math.h>
#include "superlu_zdefs.h"

void
zCreate_CompCol_Matrix_dist(SuperMatrix *A, int_t m, int_t n, int_t nnz, 
			    doublecomplex *nzval, int_t *rowind, int_t *colptr,
			    Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NCformat *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NCformat) );
    if ( !(A->Store) ) ABORT("SUPERLU_MALLOC fails for A->Store");
    Astore = (NCformat *) A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colptr = colptr;
}

void
zCreate_CompRowLoc_Matrix_dist(SuperMatrix *A, int_t m, int_t n,
			       int_t nnz_loc, int_t m_loc, int_t fst_row,
			       doublecomplex *nzval, int_t *colind, int_t *rowptr,
			       Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    NRformat_loc *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    A->Store = (void *) SUPERLU_MALLOC( sizeof(NRformat_loc) );
    if ( !(A->Store) ) ABORT("SUPERLU_MALLOC fails for A->Store");
    Astore = (NRformat_loc *) A->Store;
    Astore->nnz_loc = nnz_loc;
    Astore->fst_row = fst_row;
    Astore->m_loc = m_loc;
    Astore->nzval = nzval;
    Astore->colind = colind;
    Astore->rowptr = rowptr;
}

/*
 * Convert a row compressed storage into a column compressed storage.
 */
void
zCompRow_to_CompCol_dist(int_t m, int_t n, int_t nnz, 
                         doublecomplex *a, int_t *colind, int_t *rowptr,
                         doublecomplex **at, int_t **rowind, int_t **colptr)
{
    register int i, j, col, relpos;
    int_t *marker;

    /* Allocate storage for another copy of the matrix. */
    *at = (doublecomplex *) doublecomplexMalloc_dist(nnz);
    *rowind = intMalloc_dist(nnz);
    *colptr = intMalloc_dist(n+1);
    marker = intCalloc_dist(n);
    
    /* Get counts of each column of A, and set up column pointers */
    for (i = 0; i < m; ++i)
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
    (*colptr)[0] = 0;
    for (j = 0; j < n; ++j) {
	(*colptr)[j+1] = (*colptr)[j] + marker[j];
	marker[j] = (*colptr)[j];
    }

    /* Transfer the matrix into the compressed column storage. */
    for (i = 0; i < m; ++i) {
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
	    col = colind[j];
	    relpos = marker[col];
	    (*rowind)[relpos] = i;
	    (*at)[relpos] = a[j];
	    ++marker[col];
	}
    }

    SUPERLU_FREE(marker);
}

/* Copy matrix A into matrix B. */
void
zCopy_CompCol_Matrix_dist(SuperMatrix *A, SuperMatrix *B)
{
    NCformat *Astore, *Bstore;
    int      ncol, nnz, i;

    B->Stype = A->Stype;
    B->Dtype = A->Dtype;
    B->Mtype = A->Mtype;
    B->nrow  = A->nrow;;
    B->ncol  = ncol = A->ncol;
    Astore   = (NCformat *) A->Store;
    Bstore   = (NCformat *) B->Store;
    Bstore->nnz = nnz = Astore->nnz;
    for (i = 0; i < nnz; ++i)
	((doublecomplex *)Bstore->nzval)[i] = ((doublecomplex *)Astore->nzval)[i];
    for (i = 0; i < nnz; ++i) Bstore->rowind[i] = Astore->rowind[i];
    for (i = 0; i <= ncol; ++i) Bstore->colptr[i] = Astore->colptr[i];
}


void zPrint_CompCol_Matrix_dist(SuperMatrix *A)
{
    NCformat     *Astore;
    register int i;
    doublecomplex       *dp;
    
    printf("\nCompCol matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (NCformat *) A->Store;
    printf("nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
    if ( (dp = (doublecomplex *) Astore->nzval) != NULL ) {
        printf("nzval:\n");
        for (i = 0; i < Astore->nnz; ++i) printf("%f  ", dp[i]);
    }
    printf("\nrowind:\n");
    for (i = 0; i < Astore->nnz; ++i) printf("%d  ", Astore->rowind[i]);
    printf("\ncolptr:\n");
    for (i = 0; i <= A->ncol; ++i) printf("%d  ", Astore->colptr[i]);
    printf("\nend CompCol matrix.\n");
}

void zPrint_Dense_Matrix_dist(SuperMatrix *A)
{
    DNformat     *Astore;
    register int i;
    doublecomplex       *dp;
    
    printf("\nDense matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (DNformat *) A->Store;
    dp = (doublecomplex *) Astore->nzval;
    printf("nrow %d, ncol %d, lda %d\n", A->nrow,A->ncol,Astore->lda);
    printf("\nnzval: ");
    for (i = 0; i < A->nrow; ++i) printf("%f  ", dp[i]);
    printf("\nend Dense matrix.\n");
}

int zPrint_CompRowLoc_Matrix_dist(SuperMatrix *A)
{
    NRformat_loc     *Astore;
    int_t i, nnz_loc, m_loc;
    doublecomplex       *dp;
    
    printf("\n==== CompRowLoc matrix: ");
    printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (NRformat_loc *) A->Store;
    printf("nrow %d, ncol %d\n", A->nrow,A->ncol);
    nnz_loc = Astore->nnz_loc; m_loc = Astore->m_loc;
    printf("nnz_loc %d, m_loc %d, fst_row %d\n", nnz_loc, m_loc,
	   Astore->fst_row);
    PrintInt10("rowptr", m_loc+1, Astore->rowptr);
    PrintInt10("colind", nnz_loc, Astore->colind);
    if ( (dp = (doublecomplex *) Astore->nzval) != NULL )
        PrintDoublecomplex("nzval", nnz_loc, dp);
    printf("==== end CompRowLoc matrix\n");
}

int file_zPrint_CompRowLoc_Matrix_dist(FILE *fp, SuperMatrix *A)
{
    NRformat_loc     *Astore;
    int_t i, nnz_loc, m_loc;
    doublecomplex       *dp;
    
    fprintf(fp, "\n==== CompRowLoc matrix: ");
    fprintf(fp, "Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    Astore = (NRformat_loc *) A->Store;
    fprintf(fp, "nrow %d, ncol %d\n", A->nrow, A->ncol);
    nnz_loc = Astore->nnz_loc; m_loc = Astore->m_loc;
    fprintf(fp, "nnz_loc %d, m_loc %d, fst_row %d\n", nnz_loc, m_loc,
	    Astore->fst_row);
    file_PrintInt10(fp, "rowptr", m_loc+1, Astore->rowptr);
    file_PrintInt10(fp, "colind", nnz_loc, Astore->colind);
    if ( (dp = (doublecomplex *) Astore->nzval) != NULL )
        file_PrintDoublecomplex(fp, "nzval", nnz_loc, dp);
    fprintf(fp, "==== end CompRowLoc matrix\n");
}

void
zCreate_Dense_Matrix_dist(SuperMatrix *X, int_t m, int_t n, doublecomplex *x,
			  int_t ldx, Stype_t stype, Dtype_t dtype,
			  Mtype_t mtype)
{
    DNformat    *Xstore;
    
    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    X->Store = (void *) SUPERLU_MALLOC( sizeof(DNformat) );
    if ( !(X->Store) ) ABORT("SUPERLU_MALLOC fails for X->Store");
    Xstore = (DNformat *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = (doublecomplex *) x;
}

void
zCopy_Dense_Matrix_dist(int_t M, int_t N, doublecomplex *X, int_t ldx,
			doublecomplex *Y, int_t ldy)
{
/*
 *  Purpose
 *  =======
 *
 *  Copies a two-dimensional matrix X to another matrix Y.
 */
    int    i, j;
    
    for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i)
            Y[i + j*ldy] = X[i + j*ldx];
}

void
zCreate_SuperNode_Matrix_dist(SuperMatrix *L, int_t m, int_t n, int_t nnz, 
			      doublecomplex *nzval, int_t *nzval_colptr,
			      int_t *rowind, int_t *rowind_colptr,
			      int_t *col_to_sup, int_t *sup_to_col,
			      Stype_t stype, Dtype_t dtype, Mtype_t mtype)
{
    SCformat *Lstore;

    L->Stype = stype;
    L->Dtype = dtype;
    L->Mtype = mtype;
    L->nrow = m;
    L->ncol = n;
    L->Store = (void *) SUPERLU_MALLOC( sizeof(SCformat) );
    if ( !(L->Store) ) ABORT("SUPERLU_MALLOC fails for L->Store");
    Lstore = L->Store;
    Lstore->nnz = nnz;
    Lstore->nsuper = col_to_sup[n];
    Lstore->nzval = nzval;
    Lstore->nzval_colptr = nzval_colptr;
    Lstore->rowind = rowind;
    Lstore->rowind_colptr = rowind_colptr;
    Lstore->col_to_sup = col_to_sup;
    Lstore->sup_to_col = sup_to_col;

}

void
zGenXtrue_dist(int_t n, int_t nrhs, doublecomplex *x, int_t ldx)
{
    int  i, j;
    for (j = 0; j < nrhs; ++j)
	for (i = 0; i < n; ++i) {
	    if ( i % 2 ) x[i + j*ldx].r = 1.0;
	    else x[i + j*ldx].r = 2.0;
	    x[i + j*ldx].i = 0.0;
	}
}

/*
 * Let rhs[i] = sum of i-th row of A, so the solution vector is all 1's
 */
void
zFillRHS_dist(char *trans, int_t nrhs, doublecomplex *x, int_t ldx,
	      SuperMatrix *A, doublecomplex *rhs, int_t ldb)
{
    doublecomplex one = {1.0, 0.0};
    doublecomplex zero = {0.0, 0.0};

    sp_zgemm_dist(trans, "N", A->nrow, nrhs, A->ncol, one, A,
		  x, ldx, zero, rhs, ldb);

}

/* 
 * Fills a doublecomplex precision array with a given value.
 */
void 
zfill_dist(doublecomplex *a, int_t alen, doublecomplex dval)
{
    register int_t i;
    for (i = 0; i < alen; i++) a[i] = dval;
}



/* 
 * Check the inf-norm of the error vector 
 */
void zinf_norm_error_dist(int_t n, int_t nrhs, doublecomplex *x, int_t ldx,
			  doublecomplex *xtrue, int_t ldxtrue,
                          gridinfo_t *grid)
{
    double err, xnorm;
    doublecomplex *x_work, *xtrue_work;
    doublecomplex temp;
    int i, j;

    for (j = 0; j < nrhs; j++) {
      x_work = &x[j*ldx];
      xtrue_work = &xtrue[j*ldxtrue];
      err = xnorm = 0.0;
      for (i = 0; i < n; i++) {
        z_sub(&temp, &x_work[i], &xtrue_work[i]);
	err = SUPERLU_MAX(err, z_abs(&temp));
	xnorm = SUPERLU_MAX(xnorm, z_abs(&x_work[i]));
      }
      err = err / xnorm;
      printf("(%d) .. ||X-Xtrue||/||X|| = %e\n", grid->iam, err);
    }
}

void PrintDoublecomplex(char *name, int_t len, doublecomplex *x)
{
    register int_t i;
    
    printf("%10s:\tReal\tImag\n", name);
    for (i = 0; i < len; ++i)
	printf("\t%d\t%.4f\t%.4f\n", i, x[i].r, x[i].i);
}

int file_PrintDoublecomplex(FILE *fp, char *name, int_t len, doublecomplex *x)
{
    register int_t i;
    
    fprintf(fp, "%10s:\tReal\tImag\n", name);
    for (i = 0; i < len; ++i)
	fprintf(fp, "\t%d\t%.4f\t%.4f\n", i, x[i].r, x[i].i);
}

/* 
 * Print the blocks in the factored matrix L.
 */
void zPrintLblocks(int_t iam, int_t nsupers, gridinfo_t *grid,
		  Glu_persist_t *Glu_persist, LocalLU_t *Llu)
{
    register int_t c, extra, gb, j, lb, nsupc, nsupr, len, nb, ncb;
    register int_t k, mycol, r;
    int_t *xsup = Glu_persist->xsup;
    int_t *index;
    doublecomplex *nzval;

    printf("\n(%d) L BLOCKS IN COLUMN-MAJOR ORDER -->\n", iam);
    ncb = nsupers / grid->npcol;
    extra = nsupers % grid->npcol;
    mycol = MYCOL( iam, grid );
    if ( mycol < extra ) ++ncb;
    for (lb = 0; lb < ncb; ++lb) {
	index = Llu->Lrowind_bc_ptr[lb];
	if ( index ) { /* Not an empty column */
	    nzval = Llu->Lnzval_bc_ptr[lb];
	    nb = index[0];
	    nsupr = index[1];
	    gb = lb * grid->npcol + mycol;
	    nsupc = SuperSize( gb );
	    printf("(%d) block column (local) %d, # row blocks %d\n",
		   iam, lb, nb);
	    for (c = 0, k = BC_HEADER, r = 0; c < nb; ++c) {
		len = index[k+1];
		printf("(%d) row-block %d: block # %d\tlength %d\n", 
		       iam, c, index[k], len);
		PrintInt10("lsub", len, &index[k+LB_DESCRIPTOR]);
		for (j = 0; j < nsupc; ++j) {
		    PrintDoublecomplex("nzval", len, &nzval[r + j*nsupr]);
		}
		k += LB_DESCRIPTOR + len;
		r += len;
	    }
	}
	printf("(%d)", iam);
 	PrintInt10("ToSendR[]", grid->npcol, Llu->ToSendR[lb]);
	PrintInt10("fsendx_plist[]", grid->nprow, Llu->fsendx_plist[lb]);
    }
    printf("nfrecvx %4d\n", Llu->nfrecvx);
    k = CEILING( nsupers, grid->nprow );
    PrintInt10("fmod", k, Llu->fmod);
    
} /* ZPRINTLBLOCKS */


/* 
 * Print the blocks in the factored matrix U.
 */
void zPrintUblocks(int_t iam, int_t nsupers, gridinfo_t *grid, 
		  Glu_persist_t *Glu_persist, LocalLU_t *Llu)
{
    register int_t c, extra, jb, k, lb, len, nb, nrb, nsupc;
    register int_t myrow, r;
    int_t *xsup = Glu_persist->xsup;
    int_t *index;
    doublecomplex *nzval;

    printf("\n(%d) U BLOCKS IN ROW-MAJOR ORDER -->\n", iam);
    nrb = nsupers / grid->nprow;
    extra = nsupers % grid->nprow;
    myrow = MYROW( iam, grid );
    if ( myrow < extra ) ++nrb;
    for (lb = 0; lb < nrb; ++lb) {
	index = Llu->Ufstnz_br_ptr[lb];
	if ( index ) { /* Not an empty row */
	    nzval = Llu->Unzval_br_ptr[lb];
	    nb = index[0];
	    printf("(%d) block row (local) %d, # column blocks %d\n",
		   iam, lb, nb);
	    r  = 0;
	    for (c = 0, k = BR_HEADER; c < nb; ++c) {
		jb = index[k];
		len = index[k+1];
		printf("(%d) col-block %d: block # %d\tlength %d\n", 
		       iam, c, jb, index[k+1]);
		nsupc = SuperSize( jb );
		PrintInt10("fstnz", nsupc, &index[k+UB_DESCRIPTOR]);
		PrintDoublecomplex("nzval", len, &nzval[r]);
		k += UB_DESCRIPTOR + nsupc;
		r += len;
	    }

	    printf("(%d) ToSendD[] %d\n", iam, Llu->ToSendD[lb]);
	}
    }
} /* ZPRINTUBLOCKS */

int
zprint_gsmv_comm(FILE *fp, int_t m_loc, pzgsmv_comm_t *gsmv_comm,
                 gridinfo_t *grid)
{
  int_t procs = grid->nprow*grid->npcol;
  fprintf(fp, "TotalIndSend %d\tTotalValSend %d\n", gsmv_comm->TotalIndSend,
	  gsmv_comm->TotalValSend);
  file_PrintInt10(fp, "extern_start", m_loc, gsmv_comm->extern_start);
  file_PrintInt10(fp, "ind_tosend", gsmv_comm->TotalIndSend, gsmv_comm->ind_tosend);
  file_PrintInt10(fp, "ind_torecv", gsmv_comm->TotalValSend, gsmv_comm->ind_torecv);
  file_PrintInt10(fp, "ptr_ind_tosend", procs+1, gsmv_comm->ptr_ind_tosend);
  file_PrintInt10(fp, "ptr_ind_torecv", procs+1, gsmv_comm->ptr_ind_torecv);
  file_PrintInt10(fp, "SendCounts", procs, gsmv_comm->SendCounts);
  file_PrintInt10(fp, "RecvCounts", procs, gsmv_comm->RecvCounts);
}


/* cg5.cua
            b = A*x           y = L\b
   0        1 + 4.0000i       1.0000 + 4.0000i
   1        0 + 5.0000i	      1.3529 + 5.4118i
   2        1 + 4.0000i	      1.0000 + 4.0000i
   3        2 + 3.0000i	      2.0000 + 3.0000i
   4        1 + 4.0000i	      3.5882 + 4.3529i
   5        1 + 4.0000i	      4.1250 + 3.3202i
   6          + 5.0000i	      4.4640 + 3.8632i
   7        2 + 3.0000i	      2.0000 + 3.0000i
   8        2 + 3.0000i	      2.0000 + 3.0000i
   9        1 + 4.0000i	      1.0000 + 4.0000i
  10        1 + 4.0000i	      3.5882 + 4.3529i
  11          + 5.0000i	           0 + 5.0000i
  12        1 + 4.0000i	      5.1793 + 4.6604i
  13        2 + 3.0000i	      2.0000 + 3.0000i
  14        1 + 4.0000i	      1.0000 + 4.0000i
  15          + 5.0000i	      1.3529 + 5.4118i
  16        1 + 4.0000i	      4.0045 + 3.8950i
  17          + 5.0000i	      3.0338 + 4.6248i
  18        1 + 4.0000i	      5.4495 + 2.2703i
  19          + 5.0000i	      4.0980 + 3.7290i
  20          + 5.0000i	      4.2680 + 3.7739i
  21          + 5.0000i	      5.3514 + 2.9480i
  22        1 + 4.0000i	      4.4178 + 2.0476i
  23        1 + 4.0000i	      3.5615 + 2.8322i
  24          + 5.0000i	      4.7526 + 2.2605i
*/
