
#include <math.h>
#include "superlu_zdefs.h"

int zcreate_matrix(SuperMatrix *A, int nrhs, doublecomplex **rhs,
                   int *ldb, doublecomplex **x, int *ldx,
                   FILE *fp, gridinfo_t *grid)
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
 * ZCREATE_MATRIX read the matrix from data file in Harwell-Boeing format,
 * and distribute it to processors in a distributed compressed row format.
 * It also generate the distributed true solution X and the right-hand
 * side RHS.
 *
 *
 * Arguments   
 * =========      
 *
 * A     (output) SuperMatrix*
 *       Local matrix A in NR_loc format. 
 *
 * NRHS  (input) int_t
 *       Number of right-hand sides.
 *
 * RHS   (output) doublecomplex**
 *       The right-hand side matrix.
 *
 * LDB   (output) int*
 *       Leading dimension of the right-hand side matrix.
 *
 * X     (output) doublecomplex**
 *       The true solution matrix.
 *
 * LDX   (output) int*
 *       The leading dimension of the true solution matrix.
 *
 * FP    (input) FILE*
 *       The matrix file pointer.
 *
 * GRID  (input) gridinof_t*
 *       The 2D process mesh.
 *
 */

    SuperMatrix GA;              /* global A */
    doublecomplex   *b_global, *xtrue_global;  /* replicated on all processes */
    int_t    *rowind, *colptr;	 /* global */
    doublecomplex   *nzval;             /* global */
    doublecomplex   *nzval_loc;         /* local */
    int_t    *colind, *rowptr;	 /* local */
    int_t    m, n, nnz;
    int_t    m_loc, fst_row, nnz_loc;
    int_t    m_loc_fst; /* Record m_loc of the first p-1 processors,
			   when mod(m, p) is not zero. */ 
    int_t    iam, row, col, i, j, relpos;
    char     trans[1];
    int_t      *marker;

    iam = grid->iam;

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter zcreate_matrix()");
#endif

    if ( !iam ) {
        /* Read the matrix stored on disk in Harwell-Boeing format. */
        zreadhb_dist(iam, fp, &m, &n, &nnz, &nzval, &rowind, &colptr);

	/* Broadcast matrix A to the other PEs. */
	MPI_Bcast( &m,     1,   mpi_int_t,  0, grid->comm );
	MPI_Bcast( &n,     1,   mpi_int_t,  0, grid->comm );
	MPI_Bcast( &nnz,   1,   mpi_int_t,  0, grid->comm );
	MPI_Bcast( nzval,  nnz, SuperLU_MPI_DOUBLE_COMPLEX, 0, grid->comm );
	MPI_Bcast( rowind, nnz, mpi_int_t,  0, grid->comm );
	MPI_Bcast( colptr, n+1, mpi_int_t,  0, grid->comm );
    } else {
	/* Receive matrix A from PE 0. */
	MPI_Bcast( &m,   1,   mpi_int_t,  0, grid->comm );
	MPI_Bcast( &n,   1,   mpi_int_t,  0, grid->comm );
	MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid->comm );

	/* Allocate storage for compressed column representation. */
	zallocateA_dist(n, nnz, &nzval, &rowind, &colptr);

	MPI_Bcast( nzval,   nnz, SuperLU_MPI_DOUBLE_COMPLEX, 0, grid->comm );
	MPI_Bcast( rowind,  nnz, mpi_int_t,  0, grid->comm );
	MPI_Bcast( colptr,  n+1, mpi_int_t,  0, grid->comm );
    }

#if 0
    nzval[0].r = 0.1; nzval[0].i = 0.0;
#endif

    /* Compute the number of rows to be distributed to local process */
    m_loc = m / (grid->nprow * grid->npcol); 
    m_loc_fst = m_loc;
    /* When m / procs is not an integer */
    if ((m_loc * grid->nprow * grid->npcol) != m) {
      m_loc = m_loc+1;
      m_loc_fst = m_loc;
      if (iam == (grid->nprow * grid->npcol - 1)) 
	m_loc = m - m_loc_fst * (grid->nprow * grid->npcol - 1);
    }

    /* Create compressed column matrix for GA. */
    zCreate_CompCol_Matrix_dist(&GA, m, n, nnz, nzval, rowind, colptr,
				SLU_NC, SLU_Z, SLU_GE);

    /* Generate the exact solution and compute the right-hand side. */
    if ( !(b_global = doublecomplexMalloc_dist(m*nrhs)) )
        ABORT("Malloc fails for b[]");
    if ( !(xtrue_global = doublecomplexMalloc_dist(n*nrhs)) )
        ABORT("Malloc fails for xtrue[]");
    *trans = 'N';

    zGenXtrue_dist(n, nrhs, xtrue_global, n);
    zFillRHS_dist(trans, nrhs, xtrue_global, n, &GA, b_global, m);

    /*************************************************
     * Change GA to a local A with NR_loc format     *
     *************************************************/

    rowptr = (int_t *) intMalloc_dist(m_loc+1);
    marker = (int_t *) intCalloc_dist(n);

    /* Get counts of each row of GA */
    for (i = 0; i < n; ++i)
      for (j = colptr[i]; j < colptr[i+1]; ++j) ++marker[rowind[j]];
    /* Set up row pointers */
    rowptr[0] = 0;
    fst_row = iam * m_loc_fst;
    nnz_loc = 0;
    for (j = 0; j < m_loc; ++j) {
      row = fst_row + j;
      rowptr[j+1] = rowptr[j] + marker[row];
      marker[j] = rowptr[j];
    }
    nnz_loc = rowptr[m_loc];

    nzval_loc = (doublecomplex *) doublecomplexMalloc_dist(nnz_loc);
    colind = (int_t *) intMalloc_dist(nnz_loc);

    /* Transfer the matrix into the compressed row storage */
    for (i = 0; i < n; ++i) {
      for (j = colptr[i]; j < colptr[i+1]; ++j) {
	row = rowind[j];
	if ( (row>=fst_row) && (row<fst_row+m_loc) ) {
	  row = row - fst_row;
	  relpos = marker[row];
	  colind[relpos] = i;
	  nzval_loc[relpos] = nzval[j];
	  ++marker[row];
	}
      }
    }

#if ( DEBUGlevel>=2 )
    if ( !iam ) dPrint_CompCol_Matrix_dist(&GA);
#endif   

    /* Destroy GA */
    Destroy_CompCol_Matrix_dist(&GA);

    /******************************************************/
    /* Change GA to a local A with NR_loc format */
    /******************************************************/

    /* Set up the local A in NR_loc format */
    zCreate_CompRowLoc_Matrix_dist(A, m, n, nnz_loc, m_loc, fst_row,
				   nzval_loc, colind, rowptr,
				   SLU_NR_loc, SLU_Z, SLU_GE);
    
    /* Get the local B */
    if ( !((*rhs) = doublecomplexMalloc_dist(m_loc*nrhs)) )
        ABORT("Malloc fails for rhs[]");
    for (j =0; j < nrhs; ++j) {
	for (i = 0; i < m_loc; ++i) {
	    row = fst_row + i;
	    (*rhs)[j*m_loc+i] = b_global[j*n+row];
	}
    }
    *ldb = m_loc;

    /* Set the all true X */    
    *ldx = m_loc;
    if ( !((*x) = doublecomplexMalloc_dist(*ldx * nrhs)) )
        ABORT("Malloc fails for x_loc[]");

    /* Get the local part of xtrue_global */
    for (j = 0; j < nrhs; ++j) {
      for (i = 0; i < m_loc; ++i)
	(*x)[i + j*(*ldx)] = xtrue_global[i + fst_row + j*n];
    }

    SUPERLU_FREE(b_global);
    SUPERLU_FREE(xtrue_global);
    SUPERLU_FREE(marker);

#if ( DEBUGlevel>=1 )
    printf("sizeof(NRforamt_loc) %d\n", sizeof(NRformat_loc));
    CHECK_MALLOC(iam, "Exit zcreate_matrix()");
#endif
    return 0;
}
