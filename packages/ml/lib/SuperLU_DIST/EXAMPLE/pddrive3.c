/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include <math.h>
#include "superlu_ddefs.h"


main(int argc, char *argv[])
/*
 * Purpose
 * =======
 *
 * The driver program PDDRIVE3.
 *
 * This example illustrates how to use pdgssvx_ABglobal to solve
 * systems repeatedly with the same sparsity pattern and similar
 * numerical values of matrix A.
 * In this case, the symbolic factorization is performed only once.
 * The following data structures will be reused in the subsequent
 * call to pdgssvx_ABglobal:
 *        ScalePermstruct : DiagScale, R, C, perm_r, perm_c
 *        LUstruct        : Glu_persist, Llu
 *
 * NOTE:
 * The distributed nonzero structures of L and U remain the same,
 * although the numerical values are different. So 'Llu' is set up once
 * in the first call to pdgssvx_ABglobal, and reused in the subsequent call.
 *
 * On the Cray T3E, the program may be run by typing:
 *    mpprun -n <procs> pddrive3 -r <proc rows> -c <proc columns> <input_matrix>
 */
{
    superlu_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    gridinfo_t grid;
    double   *berr;
    double   *a, *a1, *b, *b1, *xtrue;
    int_t    *asub, *asub1, *xa, *xa1;
    int_t    i, j, m, n, nnz;
    int_t    nprow, npcol;
    int      iam, info, ldb, ldx, nrhs;
    char     trans[1];
    char     **cpp, c;
    FILE *fp, *fopen();


    nprow = 1;  /* Default process rows.      */
    npcol = 1;  /* Default process columns.   */
    nrhs = 1;   /* Number of right-hand side. */

    /* Parse command line argv[]. */
    for (cpp = argv+1; *cpp; ++cpp) {
	if ( **cpp == '-' ) {
	    c = *(*cpp+1);
	    ++cpp;
	    switch (c) {
	      case 'h':
		  printf("Options:\n");
		  printf("\t-r <int>: process rows    (default %d)\n", nprow);
		  printf("\t-c <int>: process columns (default %d)\n", npcol);
		  exit(0);
		  break;
	      case 'r': nprow = atoi(*cpp);
		        break;
	      case 'c': npcol = atoi(*cpp);
		        break;
	    }
	} else { /* Last arg is considered a filename */
	    if ( !(fp = fopen(*cpp, "r")) ) {
                fprintf(stderr, "File does not exist.");
                exit(-1);
            }
	    break;
	}
    }

    /* ------------------------------------------------------------
       INITIALIZE MPI ENVIRONMENT. 
       ------------------------------------------------------------*/
    MPI_Init( &argc, &argv );

    /* ------------------------------------------------------------
       INITIALIZE THE SUPERLU PROCESS GRID. 
       ------------------------------------------------------------*/
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, &grid);

    /* Bail out if I do not belong in the grid. */
    iam = grid.iam;
    if ( iam >= nprow * npcol )
	goto out;
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter main()");
#endif

    /* ------------------------------------------------------------
       PROCESS 0 READS THE MATRIX A, AND THEN BROADCASTS IT TO ALL
       THE OTHER PROCESSES.
       ------------------------------------------------------------*/
    if ( !iam ) {
	/* Print the CPP definitions. */
	cpp_defs();
	
	/* Read the matrix stored on disk in Harwell-Boeing format. */
	dreadhb(iam, fp, &m, &n, &nnz, &a, &asub, &xa);
	
	printf("\tDimension\t%dx%d\t # nonzeros %d\n", m, n, nnz);
	printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);

	/* Broadcast matrix A to the other PEs. */
	MPI_Bcast( &m,   1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( &n,   1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( a,    nnz, MPI_DOUBLE, 0, grid.comm );
	MPI_Bcast( asub, nnz, mpi_int_t,  0, grid.comm );
	MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid.comm );
    } else {
	/* Receive matrix A from PE 0. */
	MPI_Bcast( &m,   1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( &n,   1,   mpi_int_t,  0, grid.comm );
	MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid.comm );

	/* Allocate storage for compressed column representation. */
	dallocateA(n, nnz, &a, &asub, &xa);

	MPI_Bcast( a,    nnz, MPI_DOUBLE, 0, grid.comm );
	MPI_Bcast( asub, nnz, mpi_int_t,  0, grid.comm );
	MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid.comm );
    }
	
    /* Create compressed column matrix for A. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, NC, _D, GE);

    /* Generate the exact solution and compute the right-hand side. */
    if ( !(b = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for b[]");
    if ( !(xtrue = doubleMalloc(n * nrhs)) ) ABORT("Malloc fails for xtrue[]");
    *trans = 'N';
    ldx = n;
    ldb = m;
    dGenXtrue(n, nrhs, xtrue, ldx);
    dFillRHS(trans, nrhs, xtrue, ldx, &A, b, ldb);

    /* Save a copy of the right-hand side. */  
    if ( !(b1 = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for b1[]");
    for (j = 0; j < nrhs; ++j)
	for (i = 0; i < m; ++i) b1[i+j*ldb] = b[i+j*ldb];
    
    if ( !(berr = doubleMalloc(nrhs)) )
	ABORT("Malloc fails for berr[].");

    /* Save a copy of the matrix A. */
    dallocateA(n, nnz, &a1, &asub1, &xa1);
    for (i = 0; i < nnz; ++i) { a1[i] = a[i]; asub1[i] = asub[i]; }
    for (i = 0; i < n+1; ++i) xa1[i] = xa[i];


    /* ------------------------------------------------------------
       WE SOLVE THE LINEAR SYSTEM FOR THE FIRST TIME.
       ------------------------------------------------------------*/

    /* Set the default input options. */
    set_default_options(&options);

    /* Initialize ScalePermstruct and LUstruct. */
    ScalePermstructInit(m, n, &ScalePermstruct);
    LUstructInit(m, n, &LUstruct);

    /* Initialize the statistics variables. */
    PStatInit(&stat);

    /* Call the linear equation solver: factorize and solve. */
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid,
		     &LUstruct, berr, &stat, &info);

    /* Check the accuracy of the solution. */
    if ( !iam ) {
	dinf_norm_error(n, nrhs, b, ldb, xtrue, ldx);
    }
    
    
    PStatPrint(&stat, &grid);        /* Print the statistics. */
    PStatFree(&stat);
    Destroy_CompCol_Matrix(&A);      /* Deallocate storage of matrix A.     */
    SUPERLU_FREE(b);                 /* Free storage of right-hand side.    */


    /* ------------------------------------------------------------
       NOW WE SOLVE ANOTHER LINEAR SYSTEM.
       THE MATRIX A HAS THE SAME SPARSITY PATTERN AND THE SIMILAR
       NUMERICAL VALUES AS THAT IN A PREVIOUS SYSTEM.
       ------------------------------------------------------------*/
    options.Fact = SamePattern_SameRowPerm;
    PStatInit(&stat); /* Initialize the statistics variables. */

    /* Create compressed column matrix for A. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, a1, asub1, xa1, NC, _D, GE);

    /* Solve the linear system. */
    pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b1, ldb, nrhs, &grid,
		     &LUstruct, berr, &stat, &info);

    /* Check the accuracy of the solution. */
    if ( !iam ) {
	printf("Solve a system with the same sparsity pattern and similar values.\n");
	dinf_norm_error(n, nrhs, b1, ldb, xtrue, ldx);
    }

    /* Print the statistics. */
    PStatPrint(&stat, &grid);

    /* ------------------------------------------------------------
       DEALLOCATE STORAGE.
       ------------------------------------------------------------*/
    PStatFree(&stat);
    Destroy_CompCol_Matrix(&A);      /* Deallocate storage of matrix A.     */
    Destroy_LU(n, &grid, &LUstruct); /* Deallocate storage associated with    
					the L and U matrices.               */
    ScalePermstructFree(&ScalePermstruct);
    LUstructFree(&LUstruct);         /* Deallocate the structure of L and U.*/
    SUPERLU_FREE(b1);	             /* Free storage of right-hand side.    */
    SUPERLU_FREE(xtrue);             /* Free storage of the exact solution. */
    SUPERLU_FREE(berr);


    /* ------------------------------------------------------------
       RELEASE THE SUPERLU PROCESS GRID.
       ------------------------------------------------------------*/
out:
    superlu_gridexit(&grid);

    /* ------------------------------------------------------------
       TERMINATES THE MPI EXECUTION ENVIRONMENT.
       ------------------------------------------------------------*/
    MPI_Finalize();


#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit main()");
#endif

}


int cpp_defs()
{
    printf(".. CPP definitions:\n");
#if ( PRNTlevel>=1 )
    printf("\tPRNTlevel = %d\n", PRNTlevel);
#endif
#if ( DEBUGlevel>=1 )
    printf("\tDEBUGlevel = %d\n", DEBUGlevel);
#endif
#if ( PROFlevel>=1 )
    printf("\tPROFlevel = %d\n", PROFlevel);
#endif
#if ( StaticPivot>=1 )
    printf("\tStaticPivot = %d\n", StaticPivot);
#endif
    printf("....\n");
}
