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
 * The driver program PDDRIVE4.
 *
 * This example illustrates how to divide up the processes into
 * subgroups (multiple grids) such that each subgroup solves a linear
 * system independently from the other.
 *
 * In this example, there are 2 subgroups:
 *  1. subgroup 1 consists of processes 0 to 5 arranged as
 *     a 2-by-3 process grid.
 *  2. subgroup 2 consists of processes 6 to 9 arranged as
 *     a 2-by-2 process grid.
 *
 * On the Cray T3E, the program may be run by typing
 *    mpprun -n 10 pddrive4 <input_file>
 */
{
    superlu_options_t options;
    SuperLUStat_t stat;
    SuperMatrix A;
    ScalePermstruct_t ScalePermstruct;
    LUstruct_t LUstruct;
    gridinfo_t grid1, grid2;
    double   *berr;
    double   *a, *b, *xtrue;
    int_t    *asub, *xa;
    int_t    i, j, m, n, nnz;
    int_t    nprow, npcol, ldumap, p;
    int_t    usermap[6];
    int      iam, info, ldb, ldx, nprocs;
    int      nrhs = 1;   /* Number of right-hand side. */
    char     trans[1];
    char     **cpp, c;
    FILE *fp, *fopen();


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
                ABORT("File does not exist");
            }
	    break;
	}
    }

    /* ------------------------------------------------------------
       INITIALIZE MPI ENVIRONMENT. 
       ------------------------------------------------------------*/
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    if ( nprocs < 10 ) {
	fprintf(stderr, "Requires at least 10 processes\n");
	exit(-1);
    }

    /* ------------------------------------------------------------
       INITIALIZE THE SUPERLU PROCESS GRID 1. 
       ------------------------------------------------------------*/
    nprow = 2;
    npcol = 3;
    ldumap = 2;
    p = 0;    /* Grid 1 starts from process 0. */
    for (i = 0; i < nprow; ++i)
	for (j = 0; j < npcol; ++j) usermap[i+j*ldumap] = p++;
    superlu_gridmap(MPI_COMM_WORLD, nprow, npcol, usermap, ldumap, &grid1);

    grid1.ngrids = 2;
    grid1.mygrid = 0;


    /* ------------------------------------------------------------
       INITIALIZE THE SUPERLU PROCESS GRID 2. 
       ------------------------------------------------------------*/
    nprow = 2;
    npcol = 2;
    ldumap = 2;
    p = 6;   /* Grid 2 starts from process 6. */
    for (i = 0; i < nprow; ++i)
	for (j = 0; j < npcol; ++j) usermap[i+j*ldumap] = p++;
    superlu_gridmap(MPI_COMM_WORLD, nprow, npcol, usermap, ldumap, &grid2);

    grid2.ngrids = 2;
    grid2.mygrid = 1;

    /* Bail out if I do not belong in any of the 2 grids. */
    MPI_Comm_rank( MPI_COMM_WORLD, &iam );
    if ( iam >= 10 ) goto out;
    
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Enter main()");
#endif

    if ( iam >= 0 && iam < 6 ) { /* I am in grid 1. */
	iam = grid1.iam;  /* Get the logical number in the new grid. */

	/* ------------------------------------------------------------
	   PROCESS 0 READS THE MATRIX A, AND THEN BROADCASTS IT TO ALL
	   THE OTHER PROCESSES.
	   ------------------------------------------------------------*/
	if ( !iam ) {
	    /* Read the matrix stored on disk in Harwell-Boeing format. */
	    dreadhb(iam, fp, &m, &n, &nnz, &a, &asub, &xa);
	
	    printf("\tDimension\t%dx%d\t # nonzeros %d\n", m, n, nnz);
	    printf("\tProcess grid\t%d X %d\n", grid1.nprow, grid1.npcol);

	    /* Broadcast matrix A to the other PEs. */
	    MPI_Bcast( &m,   1,   mpi_int_t,  0, grid1.comm );
	    MPI_Bcast( &n,   1,   mpi_int_t,  0, grid1.comm );
	    MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid1.comm );
	    MPI_Bcast( a,    nnz, MPI_DOUBLE, 0, grid1.comm );
	    MPI_Bcast( asub, nnz, mpi_int_t,  0, grid1.comm );
	    MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid1.comm );
	} else {
	    /* Receive matrix A from PE 0. */
	    MPI_Bcast( &m,   1,   mpi_int_t,  0, grid1.comm );
	    MPI_Bcast( &n,   1,   mpi_int_t,  0, grid1.comm );
	    MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid1.comm );

	    /* Allocate storage for compressed column representation. */
	    dallocateA(n, nnz, &a, &asub, &xa);
	    
	    MPI_Bcast( a,    nnz, MPI_DOUBLE, 0, grid1.comm );
	    MPI_Bcast( asub, nnz, mpi_int_t,  0, grid1.comm );
	    MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid1.comm );
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

	if ( !(berr = doubleMalloc(nrhs)) )
	    ABORT("Malloc fails for berr[].");

	/* ------------------------------------------------------------
	   NOW WE SOLVE THE LINEAR SYSTEM.
	   ------------------------------------------------------------*/
	
	/* Set the default input options. */
	set_default_options(&options);

	/* Initialize ScalePermstruct and LUstruct. */
	ScalePermstructInit(m, n, &ScalePermstruct);
	LUstructInit(m, n, &LUstruct);

	/* Initialize the statistics variables. */
	PStatInit(&stat);
	
	/* Call the linear equation solver: factorize and solve. */
	pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid1,
			 &LUstruct, berr, &stat, &info);

	/* Check the accuracy of the solution. */
	if ( !iam ) {
	    dinf_norm_error(n, nrhs, b, ldb, xtrue, ldx);
	}
    
    
	/* Print the statistics. */
	PStatPrint(&stat, &grid1);

	/* ------------------------------------------------------------
	   DEALLOCATE STORAGE.
	   ------------------------------------------------------------*/
	PStatFree(&stat);
	Destroy_CompCol_Matrix(&A); 
	Destroy_LU(n, &grid1, &LUstruct);
	ScalePermstructFree(&ScalePermstruct);
	LUstructFree(&LUstruct);
	SUPERLU_FREE(b);
	SUPERLU_FREE(xtrue);
	SUPERLU_FREE(berr);

    } else { /* I am in grid 2. */
	iam = grid2.iam;  /* Get the logical number in the new grid. */

	/* ------------------------------------------------------------
	   PROCESS 0 READS THE MATRIX A, AND THEN BROADCASTS IT TO ALL
	   THE OTHER PROCESSES.
	   ------------------------------------------------------------*/
	if ( !iam ) {
	    /* Read the matrix stored on disk in Harwell-Boeing format. */
	    dreadhb(iam, fp, &m, &n, &nnz, &a, &asub, &xa);
	
	    printf("\tDimension\t%dx%d\t # nonzeros %d\n", m, n, nnz);
	    printf("\tProcess grid\t%d X %d\n", grid2.nprow, grid2.npcol);

	    /* Broadcast matrix A to the other PEs. */
	    MPI_Bcast( &m,   1,   mpi_int_t,  0, grid2.comm );
	    MPI_Bcast( &n,   1,   mpi_int_t,  0, grid2.comm );
	    MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid2.comm );
	    MPI_Bcast( a,    nnz, MPI_DOUBLE, 0, grid2.comm );
	    MPI_Bcast( asub, nnz, mpi_int_t,  0, grid2.comm );
	    MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid2.comm );
	} else {
	    /* Receive matrix A from PE 0. */
	    MPI_Bcast( &m,   1,   mpi_int_t,  0, grid2.comm );
	    MPI_Bcast( &n,   1,   mpi_int_t,  0, grid2.comm );
	    MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid2.comm );

	    /* Allocate storage for compressed column representation. */
	    dallocateA(n, nnz, &a, &asub, &xa);
	    
	    MPI_Bcast( a,    nnz, MPI_DOUBLE, 0, grid2.comm );
	    MPI_Bcast( asub, nnz, mpi_int_t,  0, grid2.comm );
	    MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid2.comm );
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

	if ( !(berr = doubleMalloc(nrhs)) )
	    ABORT("Malloc fails for berr[].");

	/* ------------------------------------------------------------
	   NOW WE SOLVE THE LINEAR SYSTEM.
	   ------------------------------------------------------------*/
	
	/* Set the default input options. */
	set_default_options(&options);
	
	/* Initialize ScalePermstruct and LUstruct. */
	ScalePermstructInit(m, n, &ScalePermstruct);
	LUstructInit(m, n, &LUstruct);

	/* Initialize the statistics variables. */
	PStatInit(&stat);
	
	/* Call the linear equation solver: factorize and solve. */
	pdgssvx_ABglobal(&options, &A, &ScalePermstruct, b, ldb, nrhs, &grid2,
			 &LUstruct, berr, &stat, &info);

	/* Check the accuracy of the solution. */
	if ( !iam ) {
	    dinf_norm_error(n, nrhs, b, ldb, xtrue, ldx);
	}
    
    
	/* Print the statistics. */
	PStatPrint(&stat, &grid2);

	/* ------------------------------------------------------------
	   DEALLOCATE STORAGE.
	   ------------------------------------------------------------*/
	PStatFree(&stat);
	Destroy_CompCol_Matrix(&A); 
	Destroy_LU(n, &grid2, &LUstruct);
	ScalePermstructFree(&ScalePermstruct);
	LUstructFree(&LUstruct);
	SUPERLU_FREE(b);
	SUPERLU_FREE(xtrue);
	SUPERLU_FREE(berr);
    }

    /* ------------------------------------------------------------
       RELEASE THE SUPERLU PROCESS GRIDS.
       ------------------------------------------------------------*/
    superlu_gridexit(&grid1);
    superlu_gridexit(&grid2);

out:
    /* ------------------------------------------------------------
       TERMINATES THE MPI EXECUTION ENVIRONMENT.
       ------------------------------------------------------------*/
    MPI_Finalize();

#if ( DEBUGlevel>=1 )
    CHECK_MALLOC(iam, "Exit main()");
#endif

}
