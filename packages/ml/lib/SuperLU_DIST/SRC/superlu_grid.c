/*
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 */

#include "superlu_ddefs.h"

/* Define global variables */
MPI_Datatype SuperLU_MPI_DOUBLE_COMPLEX;

/*
 * All processes in the MPI communicator must call this routine.
 */
void superlu_gridinit(MPI_Comm Bcomm, /* The base communicator upon which
					 the new grid is formed. */
		      int_t nprow, int_t npcol, gridinfo_t *grid)
{
    int Np = nprow * npcol;
    int_t *usermap;
    int i, j, info;

    /* Make a list of the processes in the new communicator. */
    usermap = (int_t *) SUPERLU_MALLOC(Np*sizeof(int_t));
    for (j = 0; j < npcol; ++j)
	for (i = 0; i < nprow; ++i) usermap[j*nprow+i] = i*npcol+j;
    
    /* Check MPI environment initialization. */
    MPI_Initialized( &info );
    if ( !info )
	ABORT("C main program must explicitly call MPI_Init()");

    MPI_Comm_size( Bcomm, &info );
    if ( info < Np )
	ABORT("Number of processes is smaller than NPROW * NPCOL");

    superlu_gridmap(Bcomm, nprow, npcol, usermap, nprow, grid);
    
    SUPERLU_FREE(usermap);
}


/*
 * All processes in the MPI communicator must call this routine.
 */
void superlu_gridmap(
		     MPI_Comm Bcomm, /* The base communicator upon which
					the new grid is formed. */
		     int_t nprow,
		     int_t npcol,
		     int_t usermap[], /* usermap(i,j) holds the process
					 number to be placed in {i,j} of
					 the process grid.  */
		     int_t ldumap,    /* The leading dimension of the
					 2D array usermap[].  */
		     gridinfo_t *grid)
{
    MPI_Group grp, mpi_base_group, superlu_grp;
    int Np = nprow * npcol, mycol, myrow;
    int *pranks;
    int i, j, info;
    static int first = 1;

    /* Create datatype in C for MPI complex. */
    if ( first ) {
	MPI_Type_contiguous( 2, MPI_DOUBLE, &SuperLU_MPI_DOUBLE_COMPLEX );
	MPI_Type_commit( &SuperLU_MPI_DOUBLE_COMPLEX );
	first = 0;
    }

    /* Check MPI environment initialization. */
    MPI_Initialized( &info );
    if ( !info )
	ABORT("C main program must explicitly call MPI_Init()");

    grid->nprow = nprow;
    grid->npcol = npcol;

    /* Make a list of the processes in the new communicator. */
    pranks = (int *) SUPERLU_MALLOC(Np*sizeof(int));
    for (j = 0; j < npcol; ++j)
	for (i = 0; i < nprow; ++i)
	    pranks[i*npcol+j] = usermap[j*ldumap+i];
    
    /*
     * Form MPI communicator for all.
     */
    /* Get the group underlying Bcomm. */
    MPI_Comm_group( Bcomm, &mpi_base_group );
    /* Create the new group. */
    MPI_Group_incl( mpi_base_group, Np, pranks, &superlu_grp );
    /* Create the new communicator. */
    /* NOTE: The call is to be executed by all processes in Bcomm,
       even if they do not belong in the new group -- superlu_grp. */
    MPI_Comm_create( Bcomm, superlu_grp, &grid->comm );

    /* Bail out if I am not in the group, superlu_group. */
    if ( grid->comm == MPI_COMM_NULL ) {
	grid->comm = Bcomm;
	MPI_Comm_rank( Bcomm, &i );
	grid->iam = i;
	/*grid->iam = -1;*/
	SUPERLU_FREE(pranks);
	return;
    }

    MPI_Comm_rank( grid->comm, &(grid->iam) );
    myrow = grid->iam / npcol;
    mycol = grid->iam % npcol;

    /*
     * Form MPI communicator for myrow, scope = COMM_ROW.
     */
#if 0
    for (i = 0; i < npcol; ++i) pranks[i] = myrow*npcol + i;
    MPI_Comm_group( grid->comm, &superlu_grp );          /* Find all's group */
    MPI_Group_incl( superlu_grp, npcol, pranks, &grp );  /* Form new group */
    MPI_Comm_create( grid->comm, grp, &grid->rscp.comm );/* Create new comm */
#else
    MPI_Comm_split(grid->comm, myrow, mycol, &(grid->rscp.comm));
#endif

    /*
     * Form MPI communicator for mycol, scope = COMM_COLUMN.
     */
#if 0
    for (i = 0; i < nprow; ++i) pranks[i] = i*npcol + mycol;
    MPI_Group_incl( superlu_grp, nprow, pranks, &grp );  /* Form new group */
    MPI_Comm_create( grid->comm, grp, &grid->cscp.comm );/* Create new comm */
#else
    MPI_Comm_split(grid->comm, mycol, myrow, &(grid->cscp.comm));
#endif

    grid->rscp.Np = npcol;
    grid->rscp.Iam = mycol;
    grid->cscp.Np = nprow;
    grid->cscp.Iam = myrow;

#if 0
    {
	int tag_ub;
	if ( !grid->iam ) {
	    MPI_Attr_get(Bcomm, MPI_TAG_UB, &tag_ub, &info);
	    printf("MPI_TAG_UB %d\n", tag_ub);
	    /* returns 4295677672
	       In reality it is restricted to no greater than 16384. */
	}
	exit(0);
    }
#endif

    SUPERLU_FREE(pranks);
}


void superlu_gridexit(gridinfo_t *grid)
{
    if ( grid->comm != MPI_COMM_NULL && grid->comm != MPI_COMM_WORLD ) {
	/* Marks the communicator objects for deallocation. */
	MPI_Comm_free( &grid->rscp.comm );
	MPI_Comm_free( &grid->cscp.comm );
	MPI_Comm_free( &grid->comm );
    }
    /*MPI_Type_free( &SuperLU_MPI_DOUBLE_COMPLEX );*/
}
