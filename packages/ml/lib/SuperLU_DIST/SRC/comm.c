#include "superlu_ddefs.h"


void
bcast_tree(void *buf, int count, MPI_Datatype dtype, int root, int tag,
	   gridinfo_t *grid, int scope, int *recvcnt)
/* 
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 * Purpose
 * =======
 *   Broadcast an array of *dtype* numbers. The communication pattern
 *   is a tree with number of branches equal to NBRANCHES.
 *   The process ranks are between 0 and Np-1.
 * 
 *   The following two pairs of graphs give different ways of viewing the same
 *   algorithm.  The first pair shows the trees as they should be visualized
 *   when examining the algorithm.  The second pair are isomorphic graphs of
 *   of the first, which show the actual pattern of data movement.
 *   Note that a tree broadcast with NBRANCHES = 2 is isomorphic with a
 *   hypercube broadcast (however, it does not require the nodes be a
 *   power of two to work).
 *
 *    TREE BROADCAST, NBRANCHES = 2     *    TREE BROADCAST, NBRANCHES = 3
 *       
 *     root=2
 * i=4   &______________                *
 *       |              \               *       root=2
 * i=2   &______         &______        * i=3     &______________________
 *       |      \        |      \       *         |          \           \
 * i=1   &__     &__     &__     &__    * i=1     &______     &______     &__
 *       |  \    |  \    |  \    |  \   *         |  \   \    |  \   \    |  \
 *       2   3   4   5   6   7   0   1  *         2   3   4   5   6   7   0   1
 *
 *
 *          ISOMORPHIC GRAPHS OF ABOVE, SHOWN IN MORE FAMILIAR TERMS:
 *
 *                2                                           2
 *       _________|_________                       ___________|____________
 *      /         |         \                     /           |      |     \
 *     6          4          3                   5            0      3      4
 *    / \         |                             / \           |
 *   0   7        5                            6   7          1
 *   |
 *   1
 *
 *
 * Arguments
 * =========
 * 
 * scope
 */
{
    int Iam, i, j, Np, nbranches = 2;
    int destdist; /* The distance of the destination node. */
    int mydist;   /* My distance from root. */
    superlu_scope_t *scp;
    MPI_Status status;

    if ( scope == COMM_COLUMN ) scp = &grid->cscp;
    else if ( scope == ROW ) scp = &grid->rscp;
    Np = scp->Np;
    if ( Np < 2 ) return;
    Iam = scp->Iam;
    
    if ( Iam == root ) {
	for (i = nbranches; i < Np; i *= nbranches);
	for (i /= nbranches; i > 0; i /= nbranches) {
	    for (j = 1; j < nbranches; ++j) {
		destdist = i*j;
		if ( destdist < Np )
		    MPI_Send( buf, count, dtype, (Iam+destdist)%Np, 
			     tag, scp->comm );
	    }
	}
    } else {
	mydist = (Np + Iam - root) % Np;
	for (i = nbranches; i < Np; i *= nbranches);
	for (i /= nbranches; (mydist%i); i /= nbranches);
/*	MPI_Probe( MPI_ANY_SOURCE, tag, scp->comm, &status );*/
	MPI_Recv( buf, count, dtype, MPI_ANY_SOURCE, tag, scp->comm, &status );
	MPI_Get_count( &status, dtype, recvcnt );

	/* I need to send data to others. */
	while ( (i > 1) && !(mydist%i) ) {
	    i /= nbranches;
	    for (j = 1; j < nbranches; ++j) {
		destdist = mydist + j*i;
		if ( destdist < Np )
		    MPI_Send( buf, *recvcnt, dtype, (root+destdist)%Np, 
			     tag, scp->comm );
	    }
	}
    }
} /* BCAST_TREE */







