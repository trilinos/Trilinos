// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stdio.h>
#include <mpi.h>
#include "comm.h"
#include "zoltan_mem.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int       Zoltan_Comm_Exchange_Sizes(
int      *sizes_to,		/* value I need to exchange (size of true msg) */
int      *procs_to,		/* procs I send to */
int       nsends,		/* number of messages I'll send */
int       self_msg,		/* do I copy data to myself? */
int      *sizes_from,		/* (returned) size of all my receives */
int      *procs_from,		/* procs I recv from */
int       nrecvs,		/* number of messages I receive */
int      *total_recv_size,	/* (returned) sum of all incoming sizes */
int       my_proc,		/* my processor number */
int       tag,			/* message tag I can use */
MPI_Comm  comm) {		/* communicator */

    int       self_index_to;	/* location of self in procs_to */
    MPI_Status status;		/* status of commuication operation */
    int       i;		/* loop counter */


    /* If sizes vary, then I need to communicate message lengths */

    self_index_to = -1;
    for (i = 0; i < nsends + self_msg; i++) {
	if (procs_to[i] != my_proc) 
	    MPI_Send((void *) &sizes_to[i], 1, MPI_INT, procs_to[i], tag, comm);
	else
	    self_index_to = i;
    }

    *total_recv_size = 0;
    for (i = 0; i < nrecvs + self_msg; i++) {
	if (procs_from[i] != my_proc) 
	    MPI_Recv((void *) &sizes_from[i], 1, MPI_INT, procs_from[i],
		     tag, comm, &status);
	else 
	    sizes_from[i] = sizes_to[self_index_to];

	*total_recv_size += sizes_from[i];
    }

    return (ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
