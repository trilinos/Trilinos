/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdio.h>
#include <mpi.h>
#include "comm.h"
#include "zoltan_mem.h"

/* Knowing who I send to, determine how many messages I'll receive, */
/* and their lengths.  Upon entry, the arrays "lengths_to" and "procs_to" */
/* contain list of processors I send to and the lengths of the corresponding */
/* messages. Upon exit, "lengths_from" and "procs_from" contain receive info. */

/* Note: By reinterpreting "to" and "from", this routine can do the opposite: */
/* given "from" data it computes "to" data. */

/* Data is allowed to be mapped from me to me. */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int       Zoltan_Comm_Invert_Map(
int      *lengths_to,		/* number of items I'm sending */
int      *procs_to,		/* procs I send to */
int       nsends,		/* number of messages I'll send */
int       self_msg,		/* do I copy data to myself? */
int     **plengths_from,	/* number of items I'm receiving */
int     **pprocs_from,		/* procs I recv lengths from */
int      *pnrecvs,		/* number of messages I receive */
int       my_proc,		/* my processor number */
int       nprocs,		/* total number of processors */
int       out_of_mem,		/* tell everyone I'm out of memory? */
int       tag,			/* message tag I can use */
MPI_Comm  comm)			/* communicator */
{
    int      *lengths_from;	/* lengths of my recvs */
    int      *procs_from;	/* procs I recv lengths from */
    int      *msg_count;	/* binary flag for procs I send to (nprocs) */
    int      *counts;		/* argument to Reduce_scatter */
    int       nrecvs=0;		/* number of messages I'll receive */
    int       i;		/* loop counter */
    MPI_Status status;		/* return MPI argument */

    msg_count = (int *) ZOLTAN_MALLOC(nprocs * sizeof(int));
    counts = (int *) ZOLTAN_MALLOC(nprocs * sizeof(int));

    if (msg_count == NULL || counts == NULL) out_of_mem = 1;

    MPI_Allreduce((void *) &out_of_mem, (void *) &i, 1, MPI_INT, MPI_MAX, comm);
    if (i) {
	ZOLTAN_FREE((void **) &counts);
	ZOLTAN_FREE((void **) &msg_count);
	return(ZOLTAN_MEMERR);
    }

    for (i = 0; i < nprocs; i++) {
	msg_count[i] = 0;
	counts[i] = 1;
    }
    for (i = 0; i < nsends + self_msg; i++)
	msg_count[procs_to[i]] = 1;

    MPI_Reduce_scatter((void *) msg_count, (void *) &nrecvs, counts, MPI_INT,
	MPI_SUM, comm);

    ZOLTAN_FREE((void **) &counts);
    ZOLTAN_FREE((void **) &msg_count);

    lengths_from = (int *) ZOLTAN_MALLOC(nrecvs*sizeof(int));
    procs_from = (int *) ZOLTAN_MALLOC(nrecvs*sizeof(int));

    /* Note: these mallocs should never fail as prior frees are larger. */


    /* Send the lengths of all my real messages to their receivers. */
    for (i = 0; i < nsends + self_msg; i++) {
        if (procs_to[i] != my_proc) {
	    MPI_Send((void *) &lengths_to[i], 1, MPI_INT, procs_to[i], tag, comm);
        }
	else {
	    /* Always put self stuff at end. */
	    lengths_from[nrecvs - 1] = lengths_to[i];
	    procs_from[nrecvs - 1] = my_proc;
	}
    }

    /* Now receive the lengths of all my messages from their senders. */
    /* Note that proc/length_from lists are ordered by sequence msgs arrive. */
    for (i = 0; i < nrecvs - self_msg; i++) {
	MPI_Recv((void *) &lengths_from[i], 1, MPI_INT, MPI_ANY_SOURCE, tag,
	comm, &status);
	procs_from[i] = status.MPI_SOURCE;
    }

    /* Note: I'm counting on having a unique tag or some of my incoming */
    /*       messages might get confused with others. */

    
    /* Sort recv lists to keep execution deterministic (e.g. for debugging) */

    Zoltan_Comm_Sort_Ints(procs_from, lengths_from, nrecvs);

    *plengths_from = lengths_from;
    *pprocs_from = procs_from;
    *pnrecvs = nrecvs - self_msg;    /* Only return number of true messages */

    return(ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
