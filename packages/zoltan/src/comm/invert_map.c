/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#include <stdio.h>
#include "mpi.h"
#include "lbi_const.h"
#include "all_allo_const.h"

/* Knowing who I send to, determine how many messages I'll receive, */
/* and their lengths.  Upon entry, the arrays "lengths_to" and "procs_to" */
/* contain list of processors I send to and the lengths of the corresponding */
/* messages. Upon exit, "lengths_from" and "procs_from" contain receive info. */

/* Data is allowed to be mapped from me to me. This self entry is always last */
/* in the list of recvs. */

/* List of messages to/from is ended by a negative value in procs array. */

/* Always keep self messages at end. */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Prototype definition */
static void      LB_Sort_Procs(int *, int *, int);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int       LB_Invert_Map(
int      *lengths_to,		/* lengths of my sends */
int      *procs_to,		/* procs I send to */
int       nsends,		/* number of messages I'll send */
int       self_msg,		/* do I copy data to myself? */
int     **plengths_from,	/* lengths of my recvs */
int     **pprocs_from,		/* procs I recv lengths from */
int      *pnrecvs,		/* number of messages I receive */
int       my_proc,		/* my processor number */
int       nprocs,		/* total number of processors */
int       tag,			/* message tag I can use */
int       deterministic,        /* flag indicating whether result should be
                                   deterministic (i.e., independent of the
                                   order in which messages are received). */
MPI_Comm  comm)			/* communicator */
{
    int      *lengths_from;	/* lengths of my recvs */
    int      *procs_from;	/* procs I recv lengths from (-1 ends) */
    int      *msg_count;	/* binary flag for procs I send to (nprocs) */
    int      *counts;		/* argument to Reduce_scatter */
    int       nrecvs;		/* number of messages I'll receive */
    int       i;		/* loop counter */
    MPI_Status status;		/* return MPI argument */

    msg_count = (int *) LB_MALLOC(nprocs * sizeof(int));
    counts = (int *) LB_MALLOC(nprocs * sizeof(int));

    for (i = 0; i < nprocs; i++) {
	msg_count[i] = 0;
	counts[i] = 1;
    }
    for (i = 0; i < nsends + self_msg; i++)
	msg_count[procs_to[i]] = 1;

    MPI_Reduce_scatter(msg_count, &nrecvs, counts, MPI_INT, MPI_SUM, comm);

    LB_FREE((void **) &counts);
    LB_FREE((void **) &msg_count);

    lengths_from = (int *) LB_MALLOC(nrecvs*sizeof(int));
    procs_from = (int *) LB_MALLOC(nrecvs*sizeof(int));

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

    /* Barrier to insure all my MPI_ANY_SOURCE messages are received.
       Otherwise some procs could proceed to comm_do and start sending to me.
    
    MPI_Barrier(comm);
    */
    
    /* Instead of barrier, I'm counting on having a unique tag. */

    
    if (deterministic) {
        /* 
         * If a deterministic ordering is needed (eg. for debugging), 
         * sort recvs. 
         * Note: self messages are kept at the end. 
         */
        LB_Sort_Procs(procs_from, lengths_from, nrecvs - self_msg);
    }

    *plengths_from = lengths_from;
    *pprocs_from = procs_from;
    *pnrecvs = nrecvs - self_msg;    /* Only return number of true messages */

    return(LB_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void      LB_Sort_Procs(
int      *procs,                /* procs I'll receive from */
int      *vals,                 /* value associated with each message */
int       nvals)		/* length of these two arrays */
{
/* This routine will ensure that the ordering */
/* produced by the invert_map routines is deterministic.  This should */
/* make bugs more reproducible.  This is accomplished by sorting */
/* the message lists by processor ID. */

    int       temp;		/* swapping value */
    int       i, j;             /* loop counter */

    /* Use a simple bubble sort on procs entries. */

    for (i = 1; i < nvals; i++) {
	j = i;
	while (j > 0 && procs[j] < procs[j - 1]) {
	    /* Flip j and j-1 entries. */
	    temp = procs[j];
	    procs[j] = procs[j - 1];
	    procs[j - 1] = temp;
	    temp = vals[j];
	    vals[j] = vals[j - 1];
	    vals[j - 1] = temp;
	    j--;
	}
    }
}
