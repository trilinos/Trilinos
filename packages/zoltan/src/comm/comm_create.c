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
#include "comm_const.h"
#include "all_allo_const.h"


/* From the mapping of where to send things, construct the communication */
/* data structures. */

int LB_Comm_Create(
struct Comm_Obj **cobj,		/* returned communicator object */
int       nvals,		/* number of values I currently own */
int      *assign,		/* processor assignment for all my values */
MPI_Comm  comm,			/* communicator for xfer operation */
int       tag,			/* message tag I can use */
int      *pnrecv)		/* returned # vals I own after communication */
{
    struct Comm_Obj *plan;	/* returned communication data structure pointer */
    int      *starts;		/* pointers into list of vals for procs */
    int      *lengths_to;	/* lengths of messages I'll send */
    int      *procs_to;		/* processors I'll send to */
    int      *indices_to;	/* local_id values I'll be sending */
    int      *lengths_from;	/* lengths of messages I'll receive */
    int      *procs_from;	/* processors I'll receive from */
    int       my_proc;		/* my processor tag in communicator */
    int       nprocs;		/* number of  processors in communicator */
    int       max_send_length;	/* size of longest message I send */
    int       total_recv_length;/* total size of messages I recv */
    int       no_send_buff;	/* is data nicely grouped by processor? */
    int       nactive;		/* number of values to remap */
    int       self_msg;		/* do I have data for myself? */
    int       nsends;		/* # procs I'll send to (including self) */
    int       nrecvs;		/* # procs I'll recv from (including self) */
    int       proc;		/* processor I communicate with */
    int       prev_proc;	/* processor on previous loop pass */
    int       index;		/* index into list of objects */
    int       lb_flag;		/* status flag */
    int       i, j;		/* loop counters */

    MPI_Comm_rank(comm, &my_proc);
    MPI_Comm_size(comm, &nprocs);

    /* First check to see if items are grouped by processor with no gaps. */
    /* If so, indices_to should be NULL (= identity) */

    /* Make data structures that will allow me to traverse arrays quickly. */
    /* Begin by determining number of objects I'm sending to each processor. */
    starts = (int *) LB_MALLOC((nprocs + 1) * sizeof(int));

    for (i = 0; i < nprocs; i++)
	starts[i] = 0;

    nactive = 0;
    no_send_buff = TRUE;
    prev_proc = nprocs;
    for (i = 0; i < nvals; i++) {
	proc = assign[i];
	if (no_send_buff && proc != prev_proc) { /* Checks if blocked by proc */
	    if (proc >= 0 && (starts[proc] || prev_proc < 0)) {
		no_send_buff = FALSE;
	    }
	    else {
	        prev_proc = proc;
	    }
	}
	if (proc >= 0) {
	    ++starts[proc];
	    ++nactive;
	}
    }

    self_msg = (starts[my_proc] != 0);

    if (no_send_buff) {
	/* Grouped by processor.  Array indices_to can be NULL (= identity) */
	nsends = 0;
	for (i = 0; i < nprocs; i++) {
	    if (starts[i] != 0) ++nsends;
	}
	indices_to = NULL;
	lengths_to = (int *) LB_MALLOC(nsends * sizeof(int));
	procs_to = (int *) LB_MALLOC(nsends * sizeof(int));
	index = 0;
	max_send_length = 0;
	for (i = 0; i < nsends; i++) {
	    proc = assign[index];
	    procs_to[i] = proc;
	    lengths_to[i] = starts[proc];
	    index += starts[proc];
	    if (proc != my_proc && lengths_to[i] > max_send_length) {
		max_send_length = lengths_to[i];
	    }
	}
    }

    else {	/* Not grouped by processor.  More complex data structures. */

	/* Sum starts values to be offsets into indices_to array. */
	nsends = (starts[0] != 0);
	for (i = 1; i < nprocs; i++) {
	    if (starts[i] != 0)
		++nsends;
	    starts[i] += starts[i - 1];
	}

	for (i = nprocs - 1; i; i--)
	    starts[i] = starts[i - 1];

	starts[0] = 0;

	indices_to = (int *) LB_MALLOC(nactive * sizeof(int));

	for (i = 0; i < nvals; i++) {
	    proc = assign[i];
	    if (proc >= 0) {
	        indices_to[starts[proc]] = i;
	        ++starts[proc];
	    }
	}

	/* Indices_to array now has all the surfaces in clumps for each processor. */
	/* Now reconstruct starts array to index into indices_to. */
	for (i = nprocs - 1; i; i--) {
	    starts[i] = starts[i - 1];
	}
	starts[0] = 0;
	starts[nprocs] = nactive;

	/* Construct lengths_to and procs_to arrays. */
	lengths_to = (int *) LB_MALLOC(nsends * sizeof(int));
	procs_to = (int *) LB_MALLOC(nsends * sizeof(int));
	j = 0;
	max_send_length = 0;
	for (i = 0; i < nprocs; i++) {
	    if (starts[i + 1] != starts[i]) {
		lengths_to[j] = starts[i + 1] - starts[i];
		if (i != my_proc && lengths_to[j] > max_send_length) {
		    max_send_length = lengths_to[j];
		}
		procs_to[j] = i;
		j++;
	    }
	}
    }

    /* Now change nsends to count only non-self messages */
    nsends -= self_msg;

    LB_FREE((void **) &starts);

    /* Determine how many messages & what length I'll receive. */
    lb_flag = LB_Invert_Map(lengths_to, procs_to, nsends, self_msg,
	       &lengths_from, &procs_from, &nrecvs, my_proc, nprocs, tag, comm);
    
    if (lb_flag != LB_OK && lb_flag != LB_WARN) {
	LB_FREE((void **) &indices_to);
	LB_FREE((void **) &procs_to);
	LB_FREE((void **) &lengths_to);
	LB_FREE((void **) &starts);
	return(lb_flag);
    }

    total_recv_length = 0;
    for (i = 0; i < nrecvs + self_msg; i++)
	total_recv_length += lengths_from[i];

    plan = (struct Comm_Obj *) LB_MALLOC(sizeof(struct Comm_Obj));
    plan->lengths_to = lengths_to;
    plan->procs_to = procs_to;
    plan->indices_to = indices_to;
    plan->lengths_from = lengths_from;
    plan->procs_from = procs_from;
    plan->indices_from = NULL;
    plan->nrecvs = nrecvs;
    plan->nsends = nsends;
    plan->self_msg = self_msg;
    plan->max_send_length = max_send_length;
    plan->total_recv_length = total_recv_length;
    plan->comm = comm;
    plan->request = (MPI_Request *) LB_MALLOC(plan->nrecvs * sizeof(MPI_Request));
    plan->status = (MPI_Status *) LB_MALLOC(plan->nrecvs * sizeof(MPI_Status));

    *pnrecv = total_recv_length;
    *cobj = plan;
    return (lb_flag);
}
