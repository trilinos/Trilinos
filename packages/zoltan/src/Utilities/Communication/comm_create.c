// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stdio.h>
#include "comm.h"
#include "zoltan_mem.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* From the mapping of where to send things, construct the communication */
/* data structures. */
/* Note: This routine sets up data structure assuming all objects are the */
/* same size.  If this isn't the case, a subsequent call to Comm_Resize  */
/* is needed. */

int Zoltan_Comm_Create(
ZOLTAN_COMM_OBJ **cobj,		/* returned communicator object */
int       nvals,		/* number of values I currently own */
int      *assign,		/* processor assignment for all my values */
MPI_Comm  comm,			/* communicator for xfer operation */
int       tag,			/* message tag I can use */
int      *pnvals_recv)		/* returned # vals I own after communication */
{
    ZOLTAN_COMM_OBJ *plan;	/* returned communication data structure */
    int      *starts=NULL;	/* pointers into list of vals for procs */
    int      *lengths_to=NULL;	/* lengths of messages I'll send */
    int      *procs_to=NULL;	/* processors I'll send to */
    int      *indices_to=NULL;	/* local_id values I'll be sending */
    int      *starts_to=NULL;	/* where in list my sends begin */
    int      *lengths_from=NULL;	/* lengths of messages I'll receive */
    int      *procs_from=NULL;	/* processors I'll receive from */
    int      *starts_from=NULL;	/* pointers for where to put recv data */
    int       my_proc;		/* my processor tag in communicator */
    int       nprocs;		/* number of  processors in communicator */
    int       max_send_size =0;	/* size of longest message I send */
    int       total_recv_size;  /* total size of messages I recv */
    int       no_send_buff;	/* is data nicely grouped by processor? */
    int       nactive = 0;	/* number of values to remap */
    int       self_msg = 0;	/* do I have data for myself? */
    int       nsends=0;		/* # procs I'll send to (including self) */
    int       nrecvs=0;		/* # procs I'll recv from (including self) */
    int       proc;		/* processor I communicate with */
    int       prev_proc;	/* processor on previous loop pass */
    int       index;		/* index into list of objects */
    int       out_of_mem;	/* am I out of memory? */
    int       comm_flag;	/* status flag */
    int       i, j;		/* loop counters */
    static char *yo = "Zoltan_Comm_Create";

    if (comm == MPI_COMM_NULL){
      ZOLTAN_COMM_ERROR("Invalid communicator: MPI_COMM_NULL.", yo, -1);
      return ZOLTAN_FATAL;
    }

    comm_flag = ZOLTAN_OK;
    MPI_Comm_rank(comm, &my_proc);
    MPI_Comm_size(comm, &nprocs);

    /* First check to see if items are grouped by processor with no gaps. */
    /* If so, indices_to should be NULL (= identity) */

    /* Make data structures that will allow me to traverse arrays quickly. */
    /* Begin by determining number of objects I'm sending to each processor. */
    starts = (int *) ZOLTAN_MALLOC((nprocs + 1) * sizeof(int));

    out_of_mem = FALSE;
    if (starts == NULL) { 
	out_of_mem = TRUE;
	goto Mem_Err;
    }

    /* First, use starts to count values going to each proc. */
    for (i = 0; i < nprocs; i++) {
	starts[i] = 0;
    }

    /* Note: Negative assign value means ignore item. */
    /* Non-trailing negatives mean data not packed so need send_buf. */
    /* Could (but don't) allow negatives between processor blocks w/o buf. */
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
	lengths_to = (int *) ZOLTAN_MALLOC(nsends * sizeof(int));
	starts_to = (int *) ZOLTAN_MALLOC(nsends * sizeof(int));
	procs_to = (int *) ZOLTAN_MALLOC(nsends * sizeof(int));
        if (nsends != 0 && (lengths_to == NULL || starts_to == NULL ||
			    procs_to == NULL)) {
	    out_of_mem = TRUE;
	    goto Mem_Err;
	}
	index = 0;
        /* Note that procs_to is in the order the data was passed in. */
	for (i = 0; i < nsends; i++) {
	    starts_to[i] = index;
	    proc = assign[index];
	    procs_to[i] = proc;
	    index += starts[proc];
	}

	/* Now sort the outgoing procs. */
	/* This keeps recvs deterministic if I ever invert communication */
	/* It also allows for better balance of traffic in comm_do */
	Zoltan_Comm_Sort_Ints(procs_to, starts_to, nsends);

	max_send_size = 0;
	for (i = 0; i < nsends; i++) {
	    proc = procs_to[i];
	    lengths_to[i] = starts[proc];
	    if (proc != my_proc && lengths_to[i] > max_send_size) {
		max_send_size = lengths_to[i];
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

	indices_to = (int *) ZOLTAN_MALLOC(nactive * sizeof(int));

        if (nactive != 0 && indices_to == NULL) {
	    out_of_mem = TRUE;
	    goto Mem_Err;
	}

	for (i = 0; i < nvals; i++) {
	    proc = assign[i];
	    if (proc >= 0) {
	        indices_to[starts[proc]] = i;
	        ++starts[proc];
	    }
	}

	/* Indices_to array now has the data in clumps for each processor. */
	/* Now reconstruct starts array to index into indices_to. */
	for (i = nprocs - 1; i; i--) {
	    starts[i] = starts[i - 1];
	}
	starts[0] = 0;
	starts[nprocs] = nactive;

	/* Construct lengths_to, starts_to and procs_to arrays. */
	/* Note: If indices_to is needed, procs are in increasing order */
	lengths_to = (int *) ZOLTAN_MALLOC(nsends * sizeof(int));
	starts_to = (int *) ZOLTAN_MALLOC(nsends * sizeof(int));
	procs_to = (int *) ZOLTAN_MALLOC(nsends * sizeof(int));
        if (nsends != 0 && (lengths_to == NULL || starts_to == NULL ||
			    procs_to == NULL)) {
	    out_of_mem = TRUE;
	    goto Mem_Err;
	}

	j = 0;
	max_send_size = 0;
	for (i = 0; i < nprocs; i++) {
	    if (starts[i + 1] != starts[i]) {
		starts_to[j] = starts[i];
		lengths_to[j] = starts[i + 1] - starts[i];
		if (i != my_proc && lengths_to[j] > max_send_size) {
		    max_send_size = lengths_to[j];
		}
		procs_to[j] = i;
		j++;
	    }
	}
    }

    /* Now change nsends to count only non-self messages */
    nsends -= self_msg;

Mem_Err:

    ZOLTAN_FREE(&starts);

    /* Determine how many messages & what length I'll receive. */
    comm_flag = Zoltan_Comm_Invert_Map(lengths_to, procs_to, nsends, self_msg,
	       &lengths_from, &procs_from, &nrecvs, my_proc, nprocs,
	       out_of_mem,tag, comm);

    starts_from = (int *) ZOLTAN_MALLOC((nrecvs + self_msg) * sizeof(int));
    if (starts_from == NULL && nrecvs + self_msg != 0) {
	comm_flag = ZOLTAN_MEMERR;
    }
    else {
        j = 0;
        for (i = 0; i < nrecvs + self_msg; i++) {
	    starts_from[i] = j;
	    j += lengths_from[i];
        }
    }
    
    if (comm_flag != ZOLTAN_OK) {
        if (comm_flag == ZOLTAN_MEMERR) {
	    ZOLTAN_COMM_ERROR("Out of memory", yo, my_proc);
	}
	ZOLTAN_FREE(&starts_from);
	ZOLTAN_FREE(&indices_to);
	ZOLTAN_FREE(&procs_to);
	ZOLTAN_FREE(&starts_to);
	ZOLTAN_FREE(&lengths_to);
	ZOLTAN_FREE(&lengths_from);
	ZOLTAN_FREE(&procs_from);
	return(comm_flag);
    }

    total_recv_size = 0;
    for (i = 0; i < nrecvs + self_msg; i++)
	total_recv_size += lengths_from[i];

    plan = (ZOLTAN_COMM_OBJ *) ZOLTAN_MALLOC(sizeof(ZOLTAN_COMM_OBJ));
    plan->lengths_to = lengths_to;
    plan->starts_to = starts_to;
    plan->procs_to = procs_to;
    plan->indices_to = indices_to;
    plan->lengths_from = lengths_from;
    plan->starts_from = starts_from;
    plan->procs_from = procs_from;
    plan->indices_from = NULL;


    plan->sizes = NULL;
    plan->sizes_to = NULL;
    plan->sizes_from = NULL;
    plan->starts_to_ptr = NULL;
    plan->starts_from_ptr = NULL;
    plan->indices_to_ptr = NULL;
    plan->indices_from_ptr = NULL;


    plan->nvals = nvals;
    plan->nvals_recv = total_recv_size;
    plan->nrecvs = nrecvs;
    plan->nsends = nsends;
    plan->nindices_to = nactive;
    plan->nindices_from = 0;
    plan->self_msg = self_msg;
    plan->max_send_size = max_send_size;
    plan->total_recv_size = total_recv_size;
    plan->maxed_recvs = 0;
    plan->comm = comm;
    plan->plan_reverse = NULL;

    if (MPI_RECV_LIMIT > 0){
      /* If we have a limit to the number of posted receives we are allowed,
      ** and our plan has exceeded that, then switch to an MPI_Alltoallv so
      ** that we will have fewer receives posted when we do the communication.
      */
      MPI_Allreduce(&nrecvs, &i, 1, MPI_INT, MPI_MAX, comm);
      if (i > MPI_RECV_LIMIT){
        plan->maxed_recvs = 1;
      }
    }

    if (plan->maxed_recvs){
      plan->request = NULL;
      plan->status = NULL;
    }
    else{
      plan->request = (MPI_Request *) ZOLTAN_MALLOC(plan->nrecvs * sizeof(MPI_Request));
      plan->status = (MPI_Status *) ZOLTAN_MALLOC(plan->nrecvs * sizeof(MPI_Status));
      if (plan->nrecvs && ((plan->request == NULL) || (plan->status == NULL))) 
        comm_flag = ZOLTAN_MEMERR;
    }

    *pnvals_recv = total_recv_size;
    *cobj = plan;

    return (comm_flag);
}
#define COPY_BUFFER(buf, type, num) \
  if (((num)>0) && from->buf) { \
    to->buf = (type *)ZOLTAN_MALLOC((num) * sizeof(type)); \
    if (!to->buf) { \
      ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory."); \
      Zoltan_Comm_Destroy(toptr); \
      return ZOLTAN_MEMERR; \
    } \
    memcpy(to->buf, from->buf, (num) * sizeof(type)); \
  } \
  else { \
    to->buf = NULL; \
  }
    
ZOLTAN_COMM_OBJ *Zoltan_Comm_Copy(ZOLTAN_COMM_OBJ *from)
{
  ZOLTAN_COMM_OBJ *to = NULL;

  Zoltan_Comm_Copy_To(&to, from);

  return to;
}
int Zoltan_Comm_Copy_To(ZOLTAN_COMM_OBJ **toptr, ZOLTAN_COMM_OBJ *from)
{
  static char *yo = "Zoltan_Comm_Copy_To";
  int proc = 0;
  ZOLTAN_COMM_OBJ *to= NULL;

  if (!toptr){
    return ZOLTAN_FATAL;
  }

  if (*toptr){
    Zoltan_Comm_Destroy(toptr);
  }

  if (from){

    MPI_Comm_rank(from->comm, &proc);

    to = *toptr = (ZOLTAN_COMM_OBJ *)ZOLTAN_MALLOC(sizeof(ZOLTAN_COMM_OBJ));

    *to = *from;

    MPI_Comm_dup(from->comm, &(to->comm));

    COPY_BUFFER(procs_to, int, to->nsends);
    COPY_BUFFER(procs_from, int, to->nrecvs);
    COPY_BUFFER(lengths_to, int, to->nsends);
    COPY_BUFFER(lengths_from, int, to->nrecvs);
    COPY_BUFFER(starts_to, int, to->nsends);
    COPY_BUFFER(starts_from, int, to->nrecvs + to->self_msg);
    COPY_BUFFER(indices_to, int, to->nindices_to);
    COPY_BUFFER(indices_from, int, to->nindices_from);
    COPY_BUFFER(sizes, int, to->nvals + 1);
    COPY_BUFFER(sizes_to, int, to->nsends + to->self_msg);
    COPY_BUFFER(sizes_from, int, to->nrecvs + to->self_msg);
    COPY_BUFFER(starts_to_ptr, int, to->nsends + to->self_msg);
    COPY_BUFFER(starts_from_ptr, int, to->nrecvs + to->self_msg);
    COPY_BUFFER(indices_to_ptr, int, to->nvals);
    COPY_BUFFER(indices_from_ptr, int, to->nvals);

    COPY_BUFFER(request, MPI_Request, to->nrecvs);
    COPY_BUFFER(status, MPI_Status, to->nrecvs);
  }

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
