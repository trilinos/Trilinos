/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
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
#include <memory.h>
#include <mpi.h>
#include "comm.h"
#include "zoltan_mem.h"


/* Given the communication object, perform a communication operation as
   efficiently as possible. */

/* Several cases are supported:

   I. All objects are of same size (sizes == NULL)
     A. Data is already packed by destination proc (indices_to == NULL)
        Can just send straight out of user space.
     B. Data not packed by destination proc (indices_to != NULL)
	Need to allocate send buff and copy it in before sending

   II. Objects of variable sizes (sizes != NULL)
       Same two cases above, but offsets are more complicated.
       Need stars_to_ptr/_from_ptr to address memory for (A), and
       indices_to_ptr/_from_ptr to address memory for (B).

   III. On receive, we need to do opposite of the send.
	Can receive directly into user space if  indices_from == NULL
	Otherwise need to receive in buffer and copy.
*/

/*****************************************************************************/

int       Zoltan_Comm_Do(
ZOLTAN_COMM_OBJ * plan,		/* communication data structure */
int tag,			/* message tag for communicating */
char *send_data,		/* array of data I currently own */
int nbytes,			/* multiplier for sizes */
char *recv_data)		/* array of data I'll own after comm */
{
   int status;
   
   status = Zoltan_Comm_Do_Post (plan, tag, send_data, nbytes, recv_data);
   if (status == ZOLTAN_OK)
      status = Zoltan_Comm_Do_Wait (plan, tag, send_data, nbytes, recv_data);
   return status;
}

/*****************************************************************************/
static char  *recv_buff;	/* place to receive messages */
/* 
 * KDDKDD This static variable will not work when multiple posts are issued
 * KDDKDD using different plans before the corresponding waits are issued. 
 */

/*****************************************************************************/

int       Zoltan_Comm_Do_Post(
ZOLTAN_COMM_OBJ * plan,		/* communication data structure */
int tag,			/* message tag for communicating */
char *send_data,		/* array of data I currently own */
int nbytes,			/* multiplier for sizes */
char *recv_data)		/* array of data I'll own after comm */
{
    char     *send_buff;	/* space to buffer outgoing data */
    int       my_proc;		/* processor ID */
    int       self_recv_address;/* where in recv_data self info starts */
    int       self_num;		/* where in send list my_proc appears */
    int       offset;		/* offset into array I'm copying into */
    int       self_index;	/* send offset for data I'm keeping */
    int       out_of_mem;	/* am I out of memory? */
    int       nblocks;		/* number of procs who need my data */
    int       proc_index;	/* loop counter over procs to send to */
    int       i, j, k, jj;	/* loop counters */

    static char *yo = "Zoltan_Comm_Do_Post";


    /* Check input parameters */
    if (!plan) {
        MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
	ZOLTAN_COMM_ERROR("Communication plan = NULL", yo, my_proc);
	return ZOLTAN_FATAL;
    }

    MPI_Comm_rank(plan->comm, &my_proc);

    if ((plan->nsends + plan->self_msg) && !send_data) {
	ZOLTAN_COMM_ERROR("nsends not zero, but send_data = NULL", yo, my_proc);
	return ZOLTAN_FATAL;
    }
    if ((plan->nrecvs + plan->self_msg) && !recv_data) {
	ZOLTAN_COMM_ERROR("nrecvs not zero, but recv_data = NULL", yo, my_proc);
	return ZOLTAN_FATAL;
    }
    if (nbytes < 0) {
	ZOLTAN_COMM_ERROR("Scale factor nbytes is negative", yo, my_proc);
	return ZOLTAN_FATAL;
    }


    /* Post irecvs */

    out_of_mem = 0;

    if (plan->indices_from == NULL) {
	/* Data can go directly into user space. */
	recv_buff = recv_data;
    }
    else {			/* Need to buffer receive to reorder */
	recv_buff = (char *) ZOLTAN_MALLOC(plan->total_recv_size * nbytes);
	if (recv_buff == NULL && plan->total_recv_size * nbytes != 0)
	    out_of_mem = 1;
    }

    if (!out_of_mem) {
	if (plan->sizes == NULL) {	/* All data the same size */
	    k = 0;
	    for (i = 0; i < plan->nrecvs + plan->self_msg; i++) {
		if (plan->procs_from[i] != my_proc) {
		    MPI_Irecv((void *) &recv_buff[plan->starts_from[i] * nbytes],
			      plan->lengths_from[i] * nbytes,
			      (MPI_Datatype) MPI_BYTE, plan->procs_from[i], tag,
			      plan->comm, &plan->request[k]);
		    k++;
		}
		else {
		    self_recv_address = plan->starts_from[i] * nbytes;
		}
	    }
	}

	else {			/* Data of varying sizes */
	    k = 0;
	    for (i = 0; i < plan->nrecvs + plan->self_msg; i++) {
		if (plan->procs_from[i] != my_proc) {
		    MPI_Irecv((void *) &recv_buff[plan->starts_from_ptr[i] * nbytes],
			      plan->sizes_from[i] * nbytes,
			      (MPI_Datatype) MPI_BYTE, plan->procs_from[i], tag,
			      plan->comm, &plan->request[k]);
		    k++;
		}
		else {
		    self_recv_address = plan->starts_from_ptr[i] * nbytes;
		}
	    }
	}
    }


    /* Do remaining allocation to check for any mem problems. */
    if (plan->indices_to != NULL) {	/* can't sent straight from input */
	send_buff = (char *) ZOLTAN_MALLOC(plan->max_send_size * nbytes);
	if (send_buff == 0 && plan->max_send_size * nbytes != 0)
	    out_of_mem = 1;
    }
    else
	send_buff = NULL;

    /* Barrier to ensure irecvs are posted before doing any sends. */
    /* Simultaneously see if anyone out of memory */

    MPI_Allreduce(&out_of_mem, &j, 1, MPI_INT, MPI_SUM, plan->comm);

    if (j > 0) {		/* Some proc is out of memory -> Punt */
	ZOLTAN_FREE((void **) &send_buff);
	if (plan->indices_from != NULL)
	    ZOLTAN_FREE((void **) &recv_buff);
	return (ZOLTAN_MEMERR);
    }

    /* Send out data */

    /* Scan through procs_to list to start w/ higher numbered procs */
    /* This should balance message traffic. */

    nblocks = plan->nsends + plan->self_msg;
    proc_index = 0;
    while (proc_index < nblocks && plan->procs_to[proc_index] < my_proc)
	proc_index++;
    if (proc_index == nblocks)
	proc_index = 0;

    if (plan->sizes == NULL) {	/* Data all of same size */
	if (plan->indices_to == NULL) {	/* data already blocked by processor. */
	    for (i = proc_index, j = 0; j < nblocks; j++) {
		if (plan->procs_to[i] != my_proc) {
		    MPI_Rsend((void *) &send_data[plan->starts_to[i] * nbytes],
			      plan->lengths_to[i] * nbytes,
			      (MPI_Datatype) MPI_BYTE, plan->procs_to[i], tag,
			      plan->comm);
		}
		else
		    self_num = i;
		if (++i == nblocks)
		    i = 0;
	    }

	    if (plan->self_msg) {	/* Copy data to self. */
		memcpy(&recv_buff[self_recv_address],
		       &send_data[plan->starts_to[self_num] * nbytes],
		       plan->lengths_to[self_num] * nbytes);
	    }
	}

	else {			/* Not blocked by processor.  Need to buffer. */
	    for (i = proc_index, jj = 0; jj < nblocks; jj++) {
		if (plan->procs_to[i] != my_proc) {
		    /* Need to pack message first. */
		    offset = 0;
		    j = plan->starts_to[i];
		    for (k = 0; k < plan->lengths_to[i]; k++) {
			memcpy(&send_buff[offset],
			       &send_data[plan->indices_to[j++] * nbytes], nbytes);
			offset += nbytes;
		    }
		    MPI_Rsend((void *) send_buff, plan->lengths_to[i] * nbytes,
		      (MPI_Datatype) MPI_BYTE, plan->procs_to[i], tag, plan->comm);
		}
		else {
		    self_num = i;
		    self_index = plan->starts_to[i];
		}
		if (++i == nblocks)
		    i = 0;
	    }
	    if (plan->self_msg) {	/* Copy data to self. */
		for (k = 0; k < plan->lengths_to[self_num]; k++) {
		    memcpy(&recv_buff[self_recv_address],
		      &send_data[plan->indices_to[self_index++] * nbytes], nbytes);
		    self_recv_address += nbytes;
		}
	    }

	    ZOLTAN_FREE((void **) &send_buff);
	}
    }
    else {			/* Data of differing sizes */
	if (plan->indices_to == NULL) {	/* data already blocked by processor. */
	    for (i = proc_index, j = 0; j < nblocks; j++) {

		if (plan->procs_to[i] != my_proc) {
		    MPI_Rsend((void *) &send_data[plan->starts_to_ptr[i] * nbytes],
			      plan->sizes_to[i] * nbytes,
			      (MPI_Datatype) MPI_BYTE, plan->procs_to[i], tag,
			      plan->comm);
		}
		else
		    self_num = i;
		if (++i == nblocks)
		    i = 0;
	    }

	    if (plan->self_msg) {	/* Copy data to self. */
		memcpy(&recv_buff[self_recv_address],
		       &send_data[plan->starts_to_ptr[self_num] * nbytes],
		       plan->sizes_to[self_num] * nbytes);
	    }
	}

	else {			/* Not blocked by processor.  Need to buffer. */
	    for (i = proc_index, jj = 0; jj < nblocks; jj++) {
		if (plan->procs_to[i] != my_proc) {
		    /* Need to pack message first. */
		    offset = 0;
		    j = plan->starts_to[i];
		    for (k = 0; k < plan->lengths_to[i]; k++) {
			memcpy(&send_buff[offset],
			       &send_data[plan->indices_to_ptr[j] * nbytes],
			       plan->sizes[plan->indices_to[j]] * nbytes);
			offset += plan->sizes[plan->indices_to[j]] * nbytes;
			j++;
		    }
		    MPI_Rsend((void *) send_buff, plan->sizes_to[i] * nbytes,
		      (MPI_Datatype) MPI_BYTE, plan->procs_to[i], tag, plan->comm);
		}
		else
		    self_num = i;
		if (++i == nblocks)
		    i = 0;
	    }
	    if (plan->self_msg) {	/* Copy data to self. */
		j = plan->starts_to[self_num];
		for (k = 0; k < plan->lengths_to[self_num]; k++) {
		    jj = plan->indices_to_ptr[j];
		    memcpy(&recv_buff[self_recv_address],
			   &send_data[jj * nbytes],
			   plan->sizes[plan->indices_to[j]] * nbytes);
		    self_recv_address += plan->sizes[plan->indices_to[j]] * nbytes;
		    j++;
		}
	    }

	    ZOLTAN_FREE((void **) &send_buff);
	}
    }
    return (ZOLTAN_OK);
}


/*****************************************************************************/

int       Zoltan_Comm_Do_Wait(
ZOLTAN_COMM_OBJ * plan,		/* communication data structure */
int tag,			/* message tag for communicating */
char *send_data,		/* array of data I currently own */
int nbytes,			/* multiplier for sizes */
char *recv_data)		/* array of data I'll own after comm */
{
    MPI_Status status;		/* return from Waitany */
    int       my_proc;		/* processor ID */
    int       self_num;		/* where in send list my_proc appears */
    int       i, j, k, jj;	/* loop counters */

    MPI_Comm_rank(plan->comm, &my_proc);    
    
    /* Wait for messages to arrive & unpack them if necessary. */
    /* Note: since request is in plan, could wait in later routine. */

    if (plan->indices_from == NULL) {	/* No copying required */
        if (plan->nrecvs > 0) {
	    MPI_Waitall(plan->nrecvs, plan->request, plan->status);
	}
    }

    else {			 	/* Need to copy into recv_data. */
	if (plan->self_msg) {		/* Unpack own data before waiting */
	    for (self_num = 0; self_num < plan->nrecvs + plan->self_msg; self_num++) 
		if (plan->procs_from[self_num] == my_proc) break;
	    k = plan->starts_from[self_num];
	    for (j = plan->lengths_from[self_num]; j; j--) {
		memcpy(&recv_data[plan->indices_from[k] * nbytes],
		    &recv_buff[k * nbytes], nbytes);
		k++;
	    }

	}
	else
	    self_num = plan->nrecvs;

	for (jj = 0; jj < plan->nrecvs; jj++) {

	    MPI_Waitany(plan->nrecvs, plan->request, &i, &status);

	    if (i >= self_num) i++;

	    k = plan->starts_from[i];
	    for (j = plan->lengths_from[i]; j; j--) {
		memcpy(&recv_data[plan->indices_from[k] * nbytes],
		    &recv_buff[k * nbytes], nbytes);
		k++;
	    }
	}

	ZOLTAN_FREE((void **) &recv_buff);
    }

    return (ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
