#include <stdio.h>
#include "mpi.h"
#include "comm_const.h"
#include "mem_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int       LB_Comm_Resize(
COMM_OBJ *plan,			/* communication plan object */
int      *sizes,		/* size of each item I'm sending */
int       tag)			/* message tag I can use */
{
    int       return_flag;	/* returned error code */
    int       my_proc;		/* my processor ID */
    int       nsends, nrecvs;	/* number of msgs I'll send & recv */
    int       self_msg;		/* do I have data for myself? */
    int       sum;		/* running total of item sizes */
    int      *offset = NULL;	/* address of each item in send array */
    int      *index = NULL;	/* indices into start_from_ptr */
    int      *sort_val = NULL;	/* start_from to order properly */
    int      *sizes_to = NULL;  /* size of each msg to send (if items vary) */
    int      *sizes_from = NULL;/* size of each msg to recv (if items vary) */
    int      *starts_to_ptr = NULL;	/* where in dense vector sends start */
    int      *starts_from_ptr = NULL;	/* where in dense vector recvs start */
    int      *indices_to_ptr = NULL;	/* where to find items I send */
					/* ordered like procs_to */
    int      *indices_from_ptr = NULL;	/* where to find items I recv */
					/* ordered like procs_from */
    int       i, j, k;		/* loop counters */
    static char *yo = "LB_Comm_Resize";


    /* If sizes vary, then I need to compute and communicate message lengths */

    MPI_Comm_rank(plan->comm, &my_proc);

    if (sizes != NULL && plan->indices_from != NULL) {
        /* Can't do w/o individual item sizes */
	COMM_ERROR("Non-blocked, variable-sized recvs not supported", yo, my_proc);
	return(COMM_FATAL);
    }

    LB_FREE((void *) &plan->sizes);
    LB_FREE((void *) &plan->sizes_to);
    LB_FREE((void *) &plan->sizes_from);
    LB_FREE((void *) &plan->starts_to_ptr);
    LB_FREE((void *) &plan->starts_from_ptr);
    LB_FREE((void *) &plan->indices_to_ptr);
    LB_FREE((void *) &plan->indices_from_ptr);

    nsends = plan->nsends;
    nrecvs = plan->nrecvs;
    self_msg = plan->self_msg;

    if (sizes == NULL) { /* Easy case.  Size = length */
	plan->total_recv_size = 0;
	for (i = 0; i < nrecvs + self_msg; i++) {
	    plan->total_recv_size += plan->lengths_from[i];
	}

	plan->max_send_size = 0;
	for (i = 0; i < nsends + self_msg; i++) {
	    if (plan->procs_to[i] != my_proc && 
		plan->lengths_to[i] > plan->max_send_size) {
	        plan->max_send_size = plan->lengths_to[i];
	    }
	}
	return_flag = COMM_OK;
    }

    else {		/* Need to actually compute message sizes */
	plan->sizes = (int *) LB_MALLOC((plan->nvals + 1) * sizeof(int));
	if (plan->sizes == NULL) {
	    return_flag = COMM_MEMERR;
	    goto Mem_Err;
	}
	for (i = 0; i < plan->nvals; i++) plan->sizes[i] = sizes[i];

	return_flag = COMM_OK;
	sizes_to = (int *) LB_MALLOC((nsends + self_msg) * sizeof(int));
	if (sizes_to == NULL && (nsends + self_msg) != 0) {
	    return_flag = COMM_MEMERR;
	    goto Mem_Err;
	}

	sizes_from = (int *) LB_MALLOC((nrecvs + self_msg) * sizeof(int));
	if (sizes_from == NULL && (nrecvs + self_msg) != 0) {
	    return_flag = COMM_MEMERR;
	    goto Mem_Err;
	}

	for (i = 0; i < nsends + self_msg; i++)
	    sizes_to[i] = 0;


	/* Several cases:
	   1. indices_to == NULL
		=> starts_to != NULL, need to allocate, set starts_to_ptr
	   2. indices_to != NULL (=> starts_to == NULL)
		need to allocate, set indices_to_ptr
	   3,4. mirror cases for _from
	*/
	starts_to_ptr = (int *) LB_MALLOC((nsends + self_msg) * sizeof(int));
	if (starts_to_ptr == NULL && (nsends + self_msg) != 0) {
	    return_flag = COMM_MEMERR;
	    goto Mem_Err;
	}

	if (plan->indices_to == NULL) {
	    /* Simpler case; sends already blocked by processor */

	    index = (int *) LB_MALLOC((nsends + self_msg) * sizeof(int));
	    sort_val = (int *) LB_MALLOC((nsends + self_msg) * sizeof(int));
	    if ((index == NULL || sort_val == NULL) && nsends + self_msg > 0) {
	        return_flag = COMM_MEMERR;
	        goto Mem_Err;
	    }

	    for (i = 0; i < nsends + self_msg; i++) {
	        j = plan->starts_to[i];
	        for (k = 0; k < plan->lengths_to[i]; k++) {
	 	    sizes_to[i] += sizes[j++];
	        }
	        if (sizes_to[i] > plan->max_send_size && 
		    plan->procs_to[i] != my_proc)
	            plan->max_send_size = sizes_to[i];
	    }

	    for (i = 0; i < nsends + self_msg; i++) {
		sort_val[i] = plan->starts_to[i];
		index[i] = i;
	    }
	    LB_Comm_Sort_Ints(sort_val, index, nsends + self_msg);


	    sum = 0;
	    for (i = 0; i < nsends + self_msg; i++) {
		starts_to_ptr[index[i]] = sum;
		sum += sizes_to[index[i]];
	    }

	    LB_FREE((void *) &index);
	    LB_FREE((void *) &sort_val);
	}

	else {		/* Harder case, sends not blocked */
	    offset = (int *) LB_MALLOC(plan->nvals * sizeof(int));
	    indices_to_ptr = (int *) LB_MALLOC(plan->nvals * sizeof(int));
	    if ((offset == NULL || indices_to_ptr == NULL) && plan->nvals != 0) {
	        return_flag = COMM_MEMERR;
	        goto Mem_Err;
	    }

	    /* Compute address for every item in send array */
	    sum = 0;
	    for (i = 0; i < plan->nvals; i++) {
		offset[i] = sum;
		sum += sizes[i];
	    }

	    sum = 0;
	    plan->max_send_size = 0;
	    for (i = 0; i < nsends + self_msg; i++) {
		starts_to_ptr[i] = sum;
	        j = plan->starts_to[i];
	        for (k = 0; k < plan->lengths_to[i]; k++) {
		    indices_to_ptr[j] = offset[plan->indices_to[j]];
	 	    sizes_to[i] += sizes[plan->indices_to[j++]];
	        }
	        if (sizes_to[i] > plan->max_send_size && 
		    plan->procs_to[i] != my_proc)
	            plan->max_send_size = sizes_to[i];
		sum += sizes_to[i];
	    }
	    LB_FREE((void *) &offset);
	}


	/* Note: This routine only gets message sizes, not object sizes. */
	/*	Anything requiring item sizes requires more code */
	/*      Should such functionality reside here? */

	LB_Comm_Exchange_Sizes(plan->procs_to, nsends, self_msg,
	    plan->procs_from, nrecvs, sizes_to, sizes_from,
	    &plan->total_recv_size, my_proc, tag, plan->comm);

	starts_from_ptr = (int *) LB_MALLOC((nrecvs + self_msg) * sizeof(int));
	if (starts_from_ptr == NULL && (nrecvs + self_msg) != 0) {
	    return_flag = COMM_MEMERR;
	    goto Mem_Err;
	}

	if (plan->indices_from == NULL) {
	    /* Simpler case; recvs already blocked by processor */
	    /* Harder case currently excluded at top of file */

	    index = (int *) LB_MALLOC((nrecvs + self_msg) * sizeof(int));
	    sort_val = (int *) LB_MALLOC((nrecvs + self_msg) * sizeof(int));
	    if ((index == NULL || sort_val == NULL) && nrecvs + self_msg > 0) {
	        return_flag = COMM_MEMERR;
	        goto Mem_Err;
	    }

	    for (i = 0; i < nrecvs + self_msg; i++) {
		sort_val[i] = plan->starts_from[i];
		index[i] = i;
	    }
	    LB_Comm_Sort_Ints(sort_val, index, nrecvs + self_msg);

	    sum = 0;
	    for (i = 0; i < nrecvs + self_msg; i++) {
		starts_from_ptr[index[i]] = sum;
		sum += sizes_from[index[i]];
	    }

	    LB_FREE((void *) &index);
	    LB_FREE((void *) &sort_val);
	}

	/*else {*/	/* Harder case, recvs not blocked */
	    /* Not currently supported */
	    /* Can't do w/o individual item sizes */
	    /* Currently checked for at top of file */
	/*}*/
    }

Mem_Err:
    if (return_flag == COMM_MEMERR) {
	LB_FREE((void *) &index);
	LB_FREE((void *) &sort_val);
	LB_FREE((void *) &offset);
	LB_FREE((void *) &plan->sizes);
	LB_FREE((void *) &plan->sizes_to);
	LB_FREE((void *) &plan->sizes_from);
	LB_FREE((void *) &plan->starts_to_ptr);
	LB_FREE((void *) &plan->starts_from_ptr);
	LB_FREE((void *) &plan->indices_to_ptr);
	LB_FREE((void *) &plan->indices_from_ptr);
    }

    plan->sizes_to = sizes_to;
    plan->sizes_from = sizes_from;
    plan->starts_to_ptr = starts_to_ptr;
    plan->starts_from_ptr = starts_from_ptr;
    plan->indices_to_ptr = indices_to_ptr;
    plan->indices_from_ptr = indices_from_ptr;

    return(return_flag);

}
