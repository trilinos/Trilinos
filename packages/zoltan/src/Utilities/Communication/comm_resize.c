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

int       Zoltan_Comm_Resize(
ZOLTAN_COMM_OBJ *plan,			/* communication plan object */
int      *sizes,		/* size of each item I'm sending */
int       tag,			/* message tag I can use */
int      *sum_recv_sizes)       /* sum of the sizes of the items I'll receive */
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
    int       var_sizes;        /* items have variable sizes? */
    int       i, j, k;		/* loop counters */
    static char *yo = "Zoltan_Comm_Resize";


    /* If sizes vary, then I need to compute and communicate message lengths */
    /* First check if sizes array is NULL on all procs. */

    MPI_Comm_rank(plan->comm, &my_proc);
    i = (sizes != NULL);
    MPI_Allreduce(&i, &var_sizes, 1, MPI_INT, MPI_LOR, plan->comm);
    

    if (var_sizes && plan->indices_from != NULL) {
        /* Can't do w/o individual item sizes */
	ZOLTAN_COMM_ERROR("Non-blocked, variable-sized recvs not supported", yo, my_proc);
	return(ZOLTAN_FATAL);
    }

    ZOLTAN_FREE(&plan->sizes);
    ZOLTAN_FREE(&plan->sizes_to);
    ZOLTAN_FREE(&plan->sizes_from);
    ZOLTAN_FREE(&plan->starts_to_ptr);
    ZOLTAN_FREE(&plan->starts_from_ptr);
    ZOLTAN_FREE(&plan->indices_to_ptr);
    ZOLTAN_FREE(&plan->indices_from_ptr);

    nsends = plan->nsends;
    nrecvs = plan->nrecvs;
    self_msg = plan->self_msg;

    if (!var_sizes) { /* Easy case.  Size = length */
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
	return_flag = ZOLTAN_OK;
    }

    else {		/* Need to actually compute message sizes */
	plan->sizes = (int *) ZOLTAN_MALLOC((plan->nvals + 1) * sizeof(int));
	if (plan->sizes == NULL) {
	    return_flag = ZOLTAN_MEMERR;
	    goto Mem_Err;
	}
	for (i = 0; i < plan->nvals; i++) plan->sizes[i] = sizes[i];

	return_flag = ZOLTAN_OK;
	sizes_to = (int *) ZOLTAN_MALLOC((nsends + self_msg) * sizeof(int));
	if (sizes_to == NULL && (nsends + self_msg) != 0) {
	    return_flag = ZOLTAN_MEMERR;
	    goto Mem_Err;
	}

	sizes_from = (int *) ZOLTAN_MALLOC((nrecvs + self_msg) * sizeof(int));
	if (sizes_from == NULL && (nrecvs + self_msg) != 0) {
	    return_flag = ZOLTAN_MEMERR;
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
	starts_to_ptr = (int *) ZOLTAN_MALLOC((nsends + self_msg) * sizeof(int));
	if (starts_to_ptr == NULL && (nsends + self_msg) != 0) {
	    return_flag = ZOLTAN_MEMERR;
	    goto Mem_Err;
	}

	if (plan->indices_to == NULL) {
	    /* Simpler case; sends already blocked by processor */

	    index = (int *) ZOLTAN_MALLOC((nsends + self_msg) * sizeof(int));
	    sort_val = (int *) ZOLTAN_MALLOC((nsends + self_msg) * sizeof(int));
	    if ((index == NULL || sort_val == NULL) && nsends + self_msg > 0) {
	        return_flag = ZOLTAN_MEMERR;
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
	    Zoltan_Comm_Sort_Ints(sort_val, index, nsends + self_msg);


	    sum = 0;
	    for (i = 0; i < nsends + self_msg; i++) {
		starts_to_ptr[index[i]] = sum;
		sum += sizes_to[index[i]];
	    }

	    ZOLTAN_FREE(&index);
	    ZOLTAN_FREE(&sort_val);
	}

	else {		/* Harder case, sends not blocked */
	    offset = (int *) ZOLTAN_MALLOC(plan->nvals * sizeof(int));
	    indices_to_ptr = (int *) ZOLTAN_MALLOC(plan->nvals * sizeof(int));
	    if ((offset == NULL || indices_to_ptr == NULL) && plan->nvals != 0) {
	        return_flag = ZOLTAN_MEMERR;
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
	    ZOLTAN_FREE(&offset);
	}


	/* Note: This routine only gets message sizes, not object sizes. */
	/*	Anything requiring item sizes requires more code */
	/*      Should such functionality reside here? */

	Zoltan_Comm_Exchange_Sizes(sizes_to, plan->procs_to, nsends, self_msg,
	    sizes_from, plan->procs_from, nrecvs, 
	    &plan->total_recv_size, my_proc, tag, plan->comm);

	starts_from_ptr = (int *) ZOLTAN_MALLOC((nrecvs + self_msg) * sizeof(int));
	if (starts_from_ptr == NULL && (nrecvs + self_msg) != 0) {
	    return_flag = ZOLTAN_MEMERR;
	    goto Mem_Err;
	}

	if (plan->indices_from == NULL) {
	    /* Simpler case; recvs already blocked by processor */
	    /* Harder case currently excluded at top of file */

	    index = (int *) ZOLTAN_MALLOC((nrecvs + self_msg) * sizeof(int));
	    sort_val = (int *) ZOLTAN_MALLOC((nrecvs + self_msg) * sizeof(int));
	    if ((index == NULL || sort_val == NULL) && nrecvs + self_msg > 0) {
	        return_flag = ZOLTAN_MEMERR;
	        goto Mem_Err;
	    }

	    for (i = 0; i < nrecvs + self_msg; i++) {
		sort_val[i] = plan->starts_from[i];
		index[i] = i;
	    }
	    Zoltan_Comm_Sort_Ints(sort_val, index, nrecvs + self_msg);

	    sum = 0;
	    for (i = 0; i < nrecvs + self_msg; i++) {
		starts_from_ptr[index[i]] = sum;
		sum += sizes_from[index[i]];
	    }

	    ZOLTAN_FREE(&index);
	    ZOLTAN_FREE(&sort_val);
	}

	/*else {*/	/* Harder case, recvs not blocked */
	    /* Not currently supported */
	    /* Can't do w/o individual item sizes */
	    /* Currently checked for at top of file */
	/*}*/
    }

Mem_Err:
    if (return_flag == ZOLTAN_MEMERR) {
	ZOLTAN_FREE(&index);
	ZOLTAN_FREE(&sort_val);
	ZOLTAN_FREE(&offset);
	ZOLTAN_FREE(&plan->sizes);
	ZOLTAN_FREE(&plan->sizes_to);
	ZOLTAN_FREE(&plan->sizes_from);
	ZOLTAN_FREE(&plan->starts_to_ptr);
	ZOLTAN_FREE(&plan->starts_from_ptr);
	ZOLTAN_FREE(&plan->indices_to_ptr);
	ZOLTAN_FREE(&plan->indices_from_ptr);
    }

    plan->sizes_to = sizes_to;
    plan->sizes_from = sizes_from;
    plan->starts_to_ptr = starts_to_ptr;
    plan->starts_from_ptr = starts_from_ptr;
    plan->indices_to_ptr = indices_to_ptr;
    plan->indices_from_ptr = indices_from_ptr;

    if (sum_recv_sizes)
      *sum_recv_sizes = plan->total_recv_size;

    return(return_flag);

}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
