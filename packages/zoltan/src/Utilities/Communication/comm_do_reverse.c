/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mem_const.h"
#include "comm_const.h"


/* Perform a reverse communication operation.  Communication object describes */
/* an action, and this routine does the opposite.  Can be used to return */
/* updated data to originating processor. */

int       LB_Comm_Do_Reverse(
struct Comm_Obj *plan,		/* communication data structure */
int       tag,			/* message tag for communicating */
char     *send_data,		/* array of data I currently own */
int       nbytes,		/* # bytes per data item */
int      *sizes,		/* variable size of objects (if not NULL) */
char     *recv_data)		/* array of data I'll own after reverse comm */
{
    struct Comm_Obj *plan_reverse;	/* communication data structure */
    int       my_proc;		/* current processor ID */
    int       total_send_length;/* total message length I send in plan */
    int       max_recv_length;	/* biggest message I recv in plan */
    int       comm_flag;		/* status flag */
    int       i;		/* loop counter */
    static char *yo = "LB_Comm_Do_Reverse";

    /* Check input parameters */
    if (!plan){
      fprintf(stderr, "Zoltan error in %s: Communication plan = NULL\n", 
        yo);
      return COMM_FATAL;
    }

    MPI_Comm_rank(plan->comm, &my_proc);

    /* Let LB_Comm_Do check the remaining parameters. */

    total_send_length = 0;
    for (i = 0; i < plan->nsends + plan->self_msg; i++) {
        total_send_length += plan->lengths_to[i];
    }

    max_recv_length = 0;
    for (i = 0; i < plan->nrecvs; i++) {
	if (plan->lengths_from[i] > max_recv_length)
	    max_recv_length = plan->lengths_from[i];
    }

    plan_reverse = (struct Comm_Obj *) LB_MALLOC(sizeof(struct Comm_Obj));

    plan_reverse->nvals = plan->nvals_recv;
    plan_reverse->nvals_recv = plan->nvals;
    plan_reverse->lengths_to = plan->lengths_from;
    plan_reverse->procs_to = plan->procs_from;
    plan_reverse->indices_to = plan->indices_from;
    plan_reverse->starts_to = plan->starts_from;
    plan_reverse->lengths_from = plan->lengths_to;
    plan_reverse->procs_from = plan->procs_to;
    plan_reverse->indices_from = plan->indices_to;
    plan_reverse->starts_from = plan->starts_to;
    plan_reverse->nrecvs = plan->nsends;
    plan_reverse->nsends = plan->nrecvs;
    plan_reverse->self_msg = plan->self_msg;
    plan_reverse->max_send_size = max_recv_length;
    plan_reverse->total_recv_size = total_send_length;
    plan_reverse->comm = plan->comm;
    plan_reverse->sizes = NULL;
    plan_reverse->sizes_to = NULL;
    plan_reverse->sizes_from = NULL;
    plan_reverse->starts_to_ptr = NULL;
    plan_reverse->starts_from_ptr = NULL;
    plan_reverse->indices_to_ptr = NULL;
    plan_reverse->indices_from_ptr = NULL;


    plan_reverse->request = (MPI_Request *)
	LB_MALLOC(plan_reverse->nrecvs * sizeof(MPI_Request));
    if (plan_reverse->request == NULL && plan_reverse->nrecvs != 0) {
        LB_FREE((void **) &plan_reverse);
	return(COMM_MEMERR);
    }

    plan_reverse->status = (MPI_Status *)
	LB_MALLOC(plan_reverse->nrecvs * sizeof(MPI_Status));
    if (plan_reverse->status == NULL && plan_reverse->nrecvs != 0) {
        LB_FREE((void **) &(plan_reverse->request));
        LB_FREE((void **) &plan_reverse);
	return(COMM_MEMERR);
    }

    if (sizes != NULL) {
	comm_flag = LB_Comm_Resize(plan_reverse, sizes, tag);

        if (comm_flag != COMM_OK && comm_flag != COMM_WARN) {
            LB_FREE((void **) &(plan_reverse->status));
            LB_FREE((void **) &(plan_reverse->request));
            LB_FREE((void **) &plan_reverse);
	    return(comm_flag);
        }
    }

    comm_flag = LB_Comm_Do(plan_reverse, tag, send_data, nbytes, recv_data);

    if (sizes != NULL) {
        LB_FREE((void *) &plan_reverse->sizes);
	LB_FREE((void *) &plan_reverse->sizes_to);
	LB_FREE((void *) &plan_reverse->sizes_from);
	LB_FREE((void *) &plan_reverse->starts_to_ptr);
	LB_FREE((void *) &plan_reverse->starts_from_ptr);
	LB_FREE((void *) &plan_reverse->indices_to_ptr);
	LB_FREE((void *) &plan_reverse->indices_from_ptr);
    }
    LB_FREE((void **) &(plan_reverse->status));
    LB_FREE((void **) &(plan_reverse->request));
    LB_FREE((void **) &plan_reverse);

    return(comm_flag);
}
