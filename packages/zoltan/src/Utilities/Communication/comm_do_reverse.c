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
#include <stdlib.h>
#include <mpi.h>
#include "zoltan_mem.h"
#include "comm.h"


/******************************************************************************/

int Zoltan_Comm_Do_Reverse(
ZOLTAN_COMM_OBJ *plan,		/* communication data structure */
int       tag,			    /* message tag for communicating */
char     *send_data,		/* array of data I currently own */
int       nbytes,		    /* # bytes per data item */
int      *sizes,		    /* variable size of objects (if not NULL) */
char     *recv_data)		/* array of data I'll own after reverse comm */
{
    int status;
    
    status = Zoltan_Comm_Do_Reverse_Post (plan, tag, send_data, nbytes, sizes,
     recv_data);
    if (status == ZOLTAN_OK)
        status = Zoltan_Comm_Do_Reverse_Wait (plan, tag, send_data, nbytes,
          sizes, recv_data);
    return status;
}    
    
/******************************************************************************/
static ZOLTAN_COMM_OBJ *plan_reverse;	/* communication data structure */

/* 
 * KDDKDD This static variable will not work when multiple posts are issued
 * KDDKDD using different plans before the corresponding waits are issued. 
 */

/******************************************************************************/
/* Perform a reverse communication operation.  Communication object describes */
/* an action, and this routine does the opposite.  Can be used to return */
/* updated data to originating processor. */

int       Zoltan_Comm_Do_Reverse_Post(
ZOLTAN_COMM_OBJ *plan,		/* communication data structure */
int       tag,			/* message tag for communicating */
char     *send_data,		/* array of data I currently own */
int       nbytes,		/* # bytes per data item */
int      *sizes,		/* variable size of objects (if not NULL) */
char     *recv_data)		/* array of data I'll own after reverse comm */
{
    int       total_send_length;/* total message length I send in plan */
    int       max_recv_length;	/* biggest message I recv in plan */
    int       sum_recv_sizes;	/* sum of the item sizes I receive */
    int       comm_flag;		/* status flag */
    int       i;		/* loop counter */
    static char *yo = "Zoltan_Comm_Do_Reverse_Post";

    /* Check input parameters */
    if (!plan){
      ZOLTAN_COMM_ERROR("Communication plan = NULL.", yo, -1);
      return ZOLTAN_FATAL;
    }

    /* Let Zoltan_Comm_Do check the remaining parameters. */

    total_send_length = 0;
    for (i = 0; i < plan->nsends + plan->self_msg; i++) {
        total_send_length += plan->lengths_to[i];
    }

    max_recv_length = 0;
    for (i = 0; i < plan->nrecvs; i++) {
	if (plan->lengths_from[i] > max_recv_length)
	    max_recv_length = plan->lengths_from[i];
    }

    plan_reverse = (ZOLTAN_COMM_OBJ *) ZOLTAN_MALLOC(sizeof(ZOLTAN_COMM_OBJ));

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
	ZOLTAN_MALLOC(plan_reverse->nrecvs * sizeof(MPI_Request));
    if (plan_reverse->request == NULL && plan_reverse->nrecvs != 0) {
        ZOLTAN_FREE((void **) &plan_reverse);
	return(ZOLTAN_MEMERR);
    }

    plan_reverse->status = (MPI_Status *)
	ZOLTAN_MALLOC(plan_reverse->nrecvs * sizeof(MPI_Status));
    if (plan_reverse->status == NULL && plan_reverse->nrecvs != 0) {
        ZOLTAN_FREE((void **) &(plan_reverse->request));
        ZOLTAN_FREE((void **) &plan_reverse);
	return(ZOLTAN_MEMERR);
    }

    comm_flag = Zoltan_Comm_Resize(plan_reverse, sizes, tag, &sum_recv_sizes);

    if (comm_flag != ZOLTAN_OK && comm_flag != ZOLTAN_WARN) {
        ZOLTAN_FREE((void **) &(plan_reverse->status));
        ZOLTAN_FREE((void **) &(plan_reverse->request));
        ZOLTAN_FREE((void **) &plan_reverse);
	return(comm_flag);
    }

    if (sum_recv_sizes != plan_reverse->total_recv_size){
       /* Sanity check */
       ZOLTAN_FREE((void **) &(plan_reverse->status));
       ZOLTAN_FREE((void **) &(plan_reverse->request));
       ZOLTAN_FREE((void **) &plan_reverse);
       return(ZOLTAN_FATAL);
    }

    comm_flag = Zoltan_Comm_Do_Post(plan_reverse, tag, send_data, nbytes, recv_data);
    return (comm_flag);
}


    
/******************************************************************************/
int Zoltan_Comm_Do_Reverse_Wait(
ZOLTAN_COMM_OBJ *plan,		/* communication data structure */
int       tag,			/* message tag for communicating */
char     *send_data,		/* array of data I currently own */
int       nbytes,		/* # bytes per data item */
int      *sizes,		/* variable size of objects (if not NULL) */
char     *recv_data)		/* array of data I'll own after reverse comm */
{
    int       comm_flag;		/* status flag */
    
    comm_flag = Zoltan_Comm_Do_Wait(plan_reverse, tag, send_data, nbytes, recv_data);
        
    if (sizes != NULL) {
        ZOLTAN_FREE((void *) &plan_reverse->sizes);
	ZOLTAN_FREE((void *) &plan_reverse->sizes_to);
	ZOLTAN_FREE((void *) &plan_reverse->sizes_from);
	ZOLTAN_FREE((void *) &plan_reverse->starts_to_ptr);
	ZOLTAN_FREE((void *) &plan_reverse->starts_from_ptr);
	ZOLTAN_FREE((void *) &plan_reverse->indices_to_ptr);
	ZOLTAN_FREE((void *) &plan_reverse->indices_from_ptr);
    }
    ZOLTAN_FREE((void **) &(plan_reverse->status));
    ZOLTAN_FREE((void **) &(plan_reverse->request));
    ZOLTAN_FREE((void **) &plan_reverse);

    return(comm_flag);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
