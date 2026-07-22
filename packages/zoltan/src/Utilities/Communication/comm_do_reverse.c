// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "zoltan_mem.h"
#include "comm.h"


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

static void free_reverse_plan(ZOLTAN_COMM_OBJ *plan)
{
  if (plan && (plan->plan_reverse != NULL)) {
    if (plan->plan_reverse->sizes != NULL) {
        ZOLTAN_FREE(&plan->plan_reverse->sizes);
        ZOLTAN_FREE(&plan->plan_reverse->sizes_to);
        ZOLTAN_FREE(&plan->plan_reverse->sizes_from);
        ZOLTAN_FREE(&plan->plan_reverse->starts_to_ptr);
        ZOLTAN_FREE(&plan->plan_reverse->starts_from_ptr);
        ZOLTAN_FREE(&plan->plan_reverse->indices_to_ptr);
        ZOLTAN_FREE(&plan->plan_reverse->indices_from_ptr);
    }
    ZOLTAN_FREE(&(plan->plan_reverse->status));
    ZOLTAN_FREE(&(plan->plan_reverse->request));
    ZOLTAN_FREE(&plan->plan_reverse);
  }
}

static int create_reverse_plan( 
ZOLTAN_COMM_OBJ *plan,		/* communication data structure */
int      tag,
int      *sizes)		/* variable size of objects (if not NULL) */
{
    int       total_send_length;/* total message length I send in plan */
    int       max_recv_length;	/* biggest message I recv in plan */
    int       sum_recv_sizes;	/* sum of the item sizes I receive */
    int       comm_flag;		/* status flag */
    int       i;		/* loop counter */
    static char *yo = "Zoltan_Comm_Reverse_Plan";

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

    plan->plan_reverse =(ZOLTAN_COMM_OBJ*)ZOLTAN_MALLOC(sizeof(ZOLTAN_COMM_OBJ));

    plan->plan_reverse->nvals = plan->nvals_recv;
    plan->plan_reverse->nvals_recv = plan->nvals;
    plan->plan_reverse->lengths_to = plan->lengths_from;
    plan->plan_reverse->procs_to = plan->procs_from;
    plan->plan_reverse->indices_to = plan->indices_from;
    plan->plan_reverse->starts_to = plan->starts_from;
    plan->plan_reverse->lengths_from = plan->lengths_to;
    plan->plan_reverse->procs_from = plan->procs_to;
    plan->plan_reverse->indices_from = plan->indices_to;
    plan->plan_reverse->starts_from = plan->starts_to;
    plan->plan_reverse->nrecvs = plan->nsends;
    plan->plan_reverse->nsends = plan->nrecvs;
    plan->plan_reverse->self_msg = plan->self_msg;
    plan->plan_reverse->max_send_size = max_recv_length;
    plan->plan_reverse->total_recv_size = total_send_length;
    plan->plan_reverse->comm = plan->comm;
    plan->plan_reverse->sizes = NULL;
    plan->plan_reverse->sizes_to = NULL;
    plan->plan_reverse->sizes_from = NULL;
    plan->plan_reverse->starts_to_ptr = NULL;
    plan->plan_reverse->starts_from_ptr = NULL;
    plan->plan_reverse->indices_to_ptr = NULL;
    plan->plan_reverse->indices_from_ptr = NULL;
    plan->plan_reverse->maxed_recvs = 0;
    plan->plan_reverse->plan_reverse = NULL;
    plan->plan_reverse->recv_buff = NULL;

    if (MPI_RECV_LIMIT > 0){
      /* If we have a limit to the number of posted receives we are allowed,
      ** and our plan has exceeded that, then switch to an MPI_Alltoallv so
      ** that we will have fewer receives posted when we do the communication.
      */
      MPI_Allreduce(&plan->nsends, &i, 1, MPI_INT, MPI_MAX, plan->comm);
      if (i > MPI_RECV_LIMIT){
        plan->plan_reverse->maxed_recvs = 1;
      }
    }

    if (plan->plan_reverse->maxed_recvs){
      plan->plan_reverse->request = NULL;
      plan->plan_reverse->status = NULL;
    }
    else{
      plan->plan_reverse->request = (MPI_Request *)
  	ZOLTAN_MALLOC(plan->plan_reverse->nrecvs * sizeof(MPI_Request));
      if (plan->plan_reverse->request == NULL && plan->plan_reverse->nrecvs != 0) {
  	return(ZOLTAN_MEMERR);
      }
  
      plan->plan_reverse->status = (MPI_Status *)
  	ZOLTAN_MALLOC(plan->plan_reverse->nrecvs * sizeof(MPI_Status));
      if (plan->plan_reverse->status == NULL && plan->plan_reverse->nrecvs != 0) {
  	return(ZOLTAN_MEMERR);
      }
    }

    comm_flag =Zoltan_Comm_Resize(plan->plan_reverse,sizes,tag,&sum_recv_sizes);

    if (comm_flag != ZOLTAN_OK) {
	return(comm_flag);
    }

    if (sum_recv_sizes != plan->plan_reverse->total_recv_size){
       /* Sanity check */
       return(ZOLTAN_FATAL);
    }

    return ZOLTAN_OK;
}

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

  /* create plan->plan_reverse
   */
  status = create_reverse_plan(plan, tag, sizes);

  if (status == ZOLTAN_OK){

    if (plan->plan_reverse->maxed_recvs){
  
      /* use MPI_Alltoallv to implement plan->plan_reverse, because comm_do_post
       * would post more receives that allowed on this machine
       */
      
      status = Zoltan_Comm_Do_AlltoAll(plan->plan_reverse, send_data, nbytes, recv_data);
    }
    else{
      /* use post/wait which is faster when each sends to few
       */
      status = Zoltan_Comm_Do_Post(plan->plan_reverse, tag, send_data, nbytes, recv_data);
    
      if (status == ZOLTAN_OK){
        status = Zoltan_Comm_Do_Wait (plan->plan_reverse, tag, send_data, nbytes, recv_data);
      }
    }
  }

  free_reverse_plan(plan);

  return status;
}    
/******************************************************************************/
/* Perform a reverse communication operation.  Communication object describes */
/* an action, and this routine does the opposite.  Can be used to return */
/* updated data to originating processor. */

int       Zoltan_Comm_Do_Reverse_Post(
ZOLTAN_COMM_OBJ *plan,          /* communication data structure */
int       tag,                  /* message tag for communicating */
char     *send_data,            /* array of data I currently own */
int       nbytes,               /* # bytes per data item */
int      *sizes,                /* variable size of objects (if not NULL) */
char     *recv_data)            /* array of data I'll own after reverse comm */
{
  int       comm_flag;                /* status flag */

  /* create plan->plan_reverse
   */
  comm_flag = create_reverse_plan(plan, tag, sizes);

  if (comm_flag == ZOLTAN_OK){

    /* post the receives of plan->plan_reverse
     */
    comm_flag = Zoltan_Comm_Do_Post(plan->plan_reverse, tag, send_data, nbytes,
                recv_data);
  }

  if (comm_flag != ZOLTAN_OK){
    free_reverse_plan(plan);
  }

  return (comm_flag);
}

/******************************************************************************/
int Zoltan_Comm_Do_Reverse_Wait(
ZOLTAN_COMM_OBJ *plan,          /* communication data structure */
int       tag,                  /* message tag for communicating */
char     *send_data,            /* array of data I currently own */
int       nbytes,               /* # bytes per data item */
int      *sizes,                /* variable size of objects (if not NULL) */
char     *recv_data)            /* array of data I'll own after reverse comm */
{
    int       comm_flag;                /* status flag */

    comm_flag = Zoltan_Comm_Do_Wait(plan->plan_reverse, tag, send_data, nbytes, recv_data);

    free_reverse_plan(plan);

    return(comm_flag);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
