/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * Zoltan is distributed under the GNU Lesser General Public License 2.1.    * 
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
#include "comm_const.h"
#include "all_allo_const.h"


/* Perform a reverse communication operation.  Communication object describes */
/* an action, and this routine does the opposite.  Can be used to return */
/* updated data to originating processor. */

int       LB_Comm_Do_Reverse(
struct Comm_Obj *plan,		/* communication data structure */
int       tag,			/* message tag for communicating */
char     *send_data,		/* array of data I currently own */
int       nsize,		/* # bytes per data item */
char     *recv_data)		/* array of data I'll own after comm */
{
    struct Comm_Obj *plan_reverse;	/* communication data structure */
    int       my_proc;		/* current processor ID */
    int       total_send_length;/* total message length I send in plan */
    int       max_recv_length;	/* biggest message I recv in plan */
    int       lb_flag;		/* status flag */
    int       i;		/* loop counter */
    static char *yo = "LB_Comm_Do_Reverse";

    /* Check input parameters */
    if (!plan){
      fprintf(stderr, "Zoltan error in %s: Communication plan = NULL\n", 
        yo);
      return LB_FATAL;
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

    plan_reverse->lengths_to = plan->lengths_from;
    plan_reverse->procs_to = plan->procs_from;
    plan_reverse->indices_to = plan->indices_from;
    plan_reverse->lengths_from = plan->lengths_to;
    plan_reverse->procs_from = plan->procs_to;
    plan_reverse->indices_from = plan->indices_to;
    plan_reverse->nrecvs = plan->nsends;
    plan_reverse->nsends = plan->nrecvs;
    plan_reverse->self_msg = plan->self_msg;
    plan_reverse->max_send_length = max_recv_length;
    plan_reverse->total_recv_length = total_send_length;
    plan_reverse->comm = plan->comm;
    plan_reverse->request = (MPI_Request *)
	LB_MALLOC(plan_reverse->nrecvs * sizeof(MPI_Request));
    plan_reverse->status = (MPI_Status *)
	LB_MALLOC(plan_reverse->nrecvs * sizeof(MPI_Status));

    lb_flag = LB_Comm_Do(plan_reverse, tag, send_data, nsize, recv_data);

    LB_FREE((void **) &(plan_reverse->status));
    LB_FREE((void **) &(plan_reverse->request));
    LB_FREE((void **) &plan_reverse);

    return(lb_flag);
}
