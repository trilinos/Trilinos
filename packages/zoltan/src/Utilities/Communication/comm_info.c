// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include <stdio.h>
#include <memory.h>
#include <mpi.h>
#include "comm.h"
#include "zoltan_mem.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Given the communication plan, return information about the plan. */
/* All array memory must be allocated by the calling routines, and  */
/* is assumed to be of the correct size.                            */
/* Values are not returned when arguments are NULL.                 */


int Zoltan_Comm_Info(
  ZOLTAN_COMM_OBJ *plan,	/* communication data structure */
  int *nsends,                  /* number of sends in plan */
  int *send_procs,              /* list of procs to which messages are sent */
  int *send_lengths,            /* number of values to be sent to each proc
                                   in procs_to */
  int *send_nvals,              /* number of values to send */
  int *send_max_size,           /* max size of msg to be sent. */
  int *send_list,               /* proc assignment of each value to be sent
                                   (same as "assign" in Zoltan_Comm_Create)  */
  int *nrecvs,                  /* number of receives in plan */
  int *recv_procs,              /* list of procs from which messages are 
                                   received. */
  int *recv_lengths,            /* number of values to be received from each
                                   proc in procs_from */
  int *recv_nvals,              /* number of values to recv */
  int *recv_total_size,         /* total size of items to be received. */
  int *recv_list,               /* proc assignment of each value to be received
                                   (inverse of sent_list) */
  int *self_msg                 /* number of self-messages in plan */
)
{
static char *yo = "Zoltan_Comm_Info";
int i, j, k, my_proc;

  /* Check input parameters */
  if (!plan) {
    MPI_Comm_rank(zoltan_get_global_comm(), &my_proc);
    ZOLTAN_COMM_ERROR("Communication plan = NULL", yo, my_proc);
    return ZOLTAN_FATAL;
  }

  if (nsends)
    *nsends = plan->nsends;

  if (send_procs)
    for (i = 0; i < plan->nsends + plan->self_msg; i++)
      send_procs[i] = plan->procs_to[i];

  if (send_lengths)
    for (i = 0; i < plan->nsends + plan->self_msg; i++)
      send_lengths[i] = plan->lengths_to[i];

  if (send_nvals)
    *send_nvals = plan->nvals;

  if (send_max_size)
    *send_max_size = plan->max_send_size;

  if (send_list) {
    for (i = 0; i < plan->nvals; i++) send_list[i] = -1;
    if (plan->indices_to == NULL) {
      for (i = 0; i < plan->nsends + plan->self_msg; i++) {
        k = plan->starts_to[i];
        for (j = 0; j < plan->lengths_to[i]; j++)
          send_list[k+j] = plan->procs_to[i];
      }
    }
    else {
      for (i = 0; i < plan->nsends + plan->self_msg; i++) {
        k = plan->starts_to[i];
        for (j = 0; j < plan->lengths_to[i]; j++) 
          send_list[plan->indices_to[k+j]] = plan->procs_to[i];
      }
    }
  }

  if (nrecvs)
    *nrecvs = plan->nrecvs;

  if (recv_procs)
    for (i = 0; i < plan->nrecvs + plan->self_msg; i++)
      recv_procs[i] = plan->procs_from[i];

  if (recv_lengths)
    for (i = 0; i < plan->nrecvs + plan->self_msg; i++)
      recv_lengths[i] = plan->lengths_from[i];

  if (recv_nvals)
    *recv_nvals = plan->nvals_recv;

  if (recv_total_size)
    *recv_total_size = plan->total_recv_size;

  if (recv_list) {
    for (i = 0; i < plan->nvals_recv; i++) recv_list[i] = -1;
    if (plan->indices_from == NULL) {
      for (i = 0; i < plan->nrecvs + plan->self_msg; i++) {
        k = plan->starts_from[i];
        for (j = 0; j < plan->lengths_from[i]; j++)
          recv_list[k+j] = plan->procs_from[i];
      }
    }
    else {
      for (i = 0; i < plan->nrecvs + plan->self_msg; i++) {
        k = plan->starts_from[i];
        for (j = 0; j < plan->lengths_from[i]; j++) 
          recv_list[plan->indices_from[k+j]] = plan->procs_from[i];
      }
    }
  }

  if (self_msg)
    *self_msg = plan->self_msg;

  return (ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
