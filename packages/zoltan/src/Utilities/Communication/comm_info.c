/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


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
    MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
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
