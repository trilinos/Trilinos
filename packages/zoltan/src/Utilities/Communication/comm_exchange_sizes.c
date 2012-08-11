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

int       Zoltan_Comm_Exchange_Sizes(
int      *sizes_to,		/* value I need to exchange (size of true msg) */
int      *procs_to,		/* procs I send to */
int       nsends,		/* number of messages I'll send */
int       self_msg,		/* do I copy data to myself? */
int      *sizes_from,		/* (returned) size of all my receives */
int      *procs_from,		/* procs I recv from */
int       nrecvs,		/* number of messages I receive */
int      *total_recv_size,	/* (returned) sum of all incoming sizes */
int       my_proc,		/* my processor number */
int       tag,			/* message tag I can use */
MPI_Comm  comm) {		/* communicator */

    int       self_index_to;	/* location of self in procs_to */
    MPI_Status status;		/* status of commuication operation */
    int       i;		/* loop counter */


    /* If sizes vary, then I need to communicate message lengths */

    self_index_to = -1;
    for (i = 0; i < nsends + self_msg; i++) {
	if (procs_to[i] != my_proc) 
	    MPI_Send((void *) &sizes_to[i], 1, MPI_INT, procs_to[i], tag, comm);
	else
	    self_index_to = i;
    }

    *total_recv_size = 0;
    for (i = 0; i < nrecvs + self_msg; i++) {
	if (procs_from[i] != my_proc) 
	    MPI_Recv((void *) &sizes_from[i], 1, MPI_INT, procs_from[i],
		     tag, comm, &status);
	else 
	    sizes_from[i] = sizes_to[self_index_to];

	*total_recv_size += sizes_from[i];
    }

    return (ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
