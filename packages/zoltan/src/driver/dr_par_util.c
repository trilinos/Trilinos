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

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "dr_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Author(s):	Matthew M. St.John (SNL 9226)
 *		or as noted
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *      print_sync_start()
 *      print_sync_end()
 *      boundary_exchange()
 *+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void print_sync_start(int proc, int do_print_line)

/* Routine to allow IO between print_sync_start and print_sync_end to be
   printed by each processor entirely before the next processor begins its IO.
   The printing sequence is from proc = 0 to the last processor,
   number_of_procs = nprocs - 1.

   The last argument is a boolean variable.  If true, a line of # is printed to
   indicate the start of a print_sync I/O block.

   NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

   Author: John Shadid (1421, SNL)

*/

{
  int        from, flag;
  MPI_Status status;

  MPI_Barrier(MPI_COMM_WORLD);
  if ( proc != 0) {
    from = proc - 1;
    MPI_Recv((void *) &flag, 1, MPI_INT, from, 0, MPI_COMM_WORLD, &status);
  }
  else {
    if (do_print_line) {
      (void) printf("\n");
      for (flag = 0; flag < 36; flag++) (void) printf("#");
      (void) printf(" PRINT_SYNC_START ");
      for (flag = 0; flag < 25; flag++) (void) printf("#");
      (void) printf("\n");
    }
  }

#ifdef DEBUG_PSYNC
  (void) printf("\t\tSTART OF PRINT_SYNC SECTION, Proc = %4d\n", proc);
#endif
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_sync_end(int proc, int nprocs, int do_print_line)

/*

  Routine to allow IO between print_sync_start and print_sync_end to be printed
  by each processor entirely before the next processor begins its IO.  The
  printing sequence is from proc = 0 to the last processor, number_of_procs =
  nprocs - 1.

  The last argument is a boolean variable.  If true, a line of # is printed to
  indicate the start of a print_sync I/O block.

  NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

  Author: John Shadid (1421, SNL)

*/

{
  int         from, to, flag=0;
  MPI_Status  status;
  MPI_Request req;

#ifdef DEBUG_PSYNC
  (void) printf("\t\tEND OF PRINT_SYNC SECTION, Proc = %4d\n", proc);
#endif
  if (proc < nprocs -1)
    to = proc + 1;
  else {
    to = 0;
    if (do_print_line) {
      (void) printf("\n");
      for (flag = 0; flag < 36; flag++) (void) printf("#");
      (void) printf(" PRINT_SYNC_END__ ");
      for (flag = 0; flag < 25; flag++) (void) printf("#");
      (void) printf("\n\n");
    }
  }

  if (proc == 0) {
    from = nprocs - 1;
    MPI_Irecv((void *) &flag, 1, MPI_INT, from, 0, MPI_COMM_WORLD, &req);

#ifdef DEBUG_PSYNC
    (void) printf("\t\t\t Proc 0 posted receive from %5d\n", from);
#endif
  }

  MPI_Send((void *) &flag, 1, MPI_INT, to, 0, MPI_COMM_WORLD);

  if (proc == 0) {
    MPI_Wait(&req, &status);
  }

  /*
   * Do a final sync amongst all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from Proc
   * (Num_Proc-1)
   */

  MPI_Barrier(MPI_COMM_WORLD);

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void boundary_exchange(
  MESH_INFO_PTR mesh,
  int vec_len,           /* Length of vector for each element.               */
  int *send_vec,         /* Vector of values to be sent.                     */
  int *recv_vec          /* Vector of values to be received.                 */
)
{
int i, offset;
int msg_type = 111;

MPI_Status  *status = NULL;
MPI_Request *req = NULL;

  req = (MPI_Request *) malloc(mesh->necmap * sizeof(MPI_Request));
  status = (MPI_Status *) malloc(mesh->necmap * sizeof(MPI_Status));

  /* Post receives */
  offset = 0;
  for (i = 0; i < mesh->necmap; i++) {
    MPI_Irecv(&(recv_vec[offset]), mesh->ecmap_cnt[i]*vec_len, MPI_INT, 
              mesh->ecmap_id[i], msg_type, MPI_COMM_WORLD, &(req[i]));
    offset += mesh->ecmap_cnt[i]*vec_len;
  }

  /* Send messages */
  offset = 0;
  for (i = 0; i < mesh->necmap; i++) {
    MPI_Send(&(send_vec[offset]), mesh->ecmap_cnt[i]*vec_len, MPI_INT, 
             mesh->ecmap_id[i], msg_type, MPI_COMM_WORLD);
    offset += mesh->ecmap_cnt[i]*vec_len;
  }

  /* Receive messages */
  MPI_Waitall(mesh->necmap, req, status);

  free(req);
  free(status);
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
