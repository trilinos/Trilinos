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


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "dr_const.h"

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
  int         from, to, flag;
  MPI_Status  status;

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

  MPI_Send((void *) &flag, 1, MPI_INT, to, 0, MPI_COMM_WORLD);
  if (proc == 0) {
    from = nprocs - 1;
    MPI_Recv((void *) &flag, 1, MPI_INT, from, 0, MPI_COMM_WORLD, &status);

#ifdef DEBUG_PSYNC
    (void) printf("\t\t\t Proc 0 received message from %5d\n", from);
#endif
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
