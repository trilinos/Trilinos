/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/
#ifndef lint
static char *cvs_dr_par_util = "$Id$";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Author(s):	Matthew M. St.John (SNL 9226)
 *		or as noted
 *----------------------------------------------------------------------------
 * Functions contained in this file:
 *      print_sync_start()
 *      print_sync_end()
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
    MPI_Recv(NULL, 0, MPI_INT, from, 0, MPI_COMM_WORLD, &status);
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

  MPI_Send(NULL, 0, MPI_INT, to, 0, MPI_COMM_WORLD);
  if (proc == 0) {
    from = nprocs - 1;
    MPI_Recv(NULL, 0, MPI_INT, from, 0, MPI_COMM_WORLD, &status);

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

