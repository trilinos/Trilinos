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
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_parutil_id =
  "$Id$";
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "par_const.h"
#define PRINT_SYNC 5000

void LB_print_sync_start(int do_print_line)
{
/* 
 * Routine to allow IO between print_sync_start and print_sync_end to be printed
 * by each processor entirely before the next processor begins its IO.  The
 * printing sequence is from proc = 0 to the last processor, number_of_procs =
 * LB_Num_Proc - 1.
 *
 * The argument is a boolean variable.  If true, a line of # is printed to
 * indicate the start of a print_sync I/O block.
 *
 * NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.
 *
 * Author: John Shadid (1421, SNL)
 */

int        flag = 1, from, type;
static int offset = 0;
MPI_Status st;

  offset = (offset + 1)%100;
  type   = PRINT_SYNC + offset;

  if ( LB_Proc != 0) {
    from = LB_Proc -1;
    if (MPI_Recv((void *) &flag, 1, MPI_INT, from, type, MPI_COMM_WORLD, &st)
        != 0) {
      fprintf(stderr, "LB_print_sync_start: ERROR on node %d\n", LB_Proc);
      fprintf(stderr, "MPI_Recv failed, message type %d\n", type);
      exit (-1);
    }
  }
  else {
    if (do_print_line) {
      printf("\n");
      for (flag = 0; flag < 37; flag++) printf("#");
      printf(" PRINT_SYNC_START ");
      for (flag = 0; flag < 25; flag++) printf("#");
      printf("\n");
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_print_sync_end(int do_print_line)
{
/*
 * Routine to allow IO between print_sync_start and print_sync_end to be printed
 * by each processor entirely before the next processor begins its IO.  The
 * printing sequence is from proc = 0 to the last processor, number_of_procs =
 * LB_Num_Proc - 1.
 *
 * The argument is a boolean variable.  If true, a line of # is printed to
 * indicate the start of a print_sync I/O block.
 *
 * NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.
 *
 * Author: John Shadid (1421, SNL)
 */

int         flag = 1, from, type, to;
static int  offset = 0;
MPI_Status  st;

  fflush(stdout);

  offset = (offset + 1)%100;
  type   = PRINT_SYNC + offset;

  if (LB_Proc < LB_Num_Proc -1)
    to = LB_Proc + 1;
  else {
    to = 0;
    if (do_print_line) {
      printf("\n");
      for (flag = 0; flag < 37; flag++) printf("#");
      printf(" PRINT_SYNC_END__ ");
      for (flag = 0; flag < 25; flag++) printf("#");
      printf("\n\n");
    }
  }

  if (MPI_Send((void *) &flag, 1, MPI_INT, to, type, MPI_COMM_WORLD) != 0 ) {
    fprintf(stderr, "LB_print_sync_end: ERROR on node %d\n", LB_Proc);
    fprintf(stderr, "MPI_Send failed, message type %d\n", type);
    exit (-1);
  }
  if (LB_Proc == 0) {
    from = LB_Num_Proc -1;
    if (MPI_Recv((void *) &flag, 1, MPI_INT, from, type, MPI_COMM_WORLD, &st)
        != 0) {
      fprintf(stderr, "LB_print_sync_end: ERROR on node %d\n", LB_Proc);
      fprintf(stderr, "MPI_Recv failed, message type %d/n", type);
      exit (-1);
    }
  }

  /*
   * Do a final sync amongst all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from Proc
   * (Num_Proc-1)
   */

  MPI_Barrier(MPI_COMM_WORLD);
}
