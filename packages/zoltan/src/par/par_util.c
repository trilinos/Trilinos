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
#include "lb_const.h"
#include "par_const.h"
#define PRINT_SYNC 5000

void LB_print_sync_start(LB *lb, int do_print_line)
{
/* 
 * Routine to allow I/O between print_sync_start and print_sync_end to be 
 * printed by each processor in the lb->Communicator entirely before the next
 * processor begins its I/O.  The printing sequence is from proc = 0 to the
 * last processor, where the last processor is lb->Num_Proc - 1.
 *
 * The do_print_line argument is a boolean variable.  If true, a line of # 
 * is printed to indicate the start of a print_sync I/O block.
 *
 * NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.
 *
 * Author: John Shadid (9221, SNL)
 */

int        flag = 1, from, type;
static int offset = 0;
MPI_Status st;
char *yo = "LB_print_sync_start";

  offset = (offset + 1)%100;
  type   = PRINT_SYNC + offset;

  if (lb->Proc != 0) {
    from = lb->Proc -1;
    if (MPI_Recv((void *) &flag, 1, MPI_INT, from, type, lb->Communicator, &st)
        != 0) {
      fprintf(stderr, "%s: ERROR on processor %d\n", yo, lb->Proc);
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

void LB_print_sync_end(LB *lb, int do_print_line)
{
/*
 * Routine to allow I/O between print_sync_start and print_sync_end to be 
 * printed by each processor in the lb->Communicator entirely before the next
 * processor begins its I/O.  The printing sequence is from proc = 0 to the
 * last processor, where the last processor is lb->Num_Proc - 1.
 *
 * The do_print_line argument is a boolean variable.  If true, a line of # 
 * is printed to indicate the start of a print_sync I/O block.
 *
 * NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.
 *
 * Author: John Shadid (9221, SNL)
 */

int         flag = 1, from, type, to;
static int  offset = 0;
MPI_Status  st;
char *yo = "LB_print_sync_end";

  fflush(stdout);

  offset = (offset + 1)%100;
  type   = PRINT_SYNC + offset;

  if (lb->Proc < lb->Num_Proc -1)
    to = lb->Proc + 1;
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

  if (MPI_Send((void *) &flag, 1, MPI_INT, to, type, lb->Communicator) != 0 ) {
    fprintf(stderr, "%s: ERROR on node %d\n", yo, lb->Proc);
    fprintf(stderr, "MPI_Send failed, message type %d\n", type);
    exit (-1);
  }
  if (lb->Proc == 0) {
    from = lb->Num_Proc -1;
    if (MPI_Recv((void *) &flag, 1, MPI_INT, from, type, lb->Communicator, &st)
        != 0) {
      fprintf(stderr, "%s: ERROR on node %d\n", yo, lb->Proc);
      fprintf(stderr, "MPI_Recv failed, message type %d/n", type);
      exit (-1);
    }
  }

  /*
   * Do a final sync among all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from Proc
   * (lb->Num_Proc-1)
   */

  MPI_Barrier(lb->Communicator);
}
