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

/* This routine divides the current machine (defined by the communicator)
   into two pieces.
   For now, it simply divides the machine in half.  In the future, it will
   be a more complicated routine taking into account the architecture of
   the machine and communication network. */

#include "lb_const.h"

int LB_divide_machine(
   LB *lb,             /* The load-balancing structure (not used now, will be
                          used for pointer to machine details */
   int proc,           /* my processor number in global sense */
   MPI_Comm comm,      /* communicator for part of machine to be divided */
   int *set,           /* part that proc is in after divide (lowest global
                          numbered processor in set 0) */
   int *proclower,     /* lowest numbered processor in first part */
   int *procmid,       /* lowest numbered processor in second part */
   int *num_procs,     /* on input, number of procs in the part to be divided
                          on exit, number of procs in part that proc is in */
   double *fractionlo  /* actual division of machine */
)
{
   int ierr;           /* error flag */

   /* assume for now that processors coming in are in a contiguously
      numbered set and we will divide them into two roughly equal
      contiguously numbered sets */

   ierr = MPI_Allreduce(&proc, proclower, 1, MPI_INT, MPI_MIN, comm);

   *procmid = *proclower + (*num_procs - 1)/2 + 1;

   *fractionlo = ((double) (*procmid - *proclower))/(*num_procs);

   if (proc < *procmid) {
      *set = 0;
      *num_procs = *procmid - *proclower;
   } else {
      *set = 1;
      *num_procs = *num_procs - *procmid + *proclower;
   }

   return LB_OK;
}
