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

#include "zz_const.h"


int Zoltan_Divide_Machine(
   ZZ *zz,             /* The Zoltan structure (not used now, will be
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
/* This routine divides the current machine (defined by the communicator)
   into two pieces.
   For now, it simply divides the machine in half.  In the future, it will
   be a more complicated routine taking into account the architecture of
   the machine and communication network. */
/* If using Tflops_Special, then we have to divide the set into two
   contiguously numbered sets. */

   /* assume for now that processors coming in are in a contiguously
      numbered set and we will divide them into two roughly equal
      contiguously numbered sets */

   /* The following statement assumes that proclower is being set correctly in
      the calling routine if Tflops_Special flag is set */
   if (!zz->Tflops_Special)
      MPI_Allreduce(&proc, proclower, 1, MPI_INT, MPI_MIN, comm);

   *procmid = *proclower + (*num_procs - 1)/2 + 1;

   *fractionlo = ((double) (*procmid - *proclower))/(*num_procs);

   if (proc < *procmid) {
      *set = 0;
      *num_procs = *procmid - *proclower;
   } else {
      *set = 1;
      *num_procs = *num_procs - *procmid + *proclower;
   }

   return ZOLTAN_OK;
}
