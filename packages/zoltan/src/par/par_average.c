// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "par_median_const.h"
#include "zz_const.h"

#define TINY   1.0e-6

/* prototypes for TFLOPS_SPECIAL */
static void Zoltan_average_cuts_reduce(int, int, int, double *, double *,
   int, MPI_Datatype, MPI_Comm);


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double Zoltan_RB_Average_Cut(
  int Tflops_Special,   /* Flag indicating whether Tflops_Special handling
                           of communicators should be done (to avoid memory
                            leaks in tflops' Comm_Dup and Comm_Split).        */
  double *dots,         /* array of coordinates                              */
  int *dotmark,         /* list of which side of the median each dot is on:
                                0 - dot is in lo set
                                1 - dot is in hi set                         */
  int dotnum,           /* number of dots (length of two previous arrays     */
  int num_procs,        /* Number of procs in set (used with Tflops_Special) */
  int rank,             /* rank in partition (used with Tflops_Special)      */
  int proc,             /* this proc number (used with Tflops_Special)       */
  MPI_Comm local_comm,  /* MPI communicator on which to find median          */
  double oldvaluehalf   /* Cut computed before averaging is done.            */
)
{
/* Compute a median value that is exactly between two closest dots. 
 * Routine is called when parameter AVERAGE_CUTS == 1.
 */
double val[2] = {-DBL_MAX, DBL_MAX};
double gval[2];
double valuehalf = oldvaluehalf;
int i;

  if (!Tflops_Special || num_procs > 1) { 
    /* Don't include dot info if going thru loop only due to Tflops_Special */
    for (i = 0; i < dotnum; i++) {
/*
printf("KDDDDD %d proclower=%d num_parts=%d numlist=%d i=%d dotmark[i]=%d dots[i]=%e\n", proc, proclower, num_parts, numlist, i, dotmark[i], dots[i]);
*/
      if (dotmark[i] == 0) {            /* in lower part */
        if (dots[i] > val[0]) val[0] = dots[i];
      } 
      else {   /* in upper part */
        if (dots[i] < val[1]) val[1] = dots[i];
      }
    }
    if (!Tflops_Special) {
      MPI_Allreduce(&val[0], &gval[0], 1, MPI_DOUBLE, MPI_MAX, local_comm);
      MPI_Allreduce(&val[1], &gval[1], 1, MPI_DOUBLE, MPI_MIN, local_comm);
    }
    else
      Zoltan_average_cuts_reduce(num_procs, rank, proc, val, gval, 2, 
                                 MPI_DOUBLE, local_comm);

    valuehalf = 0.5 * (gval[0] + gval[1]);
  }
/*
printf("KDDKDD %d num_procs=%d tmp_half=%e valuehalf=%e val=(%e,%e) gval=(%e,%e)\n", proc, num_procs, tmp_half, *valuehalf, val[0], val[1], gval[0], gval[1]);
*/
  return valuehalf;
}

static void Zoltan_average_cuts_reduce(
   int nproc,             /* number of processors in partition */
   int rank,              /* rank within partition */
   int proc,              /* global processor number */
   double *in,            /* input median */
   double *inout,         /* output median */
   int len,               /* # of double in in */
   MPI_Datatype datatype, /* MPI datatype for in */
   MPI_Comm comm          /* MPI communicator */
)
{
   double tmp[2];         /* temporary to receive information */
   int to;                /* communication partner */
   int tag = 32005;       /* message tag */
   int nprocs_small;      /* largest power of 2 contained in nproc */
   int hbit;              /* 2^hbit = nproc_small */
   int mask;              /* mask to determine communication partner */
   MPI_Status status;

   /* find largest power of two that is less than or equal to number of
      processors */
   for (hbit = 0; (nproc >> hbit) != 1; hbit++);

   nprocs_small = 1 << hbit;
   if (nprocs_small * 2 == nproc) {
      nprocs_small *= 2;
      hbit++;
   }

   /* get input from the processors that are larger than nprocs_small in
      the local partition of processors */
   to = proc - rank + (rank ^ nprocs_small);
   if (rank & nprocs_small)
      MPI_Send(in, len, datatype, to, tag, comm);
   else
      if (rank + nprocs_small < nproc) {
         MPI_Recv(inout, len, datatype, to, tag, comm, &status);
         if (in[0] > inout[0]) inout[0] = in[0];
         if (in[1] < inout[1]) inout[1] = in[1];
      }
      else {
         inout[0] = in[0];
         inout[1] = in[1];
      }

   if (!(rank & nprocs_small))    /* binary exchange on nprocs_small procs */
      for (mask = nprocs_small >> 1; mask; mask >>= 1) {
         tag++;
         to = proc - rank + (rank ^ mask);
         MPI_Send(inout, len, datatype, to, tag, comm);
         MPI_Recv(tmp, len, datatype, to, tag, comm, &status);
         if (tmp[0] > inout[0]) inout[0] = tmp[0];
         if (tmp[1] < inout[1]) inout[1] = tmp[1];
      }
   else
      tag += hbit;

   /* send results to the processors that are larger than nprocs_small in
      the local partition of processors */
   tag++;
   to = proc - rank + (rank ^ nprocs_small);
   if (rank & nprocs_small)
      MPI_Recv(inout, len, datatype, to, tag, comm, &status);
   else
      if (rank + nprocs_small < nproc)
         MPI_Send(inout, len, datatype, to, tag, comm);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
