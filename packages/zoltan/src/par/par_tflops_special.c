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
 *    Revision: 1.6.2.1 $
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "par_tflops_special_const.h"
#include "zoltan_mem.h"

/* Generic Tflops_special routines to avoid certain types of 
   collective communication routines. */

void Zoltan_RB_scan_double(
   double *wtok,          /* in: local weight */
   double *wtupto,        /* out: sum of weights for processors <= rank */
   int count,             /* number of elements to scan */
   MPI_Comm local_comm,   /* MPI Communicator */
   int proc,              /* global processor number */
   int rank,              /* rank in this partition */
   int num_procs          /* number of processors in this partition */
)
{
   int to;                /* communication partner (global) */
   int tor;               /* rank of partner in this partition */
   int nprocs_large;      /* power of 2 processor that contains num_procs */
   int hbit;              /* 2^hbit = nprocs_large */
   int mask;              /* mask to determine communication partner */
   int i;                 /* loop index */
   int tag = 32108;       /* message tag */
   double *tmpin;         /* temporary double array to recieve */
   double *tmpout;        /* temporary double array to send */
   double temp[2];        /* temp buffer to avoid malloc when count==1 */
   MPI_Status status;

   /* this subroutine performs a scan operation summing doubles on a subset
      of a set of processors for Tflops_Special */

   if (count==0)
     return;
   else if (count==1){ /* avoid malloc */
     tmpin = temp;
     tmpout = &temp[1];
   }
   else /* count>1 */ {
     tmpin = (double *) ZOLTAN_MALLOC(2*count*sizeof(double));
     tmpout = tmpin + count;
   }

   /* Find next lower power of 2. */
   for (hbit = 0; (num_procs >> hbit) != 0; hbit++);

   nprocs_large = 1 << hbit;
   if (nprocs_large == 2*num_procs) nprocs_large = num_procs;

   for (i=0; i<count; i++)
      tmpout[i] = wtupto[i] = wtok[i];
   for (mask = 1; mask <= nprocs_large; mask *= 2) {
      tag++;
      tor = (rank ^ mask);
      to = proc - rank + tor;
      if (tor < num_procs) {
         MPI_Send(tmpout, count, MPI_DOUBLE, to, tag, local_comm);
         MPI_Recv(tmpin, count, MPI_DOUBLE, to, tag, local_comm, &status);
         for (i=0; i<count; i++){
            tmpout[i] += tmpin[i];
            if (to < proc) 
               wtupto[i] += tmpin[i];
         }
      }
   }

   if (count>1)
      ZOLTAN_FREE(&tmpin);
}

void Zoltan_RB_sum_double(
   double   *x,               /* doubles to be summed (answer overwrites!) */
   int      count,            /* number of elements in x */
   int      proclower,        /* smallest processor in partition */
   int      rank,             /* rank of processor in partition */
   int      nprocs,           /* number of processors in partition */
   MPI_Comm comm
)
{
   double   temp, *tmp;       /* temporary for sum */
   int      tag = 32100;      /* message tag */
   int      partner;          /* message partner in binary exchange */
   int      to;               /* message partner not in binary exchange */
   int      mask;             /* mask to determine communication partner */
   int      nprocs_small;     /* largest power of 2 contained in nprocs */
   int      hbit;             /* 2^hbit = nproc_small */
   int      i;                /* loop index */
   MPI_Status status;

   /* This routine sums doubles on a subset of processors */
 
   if (count==0)
     return;
   else if (count==1){ /* avoid malloc */
     tmp = &temp;
   }
   else /* count>1 */ {
     tmp = (double *) ZOLTAN_MALLOC(count*sizeof(double));
   }

   /* Find next lower power of 2. */
   for (hbit = 0; (nprocs >> hbit) != 1; hbit++);
 
   nprocs_small = 1 << hbit;
   if (nprocs_small * 2 == nprocs) {
      nprocs_small *= 2;
      hbit++;
   }
 
   to = proclower + (rank ^ nprocs_small);
   if (rank & nprocs_small) {  /* processors greater than largest power of 2 */
      MPI_Send(x, count, MPI_DOUBLE, to, tag, comm);
      tag += hbit + 1;
      MPI_Recv(x, count, MPI_DOUBLE, to, tag, comm, &status);
   }
   else {   /* processors within greatest power of 2 */
      if (rank + nprocs_small < nprocs) {
         MPI_Recv(tmp, count, MPI_DOUBLE, to, tag, comm, &status);
         for (i=0; i<count; i++)
            x[i] += tmp[i];
      }  
      for (mask = nprocs_small >> 1; mask; mask >>= 1) { /* binary exchange */
         tag++;
         partner = proclower + (rank ^ mask);
         MPI_Send(x, count, MPI_DOUBLE, partner, tag, comm);
         MPI_Recv(tmp, count, MPI_DOUBLE, partner, tag, comm, &status);
         for (i=0; i<count; i++)
            x[i] += tmp[i];
      }  
      tag++;
      if (rank + nprocs_small < nprocs)
         MPI_Send(x, count, MPI_DOUBLE, to, tag, comm);
   }
   if (count>1)
      ZOLTAN_FREE(&tmp);
}

void Zoltan_RB_max_double(
   double   *x,               /* maximum value (answer overwrites!) */
   int      count,            /* number of elements in x */
   int      proclower,        /* smallest processor in partition */
   int      rank,             /* rank of processor in partition */
   int      nprocs,           /* number of processors in partition */
   MPI_Comm comm
)
{
   double   temp, *tmp;       /* temporaries for min/max */
   int      tag = 32100;      /* message tag */
   int      partner;          /* message partner in binary exchange */
   int      to;               /* message partner not in binary exchange */
   int      mask;             /* mask to determine communication partner */
   int      nprocs_small;     /* largest power of 2 contained in nprocs */
   int      hbit;             /* 2^hbit = nproc_small */
   int      i;                /* loop index */
   MPI_Status status;

   /* This routine finds the max on a set of procs */

   if (count==0)
     return;
   else if (count==1){ /* avoid malloc */
     tmp = &temp;
   }
   else /* count>1 */ {
     tmp = (double *) ZOLTAN_MALLOC(count*sizeof(double));
   }

   /* Find next lower power of 2. */
   for (hbit = 0; (nprocs >> hbit) != 1; hbit++);
 
   nprocs_small = 1 << hbit;
   if (nprocs_small * 2 == nprocs) {
      nprocs_small *= 2;
      hbit++;
   }
 
   to = proclower + (rank ^ nprocs_small);
   if (rank & nprocs_small) {  /* processors greater than largest power of 2 */
      MPI_Send(x, count, MPI_DOUBLE, to, tag, comm);
      tag += hbit + 1;
      MPI_Recv(x, count, MPI_DOUBLE, to, tag, comm, &status);
   }
   else {   /* processors within greatest power of 2 */
      if (rank + nprocs_small < nprocs) {
         MPI_Recv(tmp, count, MPI_DOUBLE, to, tag, comm, &status);
         for (i=0; i<count; i++)
            if (tmp[i] > x[i]) x[i] = tmp[i];
      }
      for (mask = nprocs_small >> 1; mask; mask >>= 1) { /* binary exchange */
         tag++;
         partner = proclower + (rank ^ mask);
         MPI_Send(x, count, MPI_DOUBLE, partner, tag, comm);
         MPI_Recv(tmp, count, MPI_DOUBLE, partner, tag, comm, &status);
         for (i=0; i<count; i++)
            if (tmp[i] > x[i]) x[i] = tmp[i];
      }
      tag++;
      if (rank + nprocs_small < nprocs)
         MPI_Send(x, count, MPI_DOUBLE, to, tag, comm);
   }
   if (count>1)
      ZOLTAN_FREE(&tmp);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
