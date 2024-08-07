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


#include "zz_util_const.h"
#include "create_proc_list_const.h"

static void Zoltan_RB_Gather(ZOLTAN_GNO_TYPE *, ZOLTAN_GNO_TYPE *, int, int, int, MPI_Comm);

int Zoltan_RB_Create_Proc_List(
     ZZ       *zz,            /* Load-balancing structure. */
     int       set,           /* set that processor is in */
     int       dotnum,        /* number of dots that my processor has */
     int       outgoing,      /* number of dots that my processor is sending */
     int      *proclist,      /* processor list for my outgoing dots */
     MPI_Comm  comm,          /* communicator for partition */
     int       proclower,     /* smallest processor for Tflops_Special */
     int       numprocs       /* number of processors for Tflops_Special */
)
{
/* This routine calculates a communication pattern for the situation where
   there are two groups of processors and each processor has a number of
   items which need to be communicated to the other group.  This routine
   calculates a communication pattern which seeks to minimize the number
   of items that any one processor has after communication. */

     int  nprocs;             /* number of processors in partition */
     int  rank;               /* my processor number in partition */
     int  np_other = 0;       /* number of processors in other group */
     int  i, k;            /* loop indexes */
     int  err = ZOLTAN_OK;    /* error code */
     MPI_Datatype zoltan_gno_mpi_type;

     ZOLTAN_GNO_TYPE *send;         /* array of number of dots outgoing */
     ZOLTAN_GNO_TYPE *rem;          /* array of number of dots that remain */
     ZOLTAN_GNO_TYPE *sets;         /* set for each of the processors */
     ZOLTAN_GNO_TYPE  a;            /* number of dots that will be on each proc */
     ZOLTAN_GNO_TYPE  sum_send;     /* total number sent from my group */
     ZOLTAN_GNO_TYPE  sum_rem;      /* total number remaining in other group */
     ZOLTAN_GNO_TYPE  s, sp;        /* temporary sums */
     ZOLTAN_GNO_TYPE  num_to;       /* number of dots to send to a processor */
     ZOLTAN_GNO_TYPE *tmp_send;     /* Work vector */

     zoltan_gno_mpi_type = Zoltan_mpi_gno_type();


     /* allocate memory for arrays */
     MPI_Comm_rank(comm, &rank);
     if (zz->Tflops_Special) {
        nprocs = numprocs;
        rank -= proclower;
     }
     else
        MPI_Comm_size(comm, &nprocs);

     if ((send = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(3*nprocs*sizeof(ZOLTAN_GNO_TYPE))) == NULL)
        return ZOLTAN_MEMERR;
     if ((tmp_send = (ZOLTAN_GNO_TYPE *) ZOLTAN_MALLOC(3*nprocs*sizeof(ZOLTAN_GNO_TYPE))) == NULL) {
        ZOLTAN_FREE(&send);
        return ZOLTAN_MEMERR;
     }
     rem = &send[nprocs];
     sets = &send[2*nprocs];

     /* gather number of outgoing and remaining dots across processors */

     if (zz->Tflops_Special) {
        for (i = 3*nprocs-1; i >= 0; i--)
           send[i] = 0;
        send[rank] = outgoing;
        send[rank+nprocs] = dotnum - outgoing;
        send[rank+2*nprocs] = set;
        Zoltan_RB_Gather(send, tmp_send, proclower, rank, nprocs, comm);
     }
     else {
        for (i = 3*nprocs-1; i >= 0; i--)
           tmp_send[i] = 0;
        tmp_send[rank] = outgoing;
        tmp_send[rank+nprocs] = dotnum - outgoing;
        tmp_send[rank+2*nprocs] = set;
        i = 3*nprocs;
        MPI_Allreduce(tmp_send, send, i, zoltan_gno_mpi_type, MPI_SUM, comm);
     }

     ZOLTAN_FREE(&tmp_send);

     /* Convert processor numbers to local (for this communicator) numbering
        and determine which subset of processors a processor is in. */

     /* to determine the number of dots (a) that will be on each of the
        processors in the other group, start with the average */
     sum_send = sum_rem = 0;
     for (i = 0; i < nprocs; i++)
        if (sets[i] == set)
           sum_send += send[i];
        else {
           sum_rem += rem[i];
           np_other++;
        }

     /* Modify the value of a, which is now the number of dots that a processor
        will have after communication if the number of dots that it keeps is
        not greater than a.  Set a so there is a non-negative number left. */
     if (sum_send) {
        if (np_other == 0){
          /* This should never happen. Prevent divide-by-zero to be safe. */
          err = ZOLTAN_FATAL;
          goto End;
        }
        a = (sum_send + sum_rem)/np_other;
        sp = -1;
        k = 0;
        while (!k) {
           s = 0;
           for (i = 0; i < nprocs; i++)
              if (sets[i] != set && rem[i] < a)
                 s += a - rem[i];
           if (s == sum_send)
              k = 1;
           else if (s < sum_send) {
              if (sp > sum_send)
                 k = 1;
              else
                 a++;
           }
           else {
              a--;
              if (sp < sum_send && sp > 0)
                 k = 1;
           }
           sp = s;
        }
     } else
        a = 0;

     /* Allocate who recieves how much and if necessary give out remainder.
        The variable send is now the number that will be received by each
        processor only in the other part */
     s = 0;
     for (i = 0; i < nprocs; i++)
        if (sets[i] != set && rem[i] < a)
           s += send[i] = a - rem[i];
     while (s < sum_send)
        for (i = 0; i < nprocs && s < sum_send; i++)
           if (sets[i] != set && (send[i] || !s)) {
              send[i]++;
              s++;
           }

     /* Create list of which processor to send dots to.  Only the other half
        of send has been overwritten above.  s is the number of dots that will
        be sent by processors in my part of the partition which are less then
        my processor.  sp is the sum of what processors in the other part are
        recieving.  The processor which causes sp > s is the first processor
        that we send dots to.  We continue to send dots to processors beyond
        that one until we run out of dots to send. */
     if (outgoing) {
        if (zz->Tflops_Special) a = (ZOLTAN_GNO_TYPE)outgoing;    /* keep outgoing around in a */
        s = sp = 0;
        for (i = 0; i < rank; i++)
           if (sets[i] == set)                   /* only overwrote other half */
              s += send[i];
        for (i = 0; sp <= s; i++)
           if (sets[i] != set)
              sp += send[i];
        i--;
        num_to = (outgoing < (sp - s)) ? outgoing : (sp - s);
        for (sp = 0; sp < num_to; sp++)
           proclist[sp] = i;
        outgoing -= num_to;
        while (outgoing > 0) {
           i++;
           while (sets[i] == set)
              i++;
           num_to = (outgoing < send[i]) ? outgoing : send[i];
           outgoing -= num_to;
           for (s = 0; s < num_to; s++, sp++)
              proclist[sp] = i;
        }
        if (zz->Tflops_Special)
           for (s = 0; s < a; s++)
              proclist[s] += proclower;
     }

End:
     /* free memory and return */
     ZOLTAN_FREE(&send);


     return err;
}

static void Zoltan_RB_Gather(
   ZOLTAN_GNO_TYPE *send,                 /* input/output array */
   ZOLTAN_GNO_TYPE *tmp_send,             /* temporary array */
   int proclower,             /* smallest numbered processor in partition */
   int rank,                  /* processor number within partition */
   int nprocs,                /* number of processors in this partition */
   MPI_Comm comm              /* MPI Communicator */
)
{
   int tag = 32100;           /* message tag */
   int partner;               /* message partner in binary exchange */
   int to;                    /* message partner not in binary exchange */
   int mask;                  /* mask to determine communication partner */
   int nprocs_small;          /* largest power of 2 contained in nprocs */
   int hbit;                  /* 2^hbit = nproc_small */
   int len;                   /* message length */
   int i;                     /* loop counter */
   MPI_Status status;
   MPI_Datatype zoltan_gno_mpi_type;

   zoltan_gno_mpi_type = Zoltan_mpi_gno_type();

   /* This routine sums a vector of integers on a subset of processors */

   len = 3*nprocs;

   /* Find next lower power of 2. */
   for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

   nprocs_small = 1 << hbit;
   if (nprocs_small * 2 == nprocs) {
      nprocs_small *= 2;
      hbit++;
   }

   to = proclower + (rank ^ nprocs_small);
   if (rank & nprocs_small) {  /* processors greater than largest power of 2 */
      MPI_Send(send, len, zoltan_gno_mpi_type, to, tag, comm);
      tag += hbit + 1;
      MPI_Recv(send, len, zoltan_gno_mpi_type, to, tag, comm, &status);
   }
   else {   /* processors within greatest power of 2 */
      if (rank + nprocs_small < nprocs) {
         MPI_Recv(tmp_send, len, zoltan_gno_mpi_type, to, tag, comm, &status);
         for (i = 0; i < len; i++)
            send[i] += tmp_send[i];
      }
      for (mask = nprocs_small >> 1; mask; mask >>= 1) { /* binary exchange */
         tag++;
         partner = proclower + (rank ^ mask);
         /* Change requested by Qingyu Meng <qymeng@cs.utah.edu> to        */
         /* support mvapich 1.0 on TACC Ranger.                            */
         /* MPI_Send(send, len, zoltan_gno_mpi_type, partner, tag, comm);              */
         /* MPI_Recv(tmp_send, len, zoltan_gno_mpi_type, partner, tag, comm, &status); */
         MPI_Sendrecv(send, len, zoltan_gno_mpi_type, partner, tag,
            tmp_send, len, zoltan_gno_mpi_type, partner, tag, comm, &status);

         for (i = 0; i < len; i++)
            send[i] += tmp_send[i];
      }
      tag++;
      if (rank + nprocs_small < nprocs)
         MPI_Send(send, len, zoltan_gno_mpi_type, to, tag, comm);
   }
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
