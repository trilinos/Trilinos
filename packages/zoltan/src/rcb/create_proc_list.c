/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include "lb_const.h"
#include "comm_const.h"
#include "all_allo_const.h"
#include "create_proc_list_const.h"

int LB_Create_Proc_List(
     int       set,           /* set that processor is in */
     int       dotnum,        /* number of dots that my processor has */
     int       outgoing,      /* number of dots that my processor is sending */
     int      *proclist,      /* processor list for my outgoing dots */
     MPI_Comm  comm           /* communicator for partition */
)
{
/* This routine calculates a communication pattern for the situation where
   there are two groups of processors and each processor has a number of
   items which need to be communicated to the other group.  This routine
   calculates a communication pattern which seeks to minimize the number
   of items that any one processor has after communication. */

     int  nprocs;             /* number of processors in partition */
     int  rank;               /* my processor number in partition */
     int *send;               /* array of number of dots outgoing */
     int *rem;                /* array of number of dots that remain */
     int *sets;               /* set for each of the processors */
     int  a;                  /* number of dots that will be on each proc */
     int  sum_send;           /* total number sent from my group */
     int  sum_rem;            /* total number remaining in other group */
     int  np_other = 0;       /* number of processors in other group */
     int  s, sp;              /* temporary sums on number of dots */
     int  num_to;             /* number of dots to send to a processor */
     int  i, j, k;            /* loop indexes */

     /* allocate memory for arrays */
     MPI_Comm_rank(comm, &rank);
     MPI_Comm_size(comm, &nprocs);
     if ((send = (int *) LB_MALLOC(nprocs*sizeof(int))) == NULL)
        return LB_MEMERR;
     if ((rem = (int *) LB_MALLOC(nprocs*sizeof(int))) == NULL) {
        LB_FREE(&send);
        return LB_MEMERR;
     }
     if ((sets = (int *) LB_MALLOC(nprocs*sizeof(int))) == NULL) {
        LB_FREE(&send);
        LB_FREE(&rem);
        return LB_MEMERR;
     }

     /* gather number of outgoing and remaining dots across processors */
     rem[rank] = dotnum - outgoing;
     MPI_Allgather(&rem[rank], 1, MPI_INT, rem, 1, MPI_INT, comm);
     send[rank] = outgoing;
     MPI_Allgather(&send[rank], 1, MPI_INT, send, 1, MPI_INT, comm);
     sets[rank] = set;
     MPI_Allgather(&sets[rank], 1, MPI_INT, sets, 1, MPI_INT, comm);

     /* Convert processor numbers to local (for this communicator) numbering
        and determine which subset of processors a processor is in. */

     /* to determine the number of dots (a) that will be on each of the
        processors in the other group, start with the average */
     for (i = sum_send = sum_rem = 0; i < nprocs; i++)
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
        a = (sum_send + sum_rem)/np_other;
        sp = -1;
        k = 0;
        while (!k) {
           for (i = s = 0; i < nprocs; i++)
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
     for (i = s = 0; i < nprocs; i++)
        if (sets[i] != set && rem[i] < a)
           s += send[i] = a - rem[i];
     while (s < sum_send)
        for (i = 0; i < nprocs && s < sum_send; i++)
           if (sets[i] != set && (send[i] || !s)) {
              send[i]++;
              s++;
           }

     /* Create list of which processor to send dots to.  Only the other half
        of send has been overwritten above.  j is the number of dots that will
        be sent by processors in my part of the partition which are less then
        my processor.  k is the sum of what processors in the other part are
        recieving.  The processor which causes k > j is the first processor
        that we send dots to.  We continue to send dots to processors beyond
        that one until we run out of dots to send. */
     for (i = j = 0; i < rank; i++)
        if (sets[i] == set)                  /* only overwrote other half */
           j += send[i];
     for (i = k = 0; k <= j; i++)
        if (sets[i] != set)
           k += send[i];
     i--;
     num_to = (outgoing < (k - j)) ? outgoing : (k - j);
     for (k = 0; k < num_to; k++)
        proclist[k] = i;
     outgoing -= num_to;
     while (outgoing > 0) {
        i++;
        while (sets[i] == set)
           i++;
        num_to = (outgoing < send[i]) ? outgoing : send[i];
        outgoing -= num_to;
        for (j = 0; j < num_to; j++, k++)
           proclist[k] = i;
     }

     /* free memory and return */
     LB_FREE(&send);
     LB_FREE(&rem);
     LB_FREE(&sets);

     return LB_OK;
}
