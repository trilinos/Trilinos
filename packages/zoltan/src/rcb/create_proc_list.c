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

/* This routine calculates a communication pattern for the situation where
   there are two groups of processors and each processor has a number of
   items which need to be communicated to the other group.  This routine
   calculates a communication pattern which seeks to minimize the number
   of items that any one processor has after communication. */

#include "lb_const.h"
#include "comm_const.h"
#include "all_allo_const.h"

int LB_Create_Proc_List(
     int       proc,          /* my processor number */
     int       procmid,       /* first processor in upper part of partition */
     int       proclower,     /* first processor in partition */
     int       procupper,     /* last processor in partition */
     int       dotnum,        /* number of dots that my processor has */
     int       outgoing,      /* number of dots that my processor is sending */
     int      *proclist,      /* processor list for my outgoing dots */
     MPI_Comm  comm           /* communicator for partition */
)
{
     int  nprocs;             /* number of processors in partition */
     int  rank;               /* my processor number in partition */
     int *send;               /* array of number of dots outgoing */
     int *rem;                /* array of number of dots that remain */
     int  mb;                 /* beginning of part of the partition my
                                 processor is in */
     int  me;                 /* end of that part of the partition */
     int  ob;                 /* beginning of other part of the partition */
     int  oe;                 /* end of other part of the partition */
     int  op;                 /* number of other part of the partition */
     int  a;                  /* number of dots that will be on each proc */
     int  sum_send;           /* total number sent from my group */
     int  sum_rem;            /* total number remaining in other group */
     int  s, sp;              /* temporary sums on number of dots */
     int  num_to;             /* number of dots to send to a processor */
     int  i, j, k;            /* loop indexes */

     /* allocate memory for arrays */
     nprocs = procupper - proclower + 1;
     if ((send = (int *) LB_MALLOC(nprocs*sizeof(int))) == NULL)
        return LB_MEMERR;
     if ((rem = (int *) LB_MALLOC(nprocs*sizeof(int))) == NULL) {
        LB_FREE(&send);
        return LB_MEMERR;
     }

     /* gather number of outgoing and remaining dots across processors */
     rank = proc - proclower;
     rem[rank] = dotnum - outgoing;
     MPI_Allgather(&rem[rank], 1, MPI_INT, rem, 1, MPI_INT, comm);
     send[rank] = outgoing;
     MPI_Allgather(&send[rank], 1, MPI_INT, send, 1, MPI_INT, comm);

     /* Convert processor numbers to local (for this communicator) numbering
        and determine which subset of processors a processor is in. */
     if (proc < procmid) {
        mb = 0;
        me = procmid - proclower - 1;
        ob = procmid - proclower;
        oe = procupper - proclower;
        op = procupper - procmid + 1;
     } else {
        ob = 0;
        oe = procmid - proclower - 1;
        op = procmid - proclower;
        mb = procmid - proclower;
        me = procupper - proclower;
     }

     /* to determine the number of dots (a) that will be on each of the
        processors in the other group, start with the average */
     for (i = mb, sum_send = 0; i <= me; i++)
        sum_send += send[i];
     for (i = ob, sum_rem = 0; i <= oe; i++)
        sum_rem += rem[i];

     /* Modify the value of a, which is now the number of dots that a processor
        will have after communication if the number of dots that it keeps is
        not greater than a.  Set a so there is a non-negative number left. */
     if (sum_send) {
        a = (sum_send + sum_rem)/op;
        sp = -1;
        s = k = 0;
        while (!k) {
           for (i = ob; i <= oe; i++)
              if (rem[i] < a)
                 s += a - rem[i];
           if (s == sum_send)
              k = 1;
           else if (s < sum_send)
                   if (sp > sum_send)
                      k = 1;
                   else
                      a++;
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
     for (i = ob, s = 0; i <= oe; i++)
        if (rem[i] < a)
           s += send[i] = a - rem[i];
     while (s < sum_send)
        for (i = ob; i <= oe && s < sum_send; i++)
           if (send[i]) {
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
     for (i = mb, j = 0; i < rank; i++)      /* only overwrote other half */
        j += send[i];
     for (i = ob, k = 0; k <= j; i++)
        k += send[i];
     i--;
     num_to = (outgoing < (k - j)) ? outgoing : (k - j);
     for (k = 0; k < num_to; k++)
        proclist[k] = i;
     outgoing -= num_to;
     while (outgoing > 0) {
        i++;
        num_to = (outgoing < send[i]) ? outgoing : send[i];
        outgoing -= num_to;
        for (j = 0; j < num_to; j++, k++)
           proclist[k] = i;
     }

     /* free memory and return */
     LB_FREE(&send);
     LB_FREE(&rem);

     return LB_OK;
}
