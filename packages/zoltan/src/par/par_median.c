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
 *    Revision: 1.6.2.1 $
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "lb_const.h"
#include "par_median_const.h"
#include "mem_const.h"

#define MYHUGE 1.0e30
#define TINY   1.0e-6

/* Data structure for parallel find median routine */

struct median {          /* median cut info */
  double    totallo, totalhi;   /* weight in each half of active partition */
  double    valuelo, valuehi;   /* position of dot(s) nearest to cut */
  double    wtlo, wthi;         /* total weight of dot(s) at that position */
  int       countlo, counthi;   /* # of dots at that position */
  int       proclo, prochi;     /* unique proc who owns a nearest dot */
};


/* prototypes for TFLOPS_SPECIAL */
static void Zoltan_RB_reduce(int, int, int, struct median*, struct median*, int *,
               MPI_Datatype, MPI_Comm);
static void Zoltan_RB_scan(double *, double *, MPI_Comm, int, int, int);
static void Zoltan_RB_sum_double(double *, int, int, int, MPI_Comm);
static void Zoltan_RB_max_double(double *, int, int, int, MPI_Comm);

/************ R O U T I N E S   I N   T H I S   F I L E  **********************

       NAME                             TYPE
----------------------------------------------------------------------
	Zoltan_RB_find_median			void
	Zoltan_RB_median_merge			void

******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_find_median(
  ZZ *zz,               /* The Zoltan structure                      */
  double *dots,         /* array of coordinates                              */
  double *wgts,         /* array of weights associated with dots             */
  int *dotmark,         /* returned list of which side of the median
                           each dot is on:
                                0 - dot is < valuehalf
                                1 - dot is > valuehalf                       */
  int dotnum,           /* number of dots (length of three previous arrays   */
  int proc,             /* this proc number (rank)                           */
  double fractionlo,    /* fraction of weight that should be in bottom half  */
  MPI_Comm local_comm,  /* MPI communicator on which to find median          */
  double *valuehalf,    /* on entry - first guess at median (if first_guess set)                           on exit - the median value                        */
  int first_guess,      /* if set, use value in valuehalf as first guess     */
  int *counter,         /* returned for stats, # of median interations       */
  int nprocs,           /* Total number of processors (Tflops_Special)       */
  int num_procs,        /* Number of procs in partition (Tflops_Special)     */
  int proclower,        /* Lowest numbered proc in partition (Tflops_Special)*/
  int wgtflag,          /* True if user supplied weights */
  double valuemin,      /* minimum value in partition (input) */
  double valuemax,      /* maximum value in partition (input) */
  double weight,        /* weight of entire partition (input) */
  double *wtlo,         /* weight of lower partition (output) */
  double *wthi,         /* weight of upper partition (output) */
  int    *dotlist,      /* list of active dots */
  int rectilinear_blocks/* if set all dots with same value on same side of cut*/
)
{
/* Local declarations. */
  struct median med, medme;          /* median data */

  double  wtmax, wtsum, wtok, wtupto;/* temporary wts */
  double  tolerance;                 /* largest single weight of a dot */
  double  targetlo, targethi;        /* desired wt in lower half */
  double  weightlo, weighthi;        /* wt in lower/upper half of non-active */
  double  tmp_half;

  int     i, j, k, numlist;
  int     first_iteration;
  int     indexlo, indexhi;          /* indices of dot closest to median */
  int     breakflag;                 /* for breaking out of median iteration */
  int     markactive;                /* which side of cut is active = 0/1 */
  int     rank;                      /* rank in partition (Tflops_Special) */

  /* MPI data types and user functions */

  MPI_Op            med_op;
  MPI_Datatype      med_type;
  MPI_User_function Zoltan_RB_median_merge;


/***************************** BEGIN EXECUTION ******************************/

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(sizeof(struct median),MPI_CHAR,&med_type);
  MPI_Type_commit(&med_type);

  if (!zz->Tflops_Special)
     MPI_Op_create(&Zoltan_RB_median_merge,1,&med_op);

  /*
   * intialize the dotlist array
   * while looping through, find:
   *	wtmax		- max weight on this proc
   *
   * weight = summed weight of entire partition
   * search tolerance = largest single weight (plus epsilon)
   * targetlo = desired weight in lower half of partition
   * targethi = desired weight in upper half of partition
   */
  wtmax = 0.0;
  numlist = dotnum;
  for (i = 0; i < dotnum;i++) {
    dotlist[i] = i;

    if (wgtflag)
      if (wgts[i] > wtmax) wtmax = wgts[i];
  }

  if (zz->Tflops_Special) {
    rank = proc - proclower;
    if (wgtflag) {

      /* find tolerance (max of wtmax) */
      tolerance = wtmax;
      Zoltan_RB_max_double(&tolerance, proclower, rank, num_procs, local_comm);

      tmp_half = 0.0;    /* in case of a set with one processor return a value*/
    }
    else
      tolerance = 1.0;   /* if user did not supply weights, all are 1.0 */
  }
  else {
    if (wgtflag)
      MPI_Allreduce(&wtmax,&tolerance,1,MPI_DOUBLE,MPI_MAX,local_comm);
    else
      tolerance = 1.0;   /* if user did not supply weights, all are 1.0 */
  }

  tolerance *= 0.5 + TINY;  /* ctv - changed from        1.0 + TINY
               The larger tolerance allowed one side of the cut to be larger
               than the target weight by a node of largest weight (the other
               side would be smaller by the same amount).  In that case a
               node of largest weight could be moved from the side whose weight
               is larger than its target to the other side and they would both
               be in balance with the target weight.  A tolerance less than
               half of the largest weight would allow infinite looping as a
               node of largest weight was passed back and forth. */
  targetlo = fractionlo * weight;
  targethi = weight - targetlo;

  /* weightlo/hi = total weight in non-active parts of partition */
  weighthi = weightlo = 0.0;

  first_iteration = 1;

  /* median iteration */
  /* zoom in on bisector until correct # of dots in each half of partition */
  /* as each iteration of median-loop begins, require:
          all non-active dots are marked with 0/1 in dotmark
          valuemin <= every active dot <= valuemax
          weightlo, weighthi = total wt of non-active dots */
  /* when leave median-loop, require only:
          valuehalf = correct cut position
          all dots <= valuehalf are marked with 0 in dotmark
          all dots >= valuehalf are marked with 1 in dotmark */

  if (num_procs > 1) { /* don't need to go through if only one proc.  This
                          added for Tflops_Special */
  while (1) {

    /* choose bisector value */
    /* use old value on 1st iteration if old cut dimension is the same */
    /* on 2nd option: could push valuehalf towards geometric center 
       with "1.0-factor" to force overshoot */

      if (first_iteration && first_guess) {
        tmp_half = *valuehalf;
        if (tmp_half < valuemin || tmp_half > valuemax)
          tmp_half = 0.5 * (valuemin + valuemax);
      }
      else if (weight)
        tmp_half = valuemin + (targetlo - weightlo) /
                    (weight - weightlo - weighthi) * (valuemax - valuemin);
      else
        tmp_half = 0.5 * (valuemin + valuemax);

      first_iteration = 0;
      /* initialize local median data structure */

      medme.totallo = medme.totalhi = 0.0;
      medme.valuelo = -MYHUGE;
      medme.valuehi = MYHUGE;
      medme.wtlo = medme.wthi = 0.0;
      medme.countlo = medme.counthi = 0;
      medme.proclo = medme.prochi = proc;

      /* mark all active dots on one side or other of bisector */
      /* also set all fields in median data struct */
      /* save indices of closest dots on either side */

      for (j = 0; j < numlist; j++) {
        i = dotlist[j];
        if (dots[i] <= tmp_half) {            /* in lower part */
          medme.totallo += wgts[i];
          dotmark[i] = 0;
          if (dots[i] > medme.valuelo) {       /* my closest dot */
            medme.valuelo = dots[i];
            medme.wtlo = wgts[i];
            medme.countlo = 1;
            indexlo = i;
          }                                            /* tied for closest */
          else if (dots[i] == medme.valuelo) {
            medme.wtlo += wgts[i];
            medme.countlo++;
          }
        }
        else {                                         /* in upper part */
          medme.totalhi += wgts[i];
          dotmark[i] = 1;
          if (dots[i] < medme.valuehi) {       /* my closest dot */
            medme.valuehi = dots[i];
            medme.wthi = wgts[i];
            medme.counthi = 1;
            indexhi = i;
          }                                            /* tied for closest */
          else if (dots[i] == medme.valuehi) {
            medme.wthi += wgts[i];
            medme.counthi++;
          }
        }
      }

      /* combine median data struct across current subset of procs */
      if (counter != NULL) (*counter)++;
      if (zz->Tflops_Special) {
         i = 1;
         Zoltan_RB_reduce(num_procs, rank, proc, &medme, &med, &i, med_type,
                   local_comm);
      }
      else
         MPI_Allreduce(&medme,&med,1,med_type,med_op,local_comm);

      /* test median guess for convergence */
      /* move additional dots that are next to cut across it */

      if (weightlo + med.totallo < targetlo) {    /* lower half TOO SMALL */

        weightlo += med.totallo;
        tmp_half = med.valuehi;

        if (med.counthi == 1) {                  /* only one dot to move */
          if (weightlo + med.wthi < targetlo) {  /* move it, keep iterating */
            if (proc == med.prochi) dotmark[indexhi] = 0;
          }
          else {                                 /* only move if beneficial */
            if (weightlo + med.wthi - targetlo < targetlo - weightlo) {
              if (proc == med.prochi) dotmark[indexhi] = 0;
              weightlo += med.wthi;
            }
            weighthi = weight - weightlo;
            break;                               /* all done */
          }
        }
        else {                                   /* multiple dots to move */
          breakflag = 0;
          wtok = 0.0;
          if (medme.valuehi == med.valuehi) wtok = medme.wthi;   
          if (weightlo + med.wthi >= targetlo) {                /* all done */
            if (rectilinear_blocks) {
              if (weightlo + med.wthi - targetlo > targetlo - weightlo)
                wtok = 0.0;                      /* don't move if moving group
                                                    of dots has worse balance*/
            } else {
              if (zz->Tflops_Special)
                Zoltan_RB_scan(&wtok, &wtupto, local_comm, proc, rank, 
                               num_procs);
              else
                MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,local_comm);
              wtmax = targetlo - weightlo;
              if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
            }
            breakflag = 1;
          }                                      /* wtok = most I can move */
          for (j = 0, wtsum = 0.0; j < numlist && wtsum < wtok; j++) {
            i = dotlist[j];
            if (dots[i] == med.valuehi) { /* only move if better */
              if (wtsum + wgts[i] - wtok < wtok - wtsum) {
                dotmark[i] = 0;
                wtsum += wgts[i];  /* KDD Moved sum inside if test 1/2002 */
              }
            }
          }
          if (breakflag) {                        /* done if moved enough */
            if (zz->Tflops_Special) {
              wtok = wtsum;
              Zoltan_RB_sum_double(&wtok, proclower, rank, num_procs, local_comm);
            }
            else
              MPI_Allreduce(&wtsum, &wtok, 1, MPI_DOUBLE, MPI_SUM, local_comm);
            weightlo += wtok;
            weighthi = weight - weightlo;
            break;
          }
        }

        weightlo += med.wthi;
        if (targetlo-weightlo <= tolerance) {     /* close enough */
           weighthi = weight - weightlo;
           break;
        }

        valuemin = med.valuehi;                   /* iterate again */
        markactive = 1;
      }

      else if (weighthi + med.totalhi < targethi) {  /* upper half TOO SMALL */

        weighthi += med.totalhi;
        tmp_half = med.valuelo;

        if (med.countlo == 1) {                  /* only one dot to move */
          if (weighthi + med.wtlo < targethi) {  /* move it, keep iterating */
            if (proc == med.proclo) dotmark[indexlo] = 1;
          }
          else {                                 /* only move if beneficial */
            if (weighthi + med.wtlo - targethi < targethi - weighthi) {
              if (proc == med.proclo) dotmark[indexlo] = 1;
              weighthi += med.wtlo;
            }
            weightlo = weight - weighthi;
            break;                               /* all done */
          }
        }
        else {                                   /* multiple dots to move */
          breakflag = 0;
          wtok = 0.0;
          if (medme.valuelo == med.valuelo) wtok = medme.wtlo;   
          if (weighthi + med.wtlo >= targethi) {                /* all done */
            if (rectilinear_blocks) {
              if (weighthi + med.wtlo - targethi > targethi - weighthi)
                wtok = 0.0;                      /* don't move if moving group
                                                    of dots has worse balance*/
            } else {
              if (zz->Tflops_Special)
                Zoltan_RB_scan(&wtok, &wtupto, local_comm, proc, rank, 
                               num_procs);
              else
                MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,local_comm);
              wtmax = targethi - weighthi;
              if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
            }
            breakflag = 1;
          }                                      /* wtok = most I can move */
          for (j = 0, wtsum = 0.0; j < numlist && wtsum < wtok; j++) {
            i = dotlist[j];
            if (dots[i] == med.valuelo) { /* only move if better */
              if (wtsum + wgts[i] - wtok < wtok - wtsum) {
                dotmark[i] = 1;
                wtsum += wgts[i]; /* KDD Moved sum inside if test 1/2002 */
              }
            }
          }
          if (breakflag) {                        /* done if moved enough */
            if (zz->Tflops_Special) {
              wtok = wtsum;
              Zoltan_RB_sum_double(&wtok, proclower, rank, num_procs, local_comm);
            }
            else
              MPI_Allreduce(&wtsum, &wtok, 1, MPI_DOUBLE, MPI_SUM, local_comm);
            weighthi += wtok;
            weightlo = weight - weighthi;
            break;
          }
        }

        weighthi += med.wtlo;
        if (targethi-weighthi <= tolerance) {     /* close enough */
          weightlo = weight - weighthi;
          break;
        }

        valuemax = med.valuelo;                   /* iterate again */
        markactive = 0;
      }

      else {                /* Goldilocks result: both partitions JUST RIGHT */
        weightlo += med.totallo;
        weighthi += med.totalhi;
        break;
      }

      /* shrink the active list */
      
      k = 0;
      for (j = 0; j < numlist; j++) {
        i = dotlist[j];
        if (dotmark[i] == markactive) dotlist[k++] = i;
      }
      numlist = k;

  }
  }
  else { /* if one processor set all dots to 0 (Tflops_Special) */
    for (i = 0; i < numlist; i++)
       dotmark[i] = 0;
    weightlo = weight;
  }

  /* found median */
  *valuehalf = tmp_half;
  *wtlo = weightlo;
  *wthi = weighthi;

  MPI_Type_free(&med_type);
  if (!zz->Tflops_Special)
     MPI_Op_free(&med_op);

  return 1;

}

/* merge median data structure */
/* on input:
   in,inout->totallo, totalhi = weight in both partitions on this proc
             valuelo, valuehi = pos of nearest dot(s) to cut on this proc
             wtlo, wthi       = total wt of dot(s) at that pos on this proc
             countlo, counthi = # of dot(s) nearest to cut on this proc
             proclo, prochi = not used
   on exit:
   inout->   totallo, totalhi = total # of active dots in both partitions
             valuelo, valuehi = pos of nearest dot(s) to cut
             wtlo, wthi       = total wt of dot(s) at that position
             countlo, counthi = total # of dot(s) nearest to cut
             proclo, prochi = one unique proc who owns a nearest dot
                              all procs must get same proclo,prochi

*/
void Zoltan_RB_median_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)
{
  struct median *med1,*med2;

  med1 = (struct median *) in;
  med2 = (struct median *) inout;
 
  med2->totallo += med1->totallo;
  if (med1->valuelo > med2->valuelo) {
    med2->valuelo = med1->valuelo;
    med2->wtlo = med1->wtlo;
    med2->countlo = med1->countlo;
    med2->proclo = med1->proclo;
  }
  else if (med1->valuelo == med2->valuelo) {
    med2->wtlo += med1->wtlo;
    med2->countlo += med1->countlo;
    if (med1->proclo < med2->proclo) med2->proclo = med1->proclo;
  }

  med2->totalhi += med1->totalhi;
  if (med1->valuehi < med2->valuehi) {
    med2->valuehi = med1->valuehi;
    med2->wthi = med1->wthi;
    med2->counthi = med1->counthi;
    med2->prochi = med1->prochi;
  }
  else if (med1->valuehi == med2->valuehi) {
    med2->wthi += med1->wthi;
    med2->counthi += med1->counthi;
    if (med1->prochi < med2->prochi) med2->prochi = med1->prochi;
  }
}

static void Zoltan_RB_reduce(
   int nproc,             /* number of processors in partition */
   int rank,              /* rank within partition */
   int proc,              /* global processor number */
   struct median *in,     /* input median */
   struct median *inout,  /* output median */
   int *len,              /* length to pass to Zoltan_RB_median_merge */
   MPI_Datatype datatype, /* MPI datatype for median */
   MPI_Comm comm          /* MPI communicator */
)
{
   struct median tmp;     /* temporary to recieve information */
   int to;                /* communication partner */
   int tag = 32109;       /* message tag */
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
      MPI_Send(in, 1, datatype, to, tag, comm);
   else
      if (rank + nprocs_small < nproc) {
         MPI_Recv(inout, 1, datatype, to, tag, comm, &status);
         Zoltan_RB_median_merge(in, inout, len, &datatype);
      }
      else
         *inout = *in;

   if (!(rank & nprocs_small))    /* binary exchange on nprocs_small procs */
      for (mask = nprocs_small >> 1; mask; mask >>= 1) {
         tag++;
         to = proc - rank + (rank ^ mask);
         MPI_Send(inout, 1, datatype, to, tag, comm);
         MPI_Recv(&tmp, 1, datatype, to, tag, comm, &status);
         Zoltan_RB_median_merge(&tmp, inout, len, &datatype);
      }
   else
      tag += hbit;

   /* send results to the processors that are larger than nprocs_small in
      the local partition of processors */
   tag++;
   to = proc - rank + (rank ^ nprocs_small);
   if (rank & nprocs_small)
      MPI_Recv(inout, 1, datatype, to, tag, comm, &status);
   else
      if (rank + nprocs_small < nproc)
         MPI_Send(inout, 1, datatype, to, tag, comm);
}


static void Zoltan_RB_scan(
   double *wtok,          /* local weight */
   double *wtupto,        /* sum of weights for prcessors <= rank */
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
   int tag = 32108;       /* message tag */
   double tmp;            /* temporary double to recieve */
   double sendout;        /* temporary double to send */
   MPI_Status status;

   /* this subroutine performs a scan operation summing doubles on a subset
      of a set of processors for Tflops_Special */

   /* Find next lower power of 2. */
   for (hbit = 0; (num_procs >> hbit) != 0; hbit++);

   nprocs_large = 1 << hbit;
   if (nprocs_large == 2*num_procs) nprocs_large = num_procs;

   sendout = *wtupto = *wtok;
   for (mask = 1; mask <= nprocs_large; mask *= 2) {
      tag++;
      tor = (rank ^ mask);
      to = proc - rank + tor;
      if (tor < num_procs) {
         MPI_Send(&sendout, 1, MPI_DOUBLE, to, tag, local_comm);
         MPI_Recv(&tmp, 1, MPI_DOUBLE, to, tag, local_comm, &status);
         sendout += tmp;
         if (to < proc) *wtupto += tmp;
      }
   }
}

static void Zoltan_RB_sum_double(
   double   *x,               /* double to be summed */
   int      proclower,        /* smallest processor in partition */
   int      rank,             /* rank of processor in partition */
   int      nprocs,           /* number of processors in partition */
   MPI_Comm comm
)
{
   double   tmp;              /* temporary for sum */
   int      tag = 32100;      /* message tag */
   int      partner;          /* message partner in binary exchange */
   int      to;               /* message partner not in binary exchange */
   int      mask;             /* mask to determine communication partner */
   int      nprocs_small;     /* largest power of 2 contained in nprocs */
   int      hbit;             /* 2^hbit = nproc_small */
   MPI_Status status;

   /* This routine sums doubles on a subset of processors */
 
   /* Find next lower power of 2. */
   for (hbit = 0; (nprocs >> hbit) != 1; hbit++);
 
   nprocs_small = 1 << hbit;
   if (nprocs_small * 2 == nprocs) {
      nprocs_small *= 2;
      hbit++;
   }
 
   to = proclower + (rank ^ nprocs_small);
   if (rank & nprocs_small) {  /* processors greater than largest power of 2 */
      MPI_Send(x, 1, MPI_DOUBLE, to, tag, comm);
      tag += hbit + 1;
      MPI_Recv(x, 1, MPI_DOUBLE, to, tag, comm, &status);
   }
   else {   /* processors within greatest power of 2 */
      if (rank + nprocs_small < nprocs) {
         MPI_Recv(&tmp, 1, MPI_DOUBLE, to, tag, comm, &status);
         *x += tmp;
      }  
      for (mask = nprocs_small >> 1; mask; mask >>= 1) { /* binary exchange */
         tag++;
         partner = proclower + (rank ^ mask);
         MPI_Send(x, 1, MPI_DOUBLE, partner, tag, comm);
         MPI_Recv(&tmp, 1, MPI_DOUBLE, partner, tag, comm, &status);
         *x += tmp;
      }  
      tag++;
      if (rank + nprocs_small < nprocs)
         MPI_Send(x, 1, MPI_DOUBLE, to, tag, comm);
   }
}

static void Zoltan_RB_max_double(
   double   *x,               /* maximum value */
   int      proclower,        /* smallest processor in partition */
   int      rank,             /* rank of processor in partition */
   int      nprocs,           /* number of processors in partition */
   MPI_Comm comm
)
{
   double   tmp;              /* temporaries for min/max */
   int      tag = 32100;      /* message tag */
   int      partner;          /* message partner in binary exchange */
   int      to;               /* message partner not in binary exchange */
   int      mask;             /* mask to determine communication partner */
   int      nprocs_small;     /* largest power of 2 contained in nprocs */
   int      hbit;             /* 2^hbit = nproc_small */
   MPI_Status status;

   /* This routine finds the global max */

   /* Find next lower power of 2. */
   for (hbit = 0; (nprocs >> hbit) != 1; hbit++);
 
   nprocs_small = 1 << hbit;
   if (nprocs_small * 2 == nprocs) {
      nprocs_small *= 2;
      hbit++;
   }
 
   to = proclower + (rank ^ nprocs_small);
   if (rank & nprocs_small) {  /* processors greater than largest power of 2 */
      MPI_Send(x, 1, MPI_DOUBLE, to, tag, comm);
      tag += hbit + 1;
      MPI_Recv(x, 1, MPI_DOUBLE, to, tag, comm, &status);
   }
   else {   /* processors within greatest power of 2 */
      if (rank + nprocs_small < nprocs) {
         MPI_Recv(&tmp, 1, MPI_DOUBLE, to, tag, comm, &status);
         if (tmp > *x) *x = tmp;
      }
      for (mask = nprocs_small >> 1; mask; mask >>= 1) { /* binary exchange */
         tag++;
         partner = proclower + (rank ^ mask);
         MPI_Send(x, 1, MPI_DOUBLE, partner, tag, comm);
         MPI_Recv(&tmp, 1, MPI_DOUBLE, partner, tag, comm, &status);
         if (tmp > *x) *x = tmp;
      }
      tag++;
      if (rank + nprocs_small < nprocs)
         MPI_Send(x, 1, MPI_DOUBLE, to, tag, comm);
   }
}
