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


/************ R O U T I N E S   I N   T H I S   F I L E  **********************

       NAME                             TYPE
----------------------------------------------------------------------
	LB_find_median			void
	LB_median_merge			void

******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_find_median(
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
  int *counter          /* returned for stats, # of median interations       */
)
{
/* Local declarations. */
  char    yo[] = "LB_find_median";

  struct median med, medme;          /* median data */

  double  valuemax, valuemin;        /* lower/upper bounds of active region */
  double  localmax, localmin;        /* lower/upper bounds on this proc */
  double  wtmax, wtsum, wtok, wtupto;/* temporary wts */
  double  weight;                    /* wt of entire partition */
  double  tolerance;                 /* largest single weight of a dot */
  double  targetlo, targethi;        /* desired wt in lower half */
  double  weightlo, weighthi;        /* wt in lower/upper half of non-active */
  double  tmp_half;

  int     i, j, k, wtflag = 0, numlist;
  int     first_iteration;
  int    *dotlist = NULL;            /* list of active dots */
  int     indexlo, indexhi;          /* indices of dot closest to median */
  int     breakflag;                 /* for breaking out of median iteration */
  int     markactive;                /* which side of cut is active = 0/1 */

  /* MPI data types and user functions */

  MPI_Op            med_op;
  MPI_Datatype      med_type;
  MPI_User_function LB_median_merge;


/***************************** BEGIN EXECUTION ******************************/

  if (dotnum > 0) {
    /* allocate memory */
    dotlist = (int *) LB_MALLOC(dotnum*sizeof(int));
    if (!dotlist) {
      LB_PRINT_ERROR(proc, yo, "Insufficient memory.");
      return 0;
    }

    /*
     * Check to see if the user supplied weights. If not, allocate
     * memory and set the weights to 1.0.
     * NOTE: it will be much more efficient if weights are allocated
     * and set before calling this routine.
     */
    if (!wgts) {
      wtflag = 1;
      wgts = (double *) LB_MALLOC(dotnum*sizeof(double));
      if (!wgts) {
        LB_PRINT_ERROR(proc, yo, "Insufficient memory.");
        LB_FREE(&dotlist);
        return 0;
      }
    }
  } /* if (dotnum > 0) */

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(sizeof(struct median),MPI_CHAR,&med_type);
  MPI_Type_commit(&med_type);

  MPI_Op_create(&LB_median_merge,1,&med_op);

  /*
   * intialize the dotlist array
   * while looping through, find:
   *	localmax	- max coordinate value on this proc
   *	localmin	- min coordinate value on this proc
   *	wtsum		- sum of weights on this proc
   *	wtmax		- max weight on this proc
   *
   * weight = summed weight of entire partition
   * search tolerance = largest single weight (plus epsilon)
   * targetlo = desired weight in lower half of partition
   * targethi = desired weight in upper half of partition
   */
  localmax = -MYHUGE;
  localmin =  MYHUGE;
  wtsum = wtmax = 0.0;
  numlist = dotnum;
  for (i = 0; i < dotnum;i++) {
    dotlist[i] = i;
    if (localmax < dots[i]) localmax = dots[i];
    if (localmin > dots[i]) localmin = dots[i];

    if (wtflag) wgts[i] = 1.0;
    wtsum += wgts[i];
    if (wgts[i] > wtmax) wtmax = wgts[i];
  }

  MPI_Allreduce(&localmax,&valuemax,1,MPI_DOUBLE,MPI_MAX,local_comm);
  MPI_Allreduce(&localmin,&valuemin,1,MPI_DOUBLE,MPI_MIN,local_comm);
  MPI_Allreduce(&wtsum,&weight,1,MPI_DOUBLE,MPI_SUM,local_comm);
  MPI_Allreduce(&wtmax,&tolerance,1,MPI_DOUBLE,MPI_MAX,local_comm);

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
            if (weightlo + med.wthi - targetlo < targetlo - weightlo)
            if (proc == med.prochi) dotmark[indexhi] = 0;
            break;                               /* all done */
          }
        }
        else {                                   /* multiple dots to move */
          breakflag = 0;
          wtok = 0.0;
          if (medme.valuehi == med.valuehi) wtok = medme.wthi;   
          if (weightlo + med.wthi >= targetlo) {                /* all done */
            MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,local_comm);
            wtmax = targetlo - weightlo;
            if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
            breakflag = 1;
          }                                      /* wtok = most I can move */
          for (j = 0, wtsum = 0.0; j < numlist && wtsum < wtok; j++) {
            i = dotlist[j];
            if (dots[i] == med.valuehi) { /* only move if better */
              if (wtsum + wgts[i] - wtok < wtok - wtsum)
                dotmark[i] = 0;
              wtsum += wgts[i];
            }
          }

          if (breakflag) break;                   /* done if moved enough */
        }

        weightlo += med.wthi;
        if (targetlo-weightlo <= tolerance) break;  /* close enough */

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
          if (weighthi + med.wtlo - targethi < targethi - weighthi)
            if (proc == med.proclo) dotmark[indexlo] = 1;
          break;                               /* all done */
        }
      }
      else {                                   /* multiple dots to move */
        breakflag = 0;
        wtok = 0.0;
        if (medme.valuelo == med.valuelo) wtok = medme.wtlo;   
        if (weighthi + med.wtlo >= targethi) {                /* all done */
          MPI_Scan(&wtok,&wtupto,1,MPI_DOUBLE,MPI_SUM,local_comm);
          wtmax = targethi - weighthi;
          if (wtupto > wtmax) wtok = wtok - (wtupto - wtmax);
          breakflag = 1;
        }                                      /* wtok = most I can move */
        for (j = 0, wtsum = 0.0; j < numlist && wtsum < wtok; j++) {
          i = dotlist[j];
          if (dots[i] == med.valuelo) { /* only move if better */
            if (wtsum + wgts[i] - wtok < wtok - wtsum) 
              dotmark[i] = 1;
            wtsum += wgts[i];
          }
        }
        if (breakflag) break;                   /* done if moved enough */
      }

      weighthi += med.wtlo;
      if (targethi-weighthi <= tolerance) break;  /* close enough */

      valuemax = med.valuelo;                   /* iterate again */
      markactive = 0;
    }

    else                  /* Goldilocks result: both partitions JUST RIGHT */
      break;

    /* shrink the active list */
      
    k = 0;
    for (j = 0; j < numlist; j++) {
      i = dotlist[j];
      if (dotmark[i] == markactive) dotlist[k++] = i;
    }
    numlist = k;

  }

  /* found median */
  *valuehalf = tmp_half;

  /* free all memory */
  LB_FREE(&dotlist);
  if (wtflag) LB_FREE(&wgts);

  MPI_Type_free(&med_type);
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
void LB_median_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)
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
