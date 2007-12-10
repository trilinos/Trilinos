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
#include <float.h>
#include "zoltan_mem.h"
#include "par_median_const.h"
#include "par_tflops_special_const.h"
#include "par_average_const.h"
#include "zoltan_timer.h"

#define TINY   1.0e-6
#define ABS(x) ( ((x)>0) ? (x) : (-(x)))

/* PIVOT_CHOICE_ORDERED: call Zoltan_RB_find_median which walks through
 *   potential medians (pivots) in order by numeric value.
 * PIVOT_CHOICE_MEDIAN_OF_RANDOM: call Zoltan_RB_find_median_randomized,
 *   which makes a somewhat random choice of potential pivot, by trying
 *   the median value of a small random selection of pivots.
 * PIVOT_CHOICE_RANDOM: call Zoltan_RB_find_median_randomized,
 *   which makes a very random choice of potential pivot, by trying
 *   the somewhat random value from a small random selection of pivots.
 *   Choice of potential pivot is faster than "MEDIAN_OF_RANDOM" because
 *   we don't find the median of the small random selection of pivots.
 */
#define PIVOT_CHOICE_ORDERED 1
#define PIVOT_CHOICE_MEDIAN_OF_RANDOM 2
#define PIVOT_CHOICE_RANDOM 3


/* Data structure for parallel find median routine */

struct median {
  double    totallo;    /* weight of elements less than median */
  double    totalhi;    /* weight of elements greater than median */
  double    totalmed;   /* weight of elements equal to median */
};

static struct Zoltan_Timer *timer;
static int timerNum;
static int myProc=-1;
static char debugText[64];
static int loopCount;


/*
** First change: tmp_half is the median of a random value found
** on each process with active dots, or just the median of 3 of them (faster).
** Also, allreduce is smaller (but shouldn't make a big difference).
*/

static int med3(double *v, int a, int b, int c);
static int reorder_list(double *l, int len, int pivotIdx, int *i);
static double serial_find_median(double *dots, double *wgts, int dotnum);


/* Given a random dot from each process that still has active dots,
 * get the median of those random dots.
 * We could include weights in the median calculation, but I'm not
 * sure the extra effort would pay off.
 * 
 * Maybe this is overkill and we should just call random_candidate() below.
 */

static double random_median_candidate(MPI_Comm comm, double dot, double invalidDot)
{
int rank, size, i, ndots;
double candidate;
double *values=NULL;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  values = (double *)ZOLTAN_MALLOC(size * sizeof(double));

  if (rank==0){

    MPI_Gather(&dot, 1, MPI_DOUBLE, values, 1, MPI_DOUBLE, 0, comm);

    for (i=0, ndots=0; i<size; i++){
      if (values[i] != invalidDot){
        values[ndots++] = values[i];
      }
    }

    candidate = serial_find_median(values, NULL, ndots);
  }
  else{
    MPI_Gather(&dot, 1, MPI_DOUBLE, values, 1, MPI_DOUBLE, 0, comm);
  }

  ZOLTAN_FREE(&values);

  MPI_Bcast(&candidate, 1, MPI_DOUBLE, 0, comm);

  return candidate;
}
/* Overkill? If we just want a random value from the remaining,
 * why not pick any one?
 */
static double serial_find_median(double *dots, double *wgts, int dotnum)
{
int lb, ub, idx, i;
int *widx=NULL;
double median=-1, lowBal=0.0, upBal=0.0;
double pivotBal, w1, w2;

  if (wgts){
    widx = malloc(dotnum * sizeof(int));
    for (i=0; i<dotnum; i++){
      widx[i] = i;  /* map from dot to its weight */
    }
  }

  lb = 0;
  ub = dotnum-1;

  while (lb < ub){
    /* choose a random pivot value */
    idx = med3(dots, lb, ub, (lb+ub) >> 1);

    /* rearrange list around pivot, get index of pivot */
    if (widx){
      idx = reorder_list(dots+lb, ub-lb+1, idx-lb, widx+lb);
    }
    else{
      idx = reorder_list(dots+lb, ub-lb+1, idx-lb, NULL);
    }

    idx += lb;   /* relative to dots[0] */

    /* sum the weights in both halves */

    if (wgts){
      w1 = lowBal;
      w2 = upBal;
      for (i=lb; i<idx; i++){
        w1 += wgts[widx[i]];
      }
      for (i=idx+1; i<=ub; i++){
        w2 += wgts[widx[i]];
      }
      pivotBal = wgts[widx[idx]];
    }
    else{
      /* weights are all 1.0 */
      w1 = idx;
      pivotBal = 1.0;
      w2 = dotnum - idx - 1;
    }

    if (w1 >= (pivotBal + w2)){
      ub = idx - 1;
      if (wgts){
        upBal = pivotBal + w2;
      }
    }
    else if ((w1 + pivotBal) >= w2 ){
      median = dots[idx];
      break;
    }
    else{
      lb = idx + 1;
      if (wgts){
        lowBal = w1 + pivotBal;
      }
    }
  }
  if (widx) free(widx);

  if (lb == ub){
    median = dots[lb];
  }

  return median;
}
#define exchange(a, b) {   \
  tempDouble = vals[a];    \
  vals[a] = vals[b];       \
  vals[b] = tempDouble;    \
  if (idxList){            \
  tempInt = idxList[a];    \
  idxList[a] = idxList[b]; \
  idxList[b] = tempInt; } }

static int reorder_list(double *vals, int len, int pivotIdx, int *idxList)
{
int l, r, i, j;
double p = vals[pivotIdx];
double tempDouble;
int tempInt;

  i = l = 0;
  j = r = len-1;

  /* error here */
  exchange(l, pivotIdx);
  if (vals[r] >= p){
    exchange(r, l);
  }
  while (i < j){
    exchange(i, j);
    while (vals[++i] < p);
    while ( (j > l) && (vals[--j] >= p));
  }

  if (vals[l] == p){
    exchange(l, j);
  }
  else{
    j++;
    exchange(j, r);
  }

  return (j-l);
}
/* If random_median_candidate is overkill, just... */
static int candidate_choice=0;
static double random_candidate(MPI_Comm comm, double dot, double invalidDot)
{
int rank, size, i, ndots, offset, c2, c3;
double candidate;
double *values=NULL;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  values = (double *)ZOLTAN_MALLOC(size * sizeof(double));

  if (rank==0){

    MPI_Gather(&dot, 1, MPI_DOUBLE, values, 1, MPI_DOUBLE, 0, comm);

    for (i=0, ndots=0; i<size; i++){
      if (values[i] != invalidDot){
        values[ndots++] = values[i];
      }
    }

    candidate_choice--;
    if ((candidate_choice < 0) || (candidate_choice >= ndots)){
      candidate_choice = ndots-1;
    }

    if (ndots >= 3){
      offset = ndots / 3;
      c2 = (candidate_choice+offset)%ndots;
      c3 = (c2+offset)%ndots;
      i = med3(values, candidate_choice, c2, c3);
      candidate = values[i];
    }
    else{
      candidate = values[candidate_choice];
    }
  }
  else{
    MPI_Gather(&dot, 1, MPI_DOUBLE, values, 1, MPI_DOUBLE, 0, comm);
  }

  ZOLTAN_FREE(&values);

  MPI_Bcast(&candidate, 1, MPI_DOUBLE, 0, comm);

  return candidate;
}
/*
** index of the array element with median value of three values
*/
static int med3(double *v, int a, int b, int c)
{
double v1 = v[a];
double v2 = v[b];
double v3 = v[c];

  if (((v1 <= v2) && (v1 >= v3)) ||
      ((v1 >= v2) && (v1 <= v3))){
    return a;
  }
  else if (((v2 <= v1) && (v2 >= v3)) ||
      ((v2 >= v1) && (v2 <= v3))){
    return b;
  }
  else{
    return c;
  }
}

static int mark_median(int *dotlist, int *dotmark, 
        int start, int nmeds, double *dots, double med, int mark)
{
int i, j;
int count=0;

  for (j = start; ; j++) {
    i = dotlist[j];
    if (dots[i] == med) { 
      dotmark[i] = mark;
      if (++count == nmeds) break;
    }
  }
  return j+1;
}


/************ R O U T I N E S   I N   T H I S   F I L E  **********************

       NAME                             TYPE
----------------------------------------------------------------------
	Zoltan_RB_find_median_randomized			void

******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RB_find_median_randomized(
  int Tflops_Special,   /* Flag indicating whether Tflops_Special handling 
                           of communicators should be done (to avoid memory
                           leaks in tflops' Comm_Dup and Comm_Split).        */
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
  int num_procs,        /* Number of procs in set (Tflops_Special)     */
  int proclower,        /* Lowest numbered proc in set (Tflops_Special)*/
  int num_parts,        /* Number of partitions in set (Tflops_Special) */
  int wgtflag,          /* True if user supplied weights */
  double valuemin,      /* minimum value in partition (input) */
  double valuemax,      /* maximum value in partition (input) */
  double weight,        /* weight of entire partition (input) */
  double *wgtlo,        /* weight of lower partition (output) */
  double *wgthi,        /* weight of upper partition (output) */
  int    *dotlist,      /* list of active dots */
  int rectilinear_blocks,/*if set all dots with same value on same side of cut*/
  int average_cuts,      /* force cut to be halfway between two closest dots. */
  int pivot_choice      
)
{
/* Local declarations. */
  struct median med, medme;          /* median data */

  double  wtmax, wtsum, wtupto;
  double  tolerance;                 /* largest single weight of a dot */
  double  targetlo, targethi;        /* desired wt in lower half */
  double  weightlo, weighthi;        /* wt in lower/upper half of non-active */
  double  diff1, diff2;
  double  dot, tmp_half = 0.0;
  double  invalidDot = valuemax + 1.0;

  int     i, j, k, numlist, ndots;
  int     markactive;                /* which side of cut is active = 0/1 */
  int     rank=0;                    /* rank in partition (Tflops_Special) */

  int     indexmed;        /* index of my first element equal to the median */
  int     countmed;        /* how many elements I have equal to the median */
  int     left=0, middle=2, right=1;   
  int     leftTotal, rightTotal;

  /* MPI data types and user functions */

  MPI_Op            med_op;
  MPI_Datatype      med_type;
  MPI_User_function Zoltan_RB_median_merge2;

  if (myProc < 0) myProc = proc;
  sprintf(debugText,"(%d - %d)",proclower,proclower+num_procs-1);
  timer = Zoltan_Timer_Create(ZOLTAN_TIME_WALL);
  timerNum = Zoltan_Timer_Init(timer, 0, debugText);
  Zoltan_Timer_Start(timer, timerNum, local_comm, __FILE__, __LINE__);
  loopCount = 0;

  /**************************** BEGIN EXECUTION ******************************/

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(sizeof(struct median),MPI_CHAR,&med_type);
  MPI_Type_commit(&med_type);

  if (!Tflops_Special)
     MPI_Op_create(&Zoltan_RB_median_merge2,1,&med_op);

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

  if (Tflops_Special) {
    rank = proc - proclower;
    if (wgtflag) {

      /* find tolerance (max of wtmax) */
      tolerance = wtmax;
      Zoltan_RB_max_double(&tolerance, 1, proclower, rank, num_procs, local_comm);
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

  if (!Tflops_Special || num_procs > 1) { /* don't need to go thru if only
                                             one proc with Tflops_Special. 
                                             Input argument Tflops_Special 
                                             should be 0 for
                                             serial partitioning. */
    while (1) {

      loopCount++;

      /* globally choose a quick random bisector value from active dots */
    
      if (first_guess){
        tmp_half = *valuehalf;
        first_guess = 0;
      }
      else {
        if (numlist > 0){
          dot = dots[dotlist[numlist >> 1]];
        }
        else{
          dot = invalidDot;
        }
        if (pivot_choice == PIVOT_CHOICE_MEDIAN_OF_RANDOM)
          tmp_half = random_median_candidate(local_comm, dot, invalidDot);
        else 
          tmp_half = random_candidate(local_comm, dot, invalidDot);
      }

      /* initialize local median data structure */

      medme.totallo = medme.totalmed = medme.totalhi = 0.0;
      countmed = 0;
      indexmed = -1;

      /* mark all active dots on one side or other of bisector or at bisector */
      /* count number of dots equal to bisector */

      for (j = 0; j < numlist; j++) {
        i = dotlist[j];
        if (dots[i] < tmp_half) { 
          medme.totallo += wgts[i];
          dotmark[i] = left;
        }
        else if (dots[i] == tmp_half) {
          medme.totalmed += wgts[i];
          countmed++;
          dotmark[i] = middle;
          if (indexmed < 0) indexmed = j;
        }
        else{
          medme.totalhi += wgts[i];
          dotmark[i] = right;
        }
      }

      med.totallo = med.totalmed = med.totalhi = 0.0;

      /* combine median data struct across current subset of procs */
      if (counter != NULL) (*counter)++;
      if (Tflops_Special) {
         i = 1;
         Zoltan_RB_reduce(num_procs, rank, proc, (void *) &medme, (void *) &med,
                          sizeof(medme), &i, med_type, local_comm, 
                          Zoltan_RB_median_merge2);
      }
      else {
        
         MPI_Allreduce(&medme,&med,1,med_type,med_op,local_comm);
      }

      leftTotal = weightlo + med.totallo;
      rightTotal = weighthi + med.totalhi;

      if (leftTotal + med.totalmed < targetlo){  /* left half too small */
        weightlo = leftTotal + med.totalmed;
        if (indexmed >= 0){
          /* tmp_half elements go in the left half */
          mark_median(dotlist, dotmark, indexmed, countmed, dots, tmp_half, left);
        }
 
        if (targetlo - weightlo <= tolerance){  /* close enough */
          weighthi = weight - weightlo;
          break;
        }
        /* median value is in the right half */
        markactive = right;
      }
      else if (leftTotal > targetlo){          /* left half is too large */
        weighthi = rightTotal + med.totalmed;
        if (indexmed >= 0){
          /* tmp_half elements go in the right half */
          mark_median(dotlist, dotmark, indexmed, countmed, dots, tmp_half, right);
        }
        if (leftTotal - targetlo <= tolerance){  /* close enough */
          weightlo = weight - weighthi;
          break;
        }
        /* median value is in the left half */
        markactive = left;
      }
      else{                                     /* median is tmp_half */
        weightlo = leftTotal;
        weighthi = rightTotal;

        diff1 = targetlo - (leftTotal + med.totalmed);
        diff2 = targetlo - leftTotal;

        MPI_Allreduce(&countmed, &ndots, 1, MPI_INT, MPI_SUM, local_comm);

        if ((ndots == 1) ||         /* there's only one element with median value */
            (rectilinear_blocks)){  /* all median elements have to stay together */

          if (ABS(diff1) < ABS(diff2)){
            if (indexmed >= 0){
              mark_median(dotlist, dotmark, indexmed, countmed, dots, tmp_half, left);
            }
            weightlo += med.totalmed;
          }
          else{
            if (indexmed >= 0){
              mark_median(dotlist, dotmark, indexmed, countmed, dots, tmp_half, right);
            }
            weighthi += med.totalmed;
          }
        }
        else{ /* divide median elements between left & right for best balance */

          if (Tflops_Special){
            Zoltan_RB_scan_double(&medme.totalmed, &wtupto, 1, local_comm,
                                  proc, rank, num_procs);
          }
          else{
            MPI_Scan(&medme.totalmed, &wtupto, 1, MPI_DOUBLE, MPI_SUM, local_comm);
          }
          medme.totallo = medme.totalhi = 0;
          med.totallo = med.totalhi = 0;

          if (indexmed >= 0){

            if (leftTotal + wtupto - medme.totalmed >= targetlo - tolerance){
              /* all my median elements can go on the right side */
              mark_median(dotlist, dotmark, indexmed, countmed, dots, tmp_half, right);
              medme.totalhi = medme.totalmed;
            }
            else if (leftTotal + wtupto <= targetlo + tolerance){
              /* all my median elements can go on the left side */
              mark_median(dotlist, dotmark, indexmed, countmed, dots, tmp_half, left);
              medme.totallo = medme.totalmed;
            }
            else {
              /* my median elements are split between left and right sides */
              j = indexmed;
              wtsum = leftTotal + wtupto - medme.totalmed;
              k = 0;
              for (i=0; i<countmed; i++){
                j = mark_median(dotlist, dotmark, j, 1, dots, tmp_half, left);
                medme.totallo += wgts[dotlist[j-1]];
                k++;
                if (wtsum + medme.totallo >= targetlo - tolerance){
                  break;
                }
              }
              for (i=0; i < countmed-k; i++){
                j = mark_median(dotlist, dotmark, j, 1, dots, tmp_half, right);
                medme.totalhi += wgts[dotlist[j-1]];
              }
            }
          }
          if (Tflops_Special) {
             i = 1;
             Zoltan_RB_reduce(num_procs, rank, proc, (void *) &medme, (void *) &med,
                              sizeof(medme), &i, med_type, local_comm, 
                              Zoltan_RB_median_merge2);
          }
          else {
             MPI_Allreduce(&medme,&med,1,med_type,med_op,local_comm);
          }
          weightlo += med.totallo;
          weighthi += med.totalhi;
        }
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
    weighthi = 0.;
    tmp_half = valuemax;
  }

  /* found median */
  *valuehalf = tmp_half;

  if (average_cuts) 
    *valuehalf = Zoltan_RB_Average_Cut(Tflops_Special, dots, dotmark, dotnum,
                                       num_procs, rank, proc, local_comm,
                                       *valuehalf);

  *wgtlo = weightlo;
  *wgthi = weighthi;

  MPI_Type_free(&med_type);
  if (!Tflops_Special)
     MPI_Op_free(&med_op);

  Zoltan_Timer_Stop(timer, timerNum, local_comm, __FILE__, __LINE__);
  Zoltan_Timer_Print(timer, timerNum, 0, local_comm, stdout);
  Zoltan_Timer_Destroy(&timer);
  if (proclower==myProc)
    printf("%s loop count %d interval length %d median (%lf - %lf) %lf\n",
    debugText, loopCount,dotnum,
    valuemin, valuemax, *valuehalf);
  fflush(stdout);
  MPI_Barrier(local_comm);

  return 1;
}
void Zoltan_RB_median_merge2(void *in, void *inout, int *len, MPI_Datatype *dptr)
{
  struct median *med1,*med2;

  med1 = (struct median *) in;
  med2 = (struct median *) inout;

  med2->totallo += med1->totallo;
  med2->totalhi += med1->totalhi;
  med2->totalmed += med1->totalmed;
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
