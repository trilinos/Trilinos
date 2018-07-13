/*
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine      kddevin@sandia.gov
 *                    Erik Boman        egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "zoltan_types.h"
#include "zz_util_const.h"
#include "shared.h"
#include "par_bisect_const.h"
#include "par_median_const.h"
#include "par_tflops_special_const.h"
#include "par_average_const.h"

/* EBEB: The following constants, structs, and function prototypes should
   probably be moved to a header file. */

/* #define DEBUG_BISECT  Turns on extra debugging in this file. */
#define FRACTION_SMALL 0.001  /* Smallest fraction of load allowed on
                                 either side of cut */
#define ALMOST_ONE 0.99       /* For scaling, should be slightly < 1.0 */
#define MAX_BISECT_ITER 100    /* Max. no. of iterations in main bisection
                                 loop. Avoids potential infinite loops. */

/* Data structure for parallel find bisector routine */

struct bisector {          /* bisector cut info */
  double    valuelo, valuehi;   /* position of dot(s) nearest to cut */
  ZOLTAN_GNO_TYPE countlo, counthi;   /* # of dots at that position */
  int       proclo, prochi;     /* unique proc who owns a nearest dot */
  int       nwgts;              /* number of weights (per dot) */
  double  totallo[MAX_BISECT_WGTS]; /* weight in lower half of active partition */
  double  totalhi[MAX_BISECT_WGTS]; /* weight in upper half of active partition */
  double  wtlo[MAX_BISECT_WGTS];    /* total weight of dot(s) on lo boundary */
  double  wthi[MAX_BISECT_WGTS];    /* total weight of dot(s) on hi boundary */
};

/*
  We sum up weights along a coordinate axis.
  Dots that have been permanently assigned to one half-space
  are called inactive, the rest are active.
  Weight sums are stored in the following way:

  inactive  active  | active   inactive
  ------------------+-------------------
  weightlo  totallo | totalhi  weighthi
*/



/* prototypes */
#if (RB_MAX_WGTS > 1)
static void Zoltan_reduce_bisector(int, int, int, int, struct bisector*, struct bisector*, int *, MPI_Datatype *, MPI_Comm);
static void Zoltan_bisector_copy(struct bisector*, struct bisector*);
static double Zoltan_norm(int mcnorm, int n, double *x, double *scal);
static void Zoltan_daxpy(int n, double a, double *x, double *y, double *z);
static double eval_cut_quality(int, double *, double *, double *, double *,
              int);
static void compute_weight_sums( int, int *, int, double *, double,
            double *, double *, MPI_Comm, int, int, int, int);
#endif /* RB_MAX_WGTS > 1 */

/*****************************************************************************/
/***  Main routine:  Zoltan_RB_find_bisector()                             ***/
/***                                                                       ***/
/***  Finds a bisector (cut) that partitions the set of "dots"             ***/
/***  such that the norms of the weights of the dots in the two            ***/
/***  halves are about the same. The norm function will scale              ***/
/***  each weight dimension by the fractionlo input vector,                ***/
/***  to allow for different balancing conditions for each weight.         ***/
/***  Several different norms (1,2, inf) are supported.                    ***/
/***                                                                       ***/
/***  Note: Zoltan_RB_find_bisector generalizes Zoltan_RB_find_median().   ***/
/*****************************************************************************/

int Zoltan_RB_find_bisector(
  ZZ *zz,               /* Zoltan struct, contains lots of fields. */
  int Tflops_Special,   /* usually same as zz->Tflops_Special */
  double *dots,         /* array of coordinates                              */
  double *wgts,         /* array of (multidimensional) weights associated with dots  */
  double uniformWeight, /* weight of each dot if wgts==NULL */
  int *dotmark,         /* returned list of which side of the bisector
                           each dot is on:
                                0 - dot is < valuehalf
                                1 - dot is > valuehalf                       */
  int dotnum,           /* number of dots (length of three previous arrays   */
  int nwgts,            /* number of weights (per dot) in wgts array         */
  int mcnorm,           /* norm to be used for multiweights: 1,2, or 3       */
  double *fraclo,       /* fraction of weight that should be in bottom half  */
  MPI_Comm local_comm,  /* MPI communicator on which to find bisector        */
  double *valuehalf,    /* on entry - first guess at median (if first_guess set)                           on exit - the median value                        */
  int first_guess,      /* if set, use value in valuehalf as first guess     */
  int num_procs,        /* Number of procs in partition (Tflops_Special)     */
  int proclower,        /* Lowest numbered proc in partition (Tflops_Special)*/
  int num_parts,        /* Number of partitions in set (Tflops_Special)      */
  double valuemin,      /* minimum value in partition (input) */
  double valuemax,      /* maximum value in partition (input) */
  double *weight,       /* weight of entire partition (input) */
  double *weightlo,     /* weight of lower partition (output) */
  double *weighthi,     /* weight of upper partition (output) */
  double *quality,      /* quality of the cut (output) */
  int    *dotlist,      /* list of active dots. */
  int rectilinear,      /* if 1, all dots with same value on same side of cut*/
  int average_cuts      /* force cut to be halfway between two closest dots. */
)
{
/* Local declarations. */
  char    yo[] = "Zoltan_find_bisector";
  int     proc   = zz->Proc;         /* My proc rank. */

#if (RB_MAX_WGTS <= 1)

  ZOLTAN_PRINT_ERROR(proc, yo, "Not applicable when RB_MAX_WGTS <= 1.");
  return(ZOLTAN_FATAL);

#else /* RB_MAX_WGTS > 1 */

  int     rank = proc - proclower;   /* rank in partition (Tflops_Special) */
  int     ierr = ZOLTAN_OK;          /* error code */
  struct bisector *med = NULL;       /* bisector data */
  struct bisector *medme = NULL;     /* bisector data */
  double  localmax, localmin;        /* lower/upper bounds on this proc */
  double  valuemax2, valuemin2;      /* test valuemin, valuemax */
  double  *localsum = NULL;          /* temporary sum of wts */
  double  *wtsum = NULL;             /* temporary sum of wts */
  double  *wtupto = NULL;            /* temporary sum of wts */
  double  tmp_half=0.;               /* guess for new bisection */
  double  *tmp = NULL;               /* temp array for Tflops_Special */
  double  *tmplo = NULL;             /* temp arrays for norm calculations */
  double  *tmphi = NULL;             /* temp arrays for norm calculations */
  double  *scale = NULL;             /* scaling for weights */
  double  *scalelo = NULL;           /* scaling for weights in lo */
  double  *scalehi = NULL;           /* scaling for weights in hi */
  double  normlo=0.0, normhi=0.0;    /* norms of weight vectors */
  double  oldnorm;                   /* temp norm */
  double  eps;                       /* abs. tolerance for imbalance */
  double  temp;                      /* temp variable */
  int     nprocs = zz->Num_Proc;     /* Total number of processors */
  int     indexlo=0, indexhi=0;      /* indices of dots closest to bisector */
  int     breakflag=0;               /* for breaking out of bisector iteration */
  int     markactive;                /* which side of cut is active = 0/1 */
  int     iteration = 0;             /* bisection iteration no. */
  int     i, j, k, numlist;
  int     tfs_early_exit = 0;        /* Flag used only with Tflops_Special,
                                        indicating early exit when all weight
                                        is to be put in one partition. */
  char  msg[256];                    /* for error messages */

  /* MPI data types and user functions */

  int               med_type_defined = 0;
  MPI_Op            med_op;
  MPI_Datatype      med_type;
  MPI_User_function Zoltan_bisector_merge;


/***************************** BEGIN EXECUTION ******************************/

  ZOLTAN_TRACE_ENTER(zz, yo);

#ifdef DEBUG_BISECT
  printf("[%2d] Debug: Entering Zoltan_find_bisection, "
         "nwgts=%2d, fraclo[0] = %lf \n", proc, nwgts, fraclo[0]);
  printf("[%2d] Debug: Tflops_Special=%1d, num_procs=%1d, num_parts=%1d\n",
          proc, Tflops_Special, num_procs, num_parts);
  printf("[%2d] Debug: %d dots on this proc\n", proc, dotnum);
  printf("[%2d] Debug: Coordinates = (", proc);
  for (i=0; i<dotnum; i++)
    printf("%lf  ", dots[i]);
  printf(")\n");
  printf("[%2d] Debug: valuemin=%lf, valuemax=%lf\n", proc, valuemin, valuemax);
  printf("[%2d] Debug: total weight = (%lf %lf)\n", proc, weight[0], weight[1]);
#endif


  /* Initialize scale vectors. */
  scale = (double *) ZOLTAN_MALLOC(3*nwgts*sizeof(double));
  if (!scale){
    ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  scalelo = &scale[nwgts];
  scalehi = &scale[2*nwgts];
  for (j=0; j<nwgts; j++){
    /* Make sure we avoid divide by zero. */
    temp = fraclo[j];
    if (temp < FRACTION_SMALL)           temp = FRACTION_SMALL;
    else if (temp > 1. - FRACTION_SMALL) temp = 1. - FRACTION_SMALL;
    /* Scale to make weight-norm larger where the tolerance is low. */
    /* Scale so weights are unchanged when Tol=1.1 (arbitrary choice) */
    scale[j] = 0.1/(zz->LB.Imbalance_Tol[j]-ALMOST_ONE);
    scalelo[j] = scale[j]/temp;
    scalehi[j] = scale[j]/(1.0-temp);
  }

  /* Check fraclo for incorrect and  trivial cases. */
  k = 0;
  for (j=0; j<nwgts; j++){
    if ((fraclo[j] >= 0.) && (fraclo[j] <= 1.)){
      if (fraclo[j] < FRACTION_SMALL)      k--;
      if (fraclo[j] > 1. - FRACTION_SMALL) k++;
    }
    else {
      sprintf(msg, "Invalid value for input parameter fraclo[%1d], %f",
              j, fraclo[j]);
      ZOLTAN_PRINT_ERROR(proc, yo, msg);
      ierr = ZOLTAN_FATAL;
      goto End;
    }
  }

  /* Check for early exit (unless Tflops_Special is set) */
  if (k == -nwgts){
    /* Put all dots in upper half */
    *valuehalf = valuemin;
    for (i = 0; i < dotnum; i++)
       dotmark[i] = 1;
    for (j=0; j<nwgts; j++){
      weighthi[j] = weight[j];
      weightlo[j] = 0.0;
    }
    ierr = ZOLTAN_OK;
    if (Tflops_Special)
      tfs_early_exit = 1;
    else
      goto End;
  }
  else if (k == nwgts){
    /* Put all dots in lower half */
    *valuehalf = valuemax;
    for (i = 0; i < dotnum; i++)
       dotmark[i] = 0;
    for (j=0; j<nwgts; j++){
      weightlo[j] = weight[j];
      weighthi[j] = 0.0;
    }
    ierr = ZOLTAN_OK;
    if (Tflops_Special)
      tfs_early_exit = 1;
    else
      goto End;
  }

  if (dotnum > 0) {
    /* check for illegal NULL pointers */
    if ((!dots) || (!dotmark) || (!dotlist)){
      ZOLTAN_PRINT_ERROR(proc, yo, "Required input is NULL.");
      ierr = ZOLTAN_FATAL;
      goto End;
    }

  } /* if (dotnum > 0) */

  /* Allocate space for bisector structs */
  med   = (struct bisector *) ZOLTAN_CALLOC(1, sizeof(struct bisector));
  medme = (struct bisector *) ZOLTAN_CALLOC(1, sizeof(struct bisector));
  if ((!med) || (!medme)){
    ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }



  /* Allocate space for various weight arrays in a single malloc */
  localsum = (double *) ZOLTAN_CALLOC(5*nwgts, sizeof(double));
  if (!localsum){
    ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto End;
  }
  /* Set pointers to point to distinct sections of the allocated space */
  wtsum  = &(localsum[nwgts]);
  wtupto = &(localsum[2*nwgts]);
  tmplo = &(localsum[3*nwgts]);
  tmphi = &(localsum[4*nwgts]);

  /* create MPI data and function types for bisector */
  {
    /* Describe struct bisector to MPI. Add MPI_UB at the end just to be safe.*/
    int lengths[5] = {2,2,3,4*MAX_BISECT_WGTS,1};
    MPI_Aint ind[5], offset;
    MPI_Datatype types[5];

    types[0] = types[3] = MPI_DOUBLE;
    types[1] = Zoltan_mpi_gno_type();
    types[2] = MPI_INT;
#if MPI_VERSION >= 2
    /* Since we do not use arrays of our defined type med_type, we don't need */
    /* MPI_UB or call to MPI_Type_create_resized to mark end of med_type.     */
#else /* MPI 1.x */
    /* MPI_UB is deprecated in MPI 2.0 */
    types[4] = MPI_UB;
#endif /* MPI_VERSION >= 2 */

#if MPI_VERSION >= 2
    MPI_Get_address(med, &offset);
    MPI_Get_address(&(med->countlo), &(ind[1]));
    MPI_Get_address(&(med->proclo), &(ind[2]));
    MPI_Get_address(&(med->totallo[0]), &(ind[3]));
#else /* MPI 1.x */
    /* MPI_Address is deprecated in MPI 2.0 */
    MPI_Address(med, &offset); 
    MPI_Address(&(med->countlo), &(ind[1])); 
    MPI_Address(&(med->proclo), &(ind[2])); 
    MPI_Address(&(med->totallo[0]), &(ind[3])); 
#endif /* MPI_VERSION >= 2 */

    ind[0] = 0;
    ind[1] -= offset;
    ind[2] -= offset;
    ind[3] -= offset;
    ind[4] = sizeof(struct bisector);

#if MPI_VERSION >= 2
    MPI_Type_create_struct(4, lengths, ind, types, &med_type);
    /* Since we do not use arrays of med_type, we don't need     */
    /* to call MPI_Type_create_resized to mark end of med_type.  */
#else /* MPI 1.x */
    /* MPI_Type_struct is deprecated in MPI 2.0 */
    MPI_Type_struct(5, lengths, ind, types, &med_type);
#endif /* MPI_VERSION >= 2 */

    MPI_Type_commit(&med_type);

    MPI_Op_create(&Zoltan_bisector_merge, 1, &med_op);
    med_type_defined = 1;
  }

/* EBEB: Some of the quantities computed below are given on input,
   so we do not need to recompute them. Currently there are
   a few inconsistencies so we recompute everything to be safe. */

  /*
   * intialize the dotlist array
   * while looping through, find:
   *    localmax        - max coordinate value on this proc
   *    localmin        - min coordinate value on this proc
   *    localsum        - sum of weights on this proc
   *
   * wtsum = summed weight of entire partition
   */
  localmax = -DBL_MAX;
  localmin =  DBL_MAX;
  for (j=0; j<nwgts; j++)
    localsum[j] =  0.0;
  numlist = dotnum;
  for (i = 0; i < dotnum;i++) {
    dotlist[i] = i;
    if (localmax < dots[i]) localmax = dots[i];
    if (localmin > dots[i]) localmin = dots[i];

#if 0
    if (wtflag) wgts[i] = 1.0;
#endif

    if (wgts){
      for (j=0; j<nwgts; j++)
        localsum[j] += wgts[i*nwgts+j];
    }
  }

  if (!wgts){
    localsum[0] = uniformWeight * dotnum;
  }
  if (Tflops_Special) {
     tmp = (double *) ZOLTAN_MALLOC(nprocs*sizeof(double));
     if (!tmp) {
        ZOLTAN_PRINT_ERROR(proc, yo, "Insufficient memory.");
        ierr = ZOLTAN_MEMERR;
        goto End;
     }

     /* find valuemax */
     tmp[proc] = valuemax2 = localmax;
     MPI_Allgather(&localmax, 1, MPI_DOUBLE, tmp, 1, MPI_DOUBLE, local_comm);
     for (i = proclower; i < proclower + num_procs; i++)
        if (tmp[i] > valuemax2) valuemax2 = tmp[i];

     /* find valuemin */
     tmp[proc] = valuemin2 = localmin;
     MPI_Allgather(&localmin, 1, MPI_DOUBLE, tmp, 1, MPI_DOUBLE, local_comm);
     for (i = proclower; i < proclower + num_procs; i++)
        if (tmp[i] < valuemin2) valuemin2 = tmp[i];

     /* find sum of weights - should use Zoltan_RB_sum_double? */
     for (j=0; j<nwgts; j++){
       tmp[proc] = wtsum[j] = localsum[j];
       MPI_Allgather(localsum + j, 1, MPI_DOUBLE, tmp, 1, MPI_DOUBLE, local_comm);
       for (wtsum[j] = 0.0, i = proclower; i < proclower + num_procs; i++)
          wtsum[j] += tmp[i];
     }

     ZOLTAN_FREE(&tmp);
     tmp_half = 0.0;    /* in case of a set with one processor return a value*/
  }
  else {
     MPI_Allreduce(&localmax, &valuemax2, 1, MPI_DOUBLE, MPI_MAX, local_comm);
     MPI_Allreduce(&localmin, &valuemin2, 1, MPI_DOUBLE, MPI_MIN, local_comm);
     MPI_Allreduce(localsum, wtsum, nwgts, MPI_DOUBLE, MPI_SUM, local_comm);
  }

  if (tfs_early_exit) goto End;

  /* The actual coordinate interval is [valuemin2,valuemax2].
   * The given input [[valuemin,valuemax] may be a slightly larger interval,
   * which is OK. We don't really need to compute [valuemin2,valuemax2].
   */
  if (valuemin2 < valuemin){
    if (zz->Debug_Level > ZOLTAN_DEBUG_NONE)
      printf("[%2d] Warning in %s: computed valuemin %f < input %f\n",
             proc, yo, valuemin2, valuemin);
    /* Safest bet is to use smaller value. */
    valuemin = valuemin2;
  }
  if (valuemax2 > valuemax){
    if (zz->Debug_Level > ZOLTAN_DEBUG_NONE)
      printf("[%2d] Warning in %s: computed valuemax %f > input %f\n",
             proc, yo, valuemax2, valuemax);
    /* Safest bet is to use larger value. */
    valuemax = valuemax2;
  }

  /* Verify that input 'weight' equals computed 'wtsum'. */
  for (j=0; j<nwgts; j++){
#ifdef DEBUG_BISECT
      printf("[%2d] Debug: computed wtsum[%1d] = %f, input weight = %f\n",
             proc, j, wtsum[j], weight[j]);
#endif
    /* Disable sanity check (for now) because 'weight' is sometimes incorrect.*/
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      if (wtsum[j] != weight[j])
        printf("[%2d] Warning: computed wtsum[%1d] %f does not match "
               "input %f\n", proc, j, wtsum[j], weight[j]);
    /* Use the weight sum we computed, do not trust input 'weight'. */
    weight[j] = wtsum[j];
  }

  /* weightlo/hi = total weight in non-active parts of partition */
  for (j=0; j<nwgts; j++)
    weighthi[j] = weightlo[j] = 0.0;

  /* zoom in on bisector until correct # of dots in each half of partition */
  /* as each iteration of bisector-loop begins, require:
          all non-active dots are marked with 0/1 in dotmark
          valuemin <= every active dot <= valuemax
          weightlo, weighthi = total wt of non-active dots */
  /* when leave bisector-loop, require only:
          valuehalf = correct cut position
          all dots <= valuehalf are marked with 0 in dotmark
          all dots >= valuehalf are marked with 1 in dotmark */

  if (!Tflops_Special || num_procs > 1) { /* don't need to go thru if only
                                             one proc with Tflops_Special.
                                             Input argument Tflops_Special
                                             should be 0 for
                                             serial partitioning. */

#ifdef DEBUG_BISECT
    printf("[%2d] Debug: Tflops_Special=%1d, num_procs=%1d, num_parts=%1d\n",
           proc, Tflops_Special, num_procs, num_parts);
#endif

    /* Set tolerance for each cut to min(imbal_tol)/log(p) */
    /* The imbalance tol vector is used implicitly through scaling. */
    /* EBEB Later we should implement dynamic balance tolerances. */
    temp = zz->LB.Imbalance_Tol[0];
    for (j=1; j<nwgts; j++)
      if (zz->LB.Imbalance_Tol[j]<temp)
        temp = zz->LB.Imbalance_Tol[j];
    eps = (temp-1.) / (log(num_parts)/log(2.0))
          * 0.5*Zoltan_norm(mcnorm, nwgts, wtsum, scale);

    /* bisector iteration */
    iteration = 0;
    while (++iteration <= MAX_BISECT_ITER){

      /* choose bisector value */
      /* use old value on 1st iteration if old cut dimension is the same */

      if ((iteration == 1)  && first_guess) {
        tmp_half = *valuehalf;
        if (tmp_half < valuemin || tmp_half > valuemax)
          tmp_half = 0.5 * (valuemin + valuemax);
      }
      else
        /* standard bisection */
        /* Note: could do interpolated search to improve convergence */
        tmp_half = 0.5 * (valuemin + valuemax);

#ifdef DEBUG_BISECT
      printf("[%2d] Debug: Iteration %d,  tmp_half = %f\n", 
             proc, iteration, tmp_half);
#endif

      /* initialize local bisector data structure */
      medme->valuelo = valuemin2;
      medme->valuehi = valuemax2;
      medme->countlo = medme->counthi = 0;
      medme->proclo = medme->prochi = proc;
      medme->nwgts = nwgts;
      for (j=0; j<nwgts; j++){
        medme->totallo[j] = 0.;
        medme->totalhi[j] = 0.;
        medme->wtlo[j] = 0.;
        medme->wthi[j] = 0.;
      }

      /* mark all active dots on one side or other of bisector */
      /* also set all fields in bisector data struct */
      /* save indices of closest dots on either side */
#ifdef DEBUG_BISECT
      {int nlo=0, nhi=0;
#endif

      for (j = 0; j < numlist; j++) {
        i = dotlist[j];
        if (dots[i] <= tmp_half) {            /* in lower part */
#ifdef DEBUG_BISECT
          nlo++;
#endif
          if (wgts){
            for (k=0; k<nwgts; k++)
              medme->totallo[k] += wgts[i*nwgts+k];
          }
          else{
            medme->totallo[0] += uniformWeight;
          }

          dotmark[i] = 0;
          if (dots[i] > medme->valuelo) {       /* my closest dot */
            medme->valuelo = dots[i];
            medme->countlo = 1;
            indexlo = i;
            if (wgts){
              for (k=0; k<nwgts; k++)
                medme->wtlo[k] = wgts[i*nwgts+k];
            }
            else{
              medme->wtlo[0] = uniformWeight;
            }
          }                                            /* tied for closest */
          else if (dots[i] == medme->valuelo) {
            if (wgts){
              for (k=0; k<nwgts; k++)
                medme->wtlo[k] += wgts[i*nwgts+k];
            }
            else{
              medme->wtlo[0] += uniformWeight;
            }
            medme->countlo++;
          }
        }
        else {                                         /* in upper part */
#ifdef DEBUG_BISECT
          nhi++;
#endif
          if (wgts){
            for (k=0; k<nwgts; k++)
              medme->totalhi[k] += wgts[i*nwgts+k];
          }
          else{
            medme->totalhi[0] += uniformWeight;
          }
          dotmark[i] = 1;
          if (dots[i] < medme->valuehi) {       /* my closest dot */
            medme->valuehi = dots[i];
            medme->counthi = 1;
            indexhi = i;
            if (wgts){
              for (k=0; k<nwgts; k++)
                medme->wthi[k] = wgts[i*nwgts+k];
            }
            else{
              medme->wthi[0] = uniformWeight;
            }
          }                                            /* tied for closest */
          else if (dots[i] == medme->valuehi) {
            if (wgts){
              for (k=0; k<nwgts; k++)
                medme->wthi[k] += wgts[i*nwgts+k];
            }
            else{
              medme->wthi[0] += uniformWeight;
            }
            medme->counthi++;
          }
        }
      }
#ifdef DEBUG_BISECT
      printf("[%2d] Debug: %d active dots, %d in lower, %d in upper\n",
      proc, nlo+nhi, nlo, nhi);
      }
#endif

      /* combine bisector data struct across current subset of procs */

      med->valuelo = -DBL_MAX;
      med->valuehi = DBL_MAX;
      med->countlo = med->counthi = 0;
      med->proclo = med->prochi = proc;
      med->nwgts = nwgts;
      for (j=0; j<nwgts; j++){
        med->totallo[j] = 0.;
        med->totalhi[j] = 0.;
        med->wtlo[j] = 0.;
        med->wthi[j] = 0.;
      }

      if (Tflops_Special) {
         i = 1;
         Zoltan_reduce_bisector(num_procs, rank, proc, 1, medme, med, &i,
                &med_type, local_comm);
      }
      else
         MPI_Allreduce(medme, med, 1, med_type, med_op, local_comm);

      /* test bisector guess for convergence */
      /* move additional dots that are next to cut across it */

      /* add up sums in active and inactive sectors */
      /* tmpmlo = sum of weights in lower half */
      /* tmpmhi = sum of weights in upper half */
      /* normlo = norm of weights in lower half */
      /* normhi = norm of weights in upper half */
      Zoltan_daxpy(nwgts, 1., weightlo, med->totallo, tmplo);
      normlo = Zoltan_norm(mcnorm, nwgts, tmplo, scalelo);
      Zoltan_daxpy(nwgts, 1., weighthi, med->totalhi, tmphi);
      normhi = Zoltan_norm(mcnorm, nwgts, tmphi, scalehi);

#ifdef DEBUG_BISECT
      printf("[%2d] Debug: med->valuelo = %f, med->valuehi = %f\n",
              proc, med->valuelo, med->valuehi);
      printf("[%2d] Debug: medme->totallo = (%f, %f), "
             "medme->totalhi = (%f, %f)\n",
              proc, medme->totallo[0], medme->totallo[1],
              medme->totalhi[0], medme->totalhi[1]);
      printf("[%2d] Debug: med->totallo = (%f, %f), med->totalhi = (%f, %f)\n",
              proc, med->totallo[0], med->totallo[1], 
              med->totalhi[0], med->totalhi[1]);
      printf("[%2d] Debug: med->wtlo = (%f, %f), med->wthi = (%f, %f)\n",
              proc, med->wtlo[0], med->wtlo[1], med->wthi[0], med->wthi[1]);
      printf("[%2d] Debug: weightlo = (%f, %f), weighthi = (%f, %f)\n",
              proc, weightlo[0], weightlo[1], weighthi[0], weighthi[1]);
      printf("[%2d] Debug: normlo = %f, normhi = %f, eps = %f\n",
              proc, normlo, normhi, eps);
#endif

      if (normlo < normhi - eps) {                      /* lower half TOO SMALL */

        /* move bisector to closest dot in upper half */
        tmp_half = med->valuehi;
        /* update weightlo to include weights in the active lower half */
        for (k=0; k<nwgts; k++)
          weightlo[k] += med->totallo[k];

        /* tmplo = tmplo + med->wthi */
        Zoltan_daxpy(nwgts, 1., med->wthi, tmplo, tmplo);
        /* tmphi = tmphi - med->wthi */
        Zoltan_daxpy(nwgts, -1., med->wthi, tmphi, tmphi);
        /* rectilinear case: treat several dots as a single dot. */
        if (rectilinear || (med->counthi == 1)) { /* only one dot to move */
          if (Zoltan_norm(mcnorm, nwgts, tmplo, scalelo) <
              Zoltan_norm(mcnorm, nwgts, tmphi, scalehi) ){
                                                 /* move it, keep iterating */
            if (med->counthi == 1){    /* single dot */
              if (proc == med->prochi)
                dotmark[indexhi] = 0;  /* weightlo will be updated later */
            }
            else{
              /* multiple dots on a line, move them all. */
              for (j = 0; j < numlist; j++) {
                i = dotlist[j];
                if (dots[i] == med->valuehi)
                  dotmark[i] = 0;        /* weightlo will be updated later */
              }
            }
          }
          else {                                 /* only move if beneficial */
            if (Zoltan_norm(mcnorm, nwgts, tmplo, scalelo) < normhi){
              if (med->counthi == 1){    /* single dot */
                if (proc == med->prochi)
                  dotmark[indexhi] = 0;
              }
              else{
                /* multiple dots on a line, move them all. */
                for (j = 0; j < numlist; j++) {
                  i = dotlist[j];
                  if (dots[i] == med->valuehi)
                    dotmark[i] = 0;
                }
              }
              for (k=0; k<nwgts; k++)
                weightlo[k] = tmplo[k];
            }
            break;                               /* all done */
          }
        }
        else {                                   /* multiple dots to move */
          breakflag = 0;                         /* must decide which ones */
          for (k=0; k<nwgts; k++){
            localsum[k] = 0.0;
            wtupto[k] = 0.0;
          }
          if (medme->valuehi == med->valuehi) {
            for (k=0; k<nwgts; k++)
              localsum[k] = medme->wthi[k];
          }
          if (Zoltan_norm(mcnorm, nwgts, tmplo, scalelo) >=
              Zoltan_norm(mcnorm, nwgts, tmphi, scalehi)){
            /* move some dots and all done */
            /* scan to figure out how many dots to move */
            /* wtupto will contain cumulative sum up to current proc */
            if (Tflops_Special){
              Zoltan_RB_scan_double(localsum, wtupto, nwgts, local_comm,
                proc, rank, num_procs);
            }
            else
              MPI_Scan(localsum, wtupto, nwgts, MPI_DOUBLE, MPI_SUM, local_comm);
            /* MPI_Scan is inclusive, we want to exclude my local weight */
            Zoltan_daxpy(nwgts, -1., localsum, wtupto, wtupto);
            breakflag = 1;
#ifdef DEBUG_BISECT
            printf("[%2d] Debug: breakflag = %d, moving some dots on "
                   "boundary across\n", proc, breakflag);
#endif
          }
#ifdef DEBUG_BISECT
          else
            printf("[%2d] Debug: breakflag = %d, moving all dots on "
                   "boundary across\n", proc, breakflag);
#endif
          oldnorm = normhi;
          /* reset tmplo and tmphi to undo med->wthi correction */
          /* then add wtupto contribution (different on each proc!) */
          /* tmplo -= med->wthi */
          Zoltan_daxpy(nwgts, -1.0, med->wthi, tmplo, tmplo);
          /* tmplo += wtupto */
          Zoltan_daxpy(nwgts, 1.0, wtupto, tmplo, tmplo);
          /* tmphi += med->wthi */
          Zoltan_daxpy(nwgts, 1.0, med->wthi, tmphi, tmphi);
          /* tmphi -= wtupto */
          Zoltan_daxpy(nwgts, -1.0, wtupto, tmphi, tmphi);
          /* wtsum = local weight moved across bisector */
          for (k=0; k<nwgts; k++)
            wtsum[k] = 0.0;
          /* loop thru dots; move from hi to lo if better */
          for (j = 0; j < numlist; j++) {
            i = dotlist[j];
            if (dots[i] == med->valuehi){
              if (breakflag){              /* only move if better */
                /* tmplo += wgts[i] */
                if (wgts)
                  Zoltan_daxpy(nwgts, 1., &wgts[i*nwgts], tmplo, tmplo);
                else
                  Zoltan_daxpy(nwgts, 1., &uniformWeight, tmplo, tmplo);
#ifdef DEBUG_BISECT
                printf("[%2d] Examining dot %2d = %f, norm= %f, oldnorm= %f\n",
                       proc, i, dots[i],
                       Zoltan_norm(mcnorm, nwgts, tmplo, scalelo), oldnorm);
                printf("[%2d] tmplo = (%f, %f)\n", proc, tmplo[0], tmplo[1]);
                printf("[%2d] tmphi = (%f, %f)\n", proc, tmphi[0], tmphi[1]);
#endif
                if (Zoltan_norm(mcnorm, nwgts, tmplo, scalelo) < oldnorm){
                  dotmark[i] = 0;  /* weightlo will be updated later */
                  if (wgts)
                    Zoltan_daxpy(nwgts, 1., &wgts[i*nwgts], wtsum, wtsum);
                  else
                    Zoltan_daxpy(nwgts, 1., &uniformWeight, wtsum, wtsum);

#ifdef DEBUG_BISECT
                  printf("[%2d] Debug: moving dot %d to other half, "
                         "norm(tmplo) = %g, norm(tmphi) = %g\n", 
                         proc, i, Zoltan_norm(mcnorm, nwgts, tmplo, scalelo), 
                         Zoltan_norm(mcnorm, nwgts, tmphi, scalehi));
#endif
                }
                /* tmphi -= wgts[i] */
                if (wgts)
                  Zoltan_daxpy(nwgts, -1., &wgts[i*nwgts], tmphi, tmphi);
                else
                  Zoltan_daxpy(nwgts, -1., &uniformWeight, tmphi, tmphi);
                oldnorm = Zoltan_norm(mcnorm, nwgts, tmphi, scalehi);
              }
              else                        /* move all */
                dotmark[i] = 0;  /* weightlo will be updated later */
            }
          }

#ifdef DEBUG_BISECT
          printf("[%2d] Debug: bisect value too small, breakflag = %d\n", 
                 proc, breakflag);
#endif

          if (breakflag){                  /* done if moved enough */
            /* update weightlo; add weights on cut that we moved */
            /* copy wtsum into wtupto, then sum across procs */
            for (k=0; k<nwgts; k++)
              wtupto[k] = wtsum[k];
            if (Tflops_Special){
              Zoltan_RB_sum_double(wtsum, nwgts, proclower, rank, num_procs,
                local_comm);
            }
            else
              MPI_Allreduce(wtupto, wtsum, nwgts, MPI_DOUBLE, MPI_SUM,
                local_comm);
            Zoltan_daxpy(nwgts, 1., wtsum, weightlo, weightlo);
            break;
          }
        }

        /* Didn't break out, so must have moved all closest dots across */
        for (k=0; k<nwgts; k++)
          weightlo[k] += med->wthi[k];

        valuemin = med->valuehi;                   /* iterate again */
        markactive = 1;
      }

      else if (normlo > normhi + eps) {  /* lower half TOO BIG */

        /* move bisector to closest dot in lower half */
        tmp_half = med->valuelo;
        /* update weighthi to include weights in the active upper half */
        for (k=0; k<nwgts; k++)
          weighthi[k] += med->totalhi[k];
#ifdef DEBUG_BISECT
        printf("[%2d] Debug: new weighthi = (%f, %f)\n",
               proc, weighthi[0], weighthi[1]);
#endif

        /* Update tmplo, tmphi such that dots on cut are in upper half. */
        /* tmphi = tmphi + med->wtlo */
        Zoltan_daxpy(nwgts, 1., med->wtlo, tmphi, tmphi);
        /* tmplo = tmplo - med->wtlo */
        Zoltan_daxpy(nwgts, -1., med->wtlo, tmplo, tmplo);
        /* rectilinear case: treat several dots as a single dot. */
        if (rectilinear || (med->countlo == 1)) {  /* only one dot to move */
          if (Zoltan_norm(mcnorm, nwgts, tmphi, scalehi) <
              Zoltan_norm(mcnorm, nwgts, tmplo, scalelo)){
                                                 /* move it, keep iterating */
            if (med->counthi == 1){    /* single dot */
              if (proc == med->proclo)
                dotmark[indexlo] = 1;  /* weighthi will be updated later */
            }
            else{
              /* multiple tied dots, move them all. */
              for (j = 0; j < numlist; j++) {
                i = dotlist[j];
                if (dots[i] == med->valuelo)
                  dotmark[i] = 1;        /* weighthi will be updated later */
              }
            }
          }
          else {                                 /* only move if beneficial */
            if (Zoltan_norm(mcnorm, nwgts, tmphi, scalehi) < normlo){
              if (med->counthi == 1){    /* single dot */
                if (proc == med->proclo)
                  dotmark[indexlo] = 1;
              }
              else{
                /* multiple tied dots, move them all. */
                for (j = 0; j < numlist; j++) {
                  i = dotlist[j];
                  if (dots[i] == med->valuelo)
                    dotmark[i] = 1;
                }
              }
              for (k=0; k<nwgts; k++)
                weighthi[k] = tmphi[k];
            }
            break;                               /* all done */
          }
        }
        else {                                   /* multiple dots to move */
          breakflag = 0;                         /* must decide which ones */
          for (k=0; k<nwgts; k++){
            localsum[k] = 0.0;
            wtupto[k] = 0.0;
          }
          if (medme->valuelo == med->valuelo) {
            for (k=0; k<nwgts; k++)
              localsum[k] = medme->wtlo[k];
          }
#ifdef DEBUG_BISECT
          printf("[%2d] Debug: tmplo = (%f, %f)\n", proc, tmplo[0], tmplo[1]);
          printf("[%2d] Debug: tmphi = (%f, %f)\n", proc, tmphi[0], tmphi[1]);
#endif
          if (Zoltan_norm(mcnorm, nwgts, tmphi, scalehi) >=
              Zoltan_norm(mcnorm, nwgts, tmplo, scalelo) ){
            /* move some dots and all done */
            /* scan to figure out how many dots to move */
            /* wtupto will contain cumulative sum up to current proc */
            if (Tflops_Special){
               Zoltan_RB_scan_double(localsum, wtupto, nwgts, local_comm,
                 proc, rank, num_procs);
           }
            else
               MPI_Scan(localsum, wtupto, nwgts, MPI_DOUBLE, MPI_SUM, local_comm);
            /* MPI_Scan is inclusive, we want to exclude my local weight */
            Zoltan_daxpy(nwgts, -1., localsum, wtupto, wtupto);
            breakflag = 1;
#ifdef DEBUG_BISECT
            printf("[%2d] Debug: breakflag = %d, moving some dots on "
                   "boundary across\n", proc, breakflag);
#endif
          }
#ifdef DEBUG_BISECT
          else
            printf("[%2d] Debug: breakflag = %d, moving all dots on "
                   "boundary across\n", proc, breakflag);
#endif
          oldnorm = normlo;
          /* reset tmplo and tmphi to undo earlier med->wtxx correction, */
          /* then add wtupto contribution (different on each proc!) */
          /* tmphi -= med->wtlo */
          Zoltan_daxpy(nwgts, -1.0, med->wtlo, tmphi, tmphi);
          /* tmphi += wtupto */
          Zoltan_daxpy(nwgts, 1.0, wtupto, tmphi, tmphi);
          /* tmplo += med->wtlo */
          Zoltan_daxpy(nwgts, 1.0, med->wtlo, tmplo, tmplo);
          /* tmplo -= wtupto */
          Zoltan_daxpy(nwgts, -1.0, wtupto, tmplo, tmplo);
          /* wtsum = local weight moved across bisector */
          for (k=0; k<nwgts; k++)
            wtsum[k] = 0.0;
          /* loop thru dots; move from lo to hi if better */
          for (j = 0; j < numlist; j++) {
            i = dotlist[j];
            if (dots[i] == med->valuelo) {
              if (breakflag){              /* only move if better */
                /* tmphi += wgts[i] */
                if (wgts)
                  Zoltan_daxpy(nwgts, 1., &wgts[i*nwgts], tmphi, tmphi);
                else
                  Zoltan_daxpy(nwgts, 1., &uniformWeight, tmphi, tmphi);
#ifdef DEBUG_BISECT
                printf("[%2d] Examining dot %2d = %f, norm= %f, oldnorm= %f\n",
                  proc, i, dots[i], 
                  Zoltan_norm(mcnorm, nwgts, tmphi, scalehi), oldnorm);
                printf("[%2d] tmplo = (%f, %f)\n", proc, tmplo[0], tmplo[1]);
                printf("[%2d] tmphi = (%f, %f)\n", proc, tmphi[0], tmphi[1]);
#endif
                if (Zoltan_norm(mcnorm, nwgts, tmphi, scalehi) < oldnorm){
                  dotmark[i] = 1;  /* weighthi will be updated later */
                  if (wgts)
                    Zoltan_daxpy(nwgts, 1., &wgts[i*nwgts], wtsum, wtsum);
                  else
                    Zoltan_daxpy(nwgts, 1., &uniformWeight, wtsum, wtsum);
#ifdef DEBUG_BISECT
                  printf("[%2d] Debug: moving dot %d to other half, "
                         "norm(tmplo) = %g, norm(tmphi) = %g\n",
                         proc, i, Zoltan_norm(mcnorm, nwgts, tmplo, scalelo),
                         Zoltan_norm(mcnorm, nwgts, tmphi, scalehi));
#endif
                }
                /* tmplo -= wgts[i] */
                if (wgts)
                  Zoltan_daxpy(nwgts, -1., &wgts[i*nwgts], tmplo, tmplo);
                else
                  Zoltan_daxpy(nwgts, -1., &uniformWeight, tmplo, tmplo);
                oldnorm = Zoltan_norm(mcnorm, nwgts, tmplo, scalelo);
              }
              else                        /* move all */
                dotmark[i] = 1;  /* weighthi will be updated later */
            }
          }
#ifdef DEBUG_BISECT
          printf("[%2d] Debug: bisect value too big, breakflag = %d\n", 
                 proc, breakflag);
#endif
          if (breakflag){                  /* done if moved enough */
            /* update weighthi; add weights on cut that we moved */
            /* copy wtsum into wtupto, then sum across procs */
            for (k=0; k<nwgts; k++)
              wtupto[k] = wtsum[k];
            if (Tflops_Special){
              Zoltan_RB_sum_double(wtsum, nwgts, proclower, rank, num_procs,
                local_comm);
            }
            else
              MPI_Allreduce(wtupto, wtsum, nwgts, MPI_DOUBLE, MPI_SUM,
                local_comm);
            Zoltan_daxpy(nwgts, 1., wtsum, weighthi, weighthi);
            break;
          }
        }

#ifdef DEBUG_BISECT
        printf("[%2d] Debug: A weighthi = (%f, %f)\n",
               proc, weighthi[0], weighthi[1]);
#endif
        /* Didn't break out, so must have moved all closest dots across */
        for (k=0; k<nwgts; k++)
          weighthi[k] += med->wtlo[k];
#ifdef DEBUG_BISECT
        printf("[%2d] Debug: B weighthi = (%f, %f)\n",
                proc, weighthi[0], weighthi[1]);
#endif

        valuemax = med->valuelo;                   /* iterate again */
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
    if (iteration > MAX_BISECT_ITER){
#undef KDD_NOISY
#ifdef KDD_NOISY
      char fname[37];
      FILE *fpg, *fpc, *fpa;
      int kdd;
      sprintf(fname, "kddout.graph.%04d.%04d", proc, rank);
      fpg = fopen(fname, "w");
      sprintf(fname, "kddout.coords.%04d.%04d", proc, rank);
      fpc = fopen(fname, "w");
      sprintf(fname, "kddout.assign.%04d.%04d", proc, rank);
      fpa = fopen(fname, "w");
      for (kdd = 0; kdd < dotnum; kdd++) {
        int nn;
        fprintf(fpc, "%.20f\n", dots[kdd]);
        for (nn = 0; nn < nwgts; nn++)
          fprintf(fpg, "%.20f ", wgts[kdd*nwgts+nn]);
        fprintf(fpg, "\n");
        fprintf(fpa, "%d\n", rank);
      }
      fclose(fpg);
      fclose(fpc);
      fclose(fpa);
#endif
      ierr = ZOLTAN_WARN;
      ZOLTAN_PRINT_WARN(proc, yo,
                        "MAX_BISECT_ITER reached in Zoltan/RCB.");
    }
  }
  else { /* if one processor set all dots to 0 (Tflops_Special) */
    tmp_half = valuemax;
    for (j = 0; j < dotnum; j++)
      dotmark[j] = 0;
    for (j=0; j<nwgts; j++) {
      weightlo[j] = weight[j];
      weighthi[j] = 0.0;
    }
  }

  /* found bisector */
  *valuehalf = tmp_half;
#ifdef DEBUG_BISECT
  printf("[%2d] Final bisector valuehalf = %f\n", proc, *valuehalf);
#endif

  if (average_cuts)
    *valuehalf = Zoltan_RB_Average_Cut(Tflops_Special, dots, dotmark, dotnum,
                                       num_procs, rank, proc, local_comm,
                                       *valuehalf);

End:
  /* Recompute weightlo/hi. This is only necessary because
     there is a bug so weightlo/hi aren't always correct.
     Remove this step later! */

  if (ierr == ZOLTAN_OK){
    compute_weight_sums(dotnum, dotmark, nwgts, wgts, uniformWeight, weightlo, weighthi,
      local_comm, Tflops_Special, proclower, rank, num_procs);

    /* Evaluate cut quality. */
    *quality = eval_cut_quality(nwgts, scalelo, scalehi, weightlo, weighthi,
               mcnorm);
  }

  /* free all memory */
  ZOLTAN_FREE(&med);
  ZOLTAN_FREE(&medme);
  ZOLTAN_FREE(&localsum);
  ZOLTAN_FREE(&scale);
#if 0
  if (wtflag) ZOLTAN_FREE(&wgts);
#endif

  if (med_type_defined) {
    MPI_Type_free(&med_type);
    MPI_Op_free(&med_op);
  }

  if (ierr == ZOLTAN_OK){
    par_median_accumulate_counts(nprocs, num_procs, rank, iteration);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;

#endif /* RB_MAX_WGTS > 1 */
}

#if (RB_MAX_WGTS > 1)

/* Compute weight sums in lower/upper part. */
static void compute_weight_sums(
  int dotnum,           /* Number of dots (local) */
  int *dotmark,         /* 0 for lo, 1 for hi. */
  int nwgts,            /* Number of weights (per dot). */
  double *wgts,         /* Weights array */
  double uniformWeight,
  double *weightlo,     /* Sum of weights in lower (output) */
  double *weighthi,     /* Sum of weights in upper (output) */
  MPI_Comm local_comm,  /* MPI communicator */
  int ts,               /* Tflops special */
  int proclower,
  int rank,
  int num_procs)
{
  double sumlo[MAX_BISECT_WGTS];
  double sumhi[MAX_BISECT_WGTS];
  int i,j;

  for (j=0; j<nwgts; j++){
    sumlo[j] = 0.0;
    sumhi[j] = 0.0;
  }

  if (wgts){
    for (i=0; i<dotnum; i++){
      for (j=0; j<nwgts; j++){
        if (dotmark[i]==0)
          sumlo[j] += wgts[i*nwgts+j];
        else
          sumhi[j] += wgts[i*nwgts+j];
      }
    }
  }
  else{
    for (i=0; i<dotnum; i++){
      if (dotmark[i]==0)
        sumlo[j] += uniformWeight;
      else
        sumhi[j] += uniformWeight;
    }
  }

  /* Sum local sums into global sums. */
  if (ts){
    Zoltan_RB_sum_double(sumlo, nwgts, proclower, rank, num_procs, local_comm);
    Zoltan_RB_sum_double(sumhi, nwgts, proclower, rank, num_procs, local_comm);
    for (j=0; j<nwgts; j++){
      weightlo[j] = sumlo[j];
      weighthi[j] = sumhi[j];
    }
  }
  else{
    MPI_Allreduce(sumlo, weightlo, nwgts, MPI_DOUBLE, MPI_SUM, local_comm);
    MPI_Allreduce(sumhi, weighthi, nwgts, MPI_DOUBLE, MPI_SUM, local_comm);
  }

  return;
}

/* Evaluate cut quality. */
static double eval_cut_quality(
  int nwgts,          /* Number of weights per dot. */
  double *scalelo,    /* Scale vector for lower part. */
  double *scalehi,    /* Scale vector for upper part. */
  double *weightlo,   /* Sum of weights in lower part. */
  double *weighthi,   /* Sum of weights in upper part. */
  int norm)           /* Type of norm (1,2,3) */
{
  double temp[MAX_BISECT_WGTS];
  int i;

  for (i=0; i<nwgts; i++)
     temp[i] = MAX(scalelo[i]*weightlo[i], scalehi[i]*weighthi[i]);

  return Zoltan_norm(norm, nwgts, temp, NULL);
}


/* merge bisector data structure */
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
void Zoltan_bisector_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)
{
  struct bisector *med1, *med2;
  int i, nwgts;
  char *yo="Zoltan_bisector_merge";

  med1 = (struct bisector *) in;
  med2 = (struct bisector *) inout;

  /* make sure both bisectors use the same # of weights */
  nwgts = med1->nwgts;
  if (med2->nwgts != nwgts){
    ZOLTAN_PRINT_ERROR(-1, yo, 
                       "Inconsistent number of weights in bisector structs!");
    return;
  }

  /* sum up total weights in low half */
  for (i=0; i<nwgts; i++)
    med2->totallo[i] += med1->totallo[i];

  if (med1->valuelo > med2->valuelo) {
    med2->valuelo = med1->valuelo;
    med2->countlo = med1->countlo;
    med2->proclo = med1->proclo;
    for (i=0; i<nwgts; i++)
      med2->wtlo[i] = med1->wtlo[i];
  }
  else if (med1->valuelo == med2->valuelo) {
    med2->countlo += med1->countlo;
    /* choose lowest rank processor as representative when tied (arbitrary) */
    if (med1->proclo < med2->proclo) med2->proclo = med1->proclo;
    for (i=0; i<nwgts; i++)
      med2->wtlo[i] += med1->wtlo[i];
  }

  /* sum up total weights in high half */
  for (i=0; i<nwgts; i++)
    med2->totalhi[i] += med1->totalhi[i];

  if (med1->valuehi < med2->valuehi) {
    med2->valuehi = med1->valuehi;
    med2->counthi = med1->counthi;
    med2->prochi = med1->prochi;
    for (i=0; i<nwgts; i++)
      med2->wthi[i] = med1->wthi[i];
  }
  else if (med1->valuehi == med2->valuehi) {
    med2->counthi += med1->counthi;
    /* choose lowest rank processor as representative when tied (arbitrary) */
    if (med1->prochi < med2->prochi) med2->prochi = med1->prochi;
    for (i=0; i<nwgts; i++)
      med2->wthi[i] += med1->wthi[i];
  }

}

static void Zoltan_reduce_bisector(int nproc, int rank, int proc, int n,
               struct bisector *in, struct bisector *inout, int *len,
               MPI_Datatype *datatype, MPI_Comm comm)
{
   struct bisector *tmp = NULL;
   int m, to, tag = 32109;
   MPI_Status status;

   /* this is a recursive function for Tflops_Special that takes a structure
      in in and returns the merged structure inout for a subset of processors
      of the entire number of processors.  rank is a processors rank within
      its partition of size nproc while proc is the rank of the processor with
      the entire number of processors being used. */

   m = 2*n;
   if (rank%m) {
      to = proc - n;
      MPI_Send(in, 1, *datatype, to, tag, comm);
      MPI_Recv(inout, 1, *datatype, to, tag, comm, &status);
   }
   else
      if (rank + n < nproc) {
         to = proc + n;
         MPI_Recv(inout, 1, *datatype, to, tag, comm, &status);
         Zoltan_bisector_merge(in, inout, len, datatype);
         tmp = (struct bisector *) ZOLTAN_MALLOC(sizeof(struct bisector) + 4*(in->nwgts)*sizeof(double));
         Zoltan_bisector_copy(inout, tmp);
         if (m < nproc)
            Zoltan_reduce_bisector(nproc, rank, proc, m, tmp, inout,
              len, datatype, comm);
         MPI_Send(inout, 1, *datatype, to, tag, comm);
         ZOLTAN_FREE(&tmp);
      }
      else
         Zoltan_reduce_bisector(nproc, rank, proc, m, in, inout, len,
           datatype, comm);
}

static void Zoltan_bisector_copy(struct bisector *from, struct bisector *to)
{
  /* Copy a bisector structure */
  int i;
  to->valuelo = from->valuelo;
  to->valuehi = from->valuehi;
  to->countlo = from->countlo;
  to->counthi = from->counthi;
  to->proclo  = from->proclo;
  to->prochi  = from->prochi;
  to->nwgts   = from->nwgts;
  for (i=0; i<from->nwgts; i++){
    to->totallo[i] = from->totallo[i];
    to->totalhi[i] = from->totalhi[i];
    to->wtlo[i] = from->wtlo[i];
    to->wthi[i] = from->wthi[i];
  }
}

/* Compute the appropriate norm of a vector.
    1 - 1-norm = Manhattan norm
    2 - 2-norm = sqrt of sum of squares
    3 - inf-norm = maximum norm

   If two input vectors x and scal are given, compute norm(scal.*x).
*/
static double Zoltan_norm(int mcnorm, int n, double *x, double *scal)
{
  int i;
  double tmp, result = 0.0;

  for (i=0; i<n; i++){
    if (scal == NULL)
      tmp = x[i];
    else
      tmp = x[i] * scal[i];

    if (tmp < 0.) tmp = -tmp;

    if (mcnorm==1)
      result += tmp;
    else if (mcnorm==2)
      result += tmp*tmp;
    else if (mcnorm>2) /* use infinity norm */
      if (tmp > result) result = tmp;
  }

  if (mcnorm==2) result = sqrt(result);
  return result;
}

/*
 * Compute z = alpha*x + y, where z may be the same array as x or y.
 */
static void Zoltan_daxpy(int n, double alpha, double *x, double *y, double *z)
{
  int i;
  for (i=0; i<n; i++)
    z[i] = alpha*x[i]+y[i];
}

#endif /* RB_MAX_WGTS > 1 */

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

