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
#ifndef lint
static char *cvs_rcbc_id = "$Id$";
#endif

/* Recursive coordinate bisectioning (RCB) routine
   operates on "dots" as defined in rcb.h
*/

/* Steve Plimpton, Sandia National Labs, ABQ, NM  87185
   Dept 9221, MS 1111
   (505) 845-7873
   sjplimp@cs.sandia.gov
*/

/* Notes:
   dots are balanced across procs by weight (if used)
   on return, proc owns dotnum "dots" in dense array of max-length dotmax
   all dots will be inside (or on surface of) 3-d box defined by rcbbox
   input weights (if used) are real numbers > 0.0
   can extend "rcb_dot" data structure in calling program, see rcb.h
   returned RCB tree only contains one cut on each proc,
     need to do MPI_Allgather if wish to collect it on all procs
   set RCB_CHECK and RCB_STATS at top of rcb.c as desired
*/

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "lb_const.h"
#include "rcb_const.h"
#include "all_allo_const.h"

#define MYHUGE 1.0e30
#define TINY   1.0e-6

#define RCB_DEFAULT_OVERALLOC 1.0
#define RCB_DEFAULT_REUSE FALSE

/*  RCB_CHECK = 0  No consistency check on input or results */
/*  RCB_CHECK = 1  Check input weights and final results for consistency */
static int RCB_CHECK = 1;

/*  RCB_STATS = 0  No statistics logging */
/*  RCB_STATS = 1  Log times and counts, print summary */
/*  RCB_STATS = 2  Log times and counts, print for each proc */
static int RCB_STATS = 1;

void lb_rcb(
  LB *lb,                     /* The load-balancing structure with info for
                                 the RCB balancer.                           */
  int *num_import,            /* Number of non-local objects assigned to this
                                 processor in the new decomposition.         */
  LB_GID **import_global_ids, /* Returned value:  array of global IDs for
                                 non-local objects in this processor's new
                                 decomposition.                              */
  LB_LID **import_local_ids,  /* Returned value:  array of local IDs for
                                 non-local objects in this processor's new
                                 decomposition.                              */
  int **import_procs          /* Returned value:  array of processor IDs for
                                 processors owning the non-local objects in
                                 this processor's new decomposition.         */
)
{
  int    proc,nprocs;              /* my proc id, total # of procs */
  struct rcb_dot *dotbuf, *dotpt;  /* local dot arrays */
  struct rcb_median med, medme;    /* median data */
  struct rcb_box boxtmp;           /* tmp rcb box */
  int    keep, outgoing;           /* message exchange counters */
  int    incoming, incoming2;      /* message exchange counters */
  int    pdotnum;                  /* # of dots - decomposition changes it */
  int    pdottop;                  /* dots >= this index are new */
  int   *dotmark;                  /* which side of median for each dot */
  int   *dotlist;                  /* list of active dots */
  int    dotnum;                   /* number of dots */
  int    dotmax = 0;               /* max # of dots arrays can hold */
  int    dottop;                   /* dots >= this index are new */
  int    dotnew;                   /* # of new dots after send/recv */
  int    numlist;                  /* number of active dots I own */
  int    proclower, procupper;     /* lower/upper proc in partition */
  int    procmid;                  /* 1st proc in upper half of part */
  int    procpartner, procpartner2; /* proc(s) to exchange with */
  int    readnumber;               /* # of proc partner(s) to read from */
  int    markactive;               /* which side of cut is active = 0/1 */
  int    dim;                      /* which of 3 axes median cut is on */
  int    indexlo, indexhi;         /* indices of dot closest to median */
  double valuemin, valuemax;       /* lower/upper bounds of active region */
  double valuehalf;                /* median cut position */
  double targetlo, targethi;       /* desired wt in lower/upper half */
  double weight;                   /* wt of entire partition */
  double weightlo, weighthi;       /* wt in lower/upper half of non-active */
  double wtsum,wtok,wtupto,wtmax;  /* temporary wts */
  double tolerance;                /* largest single weight of a dot */
  int    first_iteration;          /* flag for 1st time thru median search */
  int    allocflag;                /* have to re-allocate space */
  int    breakflag;                /* for breaking out of median iteration */
  int    length;                   /* message length */
  double time1,time2,time3,time4;  /* timers */
  double timestart,timestop;       /* timers */
  double timers[4];                /* diagnostic timers 
			              0 = start-up time before recursion
				      1 = time before median iterations
				      2 = time in median iterations
				      3 = communication time */
  int    counters[7];              /* diagnostic counts
			              0 = # of median iterations
				      1 = # of dots sent
				      2 = # of dots received
				      3 = most dots this proc ever owns
				      4 = most dot memory this proc ever allocs
				      5 = # of times a previous cut is re-used
				      6 = # of reallocs of dot array */
  int    i,ii,j,k;                 /* local variables */

  RCB_STRUCT *rcb;                 /* Pointer to data structures for RCB.  */
  int wtflag;                      /* (0) do not (1) do use weights.  */
  double overalloc;                /* amount to overallocate by when realloc
                                      of dot array must be done.     
                                      1.0 = no extra; 1.5 = 50% extra; etc. */
  int reuse;                       /* (0) don't use (1) use previous cuts
                                      stored in treept at initial guesses.  */
  struct rcb_box *rcbbox;          /* bounding box of final RCB sub-domain */
  struct rcb_tree *treept;         /* tree of RCB cuts - only single cut on 
                                      exit */

  double LB_start_time, LB_end_time;
  double LB_time[2], LB_max_time[2];

  /* function prototypes */

  void rcb_error(int);
  void rcb_check(struct rcb_dot *, int, int, struct rcb_box *);
  void rcb_stats(double, struct rcb_dot *,int, double *, 
		 int *, struct rcb_box *, int);

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;
  MPI_Request request, request2;
  MPI_Status status;
  MPI_Op box_op, med_op;
  MPI_Datatype box_type, med_type;
  MPI_User_function rcb_box_merge, rcb_median_merge;

  if (RCB_STATS) {
    MPI_Barrier(MPI_COMM_WORLD);
    timestart = time1 = MPI_Wtime();
  }

  /* setup for parallel */

  proc = LB_Proc;
  nprocs = LB_Num_Proc;

  /*
   *  Build the RCB Data structure and 
   *  set pointers to information in it.
   */

  LB_start_time = MPI_Wtime();
  rcb_build_data_structure(lb, &pdotnum, &dotmax);

  rcb = (RCB_STRUCT *) (lb->Data_Structure);

  dotpt  = rcb->Dots; 
  rcbbox = rcb->Box;
  treept = rcb->Tree_Ptr;
  wtflag = (lb->Get_Obj_Weight != NULL);   /* Use weights if application 
                                             specified a weight function  */
  LB_end_time = MPI_Wtime();
  LB_time[0] = LB_end_time - LB_start_time;
  LB_start_time = LB_end_time;

  if (lb->Params == NULL) {
    /* 
     *  No application-specified parameters; use defaults.
     */
    overalloc = RCB_DEFAULT_OVERALLOC;
    reuse = RCB_DEFAULT_REUSE;
  }
  else {
    if (lb->Params[0] == LB_PARAMS_INIT_VALUE)
      overalloc = RCB_DEFAULT_OVERALLOC;
    else 
      overalloc = lb->Params[0];

    if (lb->Params[1] == LB_PARAMS_INIT_VALUE)
      reuse = RCB_DEFAULT_REUSE;
    else
      reuse = lb->Params[1];
  }
  
  /* local copies of calling parameters */

  dottop = dotnum = pdotnum;

  /* initialize timers and counters */

  if (RCB_STATS) {
    counters[0] = 0;
    counters[1] = 0;
    counters[2] = 0;
    counters[3] = dotnum;
    counters[4] = dotmax;
    counters[5] = 0;
    counters[6] = 0;
    timers[1] = 0.0;
    timers[2] = 0.0;
    timers[3] = 0.0;
  }

  /* create mark and list arrays for dots */

  allocflag = 0;
  if (dotmax > 0) {
    dotmark = (int *) malloc((unsigned) dotmax * sizeof(int));
    dotlist = (int *) malloc((unsigned) dotmax * sizeof(int));
    if (dotmark == NULL || dotlist == NULL) rcb_error(dotmax*sizeof(int));
  }
  else {
    dotmark = NULL;
    dotlist = NULL;
  }

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(6,MPI_DOUBLE,&box_type);
  MPI_Type_commit(&box_type);
  MPI_Type_contiguous(sizeof(struct rcb_median),MPI_CHAR,&med_type);
  MPI_Type_commit(&med_type);

  MPI_Op_create(&rcb_box_merge,1,&box_op);
  MPI_Op_create(&rcb_median_merge,1,&med_op);

  /* set dot weights = 1.0 if user didn't */

  if (!wtflag)
    for (i = 0; i < dotnum; i++) dotpt[i].Weight = 1.0;

  /* check that all weights > 0 */

  if (RCB_CHECK) {
    j = 0;
    for (i = 0; i < dotnum; i++) if (dotpt[i].Weight <= 0.0) j++;
    MPI_Allreduce(&j,&k,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (k > 0) {
      if (proc == 0) printf("RCB ERROR: %d dot weights are <= 0\n",k);
      return;
    }
  }

  /* initialize sub-domain bounding box to entire domain */

  boxtmp.lo[0] = boxtmp.lo[1] = boxtmp.lo[2] = MYHUGE;
  boxtmp.hi[0] = boxtmp.hi[1] = boxtmp.hi[2] = -MYHUGE;
  
  for (i = 0; i < dotnum; i++) {
    for (j = 0; j < 3; j++) {
      if (dotpt[i].X[j] < boxtmp.lo[j])
	boxtmp.lo[j] = dotpt[i].X[j];
      if (dotpt[i].X[j] > boxtmp.hi[j])
	boxtmp.hi[j] = dotpt[i].X[j];
    }
  }

  MPI_Allreduce(&boxtmp,rcbbox,1,box_type,box_op,MPI_COMM_WORLD);

  /* create local communicator for use in recursion */

  MPI_Comm_dup(MPI_COMM_WORLD,&local_comm);

  if (RCB_STATS) {
    time2 = MPI_Wtime();
    timers[0] = time2 - time1;
  }

  /* recursively halve until just one proc in partition */
  
  proclower = 0;
  procupper = nprocs - 1;

  while (proclower != procupper) {

    if (RCB_STATS) time1 = MPI_Wtime();
    
    /* procmid = 1st proc in upper half of partition */
    /* if odd # of procs, lower partition gets extra one */

    procmid = proclower + (procupper - proclower) / 2 + 1;

    /* determine communication partner(s) */

    if (proc < procmid)
      procpartner = proc + (procmid - proclower);
    else
      procpartner = proc - (procmid - proclower);
    
    readnumber = 1;
    if (procpartner > procupper) {
      readnumber = 0;
      procpartner--;
    }
    if (proc == procupper && procpartner != procmid - 1) {
      readnumber = 2;
      procpartner2 = procpartner + 1;
    }
    
    /* weight = summed weight of entire partition */
    /* search tolerance = largest single weight (plus epsilon) */
    /* targetlo = desired weight in lower half of partition */
    /* targethi = desired weight in upper half of partition */

    wtmax = wtsum = 0.0;
    for (i = 0; i < dotnum; i++) {
      wtsum += dotpt[i].Weight;
      if (dotpt[i].Weight > wtmax) wtmax = dotpt[i].Weight;
    }

    MPI_Allreduce(&wtsum,&weight,1,MPI_DOUBLE,MPI_SUM,local_comm);
    MPI_Allreduce(&wtmax,&tolerance,1,MPI_DOUBLE,MPI_MAX,local_comm);

    tolerance *= 1.0 + TINY;
    targetlo = ((double) (procmid - proclower)) / 
      (procupper + 1 - proclower) * weight;
    targethi = weight - targetlo;

    /* dim = dimension (xyz = 012) to bisect on */

    dim = 0;
    if (rcbbox->hi[1] - rcbbox->lo[1] >
	rcbbox->hi[0] - rcbbox->lo[0])
      dim = 1;
    if (dim == 0 && rcbbox->hi[2] - rcbbox->lo[2] >
	rcbbox->hi[0] - rcbbox->lo[0])
      dim = 2;
    if (dim == 1 && rcbbox->hi[2] - rcbbox->lo[2] >
	rcbbox->hi[1] - rcbbox->lo[1])
      dim = 2;
    
    /* create mark array and active list for dots */

    if (allocflag) {
      allocflag = 0;
      free(dotlist);
      free(dotmark);
      dotmark = (int *) malloc((unsigned) dotmax * sizeof(int));
      dotlist = (int *) malloc((unsigned) dotmax * sizeof(int));
      if (dotmark == NULL || dotlist == NULL) rcb_error(dotmax*sizeof(int));
    }

    /* initialize active list to all dots */

    numlist = dotnum;
    for (i = 0; i < dotnum; i++) dotlist[i] = i;

    /* weightlo/hi = total weight in non-active parts of partition */

    weightlo = weighthi = 0.0;
    valuemin = rcbbox->lo[dim];
    valuemax = rcbbox->hi[dim];
    first_iteration = 1;

    if (RCB_STATS) time2 = MPI_Wtime();

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

      if (first_iteration && reuse && dim == treept[procmid].dim) {
	if (RCB_STATS) counters[5]++;
	valuehalf = treept[procmid].cut;
	if (valuehalf < valuemin || valuehalf > valuemax)
	  valuehalf = 0.5 * (valuemin + valuemax);	  
      }
      else if (weight)
	valuehalf = valuemin + (targetlo - weightlo) /
	  (weight - weightlo - weighthi) * (valuemax - valuemin);
      else
	valuehalf = 0.5 * (valuemin + valuemax);

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
	if (dotpt[i].X[dim] <= valuehalf) {            /* in lower part */
	  medme.totallo += dotpt[i].Weight;
	  dotmark[i] = 0;
	  if (dotpt[i].X[dim] > medme.valuelo) {       /* my closest dot */
	    medme.valuelo = dotpt[i].X[dim];
	    medme.wtlo = dotpt[i].Weight;
	    medme.countlo = 1;
	    indexlo = i;
	  }                                            /* tied for closest */
	  else if (dotpt[i].X[dim] == medme.valuelo) {
	    medme.wtlo += dotpt[i].Weight;
	    medme.countlo++;
	  }
	}
	else {                                         /* in upper part */
	  medme.totalhi += dotpt[i].Weight;
	  dotmark[i] = 1;
	  if (dotpt[i].X[dim] < medme.valuehi) {       /* my closest dot */
	    medme.valuehi = dotpt[i].X[dim];
	    medme.wthi = dotpt[i].Weight;
	    medme.counthi = 1;
	    indexhi = i;
	  }                                            /* tied for closest */
	  else if (dotpt[i].X[dim] == medme.valuehi) {
	    medme.wthi += dotpt[i].Weight;
	    medme.counthi++;
	  }
	}
      }

      /* combine median data struct across current subset of procs */

      if (RCB_STATS) counters[0]++;
      MPI_Allreduce(&medme,&med,1,med_type,med_op,local_comm);

      /* test median guess for convergence */
      /* move additional dots that are next to cut across it */

      if (weightlo + med.totallo < targetlo) {    /* lower half TOO SMALL */

	weightlo += med.totallo;
	valuehalf = med.valuehi;

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
	    if (dotpt[i].X[dim] == med.valuehi) { /* only move if better */
	      if (wtsum + dotpt[i].Weight - wtok < wtok - wtsum)
		dotmark[i] = 0;
	      wtsum += dotpt[i].Weight;
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
	valuehalf = med.valuelo;

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
	    if (dotpt[i].X[dim] == med.valuelo) { /* only move if better */
	      if (wtsum + dotpt[i].Weight - wtok < wtok - wtsum) 
		dotmark[i] = 1;
	      wtsum += dotpt[i].Weight;
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

    if (RCB_STATS) time3 = MPI_Wtime();

    /* store cut info in tree only if I am procmid */

    if (proc == procmid) {
      treept[proc].dim = dim;
      treept[proc].cut = valuehalf;
    }
    
    /* use cut to shrink RCB domain bounding box */

    if (proc < procmid)
      rcbbox->hi[dim] = valuehalf;
    else
      rcbbox->lo[dim] = valuehalf;

    /* outgoing = number of dots to ship to partner */
    /* dottop = number of dots that have never migrated */

    markactive = (proc < procpartner);
    for (i = 0, keep = 0, outgoing = 0; i < dotnum; i++)
      if (dotmark[i] == markactive)
	outgoing++;
      else if (i < dottop)
	keep++;
    dottop = keep;
    
    /* alert partner how many dots I'll send, read how many I'll recv */

    MPI_Send(&outgoing,1,MPI_INT,procpartner,0,MPI_COMM_WORLD);
    incoming = 0;
    if (readnumber) {
      MPI_Recv(&incoming,1,MPI_INT,procpartner,0,MPI_COMM_WORLD,&status);
      if (readnumber == 2) {
	MPI_Recv(&incoming2,1,MPI_INT,procpartner2,0,MPI_COMM_WORLD,&status);
	incoming += incoming2;
      }
    }

    /* check if need to malloc more space */

    dotnew = dotnum - outgoing + incoming;

    if (dotnew > dotmax) {
      allocflag = 1;
      dotmax = overalloc * dotnew;
      if (dotmax < dotnew) dotmax = dotnew;
      dotpt = (struct rcb_dot *) 
	realloc(dotpt,(unsigned) dotmax * sizeof(struct rcb_dot));
      if (dotpt == NULL) rcb_error(dotmax*sizeof(struct rcb_dot));
      rcb->Dots = dotpt;
      if (RCB_STATS) counters[6]++;
    }

    if (RCB_STATS) {
      counters[1] += outgoing;
      counters[2] += incoming;
      if (dotnew > counters[3]) counters[3] = dotnew;
      if (dotmax > counters[4]) counters[4] = dotmax;
    }
    
    /* malloc comm send buffer */

    if (outgoing > 0) {
      dotbuf = (struct rcb_dot *)
        malloc((unsigned) outgoing * sizeof(struct rcb_dot));
      if (dotbuf == NULL) rcb_error(outgoing*sizeof(struct rcb_dot));
    }
    else 
      dotbuf = NULL;

    /* fill buffer with dots that are marked for sending */
    /* pack down the unmarked ones */
    
    keep = outgoing = 0;
    for (i = 0; i < dotnum; i++) {
      if (dotmark[i] == markactive)
	memcpy((char *) &dotbuf[outgoing++], (char *) &dotpt[i], 
               sizeof(struct rcb_dot));
      else
	memcpy((char *) &dotpt[keep++], (char *) &dotpt[i], 
               sizeof(struct rcb_dot));
    }

    /* post receives for dot data */

    if (readnumber > 0) {
      length = incoming * sizeof(struct rcb_dot);
      MPI_Irecv(&dotpt[keep],length,MPI_CHAR,
			   procpartner,1,MPI_COMM_WORLD,&request);
      if (readnumber == 2) {
	keep += incoming - incoming2;
	length = incoming2 * sizeof(struct rcb_dot);
	MPI_Irecv(&dotpt[keep],length,MPI_CHAR,
			      procpartner2,1,MPI_COMM_WORLD,&request2);
      }
    }
    
    /* handshake before sending data to insure recvs have been posted */
    
    if (readnumber > 0) {
      MPI_Send(NULL,0,MPI_INT,procpartner,0,MPI_COMM_WORLD);
      if (readnumber == 2)
	MPI_Send(NULL,0,MPI_INT,procpartner2,0,MPI_COMM_WORLD);
    }
    MPI_Recv(NULL,0,MPI_INT,procpartner,0,MPI_COMM_WORLD,&status);

    /* send dot data to partner */

    length = outgoing * sizeof(struct rcb_dot);
    MPI_Rsend(dotbuf,length,MPI_CHAR,procpartner,1,MPI_COMM_WORLD);
    free(dotbuf);
    
    dotnum = dotnew;

    /* wait until all dots are received */

    if (readnumber > 0) {
      MPI_Wait(&request,&status);
      if (readnumber == 2) MPI_Wait(&request2,&status);
    }

    /* cut partition in half, create new communicators of 1/2 size */

    if (proc < procmid) {
      procupper = procmid - 1;
      i = 0;
    }
    else {
      proclower = procmid;
      i = 1;
    }

    MPI_Comm_split(local_comm,i,proc,&tmp_comm);
    MPI_Comm_free(&local_comm);
    local_comm = tmp_comm;

    if (RCB_STATS) {
      time4 = MPI_Wtime();
      timers[1] += time2 - time1;
      timers[2] += time3 - time2;
      timers[3] += time4 - time3;
    }

  }

  /* have recursed all the way to final single sub-domain */

  /* free all memory used by RCB and MPI */

  MPI_Comm_free(&local_comm);
  MPI_Type_free(&med_type);
  MPI_Type_free(&box_type);
  MPI_Op_free(&box_op);
  MPI_Op_free(&med_op);

  free(dotlist);
  free(dotmark);

  LB_end_time = MPI_Wtime();
  LB_time[1] = LB_end_time - LB_start_time;

  if (RCB_STATS) {
    MPI_Barrier(MPI_COMM_WORLD);
    timestop = time1 = MPI_Wtime();
  }

  /* error checking and statistics */

  if (RCB_CHECK) rcb_check(dotpt,dotnum,pdotnum,rcbbox);
  if (RCB_STATS) rcb_stats(timestop-timestart,dotpt,dotnum,
			   timers,counters,rcbbox,reuse);

  /* update calling routine parameters */
  
  LB_start_time = MPI_Wtime();

  pdotnum = dotnum;
  pdottop = dottop;

  /*  build return arguments */

  *num_import = dotnum - dottop;
  if (*num_import > 0) {
    *import_global_ids = (LB_GID *) LB_array_alloc(__FILE__, __LINE__, 1,
                                                  *num_import, sizeof(LB_GID));
    *import_local_ids  = (LB_LID *) LB_array_alloc(__FILE__, __LINE__, 1,
                                                  *num_import, sizeof(LB_LID));
    *import_procs      = (int *) LB_array_alloc(__FILE__, __LINE__, 1,
                                                  *num_import, sizeof(int));

    for (i = 0; i < *num_import; i++) {
      ii = i + dottop;
      (*import_global_ids)[i] = dotpt[ii].Tag.Global_ID;
      (*import_local_ids)[i]  = dotpt[ii].Tag.Local_ID;
      (*import_procs)[i]      = dotpt[ii].Tag.Proc;
    }
  }
  LB_end_time = MPI_Wtime();
  LB_time[0] += (LB_end_time - LB_start_time);

  MPI_Allreduce(LB_time, LB_max_time, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (LB_Proc == 0) {
    printf("DLBLIB RCB Times:  \n");
    printf("DLBLIB     Build:  %f\n", LB_max_time[0]);
    printf("DLBLIB     RCB:    %f\n", LB_max_time[1]);
  }

  if (LB_Debug > 6) {
    int i;
    LB_print_sync_start(lb, TRUE);
    printf("DLBLIB RCB Proc %d  Num_Obj=%d  Num_Keep=%d  Num_Non_Local=%d\n", 
           LB_Proc, pdotnum, pdottop, *num_import);
    printf("  Assigned objects:\n");
    for (i = 0; i < pdotnum; i++) {
      printf("    Obj:  %10d      Orig: %4d\n", rcb->Dots[i].Tag.Global_ID,
             rcb->Dots[i].Tag.Proc);
    }
    printf("  Non_locals:\n");
    for (i = 0; i < *num_import; i++) {
      printf("    Obj:  %10d      Orig: %4d\n",
             (*import_global_ids)[i], (*import_procs)[i]);
    }
    LB_print_sync_end(lb, TRUE);
  }
}



/* ----------------------------------------------------------------------- */

/* error message for malloc/realloc overflow */

void rcb_error(int size)

{
  int proc;

  MPI_Comm_rank(MPI_COMM_WORLD,&proc);
  printf("RCB ERROR: proc = %d could not malloc/realloc %d bytes",proc,size);
  exit(1);
}


/* MPI user-defined reduce operations */

/* min/max merge of each component of a rcb_box */

void rcb_box_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)

{
  int i;
  struct rcb_box *box1,*box2;

  box1 = (struct rcb_box *) in;
  box2 = (struct rcb_box *) inout;

  for (i = 0; i < 3; i++) {
    if (box1->lo[i] < box2->lo[i])
      box2->lo[i] = box1->lo[i];
    if (box1->hi[i] > box2->hi[i])
      box2->hi[i] = box1->hi[i];
  }
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

void rcb_median_merge(void *in, void *inout, int *len, MPI_Datatype *dptr)

{
  struct rcb_median *med1,*med2;

  med1 = (struct rcb_median *) in;
  med2 = (struct rcb_median *) inout;
  
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


/* consistency checks on RCB results */

void rcb_check(struct rcb_dot *dotpt, int dotnum, int dotorig,
	       struct rcb_box *rcbbox)

{
  int i,iflag,proc,nprocs,total1,total2;
  double weight,wtmax,wtmin,wtone,tolerance;

  MPI_Comm_rank(MPI_COMM_WORLD,&proc);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

  /* check that total # of dots remained the same */

  MPI_Allreduce(&dotorig,&total1,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&dotnum,&total2,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (total1 != total2) {
    if (proc == 0) 
      printf("ERROR: Points before RCB = %d, Points after RCB = %d\n",
	     total1,total2);
  }
  
  /* check that result is load-balanced within log2(P)*max-wt */

  weight = wtone = 0.0;
  for (i = 0; i < dotnum; i++) {
    weight += dotpt[i].Weight;
    if (dotpt[i].Weight > wtone) wtone = dotpt[i].Weight;
  }

  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&wtone,&tolerance,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  /* i = smallest power-of-2 >= nprocs */
  /* tolerance = largest-single-weight*log2(nprocs) */

  for (i = 0; (nprocs >> i) != 0; i++);
  tolerance = tolerance * i * (1.0 + TINY);

  if (wtmax - wtmin > tolerance) {
    if (proc == 0) 
      printf("ERROR: Load-imbalance > tolerance of %g\n",tolerance);
    MPI_Barrier(MPI_COMM_WORLD);
    if (weight == wtmin) printf("  Proc %d has weight = %g\n",proc,weight);
    if (weight == wtmax) printf("  Proc %d has weight = %g\n",proc,weight);
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  /* check that final set of points is inside RCB box of each proc */
  
  iflag = 0;
  for (i = 0; i < dotnum; i++) {
    if (dotpt[i].X[0] < rcbbox->lo[0] || dotpt[i].X[0] > rcbbox->hi[0] ||
	dotpt[i].X[1] < rcbbox->lo[1] || dotpt[i].X[1] > rcbbox->hi[1] ||
	dotpt[i].X[2] < rcbbox->lo[2] || dotpt[i].X[2] > rcbbox->hi[2])
      iflag++;
  }
  if (iflag > 0) 
    printf("ERROR: %d points are out-of-box on proc %d\n",iflag,proc);
  
  MPI_Barrier(MPI_COMM_WORLD);
    
}


/* RCB statistics */

void rcb_stats(double timetotal, struct rcb_dot *dotpt,
	       int dotnum, double *timers, int *counters,
	       struct rcb_box *rcbbox, int reuse)

{
  int i,iflag,proc,nprocs,sum,min,max;
  double ave,rsum,rmin,rmax;
  double weight,wttot,wtmin,wtmax;

  MPI_Comm_rank(MPI_COMM_WORLD,&proc);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  
  if (proc == 0) printf("RCB total time: %g (secs)\n",timetotal);

  if (proc == 0) printf("RCB Statistics:\n");

  MPI_Barrier(MPI_COMM_WORLD);

  /* distribution info */
  
  for (i = 0, weight = 0.0; i < dotnum; i++) weight += dotpt[i].Weight;
  MPI_Allreduce(&weight,&wttot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (proc == 0) {
    printf(" Total weight of dots = %g\n",wttot);
    printf(" Weight on each proc: ave = %g, max = %g, min = %g\n",
	   wttot/nprocs,wtmax,wtmin);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2) printf("    Proc %d has weight = %g\n",proc,weight);

  for (i = 0, weight = 0.0; i < dotnum; i++) 
    if (dotpt[i].Weight > weight) weight = dotpt[i].Weight;
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  
  if (proc == 0) printf(" Maximum weight of single dot = %g\n",wtmax);

  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2) printf("    Proc %d max weight = %g\n",proc,weight);

  /* counter info */

  MPI_Allreduce(&counters[0],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[0],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[0],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Median iter: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2) 
    printf("    Proc %d median count = %d\n",proc,counters[0]);

  MPI_Allreduce(&counters[1],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[1],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[1],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Send count: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d send count = %d\n",proc,counters[1]);
  
  MPI_Allreduce(&counters[2],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[2],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[2],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Recv count: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d recv count = %d\n",proc,counters[2]);
  
  MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Max dots: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d max dots = %d\n",proc,counters[3]);
  
  MPI_Allreduce(&counters[4],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[4],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[4],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Max memory: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d max memory = %d\n",proc,counters[4]);
  
  if (reuse) {
    MPI_Allreduce(&counters[5],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&counters[5],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&counters[5],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    ave = ((double) sum)/nprocs;
    if (proc == 0) 
      printf(" # of Reuse: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(MPI_COMM_WORLD);
    if (RCB_STATS == 2)
      printf("    Proc %d # of Reuse = %d\n",proc,counters[5]);
  }
  
  MPI_Allreduce(&counters[6],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[6],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[6],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" # of OverAlloc: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d # of OverAlloc = %d\n",proc,counters[6]);

  /* timer info */
  
  MPI_Allreduce(&timers[0],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[0],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[0],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Start-up time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d start-up time = %g\n",proc,timers[0]);
  
  MPI_Allreduce(&timers[1],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[1],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[1],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Pre-median time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d pre-median time = %g\n",proc,timers[1]);
  
  MPI_Allreduce(&timers[2],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[2],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[2],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Median time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d median time = %g\n",proc,timers[2]);
  
  MPI_Allreduce(&timers[3],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[3],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[3],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Comm time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (RCB_STATS == 2)
    printf("    Proc %d comm time = %g\n",proc,timers[3]);
  
  /* RCB boxes for each proc */
  
  if (RCB_STATS == 2) {
    if (proc == 0) printf(" RCB sub-domain boxes:\n");
    for (i = 0; i < 3; i++) {
      MPI_Barrier(MPI_COMM_WORLD);
      if (proc == 0) printf("    Dimension %d\n",i+1);
      MPI_Barrier(MPI_COMM_WORLD);
      printf("      Proc = %d: Box = %g %g\n",
	     proc,rcbbox->lo[i],rcbbox->hi[i]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

}
