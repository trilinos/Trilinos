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
#include "par_const.h"
#include "params_const.h"

#define MYHUGE 1.0e30
#define TINY   1.0e-6

#define RCB_DEFAULT_OVERALLOC 1.0
#define RCB_DEFAULT_REUSE FALSE

static int rcb(LB *, int *, LB_GID **, LB_LID **, int **, double, int, int);

/*  RCB_CHECK = 0  No consistency check on input or results */
/*  RCB_CHECK = 1  Check input weights and final results for consistency */
static int RCB_CHECK = 1;

/*  RCB_STATS = 0  No statistics logging */
/*  RCB_STATS = 1  Log times and counts, print summary */
/*  RCB_STATS = 2  Log times and counts, print for each proc */
static int RCB_STATS = 1;

/*---------------------------------------------------------------------------*/

int LB_RCB_Set_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */
    PARAM_VARS rcb_params[] = {
	{ "RCB_OVERALLOC", NULL, "DOUBLE" },
	{ "RCB_REUSE", NULL, "INT" },
	{ "RCB_WGTFLAG", NULL, "INT" },
	{ NULL, NULL, NULL } };

    status = LB_Check_Param(name, val, rcb_params, &result, &index);

    return(status);
}

/*---------------------------------------------------------------------------*/

int LB_rcb(
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
    /* Wrapper routine to set parameter values and call the real rcb. */
    double overalloc;         /* amount to overallocate by when realloc
                                 of dot array must be done.     
                                 1.0 = no extra; 1.5 = 50% extra; etc. */
    int reuse;                /* (0) don't use (1) use previous cuts
                                 stored in treept at initial guesses.  */
    int wgtflag;              /* (0) do not (1) do use weights.
                                 Multidimensional weights not supported */
    PARAM_VARS rcb_params[] = {
	{ "RCB_OVERALLOC", NULL, "DOUBLE" },
	{ "RCB_REUSE", NULL, "INT" },
	{ "RCB_WGTFLAG", NULL, "INT" },
	{ NULL, NULL, NULL } };
    
    rcb_params[0].ptr = (void *) &overalloc;
    rcb_params[1].ptr = (void *) &reuse;
    rcb_params[2].ptr = (void *) &wgtflag;

    overalloc = RCB_DEFAULT_OVERALLOC;
    reuse = RCB_DEFAULT_REUSE;
    wgtflag = 0;

    LB_Assign_Param_Vals(lb->Params, rcb_params);

    return(rcb(lb, num_import, import_global_ids, import_local_ids,
		 import_procs, overalloc, reuse, wgtflag));
}

/*---------------------------------------------------------------------------*/

static int rcb(
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
  int **import_procs,         /* Returned value:  array of processor IDs for
                                 processors owning the non-local objects in
                                 this processor's new decomposition.         */
  double overalloc,           /* amount to overallocate by when realloc
                                 of dot array must be done.     
                                 1.0 = no extra; 1.5 = 50% extra; etc. */
  int reuse,                  /* (0) don't use (1) use previous cuts
                                 stored in treept at initial guesses.  */
  int wgtflag                 /* (0) do not (1) do use weights.
                                      Multidimensional weights not supported */
)
{
  char    yo[] = "rcb";
  int     proc,nprocs;              /* my proc id, total # of procs */
  struct rcb_dot *dotbuf, *dotpt;   /* local dot arrays */
  struct rcb_box boxtmp;            /* tmp rcb box */
  int     keep, outgoing;           /* message exchange counters */
  int     incoming, incoming2;      /* message exchange counters */
  int     pdotnum;                  /* # of dots - decomposition changes it */
  int     pdottop;                  /* dots >= this index are new */
  int    *dotmark;                  /* which side of median for each dot */
  int     dotnum;                   /* number of dots */
  int     dotmax = 0;               /* max # of dots arrays can hold */
  int     dottop;                   /* dots >= this index are new */
  int     dotnew;                   /* # of new dots after send/recv */
  int     proclower, procupper;     /* lower/upper proc in partition */
  int     procmid;                  /* 1st proc in upper half of part */
  int     procpartner, procpartner2; /* proc(s) to exchange with */
  int     readnumber;               /* # of proc partner(s) to read from */
  int     markactive;               /* which side of cut is active = 0/1 */
  int     dim;                      /* which of 3 axes median cut is on */
  double *coord, *wgts;             /* temp arrays for median_find */
  double  valuehalf;                /* median cut position */
  double  fractionlo;               /* desired wt in lower half */
  int     first_guess;              /* flag if first guess for median search */
  int     allocflag;                /* have to re-allocate space */
  int     length;                   /* message length */
  double  time1,time2,time3,time4;  /* timers */
  double  timestart,timestop;       /* timers */
  double  timers[4];                /* diagnostic timers 
			              0 = start-up time before recursion
				      1 = time before median iterations
				      2 = time in median iterations
				      3 = communication time */
  int     counters[7];              /* diagnostic counts
			              0 = # of median iterations
				      1 = # of dots sent
				      2 = # of dots received
				      3 = most dots this proc ever owns
				      4 = most dot memory this proc ever allocs
				      5 = # of times a previous cut is re-used
				      6 = # of reallocs of dot array */
  int     i,ii,j,k;                 /* local variables */

  RCB_STRUCT *rcb;                 /* Pointer to data structures for RCB.  */
  struct rcb_box *rcbbox;          /* bounding box of final RCB sub-domain */
  struct rcb_tree *treept;         /* tree of RCB cuts - only single cut on 
                                      exit */

  double LB_start_time, LB_end_time;
  double LB_time[2], LB_max_time[2];

  /* function prototypes */

  void rcb_error(LB *, int);
  void rcb_check(LB *, struct rcb_dot *, int, int, struct rcb_box *);
  void rcb_stats(LB *, double, struct rcb_dot *,int, double *, 
		 int *, struct rcb_box *, int);

  /* MPI data types and user functions */

  MPI_Comm local_comm, tmp_comm;
  MPI_Request request, request2;
  MPI_Status status;
  MPI_Op box_op;
  MPI_Datatype box_type;
  MPI_User_function rcb_box_merge;

  if (RCB_STATS) {
    MPI_Barrier(lb->Communicator);
    timestart = time1 = MPI_Wtime();
  }

  /* setup for parallel */

  proc = lb->Proc;
  nprocs = lb->Num_Proc;

  /*
   *  Build the RCB Data structure and 
   *  set pointers to information in it.
   */

  LB_start_time = MPI_Wtime();
  LB_rcb_build_data_structure(lb, &pdotnum, &dotmax, wgtflag);

  rcb = (RCB_STRUCT *) (lb->Data_Structure);

  dotpt  = rcb->Dots; 
  rcbbox = rcb->Box;
  treept = rcb->Tree_Ptr;
  LB_end_time = MPI_Wtime();
  LB_time[0] = LB_end_time - LB_start_time;
  LB_start_time = LB_end_time;

  
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
    dotmark = (int *) LB_Array_Alloc(__FILE__, __LINE__, 1, (unsigned) dotmax,
                                     sizeof(int));
    if (dotmark == NULL)
      rcb_error(lb, dotmax*sizeof(int));
    coord = (double *) LB_Array_Alloc(__FILE__, __LINE__, 1, (unsigned) dotmax,
                                      sizeof(double));
    if (coord == NULL) {
      LB_Free((void **) &dotmark);
      rcb_error(lb, dotmax*sizeof(int));
    }
    wgts = (double *) LB_Array_Alloc(__FILE__, __LINE__, 1, (unsigned) dotmax,
                                     sizeof(double));
    if (wgts == NULL) {
      LB_Free((void **) &dotmark);
      LB_Free((void **) &coord);
      rcb_error(lb, dotmax*sizeof(int));
    }
  }
  else {
    dotmark = NULL;
    coord = NULL;
    wgts = NULL;
  }

  /* create MPI data and function types for box and median */

  MPI_Type_contiguous(6,MPI_DOUBLE,&box_type);
  MPI_Type_commit(&box_type);

  MPI_Op_create(&rcb_box_merge,1,&box_op);

  /* set dot weights = 1.0 if user didn't */

  if (!wgtflag)
    for (i = 0; i < dotnum; i++) dotpt[i].Weight = 1.0;

  /* check that all weights > 0 */

  if (RCB_CHECK) {
    j = 0;
    for (i = 0; i < dotnum; i++) if (dotpt[i].Weight <= 0.0) j++;
    MPI_Allreduce(&j,&k,1,MPI_INT,MPI_SUM,lb->Communicator);
    if (k > 0) {
      if (proc == 0) printf("RCB ERROR: %d dot weights are <= 0\n",k);
      return LB_FATAL;
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

  MPI_Allreduce(&boxtmp,rcbbox,1,box_type,box_op,lb->Communicator);

  /* create local communicator for use in recursion */

  MPI_Comm_dup(lb->Communicator,&local_comm);

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
    
    /* fractionlo = desired fraction of weight in lower half of partition */

    fractionlo = ((double) (procmid - proclower)) / 
      (procupper + 1 - proclower);

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
      LB_Free((void **) &dotmark);
      LB_Free((void **) &coord);
      LB_Free((void **) &wgts);
      dotmark = (int *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                       (unsigned) dotmax, sizeof(int));
      if (dotmark == NULL)
        rcb_error(lb, dotmax*sizeof(int));
      coord = (double *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                        (unsigned) dotmax, sizeof(double));
      if (coord == NULL) {
        LB_Free((void **) &dotmark);
        rcb_error(lb, dotmax*sizeof(int));
      }
      wgts = (double *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                       (unsigned) dotmax, sizeof(double));
      if (wgts == NULL) {
        LB_Free((void **) &dotmark);
        LB_Free((void **) &coord);
        rcb_error(lb, dotmax*sizeof(int));
      }
    }

    /* copy correct coordinate value into the temporary array */
    for (i = 0; i < dotnum; i++) {
      coord[i] = dotpt[i].X[dim];
      wgts[i] = dotpt[i].Weight;
    }

    /* determine if there is a first guess to use */
    if (reuse && dim == treept[procmid].dim) {
      if (RCB_STATS) counters[5]++;
      valuehalf = treept[procmid].cut;
      first_guess = 1;
    }
    else first_guess = 0;

    if (RCB_STATS) time2 = MPI_Wtime();

    if (!LB_find_median(coord, wgts, dotmark, dotnum, proc, fractionlo,
                        local_comm, &valuehalf, first_guess, &(counters[0]))) {
      fprintf(stderr, "[%d] %s: Error returned from find_median\n", proc, yo);
      return LB_FATAL;
    }

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

    MPI_Send(&outgoing,1,MPI_INT,procpartner,0,lb->Communicator);
    incoming = 0;
    if (readnumber) {
      MPI_Recv(&incoming,1,MPI_INT,procpartner,0,lb->Communicator,&status);
      if (readnumber == 2) {
	MPI_Recv(&incoming2,1,MPI_INT,procpartner2,0,lb->Communicator,&status);
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
	LB_REALLOC(dotpt,(unsigned) dotmax * sizeof(struct rcb_dot));
      if (dotpt == NULL) rcb_error(lb, dotmax*sizeof(struct rcb_dot));
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
      dotbuf = (struct rcb_dot *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                                 (unsigned) outgoing,
                                                 sizeof(struct rcb_dot));
      if (dotbuf == NULL) {
        rcb_error(lb, outgoing*sizeof(struct rcb_dot));
      }
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
			   procpartner,1,lb->Communicator,&request);
      if (readnumber == 2) {
	keep += incoming - incoming2;
	length = incoming2 * sizeof(struct rcb_dot);
	MPI_Irecv(&dotpt[keep],length,MPI_CHAR,
			      procpartner2,1,lb->Communicator,&request2);
      }
    }
    
    /* handshake before sending data to insure recvs have been posted */
    
    if (readnumber > 0) {
      MPI_Send(NULL,0,MPI_INT,procpartner,0,lb->Communicator);
      if (readnumber == 2)
	MPI_Send(NULL,0,MPI_INT,procpartner2,0,lb->Communicator);
    }
    MPI_Recv(NULL,0,MPI_INT,procpartner,0,lb->Communicator,&status);

    /* send dot data to partner */

    length = outgoing * sizeof(struct rcb_dot);
    MPI_Rsend(dotbuf,length,MPI_CHAR,procpartner,1,lb->Communicator);
    LB_Free((void **) &dotbuf);
    
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
  MPI_Type_free(&box_type);
  MPI_Op_free(&box_op);

  LB_Free((void **) &coord);
  LB_Free((void **) &wgts);
  LB_Free((void **) &dotmark);

  LB_end_time = MPI_Wtime();
  LB_time[1] = LB_end_time - LB_start_time;

  if (RCB_STATS) {
    MPI_Barrier(lb->Communicator);
    timestop = time1 = MPI_Wtime();
  }

  /* error checking and statistics */

  if (RCB_CHECK) rcb_check(lb, dotpt,dotnum,pdotnum,rcbbox);
  if (RCB_STATS) rcb_stats(lb, timestop-timestart,dotpt,dotnum,
			   timers,counters,rcbbox,reuse);

  /* update calling routine parameters */
  
  LB_start_time = MPI_Wtime();

  pdotnum = dotnum;
  pdottop = dottop;

  /*  build return arguments */

  *num_import = dotnum - dottop;
  if (*num_import > 0) {
    *import_global_ids = (LB_GID *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                                  *num_import, sizeof(LB_GID));
    if (!(*import_global_ids))
      rcb_error(lb, *num_import*sizeof(LB_GID));
    *import_local_ids  = (LB_LID *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                                  *num_import, sizeof(LB_LID));
    if (!(*import_local_ids)) {
      LB_Free((void **) import_global_ids);
      rcb_error(lb, *num_import*sizeof(LB_LID));
    }
    *import_procs      = (int *) LB_Array_Alloc(__FILE__, __LINE__, 1,
                                                  *num_import, sizeof(int));
    if (!(*import_procs)) {
      LB_Free((void **) import_global_ids);
      LB_Free((void **) import_local_ids);
      rcb_error(lb, *num_import*sizeof(int));
    }


    for (i = 0; i < *num_import; i++) {
      ii = i + dottop;
      LB_SET_GID((*import_global_ids)[i], dotpt[ii].Tag.Global_ID);
      LB_SET_LID((*import_local_ids)[i], dotpt[ii].Tag.Local_ID);
      (*import_procs)[i]      = dotpt[ii].Tag.Proc;
    }
  }
  LB_end_time = MPI_Wtime();
  LB_time[0] += (LB_end_time - LB_start_time);

  MPI_Allreduce(LB_time, LB_max_time, 2, MPI_DOUBLE, MPI_MAX,
                lb->Communicator);
  if (lb->Proc == 0) {
    printf("LBLIB RCB Times:  \n");
    printf("LBLIB     Build:  %f\n", LB_max_time[0]);
    printf("LBLIB     RCB:    %f\n", LB_max_time[1]);
  }

  if (lb->Debug > 6) {
    int i;
    LB_print_sync_start(lb, TRUE);
    printf("LBLIB RCB Proc %d  Num_Obj=%d  Num_Keep=%d  Num_Non_Local=%d\n", 
           lb->Proc, pdotnum, pdottop, *num_import);
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

  /* Temporary return value until error codes are fully implemented */
  return(LB_OK);  
}



/* ----------------------------------------------------------------------- */

/* error message for malloc/realloc overflow */

void rcb_error(LB *lb, int size)

{
  int proc;

  MPI_Comm_rank(lb->Communicator,&proc);
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

/* consistency checks on RCB results */

void rcb_check(LB *lb, struct rcb_dot *dotpt, int dotnum, int dotorig,
	       struct rcb_box *rcbbox)

{
  int i,iflag,proc,nprocs,total1,total2;
  double weight,wtmax,wtmin,wtone,tolerance;

  MPI_Comm_rank(lb->Communicator,&proc);
  MPI_Comm_size(lb->Communicator,&nprocs);

  /* check that total # of dots remained the same */

  MPI_Allreduce(&dotorig,&total1,1,MPI_INT,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&dotnum,&total2,1,MPI_INT,MPI_SUM,lb->Communicator);
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

  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  MPI_Allreduce(&wtone,&tolerance,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);

  /* i = smallest power-of-2 >= nprocs */
  /* tolerance = largest-single-weight*log2(nprocs) */

  for (i = 0; (nprocs >> i) != 0; i++);
  tolerance = tolerance * i * (1.0 + TINY);

  if (wtmax - wtmin > tolerance) {
    if (proc == 0) 
      printf("ERROR: Load-imbalance > tolerance of %g\n",tolerance);
    MPI_Barrier(lb->Communicator);
    if (weight == wtmin) printf("  Proc %d has weight = %g\n",proc,weight);
    if (weight == wtmax) printf("  Proc %d has weight = %g\n",proc,weight);
  }
  
  MPI_Barrier(lb->Communicator);
  
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
  
  MPI_Barrier(lb->Communicator);
    
}


/* RCB statistics */

void rcb_stats(LB *lb, double timetotal, struct rcb_dot *dotpt,
	       int dotnum, double *timers, int *counters,
	       struct rcb_box *rcbbox, int reuse)

{
  int i,proc,nprocs,sum,min,max;
  double ave,rsum,rmin,rmax;
  double weight,wttot,wtmin,wtmax;

  MPI_Comm_rank(lb->Communicator,&proc);
  MPI_Comm_size(lb->Communicator,&nprocs);
  
  if (proc == 0) printf("RCB total time: %g (secs)\n",timetotal);

  if (proc == 0) printf("RCB Statistics:\n");

  MPI_Barrier(lb->Communicator);

  /* distribution info */
  
  for (i = 0, weight = 0.0; i < dotnum; i++) weight += dotpt[i].Weight;
  MPI_Allreduce(&weight,&wttot,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);

  if (proc == 0) {
    printf(" Total weight of dots = %g\n",wttot);
    printf(" Weight on each proc: ave = %g, max = %g, min = %g\n",
	   wttot/nprocs,wtmax,wtmin);
  }

  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2) printf("    Proc %d has weight = %g\n",proc,weight);

  for (i = 0, weight = 0.0; i < dotnum; i++) 
    if (dotpt[i].Weight > weight) weight = dotpt[i].Weight;
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  
  if (proc == 0) printf(" Maximum weight of single dot = %g\n",wtmax);

  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2) printf("    Proc %d max weight = %g\n",proc,weight);

  /* counter info */

  MPI_Allreduce(&counters[0],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&counters[0],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&counters[0],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Median iter: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2) 
    printf("    Proc %d median count = %d\n",proc,counters[0]);

  MPI_Allreduce(&counters[1],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&counters[1],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&counters[1],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Send count: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d send count = %d\n",proc,counters[1]);
  
  MPI_Allreduce(&counters[2],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&counters[2],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&counters[2],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Recv count: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d recv count = %d\n",proc,counters[2]);
  
  MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Max dots: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d max dots = %d\n",proc,counters[3]);
  
  MPI_Allreduce(&counters[4],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&counters[4],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&counters[4],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Max memory: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d max memory = %d\n",proc,counters[4]);
  
  if (reuse) {
    MPI_Allreduce(&counters[5],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
    MPI_Allreduce(&counters[5],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
    MPI_Allreduce(&counters[5],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == 0) 
      printf(" # of Reuse: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(lb->Communicator);
    if (RCB_STATS == 2)
      printf("    Proc %d # of Reuse = %d\n",proc,counters[5]);
  }
  
  MPI_Allreduce(&counters[6],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&counters[6],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&counters[6],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" # of OverAlloc: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d # of OverAlloc = %d\n",proc,counters[6]);

  /* timer info */
  
  MPI_Allreduce(&timers[0],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&timers[0],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&timers[0],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Start-up time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d start-up time = %g\n",proc,timers[0]);
  
  MPI_Allreduce(&timers[1],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&timers[1],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&timers[1],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Pre-median time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d pre-median time = %g\n",proc,timers[1]);
  
  MPI_Allreduce(&timers[2],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&timers[2],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&timers[2],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Median time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d median time = %g\n",proc,timers[2]);
  
  MPI_Allreduce(&timers[3],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&timers[3],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&timers[3],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Comm time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(lb->Communicator);
  if (RCB_STATS == 2)
    printf("    Proc %d comm time = %g\n",proc,timers[3]);
  
  /* RCB boxes for each proc */
  
  if (RCB_STATS == 2) {
    if (proc == 0) printf(" RCB sub-domain boxes:\n");
    for (i = 0; i < 3; i++) {
      MPI_Barrier(lb->Communicator);
      if (proc == 0) printf("    Dimension %d\n",i+1);
      MPI_Barrier(lb->Communicator);
      printf("      Proc = %d: Box = %g %g\n",
	     proc,rcbbox->lo[i],rcbbox->hi[i]);
    }
  }

  MPI_Barrier(lb->Communicator);

}
