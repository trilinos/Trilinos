#include <stdio.h>
#include <mpi.h>
#include <lb_const.h>
#include "all_allo_const.h"

/* void print_stats()
 *
 * Prints out statistic on the octree load balancing partitioner 
 */
void print_stats(double timetotal, double *timers, int *counters, 
		 int STATS_TYPE)
{
  LB_ID *obj_ids;                          /* pointer to all the objects ids */
  int i,                                   /* index counter */
      proc,                                /* the processor number */
      nprocs,                              /* total number of processors */
      sum,                                 /* the sum of the counters */
      min,                                 /* minimum value */
      max;                                 /* maximum value */
  double ave,                              /* average of the timers */
         rsum,                             /* the sum of the timers */
         rmin,                             /* minimum timer value */
         rmax;                             /* maximum timer value */
  double weight,                           /* weight of on object */
         wttot,                            /* total weight of all objects */
         wtmin,                            /* minimum weight */
         wtmax,                            /* maximum weight */
         mweight,                          /* max weight of an object */
         tweight;                          /* total weight on this processor */

  /* get the processor number and the total number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  /* print out some headders */
  if (proc == 0) printf("Total time: %g (secs)\n",timetotal);
  if (proc == 0) printf("Statistics:\n");

  MPI_Barrier(MPI_COMM_WORLD);

#if 0
{
  int numobj;
  LB *lb;

  obj_ids = (LB_ID *) LB_array_alloc(__FILE__, __LINE__, 1, numobj,
                                     sizeof(LB_ID));
  lb->Get_All_Local_Objs(lb->Object_Type, obj_ids);
  /* need to get weights of all the objects */
  weight = mweight = tweight = 0.0;
  if(lb->Get_Obj_Weight != NULL)
    for (i = 0; i < numobj; i++) {
      weight = lb->Get_Obj_Weight(obj_ids[i], lb->Object_Type);
      if (weight > mweight)
	mweight = weight;
      tweight += weight;
    }
  else {
    tweight = numobj;                    /* uniform weight, normalized to 1 */
    mweight = 1;
  }

  /* get information on the distribution of the objects */
  MPI_Allreduce(&tweight,&wttot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&tweight,&wtmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&tweight,&wtmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (proc == 0) {
    printf(" Total weight of Objects = %g\n",wttot);
    printf(" Weight on each proc: ave = %g, max = %g, min = %g\n",
	   wttot/nprocs, wtmax, wtmin);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)                       /* print out individual weights */
    printf("    Proc %d has weight = %g\n", proc, tweight);

  MPI_Allreduce(&mweight,&wtmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (proc == 0) printf(" Maximum weight of single object = %g\n", wtmax);

  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d max weight of a sigle object = %g\n", proc, mweight);
}
#endif

  /* counter info */
  MPI_Allreduce(&counters[0],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[0],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[0],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Number of parition iters: ave = %g, min = %d, max = %d\n", ave, min, max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2) 
    printf("    Proc %d iteration count = %d\n", proc, counters[0]);
  /* ATTN: Will, all load balancing schemes be able to use these counters? */

  MPI_Allreduce(&counters[1],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[1],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[1],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Send count: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d send count = %d\n",proc,counters[1]);
  
  MPI_Allreduce(&counters[2],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[2],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[2],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Recv count: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d recv count = %d\n",proc,counters[2]);
  
  MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Max objs: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d max objs = %d\n",proc,counters[3]);

#if 0  
  MPI_Allreduce(&counters[4],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[4],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[4],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" Max memory: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d max memory = %d\n",proc,counters[4]);
  
  /* ATTN: Not sure what this "reuse" is, used in the RCB code */
  if (0) {
    MPI_Allreduce(&counters[5],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&counters[5],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&counters[5],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    ave = ((double) sum)/nprocs;
    if (proc == 0) 
      printf(" # of Reuse: ave = %g, min = %d, max = %d\n", ave, min, max);
    MPI_Barrier(MPI_COMM_WORLD);
    if (STATS_TYPE == 2)
      printf("    Proc %d # of Reuse = %d\n", proc, counters[5]);
  }
  
  if(counters[4] > counters[3])
    counters[6] = counters[4] - counters[3];

  MPI_Allreduce(&counters[6],&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[6],&min,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&counters[6],&max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  ave = ((double) sum)/nprocs;
  if (proc == 0) 
    printf(" # of OverAlloc: ave = %g, min = %d, max = %d\n", ave, min, max);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d # of OverAlloc = %d\n", proc, counters[6]);
#endif

  /* timer info */
  MPI_Allreduce(&timers[0],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[0],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[0],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Start-up time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d start-up time = %g\n", proc, timers[0]);
  
  MPI_Allreduce(&timers[1],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[1],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[1],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Partition time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d partition time = %g\n",proc, timers[1]);
  
  MPI_Allreduce(&timers[2],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[2],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[2],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Migration notice time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d migration notice time = %g\n",proc,timers[2]);


#if 0
  /* ATTN: I don't time the communication */  
  MPI_Allreduce(&timers[3],&rsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[3],&rmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&timers[3],&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  ave = rsum/nprocs;
  if (proc == 0) 
    printf(" Comm time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(MPI_COMM_WORLD);
  if (STATS_TYPE == 2)
    printf("    Proc %d comm time = %g\n",proc,timers[3]);
#endif
  
  /* RCB boxes for each proc */
  /* ATTN: Not sure this can be used generally.... */
#if 0
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
#endif

  MPI_Barrier(MPI_COMM_WORLD);

}
