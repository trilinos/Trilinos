/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include "zz_const.h"
#include "octant_const.h"
#include "octupdate_const.h"

/* void Zoltan_Oct_print_stats()
 *
 * Prints out statistic on the octree load balancing partitioner 
 */
void Zoltan_Oct_print_stats(ZZ *zz, double timetotal, double *timers, int *counters,
                        float *c, int STATS_TYPE)
{
  int proc,                                /* the processor number */
      print_proc,                          /* the processor that prints */
      nprocs,                              /* total number of processors */
      sum,                                 /* the sum of the counters */
      min,                                 /* minimum value */
      max;                                 /* maximum value */
  float sum1,
        min1,
        max1;
  double ave,                              /* average of the timers */
         rsum,                             /* the sum of the timers */
         rmin,                             /* minimum timer value */
         rmax;                             /* maximum timer value */

  /* get the processor number and the total number of processors */
  MPI_Comm_rank(zz->Communicator, &proc);
  MPI_Comm_size(zz->Communicator, &nprocs);
  print_proc = zz->Debug_Proc;
  
  /* print out some headders */
  if (proc == print_proc) printf("Total time: %g (secs)\n",timetotal);
  if (proc == print_proc) printf("Statistics:\n");

  MPI_Barrier(zz->Communicator);

  /* counter info */
  MPI_Allreduce(&counters[0],&sum,1,MPI_INT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&counters[0],&min,1,MPI_INT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&counters[0],&max,1,MPI_INT,MPI_MAX,zz->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == print_proc) 
    printf(" Number of partition iters: ave = %g, min = %d, max = %d\n", 
	   ave, min, max);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2) 
    printf("    Proc %d iteration count = %d\n", proc, counters[0]);

  MPI_Allreduce(&counters[1],&sum,1,MPI_INT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&counters[1],&min,1,MPI_INT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&counters[1],&max,1,MPI_INT,MPI_MAX,zz->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == print_proc) 
    printf(" Objs sent during gen tree: ave = %g, min = %d, max = %d\n",
	   ave,min,max);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2)
    printf("    Proc %d send count = %d\n",proc,counters[1]);
  
  MPI_Allreduce(&counters[2],&sum,1,MPI_INT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&counters[2],&min,1,MPI_INT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&counters[2],&max,1,MPI_INT,MPI_MAX,zz->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == print_proc) 
    printf(" Objs recv during gen tree: ave = %g, min = %d, max = %d\n",
	   ave,min,max);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2)
    printf("    Proc %d recv count = %d\n",proc,counters[2]);
  
  MPI_Allreduce(&counters[4],&sum,1,MPI_INT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&counters[4],&min,1,MPI_INT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&counters[4],&max,1,MPI_INT,MPI_MAX,zz->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == print_proc) 
    printf(" Objs sent during balancing: ave = %g, min = %d, max = %d\n",
	   ave,min,max);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2)
    printf("    Proc %d send count = %d\n",proc,counters[4]);
  
  MPI_Allreduce(&counters[5],&sum,1,MPI_INT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&counters[5],&min,1,MPI_INT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&counters[5],&max,1,MPI_INT,MPI_MAX,zz->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == print_proc) 
    printf(" Objs recv during balancing: ave = %g, min = %d, max = %d\n",
	   ave,min,max);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2)
    printf("    Proc %d recv count = %d\n",proc,counters[5]);
  
  MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,zz->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == print_proc) 
    printf(" Max objs: ave = %g, min = %d, max = %d\n",ave,min,max);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2)
    printf("    Proc %d max objs = %d\n",proc,counters[3]);

  counters[3] += (counters[2] - counters[1]);
  MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,zz->Communicator);
  ave = ((double) sum)/nprocs;
  if (proc == print_proc) 
    printf(" Max objs on Proc: ave = %g, min = %d, max = %d\n",ave,min,max);

  MPI_Allreduce(&c[0],&sum1,1,MPI_FLOAT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&c[0],&min1,1,MPI_FLOAT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&c[0],&max1,1,MPI_FLOAT,MPI_MAX,zz->Communicator);
  ave = ((double) sum1)/nprocs;
  if (proc == print_proc) 
    printf(" Initial Load: ave = %g, min = %f, max = %f\n",ave,min1,max1);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2)
    printf("    Proc %d intial load = %f\n",proc,c[0]);
  
  MPI_Allreduce(&c[1],&sum1,1,MPI_FLOAT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&c[1],&min1,1,MPI_FLOAT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&c[1],&max1,1,MPI_FLOAT,MPI_MAX,zz->Communicator);
  ave = ((double) sum1)/nprocs;
  if (proc == print_proc) 
    printf(" Load Before Balancing: ave = %g, min = %f, max = %f\n", 
	   ave, min1, max1);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2)
    printf("    Proc %d load before balancing = %f\n", proc, c[1]);
  
  c[3] += (c[1] - c[2]);
  MPI_Allreduce(&c[3],&sum1,1,MPI_FLOAT,MPI_SUM,zz->Communicator);
  MPI_Allreduce(&c[3],&min1,1,MPI_FLOAT,MPI_MIN,zz->Communicator);
  MPI_Allreduce(&c[3],&max1,1,MPI_FLOAT,MPI_MAX,zz->Communicator);
  ave = ((double) sum1)/nprocs;
  if (proc == print_proc) 
    printf(" Load After Balancing: ave = %g, min = %f, max = %f\n", 
	   ave, min1, max1);
  MPI_Barrier(zz->Communicator);
  if (STATS_TYPE == 2)
    printf("    Proc %d load after balancing = %f\n", proc, c[3]);

  /* timer info */
  if (timetotal>0){
    MPI_Allreduce(&timers[0],&rsum,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&timers[0],&rmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&timers[0],&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
    ave = rsum/nprocs;
    if (proc == print_proc) 
      printf(" Start-up time %%: ave = %g, min = %g, max = %g\n",
  	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
    MPI_Barrier(zz->Communicator);
    if (STATS_TYPE == 2)
      printf("    Proc %d start-up time = %g\n", proc, timers[0]);
    
    MPI_Allreduce(&timers[1],&rsum,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&timers[1],&rmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&timers[1],&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
    ave = rsum/nprocs;
    if (proc == print_proc) 
      printf(" Partition time %%: ave = %g, min = %g, max = %g\n",
  	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
    MPI_Barrier(zz->Communicator);
    if (STATS_TYPE == 2)
      printf("    Proc %d partition time = %g\n",proc, timers[1]);
    
    MPI_Allreduce(&timers[2],&rsum,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&timers[2],&rmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&timers[2],&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
    ave = rsum/nprocs;
    if (proc == print_proc) 
      printf(" Migration notice time %%: ave = %g, min = %g, max = %g\n",
  	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
    MPI_Barrier(zz->Communicator);
    if (STATS_TYPE == 2)
      printf("    Proc %d migration notice time = %g\n",proc,timers[2]);
  
#if 0
    /* ATTN: I don't time the communication */  
    MPI_Allreduce(&timers[3],&rsum,1,MPI_DOUBLE,MPI_SUM,zz->Communicator);
    MPI_Allreduce(&timers[3],&rmin,1,MPI_DOUBLE,MPI_MIN,zz->Communicator);
    MPI_Allreduce(&timers[3],&rmax,1,MPI_DOUBLE,MPI_MAX,zz->Communicator);
    ave = rsum/nprocs;
    if (proc == print_proc) 
      printf(" Comm time %%: ave = %g, min = %g, max = %g\n",
  	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
    MPI_Barrier(zz->Communicator);
    if (STATS_TYPE == 2)
      printf("    Proc %d comm time = %g\n",proc,timers[3]);
#endif
  }
  
  MPI_Barrier(zz->Communicator);

}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
