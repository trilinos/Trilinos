// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/* MPI distributions in some version of Ubuntu and FreeBSD do not handle
 * self messages correctly.  Errors then appear in Zoltan in the Comm package,
 * or in PHG's building of the hypergraph.  If this test fails, there is 
 * something wrong with MPI, so one cannot expect Zoltan to work correctly.
 * KDD 9/2018 -- workaround for Xyce with FreeBSD and bad OpenMPI installation
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int narg, char**arg)
{
  MPI_Init(&narg, &arg);
  int nprocs, my_proc;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  int nrecvs = nprocs;
  int nsends = nprocs;
  int i;
  
  int *vals_from = (int *) malloc(nrecvs*sizeof(int));
  for (i = 0; i < nrecvs; i++) vals_from[i] = -1000-i;

  int *vals_to = (int *) malloc(nsends*sizeof(int));
  for (i = 0; i < nsends; i++) vals_to[i] = (my_proc+1)*1000+i;

  int *procs_from = (int *) malloc(nrecvs*sizeof(int));
  for (i = 0; i < nrecvs; i++) procs_from[i] = -1;

  int *procs_to = (int *) malloc(nsends*sizeof(int));
  for (i = 0; i < nsends; i++) procs_to[i] = i;

  MPI_Request *req = (MPI_Request *)malloc(nrecvs*sizeof(MPI_Request));

  int tag = 30000;
  
  for (i=0; i < nrecvs; i++){
    printf("%d posting receive %d %p\n", my_proc, i, (void *)(vals_from+i));
    MPI_Irecv(vals_from + i, 1, MPI_INT, MPI_ANY_SOURCE, tag, comm, req + i);
  }
  
  for (i=0; i < nsends; i++){
    printf("%d sending to %d value %d\n", my_proc, procs_to[i], vals_to[i]);
    MPI_Send(vals_to + i, 1, MPI_INT, procs_to[i], tag, comm);
  }
  
  MPI_Status status;
  for (i=0; i < nrecvs; i++){
    MPI_Wait(req + i, &status);
    procs_from[i] = status.MPI_SOURCE;
    int ucount;
    MPI_Get_count(&status, MPI_INT, &ucount);
    printf("%d wait source %d count %d \n", 
           my_proc, status.MPI_SOURCE, ucount);
  }

  for (i = 0; i < nrecvs; i++)
    printf("%d procs_from %d vals_from %d %s \n",
           my_proc, procs_from[i], vals_from[i], 
           (vals_from[i] < 0 ? "FAIL FAIL FAIL" : " "));

  MPI_Comm_free(&comm);
  MPI_Finalize();
  return 0;
}
