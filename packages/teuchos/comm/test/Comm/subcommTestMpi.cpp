// MPI-only version of subcommTestTeuchosComm.cpp

#include <stdio.h>
#include <mpi.h>


int main(int narg, char **arg)
{
  MPI_Init(&narg, &arg);
  MPI_Comm comm = MPI_COMM_WORLD;

  int me, np;
  MPI_Comm_rank(comm, &me);
  MPI_Comm_size(comm, &np);

  int niter = 4;
  int *ids = NULL;
  ids = new int[np/2+1];
  ids[0] = me;
  for (int i = 1; i < np/2+1; i++) {
    ids[i] = (i != me ? i : 0);
  }
  for (int i = 0; i < niter; i++) {
    MPI_Comm a;
    MPI_Group cgrp;
    MPI_Comm_group(comm, &cgrp);
    MPI_Group_incl(cgrp, np/2+1, ids, &agrp)
    MPI_Comm_create(comm, agrp, a);

    int anp;
    MPI_Comm_size(a, &anp);
    printf("Iteration %d:  New comm has %d ranks\n", i, anp);

    MPI_Group_free(&agrp);
    MPI_Group_free(&cgrp);
    MPI_Comm_free(&a);
  }
  delete [] ids;
  if (me == 0)
    printf("\nPASS\n");

  return 0;
}
