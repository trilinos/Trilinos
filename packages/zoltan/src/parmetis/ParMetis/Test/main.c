/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * main.c
 * 
 * This file contains code for testing teh adaptive partitioning routines
 *
 * Started 5/19/97
 * George
 *
 * $Id$
 *
 */

#include <parmetis.h>


/*************************************************************************
* Let the game begin
**************************************************************************/
main(int argc, char *argv[])
{
  int mype, npes;
  MPI_Comm comm;

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (argc != 2) {
    if (mype == 0)
      printf("Usage: %s <graph-file>\n", argv[0]);

    MPI_Finalize();
    exit(0);
  }

  TestParMetis(argv[1], comm); 

  MPI_Comm_free(&comm);

  MPI_Finalize();
}


