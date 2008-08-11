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
 *    Revision: 1.56.2.3 $
 ****************************************************************************/


#include "mpi.h"     /* Zoltan requires MPI. */
#include <stdio.h> 
#include <stdlib.h> 
#include "zoltan.h"  /* always include this header for for Zoltan */
#include "matrix.h"  /* data structure for this example only */

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/* Simple Zoltan example.
 *
 * We set up a sparse matrix, and load-balance
 * with respect to the nonzero entries. 
 */


int main(int argc, char *argv[])
{
  /* Local declarations. */
  struct Zoltan_Struct *zz;
  float  version;
  int    k, col; 
  int    myrank, nproc;
  int    max_nnz, error;
  Matrix A; /* a distributed sparse matrix */


  /* Fire up MPI. Zoltan will do this if you don't. */
  MPI_Init(&argc, &argv);

  /* Get the MPI rank. We use this to set up our local data. */
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  /* Set up a sparse matrix structure. 
     We store it in the (i,j,v) coordinate format.
     See matrix.h for details.
     Each processor owns a single row in this example.
   */

  A.num_globalrows = nproc;
  A.num_globalcols = 20;
  /* Allocate memory for nonzeros. */
  max_nnz = 10;
  A.entries = (Matrix_entry *) malloc(max_nnz * sizeof(Matrix_entry));

  /* Insert a row of data into the matrix. */
  k = 0;
  for (col=0; col<A.num_globalcols; col++){ 
    /* Generate random nonzeros in my row. */
    if (rand()%3){
      A.entries[k].i = myrank; /* my row is my proc rank */
      A.entries[k].j = col;    /* Add an entry in column col */
      A.entries[k].val = 1.0;
      k++;
    }
    if (k==max_nnz) break; /* only allow max_nnz per row */
  }
  A.num_mynz = k;

  /* 
   * Print matrix before load balancing.
   */
  printf("Before load balancing:\n");
  for (k=0; k<A.num_mynz; k++){
    printf("[Proc %1d] %d %d %g\n", myrank,
      A.entries[k].i, A.entries[k].j, A.entries[k].val);
  }


  /*  Initialize Zoltan. It will start MPI if we haven't already. */
  /*  Do this only once. */

  if ((error = Zoltan_Initialize(argc, argv, &version)) != ZOLTAN_OK) {
    printf("fatal: Zoltan_Initialize returned error code, %d", error);
    goto End;
  }

  /*
   *  Create a Zoltan structure.
   */
  if ((zz = Zoltan_Create(MPI_COMM_WORLD)) == NULL) {
    printf("fatal:  NULL returned from Zoltan_Create()\n");
    goto End;
  }

  /*
   *  Tell Zoltan what kind of local/global IDs we will use. 
   *  In our case, each GID is two ints and there are no local ids. 
   *  One can skip this step if the IDs are just single ints.
   */
  Zoltan_Set_Param(zz, "num_gid_entries", "2");
  Zoltan_Set_Param(zz, "num_lid_entries", "0");

  /*
   *  Set up Zoltan query functions for our Matrix data structure.
   */
  if (!setup_zoltan(zz, A)) {
    printf("fatal: Error returned from setup_zoltan\n");
    goto End;
  }


  /*
   * Run Zoltan to compute a new load balance.
   * Data migration may also happen here.
   */
  if (!run_zoltan(zz, A)) {
    printf("fatal: Error returned from run_zoltan\n");
    goto End;
  }

#if 0 /* No migration yet, so nothing has changed! */
  /* 
   * Print matrix after load balancing.
   */
  printf("After load balancing:\n");
  for (k=0; k<A.num_mynz; k++){
    printf("[Proc %1d] %d %d %g\n", myrank,
      A.entries[k].i, A.entries[k].j, A.entries[k].val);
  }
#endif

End:
  /* Destroy Zoltan structure */
  Zoltan_Destroy(&zz);

  /* End MPI. */
  MPI_Finalize();

  return 0;
}

