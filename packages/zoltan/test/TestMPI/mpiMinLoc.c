// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* Zoltan's gold-standard answers were generated with MPI that returns the
 * lowest numbered rank in MPI_Allreduce with MPI_MINLOC when the compared
 * values are equal on several processors.  
 *
 * The MPI standard indicates that MPI_MINLOC should return the lowest
 * numbered rank; see, for example, 
 * https://www.mpi-forum.org/docs/mpi-1.1/mpi-11-html/node79.html
 *
 * Some implementations of MPI (e.g., Spectrum MPI) incorrectly return
 * the largest rank with the minimum value.
 *
 * Some Zoltan PHG tests will fail if the lowest rank is not returned.
 * Zoltan_PHG_CoarsePartition tries several coarse partitions and selects
 * the best wrt hyperedge cut.  When two processors compute equally good
 * but different coarse partitions, the choice of a different coarse
 * partition leads to a different final partition.  The final partition is
 * correct, but it doesn't match the gold-standard answers and, thus, causes
 * Zoltan tests to fail.  
 *
 * See, for example, https://github.com/trilinos/Trilinos/issues/8798
 *
 * This test detects the condition needed in MINLOC to match the gold-standard
 * answers.  If this test fails, Zoltan will operate correctly and produce
 * correct partitions, but the Zoltan tests may report failures.
 */

int main(int narg, char **arg)
{
  struct {
    float val;
    int rank;
  } local, global;
  int me;

  MPI_Init(&narg, &arg);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  local.val = 1.0;
  local.rank = me;

  MPI_Allreduce(&local, &global, 1, MPI_FLOAT_INT, MPI_MINLOC, MPI_COMM_WORLD);

  printf("%d MINLOC value %f on rank %d\n\n", me, global.val, global.rank);

  if (global.rank == 0) printf("PASSED\n\n");
  else printf("FAILED\n\n");

  MPI_Finalize();
  return 0;
}
