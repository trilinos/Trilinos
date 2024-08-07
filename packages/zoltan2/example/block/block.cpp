// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file block.cpp
 *  \brief An example of partitioning global ids with Block.
 */

#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

/*! \example block.cpp
 *  An example of the use of the Block algorithm to partition data.
 *  \todo error handling
 *  \todo write some examples that don't use teuchos
 */

int main(int argc, char *argv[]) {
#ifdef HAVE_ZOLTAN2_MPI
  MPI_Init(&argc, &argv);
  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int rank=0, nprocs=1;
#endif

  // For convenience, we'll use the Tpetra defaults for local/global ID types
  // Users can substitute their preferred local/global ID types
  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;
  typedef Tpetra::Details::DefaultTypes::scalar_type scalar_t;

  ///////////////////////////////////////////////////////////////////////
  // Generate some input data.

  int localCount = 40*(rank+1);
  int totalCount = 20*nprocs*(nprocs+1);
  int targetCount = totalCount / nprocs;
  globalId_t *globalIds = new globalId_t[localCount];

  if (rank==0) {
    for (int i=0, num=40; i < nprocs ; i++, num+=40) {
      std::cout << "Rank " << i << " generates " << num << " ids." << std::endl;
    }
  }

  globalId_t offset = 0;
  for (int i=1; i <= rank; i++) {
    offset += 40*i;
  }

  for (int i=0; i < localCount; i++) {
    globalIds[i] = offset++;
  }

  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 input adapter with no weights

  // TODO explain
  typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;

  // TODO explain
  typedef Zoltan2::BasicIdentifierAdapter<myTypes> inputAdapter_t;

  std::vector<const scalar_t *> noWeights;
  std::vector<int> noStrides;

  inputAdapter_t ia(localCount, globalIds, noWeights, noStrides);

  ///////////////////////////////////////////////////////////////////////
  // Create parameters for a Block problem

  Teuchos::ParameterList params("test params");
  params.set("debug_level", "basic_status");
  params.set("debug_procs", "0");
  params.set("error_check_level", "debug_mode_assertions");

  params.set("algorithm", "block");
  params.set("imbalance_tolerance", 1.1);
  params.set("num_global_parts", nprocs);

  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 partitioning problem

  Zoltan2::PartitioningProblem<inputAdapter_t> *problem = 
      new Zoltan2::PartitioningProblem<inputAdapter_t>(&ia, &params);

  ///////////////////////////////////////////////////////////////////////
  // Solve the problem - do the partitioning

  problem->solve();

  ///////////////////////////////////////////////////////////////////////
  // Check and print the solution.
  // Count number of IDs assigned to each part; compare to targetCount

  const globalId_t *ids = NULL;
  ia.getIDsView(ids);
  std::vector<int> partCounts(nprocs, 0), globalPartCounts(nprocs, 0);

  for (size_t i = 0; i < ia.getLocalNumIDs(); i++) {
    int pp = problem->getSolution().getPartListView()[i];
    std::cout << rank << " LID " << i << " GID " << ids[i]
              << " PART " << pp << std::endl;
    partCounts[pp]++;
  }

#ifdef HAVE_ZOLTAN2_MPI
  MPI_Allreduce(&(partCounts[0]), &(globalPartCounts[0]), nprocs,
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  for (int i = 0; i < nprocs; i++) globalPartCounts[i] = partCounts[i];
#endif

  if (rank == 0) {
    int ierr = 0;
    for (int i = 0; i < nprocs; i++) {
      if (globalPartCounts[i] != targetCount) {
        std::cout << "FAIL: part " << i << " has " << globalPartCounts[i]
                  << " != " << targetCount << "; " << ++ierr << " errors"
                  << std::endl;
      }
    }
    if (ierr == 0) {
      std::cout << "PASS" << std::endl;
    }
  }

  delete [] globalIds;
  delete problem;
#ifdef HAVE_ZOLTAN2_MPI
  MPI_Finalize();
#endif
}
