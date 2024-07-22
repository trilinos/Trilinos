// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file kokkosBlock.cpp
 *  \brief An example of partitioning global ids with Block.
 */

#include <Kokkos_Core.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Zoltan2_BasicKokkosIdentifierAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

using Teuchos::Comm;
using Teuchos::RCP;

/*! \example kokkosBlock.cpp
 *  An example of the use of the Block algorithm to partition data.
 *  \todo error handling
 *  \todo write some examples that don't use teuchos
 */

int main(int narg, char *arg[]) {
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  // For convenience, we'll use the Tpetra defaults for local/global ID types
  // Users can substitute their preferred local/global ID types
  typedef Tpetra::Map<> Map_t;
  typedef Map_t::local_ordinal_type localId_t;
  typedef Map_t::global_ordinal_type globalId_t;
  typedef Tpetra::Details::DefaultTypes::scalar_type scalar_t;
  typedef Tpetra::Map<>::node_type node_t;

  ///////////////////////////////////////////////////////////////////////
  // Generate some input data.

  int localCount = 40 * (rank + 1);
  int totalCount = 20 * nprocs * (nprocs + 1);
  int targetCount = totalCount / nprocs;

  Kokkos::View<globalId_t*, typename node_t::device_type>
    globalIds(Kokkos::ViewAllocateWithoutInitializing("globalIds"), localCount);
  auto host_globalIds = Kokkos::create_mirror_view(globalIds);

  if (rank == 0) {
    for (int i = 0, num = 40; i < nprocs ; i++, num += 40) {
      std::cout << "Rank " << i << " generates " << num << " ids." << std::endl;
    }
  }

  globalId_t offset = 0;
  for (int i = 1; i <= rank; i++) {
    offset += 40 * i;
  }

  for (int i = 0; i < localCount; i++) {
    host_globalIds(i) = offset++;
  }
  Kokkos::deep_copy(globalIds, host_globalIds);

  ///////////////////////////////////////////////////////////////////////
  // Create a Zoltan2 input adapter with no weights

  typedef Zoltan2::BasicUserTypes<scalar_t, localId_t, globalId_t> myTypes;
  typedef Zoltan2::BasicKokkosIdentifierAdapter<myTypes> inputAdapter_t;

  const int nWeights = 1;
  Kokkos::View<scalar_t **, typename node_t::device_type>
    weights("weights", localCount, nWeights);
  auto host_weights = Kokkos::create_mirror_view(weights);

  for (int index = 0; index < localCount; index++) {
    host_weights(index, 0) = 1; // Error check relies on uniform weights
  }

  Kokkos::deep_copy(weights, host_weights);

  inputAdapter_t ia(globalIds, weights);

  /////////////////////////////////////////////////////////////////////////
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

  Kokkos::View<const globalId_t *, typename node_t::device_type> ids;
  ia.getIDsKokkosView(ids);

  auto host_ids = Kokkos::create_mirror_view(ids);
  Kokkos::deep_copy(host_ids, ids);

  Kokkos::View<int*, Kokkos::HostSpace> partCounts("partCounts", nprocs);

  Kokkos::View<int*, Kokkos::HostSpace>
    globalPartCounts("globalPartCounts", nprocs);

  for (size_t i = 0; i < ia.getLocalNumIDs(); i++) {
    int pp = problem->getSolution().getPartListView()[i];
    std::cout << rank << " LID " << i << " GID " << host_ids(i)
              << " PART " << pp << std::endl;
    partCounts(pp)++;
  }

  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, nprocs,
      partCounts.data(), globalPartCounts.data());

  if (rank == 0) {
    int ierr = 0;
    for (int i = 0; i < nprocs; i++) {
      if (globalPartCounts(i) != targetCount) {
        std::cout << "FAIL: part " << i << " has " << globalPartCounts(i)
                  << " != " << targetCount << "; " << ++ierr << " errors"
                  << std::endl;
      }
    }
    if (ierr == 0) {
      std::cout << "PASS" << std::endl;
    }
  }

  delete problem;
}

