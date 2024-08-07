// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Basic testing of Zoltan2::BasicKokkosIdentifierAdapter 

#include <Kokkos_Core.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_BasicKokkosIdentifierAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>

int main(int narg, char *arg[]) {

  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> userTypes_t;

  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0, gfail = 0;

  // Create global identifiers with weights
  zlno_t numLocalIds = 10;
  const int nWeights = 2;

  Kokkos::View<zgno_t *, typename znode_t::device_type>
    myIds(Kokkos::ViewAllocateWithoutInitializing("myIds"), numLocalIds);
  zgno_t myFirstId = rank * numLocalIds * numLocalIds;
  Kokkos::View<zscalar_t **, typename znode_t::device_type>
    weights(Kokkos::ViewAllocateWithoutInitializing("weights"),
    numLocalIds, nWeights);

  Kokkos::parallel_for(
    Kokkos::RangePolicy<typename znode_t::execution_space,
    zlno_t> (0, numLocalIds), KOKKOS_LAMBDA (zlno_t i) {
    myIds(i) = zgno_t(myFirstId + i);
    weights(i, 0) = 1.0;
    weights(i, 1) = (nprocs - rank) / (i + 1);
  });

  Zoltan2::BasicKokkosIdentifierAdapter<userTypes_t> ia(myIds, weights);

  if (!fail && ia.getLocalNumIDs() != size_t(numLocalIds)) {
    fail = 4;
  }
  if (!fail && ia.getNumWeightsPerID() != nWeights) {
    fail = 5;
  }

  Kokkos::View<const zgno_t *, typename znode_t::device_type> globalIdsIn;
  Kokkos::View<zscalar_t **, typename znode_t::device_type> weightsIn;

  ia.getIDsKokkosView(globalIdsIn);

  ia.getWeightsKokkosView(weightsIn);

  auto host_globalIdsIn = Kokkos::create_mirror_view(globalIdsIn);
  Kokkos::deep_copy(host_globalIdsIn, globalIdsIn);
  auto host_weightsIn = Kokkos::create_mirror_view(weightsIn);
  Kokkos::deep_copy(host_weightsIn, weightsIn);
  auto host_weights = Kokkos::create_mirror_view(weights);
  Kokkos::deep_copy(host_weights, weights);

  auto host_w0 = Kokkos::subview(host_weightsIn, Kokkos::ALL, 0);
  auto host_w1 = Kokkos::subview(host_weightsIn, Kokkos::ALL, 1);

  for (zlno_t i = 0; !fail && i < numLocalIds; i++){
    if (host_globalIdsIn(i) != zgno_t(myFirstId + i)) {
      fail = 8;
    }
    if (!fail && host_w0(i) != 1.0) {
      fail = 9;
    }
    if (!fail && host_w1(i) != host_weights(i, 1)) {
      fail = 10;
    }
  }

  gfail = globalFail(*comm, fail);
  if (gfail) {
    printFailureCode(*comm, fail); // will exit(1)
  }
  if (rank == 0) {
    std::cout << "PASS" << std::endl;
  }
}
