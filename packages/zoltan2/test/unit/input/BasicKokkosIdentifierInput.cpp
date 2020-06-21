// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
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
