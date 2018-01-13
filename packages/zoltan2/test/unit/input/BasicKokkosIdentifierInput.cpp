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
//                    Andy Wantuch      (acwantu@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Basic testing of Zoltan2::BasicKokkosIdentifierAdapter 

#include <Kokkos_Core.hpp>

#include <Zoltan2_BasicKokkosIdentifierAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;

int main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0, gfail = 0;

  // Create global identifiers with weights

  zlno_t numLocalIds = 10;
  const int nWeights = 2;

  Kokkos::View<zgno_t*> myIds("myIds", numLocalIds);
  zgno_t myFirstId = rank * numLocalIds * numLocalIds;

  Kokkos::View<zscalar_t **> weights("weights", numLocalIds, nWeights);
  for (zlno_t i = 0; i < numLocalIds; i++) {
    myIds(i) = zgno_t(myFirstId + i);
    // Fill in 2D array
    weights(i, 0) = 1.0;
    weights(i, 1) = (nprocs - rank) / (i + 1);
  }

  // Create a Zoltan2::BasicKokkosIdentifierAdapter object
  // and verify that it is correct

  // These types are from /home/acwantu/UUR_git/Trilinos/packages/zoltan2/test/helpers/Zoltan2_TestHelpers.hpp
  typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> userTypes_t;

  // The new Kokkos adapter will take in 2 args instead of 4 because it doesn't need to do striding.
  // The Kokkos::View stores a 2D array of data, so the indexing method is already known.
  // With the old method, there were different ways to do striding because 1D arrays were used.
  Zoltan2::BasicKokkosIdentifierAdapter<userTypes_t> ia(myIds, weights); // TODO: Will need to use new adapter
  //Zoltan2::BasicIdentifierAdapter<userTypes_t> ia(numLocalIds, myIds, weightValues, strides);

  if (!fail && ia.getLocalNumIDs() != size_t(numLocalIds)) {
    fail = 4;
  }

  if (!fail && ia.getNumWeightsPerID() != nWeights) {
    fail = 5;
  }

  Kokkos::View<zgno_t *> globalIdsIn; // Pointer which will later point to the IDs // TODO: 1/12/18 Ask Karen, is this actually a pointer? Bug? Make it *globalIdsIn?
  Kokkos::View<zscalar_t *> weightsIn[nWeights]; // Pointer which will later point to the weights // It's an array of Views.

  // In the old implementation, Views were pointers to memory containing C arrays.
  ia.getIDsView(globalIdsIn); // Make the function mutate globalIdsIn to point to a Kokkos::View

  for (int w = 0; !fail && w < nWeights; w++) {
    ia.getWeightsView(weightsIn[w], w); // This function will need to use the correct Kokkos subview method to get a portion of the view when it's implemented in the adapter.
  }

  Kokkos::View<zscalar_t *> w0 = weightsIn[0];
  Kokkos::View<zscalar_t *> w1 = weightsIn[1];

  for (zlno_t i = 0; !fail && i < numLocalIds; i++){
    if (globalIdsIn(i) != zgno_t(myFirstId + i)) {
      fail = 8;
    }
    if (!fail && w0(i) != 1.0) {
      fail = 9;
    }
    if (!fail && w1(i) != weights(i, 1)) {
      fail = 10;
    }
  }

  gfail = globalFail(comm, fail);
  if (gfail) {
    printFailureCode(comm, fail);   // will exit(1)
  }
  if (rank == 0) {
    std::cout << "PASS" << std::endl;
  }
}

