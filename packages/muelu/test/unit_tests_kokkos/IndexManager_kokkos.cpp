// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_IndexManager_kokkos.hpp"

namespace MueLuTests {

template <class IndexManager_type, class LO>
struct testFunctor {
  IndexManager_type IndexManager;

  testFunctor(RCP<IndexManager_type> myIM)
    : IndexManager(*myIM) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO& LID) const {
    LO fineNodeTuple[3];
    IndexManager.getFineLID2FineTuple(LID, fineNodeTuple);
    printf("LID=%d --> (i,j,k)=(%d, %d, %d)\n", static_cast<int>(LID),
           static_cast<int>(fineNodeTuple[0]),
           static_cast<int>(fineNodeTuple[1]),
           static_cast<int>(fineNodeTuple[2]));
  }
};

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IndexManager_kokkos, IndexManager, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;
  RCP<Teuchos::FancyOStream> fout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fout->setShowAllFrontMatter(false).setShowProcRank(true);

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<const Teuchos::Comm<int> > comm = MueLuTests::TestHelpers_kokkos::Parameters::getDefaultComm();

  // Set global geometric data
  const int numDimensions      = 3;
  const int interpolationOrder = 0;
  Array<GO> meshData;
  Array<GO> gNodesPerDir(3);
  Array<LO> lNodesPerDir(3);
  Array<LO> coarseRate(3);
  for (int dim = 0; dim < 3; ++dim) {
    if (dim < numDimensions) {
      // Use more nodes in 1D to have a reasonable number of nodes per procs
      gNodesPerDir[dim] = 5;
      coarseRate[dim]   = 2;
    } else {
      gNodesPerDir[dim] = -1;
      coarseRate[dim]   = 1;
    }
  }

  RCP<RealValuedMultiVector> coords =
      MueLuTests::TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions,
                                                                                       gNodesPerDir,
                                                                                       lNodesPerDir,
                                                                                       meshData,
                                                                                       "Local Lexicographic");

  std::cout << "Hello!" << std::endl;
  RCP<IndexManager_kokkos> myIndexManager = rcp(new IndexManager_kokkos(numDimensions,
                                                                        interpolationOrder,
                                                                        comm->getRank(),
                                                                        lNodesPerDir,
                                                                        coarseRate));

  using execution_space = typename IndexManager_kokkos::execution_space;
  using memory_space    = typename IndexManager_kokkos::memory_space;
  using device_type     = typename IndexManager_kokkos::device_type;

  Kokkos::fence();
  std::cout << "Allocate host views before performing checks" << std::endl;

  typename Kokkos::View<const int[3], device_type>::HostMirror coarseRate_h = Kokkos::create_mirror_view(myIndexManager->getCoarseningRates());
  Kokkos::deep_copy(coarseRate_h, myIndexManager->getCoarseningRates());

  typename Kokkos::View<LO[3], device_type>::HostMirror lFineNodesPerDir_h = Kokkos::create_mirror_view(myIndexManager->getLocalFineNodesPerDir());
  Kokkos::deep_copy(lFineNodesPerDir_h, myIndexManager->getLocalFineNodesPerDir());

  typename Kokkos::View<LO[3], device_type>::HostMirror lCoarseNodesPerDir_h = Kokkos::create_mirror_view(myIndexManager->getCoarseNodesPerDir());
  Kokkos::deep_copy(lCoarseNodesPerDir_h, myIndexManager->getCoarseNodesPerDir());

  Kokkos::fence();
  std::cout << "Start checking process" << std::endl;

  int chk = 0;
  if (myIndexManager->getInterpolationOrder() != interpolationOrder) {
    chk = -1;
  }
  if (myIndexManager->getNumDimensions() != numDimensions) {
    chk = -1;
  }
  for (int dim = 0; dim < 3; ++dim) {
    if (coarseRate_h(dim) != coarseRate[dim]) {
      chk = -1;
    }
  }
  if (comm->getSize() == 1) {
    if (myIndexManager->getNumLocalFineNodes() != 125) {
      chk = -1;
    }
    if (myIndexManager->getNumCoarseNodes() != 27) {
      chk = -1;
    }
    for (int dim = 0; dim < 3; ++dim) {
      if (lFineNodesPerDir_h(dim) != lNodesPerDir[dim]) {
        chk = -1;
      }
      if (lCoarseNodesPerDir_h(dim) != 3) {
        chk = -1;
      }
    }

  } else if (comm->getSize() == 4) {
    if (comm->getRank() == 0) {
      if (myIndexManager->getNumLocalFineNodes() != 45) {
        chk = -1;
      }
      if (myIndexManager->getNumCoarseNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3] = {3, 3, 5};
      LO lCNPD[3] = {2, 2, 3};
      for (int dim = 0; dim < 3; ++dim) {
        if (lFineNodesPerDir_h(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (lCoarseNodesPerDir_h(dim) != lCNPD[dim]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 1) {
      if (myIndexManager->getNumLocalFineNodes() != 30) {
        chk = -1;
      }
      if (myIndexManager->getNumCoarseNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3] = {2, 3, 5};
      LO lCNPD[3] = {2, 2, 3};
      for (int dim = 0; dim < 3; ++dim) {
        if (lFineNodesPerDir_h(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (lCoarseNodesPerDir_h(dim) != lCNPD[dim]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 2) {
      if (myIndexManager->getNumLocalFineNodes() != 30) {
        chk = -1;
      }
      if (myIndexManager->getNumCoarseNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3] = {3, 2, 5};
      LO lCNPD[3] = {2, 2, 3};
      for (int dim = 0; dim < 3; ++dim) {
        if (lFineNodesPerDir_h(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (lCoarseNodesPerDir_h(dim) != lCNPD[dim]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 3) {
      if (myIndexManager->getNumLocalFineNodes() != 20) {
        chk = -1;
      }
      if (myIndexManager->getNumCoarseNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3] = {2, 2, 5};
      LO lCNPD[3] = {2, 2, 3};
      for (int dim = 0; dim < 3; ++dim) {
        if (lFineNodesPerDir_h(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (lCoarseNodesPerDir_h(dim) != lCNPD[dim]) {
          chk = -1;
        }
      }
    }
  }

  int gbl_chk[1] = {0};
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &chk, gbl_chk);
  TEST_EQUALITY(gbl_chk[0], 0);

  std::cout << "Testing UIM methods call from device." << std::endl;

  testFunctor<IndexManager_kokkos, LO> myTestFunctor(myIndexManager);
  Kokkos::parallel_for("test UIM::getFineLID2FineTuple()",
                       Kokkos::RangePolicy<execution_space>(0, 125),
                       myTestFunctor);
}  // Uncoupled

#ifdef HAVE_MUELU_SYCL
// Enabling this w/ SYCL yields the error "SYCL kernel cannot call a variadic function"
#define MUELU_ETI_GROUP(SC, LO, GO, NO)
#else
#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IndexManager_kokkos, IndexManager, SC, LO, GO, NO)
#endif

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
