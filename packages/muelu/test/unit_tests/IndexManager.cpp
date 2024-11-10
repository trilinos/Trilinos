// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_EReductionType.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_UncoupledIndexManager.hpp"
#include "MueLu_LocalLexicographicIndexManager.hpp"
#include "MueLu_GlobalLexicographicIndexManager.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IndexManager, CreateGlobalLexicographicIndexManager, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Teuchos::FancyOStream> fout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fout->setShowAllFrontMatter(false).setShowProcRank(true);

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<const Teuchos::Comm<int> > comm = MueLuTests::TestHelpers::Parameters::getDefaultComm();

  // Set global geometric data
  const bool coupled           = true;
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
      MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions,
                                                                                gNodesPerDir,
                                                                                lNodesPerDir,
                                                                                meshData,
                                                                                "Global Lexicographic");

  RCP<MueLu::GlobalLexicographicIndexManager<LO, GO, NO> > myIndexManager =
      rcp(new MueLu::GlobalLexicographicIndexManager<LO, GO, NO>(comm, coupled, numDimensions,
                                                                 interpolationOrder, gNodesPerDir,
                                                                 lNodesPerDir, coarseRate,
                                                                 coords->getMap()->getMinGlobalIndex()));

  Array<LO> ghostedNodeCoarseLIDs;
  Array<int> ghostedNodeCoarsePIDs;
  Array<GO> ghostedNodeCoarseGIDs;
  myIndexManager->getGhostedNodesData(coords->getMap(), ghostedNodeCoarseLIDs,
                                      ghostedNodeCoarsePIDs, ghostedNodeCoarseGIDs);

  Array<GO> coarseNodeCoarseGIDs;
  Array<GO> coarseNodeFineGIDs;
  myIndexManager->getCoarseNodesData(coords->getMap(), coarseNodeCoarseGIDs,
                                     coarseNodeFineGIDs);

  int chk = 0;
  if (myIndexManager->isAggregationCoupled() != coupled) {
    chk = -1;
  }
  if (myIndexManager->getInterpolationOrder() != interpolationOrder) {
    chk = -1;
  }
  if (myIndexManager->getNumDimensions() != numDimensions) {
    chk = -1;
  }
  if (myIndexManager->getNumGlobalFineNodes() != 125) {
    chk = -1;
  }
  if (myIndexManager->getNumGlobalCoarseNodes() != 27) {
    chk = -1;
  }
  for (int dim = 0; dim < 3; ++dim) {
    if (myIndexManager->getCoarseningRate(dim) != coarseRate[dim]) {
      chk = -1;
    }
    if (myIndexManager->getGlobalFineNodesInDir(dim) != gNodesPerDir[dim]) {
      chk = -1;
    }
    if (myIndexManager->getGlobalCoarseNodesInDir(dim) != 3) {
      chk = -1;
    }
  }
  if (comm->getSize() == 1) {
    if (myIndexManager->getNumLocalFineNodes() != 125) {
      chk = -1;
    }
    if (myIndexManager->getNumLocalCoarseNodes() != 27) {
      chk = -1;
    }
    if (myIndexManager->getNumLocalGhostedNodes() != 27) {
      chk = -1;
    }
    for (int dim = 0; dim < 3; ++dim) {
      if (myIndexManager->getLocalFineNodesInDir(dim) != lNodesPerDir[dim]) {
        chk = -1;
      }
      if (myIndexManager->getStartIndex(dim) != 0) {
        chk = -1;
      }
      if (myIndexManager->getStartIndex(dim + 3) != gNodesPerDir[dim] - 1) {
        chk = -1;
      }
      if (myIndexManager->getOffset(dim) != 0) {
        chk = -1;
      }
      if (myIndexManager->getGhostInterface(dim) != false) {
        chk = -1;
      }
      if (myIndexManager->getGhostInterface(dim + 3) != false) {
        chk = -1;
      }
      if (myIndexManager->getLocalCoarseNodesInDir(dim) != 3) {
        chk = -1;
      }
      if (myIndexManager->getGhostedNodesInDir(dim) != 3) {
        chk = -1;
      }
    }
    GO cnfGIDs[27] = {0, 2, 4, 10, 12, 14, 20, 22, 24,
                      50, 52, 54, 60, 62, 64, 70, 72, 74,
                      100, 102, 104, 110, 112, 114, 120, 122, 124};
    for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
      if (ghostedNodeCoarseLIDs[coarseIdx] != coarseIdx) {
        chk = -1;
      }
      if (ghostedNodeCoarsePIDs[coarseIdx] != 0) {
        chk = -1;
      }
      if (ghostedNodeCoarseGIDs[coarseIdx] != coarseIdx) {
        chk = -1;
      }
      if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
        chk = -1;
      }
      if (coarseNodeCoarseGIDs[coarseIdx] != coarseIdx) {
        chk = -1;
      }
    }
  } else if (comm->getSize() == 4) {
    if (comm->getRank() == 0) {
      if (myIndexManager->getNumLocalFineNodes() != 45) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 12) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3]      = {3, 3, 5};
      LO lCNPD[3]      = {2, 2, 3};
      LO GNPD[3]       = {2, 2, 3};
      GO startIndex[6] = {0, 0, 0, 2, 2, 4};
      GO offsets[3]    = {0, 0, 0};
      bool ghostInt[6] = {false, false, false, false, false, false};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim) != startIndex[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim + 3) != startIndex[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getOffset(dim) != offsets[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != ghostInt[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != ghostInt[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
      int gncPIDs[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      GO gncGIDs[12]  = {0, 1, 3, 4, 9, 10, 12, 13, 18, 19, 21, 22};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarseGIDs[ghostedIdx] != gncGIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[12] = {0, 2, 10, 12, 50, 52, 60, 62, 100, 102, 110, 112};
      GO cncGIDs[12] = {0, 1, 3, 4, 9, 10, 12, 13, 18, 19, 21, 22};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalCoarseNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
        if (coarseNodeCoarseGIDs[coarseIdx] != cncGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 1) {
      if (myIndexManager->getNumLocalFineNodes() != 30) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 6) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3]      = {2, 3, 5};
      LO lCNPD[3]      = {1, 2, 3};
      LO GNPD[3]       = {2, 2, 3};
      GO startIndex[6] = {3, 0, 0, 4, 2, 4};
      GO offsets[3]    = {1, 0, 0};
      bool ghostInt[6] = {true, false, false, false, false, false};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim) != startIndex[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim + 3) != startIndex[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getOffset(dim) != offsets[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != ghostInt[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != ghostInt[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {1, 0, 3, 1, 5, 2, 7, 3, 9, 4, 11, 5};
      int gncPIDs[12] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
      GO gncGIDs[12]  = {1, 2, 4, 5, 10, 11, 13, 14, 19, 20, 22, 23};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarseGIDs[ghostedIdx] != gncGIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[6] = {4, 14, 54, 64, 104, 114};
      GO cncGIDs[6] = {2, 5, 11, 14, 20, 23};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalCoarseNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
        if (coarseNodeCoarseGIDs[coarseIdx] != cncGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 2) {
      if (myIndexManager->getNumLocalFineNodes() != 30) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 6) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3]      = {3, 2, 5};
      LO lCNPD[3]      = {2, 1, 3};
      LO GNPD[3]       = {2, 2, 3};
      GO startIndex[6] = {0, 3, 0, 2, 4, 4};
      GO offsets[3]    = {0, 1, 0};
      bool ghostInt[6] = {false, false, true, false, false, false};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim) != startIndex[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim + 3) != startIndex[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getOffset(dim) != offsets[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != ghostInt[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != ghostInt[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {2, 3, 0, 1, 6, 7, 2, 3, 10, 11, 4, 5};
      int gncPIDs[12] = {0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 2};
      GO gncGIDs[12]  = {3, 4, 6, 7, 12, 13, 15, 16, 21, 22, 24, 25};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarseGIDs[ghostedIdx] != gncGIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[6] = {20, 22, 70, 72, 120, 122};
      GO cncGIDs[6] = {6, 7, 15, 16, 24, 25};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalCoarseNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
        if (coarseNodeCoarseGIDs[coarseIdx] != cncGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 3) {
      if (myIndexManager->getNumLocalFineNodes() != 20) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 3) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3]      = {2, 2, 5};
      LO lCNPD[3]      = {1, 1, 3};
      LO GNPD[3]       = {2, 2, 3};
      GO startIndex[6] = {3, 3, 0, 4, 4, 4};
      GO offsets[3]    = {1, 1, 0};
      bool ghostInt[6] = {true, false, true, false, false, false};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim) != startIndex[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim + 3) != startIndex[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getOffset(dim) != offsets[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != ghostInt[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != ghostInt[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {3, 1, 1, 0, 7, 3, 3, 1, 11, 5, 5, 2};
      int gncPIDs[12] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
      GO gncGIDs[12]  = {4, 5, 7, 8, 13, 14, 16, 17, 22, 23, 25, 26};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarseGIDs[ghostedIdx] != gncGIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[3] = {24, 74, 124};
      GO cncGIDs[3] = {8, 17, 26};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalCoarseNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
        if (coarseNodeCoarseGIDs[coarseIdx] != cncGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    }
  }

  int gbl_chk[1] = {0};
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &chk, gbl_chk);
  TEST_EQUALITY(gbl_chk[0], 0);

}  // CreateGlobalLexicographicIndexManager

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IndexManager, CreateLocalLexicographicIndexManager, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Teuchos::FancyOStream> fout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fout->setShowAllFrontMatter(false).setShowProcRank(true);

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<const Teuchos::Comm<int> > comm = MueLuTests::TestHelpers::Parameters::getDefaultComm();

  // Set global geometric data
  const bool coupled           = true;
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
      MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions,
                                                                                gNodesPerDir,
                                                                                lNodesPerDir,
                                                                                meshData,
                                                                                "Local Lexicographic");

  RCP<MueLu::LocalLexicographicIndexManager<LO, GO, NO> > myIndexManager =
      rcp(new MueLu::LocalLexicographicIndexManager<LO, GO, NO>(comm, coupled, numDimensions,
                                                                interpolationOrder, comm->getRank(),
                                                                comm->getSize(), gNodesPerDir,
                                                                lNodesPerDir, coarseRate, meshData));

  Array<LO> ghostedNodeCoarseLIDs;
  Array<int> ghostedNodeCoarsePIDs;
  Array<GO> ghostedNodeCoarseGIDs;
  myIndexManager->getGhostedNodesData(coords->getMap(), ghostedNodeCoarseLIDs,
                                      ghostedNodeCoarsePIDs, ghostedNodeCoarseGIDs);

  Array<GO> coarseNodeCoarseGIDs;
  Array<GO> coarseNodeFineGIDs;
  myIndexManager->getCoarseNodesData(coords->getMap(), coarseNodeCoarseGIDs,
                                     coarseNodeFineGIDs);

  int chk = 0;
  if (myIndexManager->isAggregationCoupled() != coupled) {
    chk = -1;
  }
  if (myIndexManager->getInterpolationOrder() != interpolationOrder) {
    chk = -1;
  }
  if (myIndexManager->getNumDimensions() != numDimensions) {
    chk = -1;
  }
  if (myIndexManager->getNumGlobalFineNodes() != 125) {
    chk = -1;
  }
  if (myIndexManager->getNumGlobalCoarseNodes() != 27) {
    chk = -1;
  }
  for (int dim = 0; dim < 3; ++dim) {
    if (myIndexManager->getCoarseningRate(dim) != coarseRate[dim]) {
      chk = -1;
    }
    if (myIndexManager->getGlobalFineNodesInDir(dim) != gNodesPerDir[dim]) {
      chk = -1;
    }
    if (myIndexManager->getGlobalCoarseNodesInDir(dim) != 3) {
      chk = -1;
    }
  }
  if (comm->getSize() == 1) {
    if (myIndexManager->getNumLocalFineNodes() != 125) {
      chk = -1;
    }
    if (myIndexManager->getNumLocalCoarseNodes() != 27) {
      chk = -1;
    }
    if (myIndexManager->getNumLocalGhostedNodes() != 27) {
      chk = -1;
    }
    for (int dim = 0; dim < 3; ++dim) {
      if (myIndexManager->getLocalFineNodesInDir(dim) != lNodesPerDir[dim]) {
        chk = -1;
      }
      if (myIndexManager->getStartIndex(dim) != 0) {
        chk = -1;
      }
      if (myIndexManager->getStartIndex(dim + 3) != gNodesPerDir[dim] - 1) {
        chk = -1;
      }
      if (myIndexManager->getOffset(dim) != 0) {
        chk = -1;
      }
      if (myIndexManager->getGhostInterface(dim) != false) {
        chk = -1;
      }
      if (myIndexManager->getGhostInterface(dim + 3) != false) {
        chk = -1;
      }
      if (myIndexManager->getLocalCoarseNodesInDir(dim) != 3) {
        chk = -1;
      }
      if (myIndexManager->getGhostedNodesInDir(dim) != 3) {
        chk = -1;
      }
    }
    GO cnfGIDs[27] = {0, 2, 4, 10, 12, 14, 20, 22, 24,
                      50, 52, 54, 60, 62, 64, 70, 72, 74,
                      100, 102, 104, 110, 112, 114, 120, 122, 124};
    for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
      if (ghostedNodeCoarseLIDs[coarseIdx] != coarseIdx) {
        chk = -1;
      }
      if (ghostedNodeCoarsePIDs[coarseIdx] != 0) {
        chk = -1;
      }
      if (ghostedNodeCoarseGIDs[coarseIdx] != coarseIdx) {
        chk = -1;
      }
      if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
        chk = -1;
      }
      if (coarseNodeCoarseGIDs[coarseIdx] != coarseIdx) {
        chk = -1;
      }
    }
  } else if (comm->getSize() == 4) {
    if (comm->getRank() == 0) {
      if (myIndexManager->getNumLocalFineNodes() != 45) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 12) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3]      = {3, 3, 5};
      LO lCNPD[3]      = {2, 2, 3};
      LO GNPD[3]       = {2, 2, 3};
      GO startIndex[6] = {0, 0, 0, 2, 2, 4};
      GO offsets[3]    = {0, 0, 0};
      bool ghostInt[6] = {false, false, false, false, false, false};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim) != startIndex[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim + 3) != startIndex[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getOffset(dim) != offsets[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != ghostInt[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != ghostInt[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
      int gncPIDs[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      GO gncGIDs[12]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarseGIDs[ghostedIdx] != gncGIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[12] = {0, 2, 6, 8, 18, 20, 24, 26, 36, 38, 42, 44};
      GO cncGIDs[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalCoarseNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
        if (coarseNodeCoarseGIDs[coarseIdx] != cncGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 1) {
      if (myIndexManager->getNumLocalFineNodes() != 30) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 6) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3]      = {2, 3, 5};
      LO lCNPD[3]      = {1, 2, 3};
      LO GNPD[3]       = {2, 2, 3};
      GO startIndex[6] = {3, 0, 0, 4, 2, 4};
      GO offsets[3]    = {1, 0, 0};
      bool ghostInt[6] = {true, false, false, false, false, false};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim) != startIndex[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim + 3) != startIndex[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getOffset(dim) != offsets[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != ghostInt[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != ghostInt[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {1, 0, 3, 1, 5, 2, 7, 3, 9, 4, 11, 5};
      int gncPIDs[12] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
      GO gncGIDs[12]  = {1, 12, 3, 13, 5, 14, 7, 15, 9, 16, 11, 17};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarseGIDs[ghostedIdx] != gncGIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[6] = {46, 50, 58, 62, 70, 74};
      GO cncGIDs[6] = {12, 13, 14, 15, 16, 17};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalCoarseNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
        if (coarseNodeCoarseGIDs[coarseIdx] != cncGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 2) {
      if (myIndexManager->getNumLocalFineNodes() != 30) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 6) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3]      = {3, 2, 5};
      LO lCNPD[3]      = {2, 1, 3};
      LO GNPD[3]       = {2, 2, 3};
      GO startIndex[6] = {0, 3, 0, 2, 4, 4};
      GO offsets[3]    = {0, 1, 0};
      bool ghostInt[6] = {false, false, true, false, false, false};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim) != startIndex[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim + 3) != startIndex[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getOffset(dim) != offsets[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != ghostInt[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != ghostInt[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {2, 3, 0, 1, 6, 7, 2, 3, 10, 11, 4, 5};
      int gncPIDs[12] = {0, 0, 2, 2, 0, 0, 2, 2, 0, 0, 2, 2};
      GO gncGIDs[12]  = {2, 3, 18, 19, 6, 7, 20, 21, 10, 11, 22, 23};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarseGIDs[ghostedIdx] != gncGIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[6] = {78, 80, 90, 92, 102, 104};
      GO cncGIDs[6] = {18, 19, 20, 21, 22, 23};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalCoarseNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
        if (coarseNodeCoarseGIDs[coarseIdx] != cncGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 3) {
      if (myIndexManager->getNumLocalFineNodes() != 20) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 3) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3]      = {2, 2, 5};
      LO lCNPD[3]      = {1, 1, 3};
      LO GNPD[3]       = {2, 2, 3};
      GO startIndex[6] = {3, 3, 0, 4, 4, 4};
      GO offsets[3]    = {1, 1, 0};
      bool ghostInt[6] = {true, false, true, false, false, false};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim) != startIndex[dim]) {
          chk = -1;
        }
        if (myIndexManager->getStartIndex(dim + 3) != startIndex[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getOffset(dim) != offsets[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != ghostInt[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != ghostInt[dim + 3]) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {3, 1, 1, 0, 7, 3, 3, 1, 11, 5, 5, 2};
      int gncPIDs[12] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
      GO gncGIDs[12]  = {3, 13, 19, 24, 7, 15, 21, 25, 11, 17, 23, 26};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarseGIDs[ghostedIdx] != gncGIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[3] = {108, 116, 124};
      GO cncGIDs[3] = {24, 25, 26};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalCoarseNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
        if (coarseNodeCoarseGIDs[coarseIdx] != cncGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    }
  }

  int gbl_chk[1] = {0};
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &chk, gbl_chk);
  TEST_EQUALITY(gbl_chk[0], 0);

}  // CreateLocalLexicographicIndexManager

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IndexManager, CreateUncoupledIndexManager, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Teuchos::FancyOStream> fout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fout->setShowAllFrontMatter(false).setShowProcRank(true);

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<const Teuchos::Comm<int> > comm = MueLuTests::TestHelpers::Parameters::getDefaultComm();

  // Set global geometric data
  const bool coupled           = false;
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
      MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions,
                                                                                gNodesPerDir,
                                                                                lNodesPerDir,
                                                                                meshData,
                                                                                "Local Lexicographic");

  RCP<MueLu::UncoupledIndexManager<LO, GO, NO> > myIndexManager =
      rcp(new MueLu::UncoupledIndexManager<LO, GO, NO>(comm, coupled, numDimensions,
                                                       interpolationOrder, comm->getRank(),
                                                       comm->getSize(), gNodesPerDir,
                                                       lNodesPerDir, coarseRate, false));
  myIndexManager->computeGlobalCoarseParameters();

  Array<LO> ghostedNodeCoarseLIDs;
  Array<int> ghostedNodeCoarsePIDs;
  Array<GO> ghostedNodeCoarseGIDs;
  myIndexManager->getGhostedNodesData(coords->getMap(), ghostedNodeCoarseLIDs,
                                      ghostedNodeCoarsePIDs, ghostedNodeCoarseGIDs);

  Array<GO> coarseNodeCoarseGIDs;
  Array<GO> coarseNodeFineGIDs;
  myIndexManager->getCoarseNodesData(coords->getMap(), coarseNodeCoarseGIDs,
                                     coarseNodeFineGIDs);

  int chk = 0;
  if (myIndexManager->isAggregationCoupled() != coupled) {
    chk = -1;
  }
  if (myIndexManager->getInterpolationOrder() != interpolationOrder) {
    chk = -1;
  }
  if (myIndexManager->getNumDimensions() != numDimensions) {
    chk = -1;
  }
  if (myIndexManager->getNumGlobalFineNodes() != -1) {
    chk = -1;
  }
  for (int dim = 0; dim < 3; ++dim) {
    if (myIndexManager->getCoarseningRate(dim) != coarseRate[dim]) {
      chk = -1;
    }
  }
  if (comm->getSize() == 1) {
    if (myIndexManager->getNumLocalFineNodes() != 125) {
      chk = -1;
    }
    if (myIndexManager->getNumLocalCoarseNodes() != 27) {
      chk = -1;
    }
    if (myIndexManager->getNumLocalGhostedNodes() != 27) {
      chk = -1;
    }
    if (myIndexManager->getNumGlobalCoarseNodes() != 27) {
      chk = -1;
    }
    for (int dim = 0; dim < 3; ++dim) {
      if (myIndexManager->getLocalFineNodesInDir(dim) != lNodesPerDir[dim]) {
        chk = -1;
      }
      if (myIndexManager->getGhostInterface(dim) != false) {
        chk = -1;
      }
      if (myIndexManager->getGhostInterface(dim + 3) != false) {
        chk = -1;
      }
      if (myIndexManager->getLocalCoarseNodesInDir(dim) != 3) {
        chk = -1;
      }
      if (myIndexManager->getGhostedNodesInDir(dim) != 3) {
        chk = -1;
      }
      if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
        chk = -1;
      }
    }
    GO cnfGIDs[27] = {0, 2, 4, 10, 12, 14, 20, 22, 24,
                      50, 52, 54, 60, 62, 64, 70, 72, 74,
                      100, 102, 104, 110, 112, 114, 120, 122, 124};
    for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
      if (ghostedNodeCoarseLIDs[coarseIdx] != coarseIdx) {
        chk = -1;
      }
      if (ghostedNodeCoarsePIDs[coarseIdx] != 0) {
        chk = -1;
      }
      if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
        chk = -1;
      }
    }
  } else if (comm->getSize() == 4) {
    if (comm->getRank() == 0) {
      if (myIndexManager->getNumLocalFineNodes() != 45) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 12) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3] = {3, 3, 5};
      LO lCNPD[3] = {2, 2, 3};
      LO GNPD[3]  = {2, 2, 3};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != false) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != false) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
      int gncPIDs[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[12] = {0, 2, 6, 8, 18, 20, 24, 26, 36, 38, 42, 44};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 1) {
      if (myIndexManager->getNumLocalFineNodes() != 30) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 12) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3] = {2, 3, 5};
      LO lCNPD[3] = {2, 2, 3};
      LO GNPD[3]  = {2, 2, 3};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != false) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != false) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
      int gncPIDs[12] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[12] = {45, 46, 49, 50, 57, 58, 61, 62, 69, 70, 73, 74};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 2) {
      if (myIndexManager->getNumLocalFineNodes() != 30) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 12) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3] = {3, 2, 5};
      LO lCNPD[3] = {2, 2, 3};
      LO GNPD[3]  = {2, 2, 3};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != false) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != false) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
      int gncPIDs[12] = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[12] = {75, 77, 78, 80, 87, 89, 90, 92, 99, 101, 102, 104};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 3) {
      if (myIndexManager->getNumLocalFineNodes() != 20) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 12) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 12) {
        chk = -1;
      }
      LO lFNPD[3] = {2, 2, 5};
      LO lCNPD[3] = {2, 2, 3};
      LO GNPD[3]  = {2, 2, 3};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != false) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != false) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
          chk = -1;
        }
      }
      LO gncLIDs[12]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
      int gncPIDs[12] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[12] = {105, 106, 107, 108, 113, 114, 115, 116, 121, 122, 123, 124};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    }
  }

  int gbl_chk[1] = {0};
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &chk, gbl_chk);
  TEST_EQUALITY(gbl_chk[0], 0);

}  // CreateUncoupledIndexManager

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(IndexManager, UncoupledIndexManagerSingleCoarsePoint, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Teuchos::FancyOStream> fout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  fout->setShowAllFrontMatter(false).setShowProcRank(true);

  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  RCP<const Teuchos::Comm<int> > comm = MueLuTests::TestHelpers::Parameters::getDefaultComm();

  // Set global geometric data
  const bool coupled           = false;
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
  // Change the number of nodes in the third
  // direction to trigger singleCoarsePoint
  gNodesPerDir[2] = 2;

  RCP<RealValuedMultiVector> coords =
      MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions,
                                                                                gNodesPerDir,
                                                                                lNodesPerDir,
                                                                                meshData,
                                                                                "Local Lexicographic");

  RCP<MueLu::UncoupledIndexManager<LO, GO, NO> > myIndexManager =
      rcp(new MueLu::UncoupledIndexManager<LO, GO, NO>(comm, coupled, numDimensions,
                                                       interpolationOrder, comm->getRank(),
                                                       comm->getSize(), gNodesPerDir,
                                                       lNodesPerDir, coarseRate, true));
  myIndexManager->computeGlobalCoarseParameters();

  Array<LO> ghostedNodeCoarseLIDs;
  Array<int> ghostedNodeCoarsePIDs;
  Array<GO> ghostedNodeCoarseGIDs;
  myIndexManager->getGhostedNodesData(coords->getMap(), ghostedNodeCoarseLIDs,
                                      ghostedNodeCoarsePIDs, ghostedNodeCoarseGIDs);

  Array<GO> coarseNodeCoarseGIDs;
  Array<GO> coarseNodeFineGIDs;
  myIndexManager->getCoarseNodesData(coords->getMap(), coarseNodeCoarseGIDs,
                                     coarseNodeFineGIDs);

  int chk = 0;
  if (myIndexManager->isAggregationCoupled() != coupled) {
    chk = -1;
  }
  if (myIndexManager->getInterpolationOrder() != interpolationOrder) {
    chk = -1;
  }
  if (myIndexManager->getNumDimensions() != numDimensions) {
    chk = -1;
  }
  if (myIndexManager->getNumGlobalFineNodes() != -1) {
    chk = -1;
  }
  for (int dim = 0; dim < 3; ++dim) {
    if (myIndexManager->getCoarseningRate(dim) != coarseRate[dim]) {
      chk = -1;
    }
  }
  if (comm->getSize() == 1) {
    if (myIndexManager->getNumLocalFineNodes() != 50) {
      chk = -1;
    }
    if (myIndexManager->getNumLocalCoarseNodes() != 9) {
      chk = -1;
    }
    if (myIndexManager->getNumLocalGhostedNodes() != 9) {
      chk = -1;
    }
    if (myIndexManager->getNumGlobalCoarseNodes() != 9) {
      chk = -1;
    }
    LO lFNPD[3] = {5, 5, 2};
    LO lCNPD[3] = {3, 3, 1};
    LO GNPD[3]  = {3, 3, 1};
    for (int dim = 0; dim < 3; ++dim) {
      if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
        chk = -1;
      }
      if (myIndexManager->getGhostInterface(dim) != false) {
        chk = -1;
      }
      if (myIndexManager->getGhostInterface(dim + 3) != false) {
        chk = -1;
      }
      if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
        chk = -1;
      }
      if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
        chk = -1;
      }
      if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
        chk = -1;
      }
    }
    GO cnfGIDs[9] = {0, 2, 4, 10, 12, 14, 20, 22, 24};
    for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
      if (ghostedNodeCoarseLIDs[coarseIdx] != coarseIdx) {
        chk = -1;
      }
      if (ghostedNodeCoarsePIDs[coarseIdx] != 0) {
        chk = -1;
      }
      if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
        chk = -1;
      }
    }
  } else if (comm->getSize() == 4) {
    if (comm->getRank() == 0) {
      if (myIndexManager->getNumLocalFineNodes() != 18) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 4) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 4) {
        chk = -1;
      }
      LO lFNPD[3] = {3, 3, 2};
      LO lCNPD[3] = {2, 2, 1};
      LO GNPD[3]  = {2, 2, 1};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != false) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != false) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
          chk = -1;
        }
      }
      LO gncLIDs[4]  = {0, 1, 2, 3};
      int gncPIDs[4] = {0, 0, 0, 0};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[4] = {0, 2, 6, 8};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 1) {
      if (myIndexManager->getNumLocalFineNodes() != 12) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 2) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 2) {
        chk = -1;
      }
      LO lFNPD[3] = {2, 3, 2};
      LO lCNPD[3] = {1, 2, 1};
      LO GNPD[3]  = {1, 2, 1};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != false) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != false) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
          chk = -1;
        }
      }
      LO gncLIDs[2]  = {0, 1};
      int gncPIDs[2] = {1, 1};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[2] = {18, 22};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 2) {
      if (myIndexManager->getNumLocalFineNodes() != 12) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 2) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 2) {
        chk = -1;
      }
      LO lFNPD[3] = {3, 2, 2};
      LO lCNPD[3] = {2, 1, 1};
      LO GNPD[3]  = {2, 1, 1};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != false) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != false) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
          chk = -1;
        }
      }
      LO gncLIDs[2]  = {0, 1};
      int gncPIDs[2] = {2, 2};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[2] = {30, 32};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    } else if (comm->getRank() == 3) {
      if (myIndexManager->getNumLocalFineNodes() != 8) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalCoarseNodes() != 1) {
        chk = -1;
      }
      if (myIndexManager->getNumLocalGhostedNodes() != 1) {
        chk = -1;
      }
      LO lFNPD[3] = {2, 2, 2};
      LO lCNPD[3] = {1, 1, 1};
      LO GNPD[3]  = {1, 1, 1};
      for (int dim = 0; dim < 3; ++dim) {
        if (myIndexManager->getLocalFineNodesInDir(dim) != lFNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim) != false) {
          chk = -1;
        }
        if (myIndexManager->getGhostInterface(dim + 3) != false) {
          chk = -1;
        }
        if (myIndexManager->getLocalCoarseNodesInDir(dim) != lCNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGhostedNodesInDir(dim) != GNPD[dim]) {
          chk = -1;
        }
        if (myIndexManager->getGlobalCoarseNodesInDir(dim) != -1) {
          chk = -1;
        }
      }
      LO gncLIDs[1]  = {0};
      int gncPIDs[1] = {3};
      for (int ghostedIdx = 0; ghostedIdx < myIndexManager->getNumLocalGhostedNodes(); ++ghostedIdx) {
        if (ghostedNodeCoarseLIDs[ghostedIdx] != gncLIDs[ghostedIdx]) {
          chk = -1;
        }
        if (ghostedNodeCoarsePIDs[ghostedIdx] != gncPIDs[ghostedIdx]) {
          chk = -1;
        }
      }
      GO cnfGIDs[1] = {42};
      for (int coarseIdx = 0; coarseIdx < myIndexManager->getNumLocalGhostedNodes(); ++coarseIdx) {
        if (coarseNodeFineGIDs[coarseIdx] != cnfGIDs[coarseIdx]) {
          chk = -1;
        }
      }
    }
  }

  int gbl_chk[1] = {0};
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &chk, gbl_chk);
  TEST_EQUALITY(gbl_chk[0], 0);

}  // UncoupledIndexManagerSingleCoarsePoint

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IndexManager, CreateGlobalLexicographicIndexManager, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IndexManager, CreateLocalLexicographicIndexManager, Scalar, LO, GO, Node)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IndexManager, CreateUncoupledIndexManager, Scalar, LO, GO, Node)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(IndexManager, UncoupledIndexManagerSingleCoarsePoint, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
