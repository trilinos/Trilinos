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
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_BlackBoxPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
// #include "MueLu_CoarseMapFactory.hpp"

namespace MueLuTests {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class BlackBoxPFactoryTester {
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  BlackBoxPFactoryTester() {}

  //! Destructor.
  virtual ~BlackBoxPFactoryTester() {}
  //@}

  void TestGetGeometricData(Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> >& coordinates,
                            Array<LO> coarseRate, Array<GO> gNodesPerDim, Array<LO> lNodesPerDim,
                            LO BlkSize, Array<GO>& gIndices, Array<GO>& gCoarseNodesPerDir,
                            Array<GO>& ghostGIDs, Array<GO>& coarseNodesGIDs, Array<GO>& colGIDs,
                            GO& gNumCoarseNodes, Array<LO>& myOffset,
                            Array<LO>& lCoarseNodesPerDir, Array<LO>& glCoarseNodesPerDir,
                            Array<LO>& endRate, LO& lNumCoarseNodes, Array<bool>& ghostInterface,
                            ArrayRCP<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coarseNodes, Array<int>& boundaryFlags) const {
    // Call the method to be tested.
    RCP<typename MueLu::BlackBoxPFactory<SC, LO, GO, Node>::NodesIDs> ghostedCoarseNodes = rcp(new typename MueLu::BlackBoxPFactory<SC, LO, GO, Node>::NodesIDs{});
    MueLu::BlackBoxPFactory<SC, LO, GO, Node> mybbmgPFactory;
    mybbmgPFactory.GetGeometricData(coordinates, coarseRate, gNodesPerDim, lNodesPerDim, BlkSize,
                                    gIndices, myOffset, ghostInterface, endRate,
                                    gCoarseNodesPerDir, lCoarseNodesPerDir, glCoarseNodesPerDir,
                                    ghostGIDs, coarseNodesGIDs, colGIDs, gNumCoarseNodes,
                                    lNumCoarseNodes, coarseNodes, boundaryFlags,
                                    ghostedCoarseNodes);
  };

  void TestComputeLocalEntries(const RCP<const Matrix>& Aghost, const Array<LO> coarseRate,
                               const Array<LO> endRate, const LO BlkSize, const Array<LO> elemInds,
                               const Array<LO> lCoarseElementsPerDir,
                               const LO numDimensions, const Array<LO> lFineNodesPerDir,
                               const Array<GO> gFineNodesPerDir, const Array<GO> gIndices,
                               const Array<LO> lCoarseNodesPerDir,
                               const Array<bool> ghostInterface,
                               const Array<int> elementFlags,
                               const std::string stencilType, const std::string blockStrategy,
                               const Array<LO> elementNodesPerDir, const LO numNodesInElement,
                               const Array<GO> colGIDs,
                               Teuchos::SerialDenseMatrix<LO, SC>& Pi,
                               Teuchos::SerialDenseMatrix<LO, SC>& Pf,
                               Teuchos::SerialDenseMatrix<LO, SC>& Pe,
                               Array<LO>& dofType, Array<LO>& lDofInd) const {
    MueLu::BlackBoxPFactory<SC, LO, GO, NO> mybbmgPFactory;
    mybbmgPFactory.ComputeLocalEntries(Aghost, coarseRate, endRate, BlkSize, elemInds,
                                       lCoarseElementsPerDir, numDimensions,
                                       lFineNodesPerDir, gFineNodesPerDir, gIndices,
                                       lCoarseNodesPerDir, ghostInterface, elementFlags,
                                       stencilType, blockStrategy, elementNodesPerDir,
                                       numNodesInElement, colGIDs, Pi, Pf, Pe, dofType, lDofInd);
  };

  void TestFormatStencil(const LO BlkSize, const Array<bool> ghostInterface, const LO ie,
                         const LO je, const LO ke, const ArrayView<const SC> rowValues,
                         const Array<LO> elementNodesPerDir, const int collapseFlags[3],
                         const std::string stencilType, Array<SC>& stencil) const {
    MueLu::BlackBoxPFactory<SC, LO, GO, Node> mybbmgPFactory;
    mybbmgPFactory.FormatStencil(BlkSize, ghostInterface, ie, je, ke, rowValues,
                                 elementNodesPerDir, collapseFlags, stencilType, stencil);
  };

  void TestGetNodeInfo(const LO ie, const LO je, const LO ke, const Array<LO> elementNodesPerDir,
                       int* nodeType, LO& nodeIndex) const {
    MueLu::BlackBoxPFactory<SC, LO, GO, Node> mybbmgPFactory;
    int dummy;
    mybbmgPFactory.GetNodeInfo(ie, je, ke, elementNodesPerDir, nodeType, nodeIndex, &dummy);
  };
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GetProblemData(RCP<const Teuchos::Comm<int> >& comm, const Xpetra::UnderlyingLib lib,
                    const LocalOrdinal numDimensions, const std::string stencilType,
                    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Op,
                    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> >& Coordinates,
                    RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& map,
                    Array<GlobalOrdinal>& gNodesPerDim, Array<LocalOrdinal>& lNodesPerDim) {
#include "MueLu_UseShortNames.hpp"

  GO nx         = 7;
  GO ny         = 7;
  GO nz         = (numDimensions < 3 ? 1 : 7);
  GO gNumPoints = nx * ny * nz;

  gNodesPerDim[0] = nx;
  gNodesPerDim[1] = ny;
  gNodesPerDim[2] = nz;

  GO myOffset = 0;
  if (comm->getSize() == 1) {
    myOffset        = 0;
    lNodesPerDim[0] = nx;
    lNodesPerDim[1] = ny;
    lNodesPerDim[2] = nz;
  } else if (comm->getSize() == 4) {
    if (comm->getRank() == 0) {
      if (numDimensions == 2) {
        myOffset        = 0;
        lNodesPerDim[0] = 4;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 0;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 4;
      }
    } else if (comm->getRank() == 1) {
      if (numDimensions == 2) {
        myOffset        = 4;
        lNodesPerDim[0] = 3;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 28;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 3;
        lNodesPerDim[2] = 4;
      }
    } else if (comm->getRank() == 2) {
      if (numDimensions == 2) {
        myOffset        = 28;
        lNodesPerDim[0] = 4;
        lNodesPerDim[1] = 3;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 196;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 3;
      }
    } else if (comm->getRank() == 3) {
      if (numDimensions == 2) {
        myOffset        = 32;
        lNodesPerDim[0] = 3;
        lNodesPerDim[1] = 3;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 224;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 3;
        lNodesPerDim[2] = 3;
      }
    }
  }

  GO myZoffset = 0, myYoffset = 0, myXoffset = 0;
  if (numDimensions == 2) {
    myZoffset = 0;
    myYoffset = myOffset / gNodesPerDim[0];
    myXoffset = myOffset % gNodesPerDim[0];
  } else if (numDimensions == 3) {
    myZoffset = myOffset / (gNodesPerDim[1] * gNodesPerDim[0]);
    myYoffset = (myOffset - myZoffset * gNodesPerDim[1] * gNodesPerDim[0]) / gNodesPerDim[0];
    myXoffset = (myOffset - myZoffset * gNodesPerDim[1] * gNodesPerDim[0]) % gNodesPerDim[0];
  }
  LO lNumPoints = lNodesPerDim[0] * lNodesPerDim[1] * lNodesPerDim[2];

  // Construct map and local coordinates
  Teuchos::Array<GO> myGIDs(lNumPoints);
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> myXCoords(lNumPoints);
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> myYCoords(lNumPoints);
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> myZCoords(lNumPoints);
  for (LO k = 0; k < lNodesPerDim[2]; ++k) {
    for (LO j = 0; j < lNodesPerDim[1]; ++j) {
      for (LO i = 0; i < lNodesPerDim[0]; ++i) {
        myGIDs[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] = myOffset + k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + i;
        myXCoords[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] =
            (i + myXoffset) / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[0] - 1);
        myYCoords[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] =
            (j + myYoffset) / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[1] - 1);
        myZCoords[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] =
            (k + myZoffset) / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[2] - 1);
      }
    }
  }

  Teuchos::Array<Teuchos::ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > myCoordinates(numDimensions);
  if (numDimensions == 1) {
    myCoordinates[0] = myXCoords();
  } else if (numDimensions == 2) {
    myCoordinates[0] = myXCoords();
    myCoordinates[1] = myYCoords();
  } else if (numDimensions == 3) {
    myCoordinates[0] = myXCoords();
    myCoordinates[1] = myYCoords();
    myCoordinates[2] = myZCoords();
  }

  // Create the map and store coordinates using the above array views
  map         = MapFactory::Build(lib, gNumPoints, myGIDs(), 0, comm);
  Coordinates = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(map, myCoordinates(),
                                                                                                                     numDimensions);

  // small parameter list for Galeri
  Teuchos::ParameterList problemList;
  if (numDimensions == 1) {
    problemList.set("nx", gNodesPerDim[0]);
  } else if (numDimensions == 2) {
    problemList.set("nx", gNodesPerDim[0]);
    problemList.set("ny", gNodesPerDim[1]);
  } else if (numDimensions == 3) {
    problemList.set("nx", gNodesPerDim[0]);
    problemList.set("ny", gNodesPerDim[1]);
    problemList.set("nz", gNodesPerDim[2]);
  }
  problemList.set("keepBCs", true);

  // create Poisson problem and matrix
  if (stencilType == "reduced") {
    if (numDimensions == 1) {
      Galeri::Xpetra::
          Laplace1DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>
              Problem(problemList, map);
      Op = Problem.BuildMatrix();
    } else if (numDimensions == 2) {
      Galeri::Xpetra::
          Laplace2DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>
              Problem(problemList, map);
      Op = Problem.BuildMatrix();
    } else if (numDimensions == 3) {
      Galeri::Xpetra::
          Laplace3DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>
              Problem(problemList, map);
      Op = Problem.BuildMatrix();
    }
  } else if (stencilType == "full") {
    if (numDimensions == 1) {
      Galeri::Xpetra::
          Laplace1DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>
              Problem(problemList, map);
      Op = Problem.BuildMatrix();
    } else if (numDimensions == 2) {
      Galeri::Xpetra::
          Star2DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>
              Problem(problemList, map);
      Op = Problem.BuildMatrix();
    } else if (numDimensions == 3) {
      Galeri::Xpetra::
          Brick3DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>
              Problem(problemList, map);
      Op = Problem.BuildMatrix();
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlackBoxPFactory, Constructor, Scalar, LocalOrdinal,
                                  GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<BlackBoxPFactory> bbmgPFact = rcp(new BlackBoxPFactory);
  TEST_EQUALITY(bbmgPFact != Teuchos::null, true);

}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlackBoxPFactory, BlackBoxGhosts, Scalar, LocalOrdinal,
                                  GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // used Xpetra lib (for maps and smoothers)
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  const LO numDimensions        = 2;
  const std::string stencilType = "reduced";

  RCP<Matrix> Op;
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates;
  RCP<Map> map;
  Array<GO> gNodesPerDim(3);
  Array<LO> lNodesPerDim(3);
  GetProblemData<SC, LO, GO, NO>(comm, lib, numDimensions, stencilType, Op, coordinates, map,
                                 gNodesPerDim, lNodesPerDim);

  Array<LO> coarseRate(3);
  coarseRate[0] = 2;
  coarseRate[1] = 2;
  if (numDimensions == 2) {
    coarseRate[2] = 1;
  } else if (numDimensions == 3) {
    coarseRate[2] = 2;
  }

  BlackBoxPFactoryTester<SC, LO, GO, Node> factTester;
  Array<GO> gIndices(3), gCoarseNodesPerDir(3), ghostGIDs, coarseNodesGIDs, colGIDs;
  Array<LO> myOffset(3), lCoarseNodesPerDir(3), glCoarseNodesPerDir(3), endRate(3);
  Array<bool> ghostInterface(6);
  Array<int> boundaryFlags(3);
  Array<ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > fineNodes(numDimensions);
  ArrayRCP<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coarseNodes(numDimensions);
  GO gNumCoarseNodes = 0;
  LO lNumCoarseNodes = 0;
  for (LO dim = 0; dim < numDimensions; ++dim) {
    fineNodes[dim] = coordinates->getData(dim)();
  }

  factTester.TestGetGeometricData(coordinates, coarseRate, gNodesPerDim, lNodesPerDim, 1,
                                  gIndices, gCoarseNodesPerDir, ghostGIDs, coarseNodesGIDs,
                                  colGIDs, gNumCoarseNodes, myOffset, lCoarseNodesPerDir,
                                  glCoarseNodesPerDir, endRate, lNumCoarseNodes, ghostInterface,
                                  coarseNodes, boundaryFlags);

  Array<GO> ghostGIDs_check, coarseNodesGIDs_check, colGIDs_check;
  if (map->getComm()->getSize() == 1) {
    coarseNodesGIDs_check.resize(16);
    for (LO i = 0; i < 16; ++i) {
      coarseNodesGIDs_check[i] = i;
    }
    colGIDs_check = coarseNodesGIDs_check;
    ghostGIDs_check.resize(0);
  } else if (map->getComm()->getSize() == 4) {
    coarseNodesGIDs_check.resize(4);
    if (map->getComm()->getRank() == 0) {
      coarseNodesGIDs_check[0] = 0;
      coarseNodesGIDs_check[1] = 1;
      coarseNodesGIDs_check[2] = 4;
      coarseNodesGIDs_check[3] = 5;

      ghostGIDs_check.resize(5);
      ghostGIDs_check[0] = 4;
      ghostGIDs_check[1] = 18;
      ghostGIDs_check[2] = 28;
      ghostGIDs_check[3] = 30;
      ghostGIDs_check[4] = 32;

      colGIDs_check.resize(9);
      colGIDs_check[0] = 0;
      colGIDs_check[1] = 1;
      colGIDs_check[2] = 4;
      colGIDs_check[3] = 5;
      colGIDs_check[4] = 2;
      colGIDs_check[5] = 6;
      colGIDs_check[6] = 8;
      colGIDs_check[7] = 9;
      colGIDs_check[8] = 10;
    } else if (map->getComm()->getRank() == 1) {
      coarseNodesGIDs_check[0] = 2;
      coarseNodesGIDs_check[1] = 3;
      coarseNodesGIDs_check[2] = 6;
      coarseNodesGIDs_check[3] = 7;

      ghostGIDs_check.resize(2);
      ghostGIDs_check[0] = 32;
      ghostGIDs_check[1] = 34;

      colGIDs_check.resize(6);
      colGIDs_check[0] = 2;
      colGIDs_check[1] = 3;
      colGIDs_check[2] = 6;
      colGIDs_check[3] = 7;
      colGIDs_check[4] = 10;
      colGIDs_check[5] = 11;
    } else if (map->getComm()->getRank() == 2) {
      coarseNodesGIDs_check[0] = 8;
      coarseNodesGIDs_check[1] = 9;
      coarseNodesGIDs_check[2] = 12;
      coarseNodesGIDs_check[3] = 13;

      ghostGIDs_check.resize(2);
      ghostGIDs_check[0] = 32;
      ghostGIDs_check[1] = 46;

      colGIDs_check.resize(6);
      colGIDs_check[0] = 8;
      colGIDs_check[1] = 9;
      colGIDs_check[2] = 12;
      colGIDs_check[3] = 13;
      colGIDs_check[4] = 10;
      colGIDs_check[5] = 14;
    } else if (map->getComm()->getRank() == 3) {
      coarseNodesGIDs_check[0] = 10;
      coarseNodesGIDs_check[1] = 11;
      coarseNodesGIDs_check[2] = 14;
      coarseNodesGIDs_check[3] = 15;

      ghostGIDs_check.resize(0);

      colGIDs_check.resize(4);
      colGIDs_check[0] = 10;
      colGIDs_check[1] = 11;
      colGIDs_check[2] = 14;
      colGIDs_check[3] = 15;
    }
  }
  TEST_EQUALITY(coarseNodesGIDs_check == coarseNodesGIDs, true);
  TEST_EQUALITY(ghostGIDs_check == ghostGIDs, true);
  TEST_EQUALITY(colGIDs_check == colGIDs, true);

}  // BlackBoxGhosts

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlackBoxPFactory, GetNodeInfo, Scalar, LocalOrdinal,
                                  GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  // Creater tester factory
  BlackBoxPFactoryTester<SC, LO, GO, Node> factTester;

  // Create coarse element data
  Teuchos::Array<LO> elementNodesPerDir(3);
  elementNodesPerDir[0] = 4;
  elementNodesPerDir[1] = 4;
  elementNodesPerDir[2] = 3;

  // Instantiate outputs
  int nodeType = 0;
  LO nodeIndex = 0;

  // Create check variables
  bool checkResult            = true;
  std::vector<int> nodeTypes  = {0, 1, 1, 0, 1, 2, 2, 1, 1, 2, 2, 1, 0, 1, 1, 0,
                                 1, 2, 2, 1, 2, 3, 3, 2, 2, 3, 3, 2, 1, 2, 2, 1,
                                 0, 1, 1, 0, 1, 2, 2, 1, 1, 2, 2, 1, 0, 1, 1, 0};
  std::vector<LO> nodeIndices = {0, 0, 1, 1, 2, 0, 1, 3, 4, 2, 3, 5, 2, 6, 7, 3,
                                 8, 4, 5, 9, 6, 0, 1, 7, 8, 2, 3, 9, 10, 10, 11, 11,
                                 4, 12, 13, 5, 14, 12, 13, 15, 16, 14, 15, 17, 6, 18, 19, 7};

  LO currentIndex;
  for (LO k = 0; k < elementNodesPerDir[2]; ++k) {
    for (LO j = 0; j < elementNodesPerDir[1]; ++j) {
      for (LO i = 0; i < elementNodesPerDir[0]; ++i) {
        currentIndex = k * elementNodesPerDir[1] * elementNodesPerDir[0] + j * elementNodesPerDir[0] + i;
        factTester.TestGetNodeInfo(i, j, k, elementNodesPerDir, &nodeType, nodeIndex);
        if ((nodeTypes[currentIndex] != nodeType) || (nodeIndices[currentIndex] != nodeIndex)) {
          checkResult = false;
        }
      }
    }
  }
  TEST_EQUALITY(checkResult, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlackBoxPFactory, Prolongator, Scalar, LocalOrdinal,
                                  GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // used Xpetra lib (for maps and smoothers)
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  const LO numDimensions        = 3;
  const std::string stencilType = "reduced";

  RCP<Matrix> A;
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates;
  RCP<Map> map;
  Array<GO> gNodesPerDim(3);
  Array<LO> lNodesPerDim(3);
  GetProblemData(comm, lib, numDimensions, stencilType, A, coordinates, map, gNodesPerDim,
                 lNodesPerDim);

  // Creater tester factory
  BlackBoxPFactoryTester<SC, LO, GO, Node> factTester;

  LO BlkSize = 1;
  Array<GO> gIndices(3), gCoarseNodesPerDir(3), ghostGIDs, coarseNodesGIDs, colGIDs;
  Array<LO> myOffset(3), coarseRate(3), endRate(3), lCoarseNodesPerDir(3);
  Array<LO> glCoarseNodesPerDir(3), lCoarseElementsPerDir(3);
  Array<bool> ghostInterface(6);
  Array<int> boundaryFlags(3);
  ArrayRCP<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coarseNodes(numDimensions);
  for (int i = 0; i < numDimensions; ++i) {
    coarseRate[i] = 2;
  }
  GO gNumCoarseNodes = 0;
  LO lNumCoarseNodes = 0;
  factTester.TestGetGeometricData(coordinates, coarseRate, gNodesPerDim, lNodesPerDim, 1,
                                  gIndices, gCoarseNodesPerDir, ghostGIDs, coarseNodesGIDs,
                                  colGIDs, gNumCoarseNodes, myOffset, lCoarseNodesPerDir,
                                  glCoarseNodesPerDir, endRate, lNumCoarseNodes, ghostInterface,
                                  coarseNodes, boundaryFlags);

  // From A on local rank, let us construct Aghost
  Array<GO> ghostRowGIDs, ghostColGIDs, nodeSteps(3);
  nodeSteps[0] = 1;
  nodeSteps[1] = gNodesPerDim[0];
  nodeSteps[2] = gNodesPerDim[0] * gNodesPerDim[1];
  Array<LO> range(3);
  GO startingGID = A->getRowMap()->getMinGlobalIndex();
  for (LO dim = 0; dim < 3; ++dim) {
    LO numCoarseNodes = 0;
    if (dim < numDimensions) {
      startingGID -= myOffset[dim] * nodeSteps[dim];
      numCoarseNodes = lCoarseNodesPerDir[dim];
      if (ghostInterface[2 * dim]) {
        ++numCoarseNodes;
      }
      if (ghostInterface[2 * dim + 1]) {
        ++numCoarseNodes;
      }
      if (gIndices[dim] + lNodesPerDim[dim] == gNodesPerDim[dim]) {
        range[dim] = (numCoarseNodes - 2) * coarseRate[dim] + endRate[dim] + 1;
      } else {
        range[dim] = (numCoarseNodes - 1) * coarseRate[dim] + 1;
      }
    } else {
      range[dim] = 1;
    }
  }
  ghostRowGIDs.resize(range[0] * range[1] * range[2] * BlkSize);
  for (LO k = 0; k < range[2]; ++k) {
    for (LO j = 0; j < range[1]; ++j) {
      for (LO i = 0; i < range[0]; ++i) {
        for (LO l = 0; l < BlkSize; ++l) {
          ghostRowGIDs[(k * range[1] * range[0] + j * range[0] + i) * BlkSize + l] = startingGID + (k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + i) * BlkSize + l;
        }
      }
    }
  }

  // Looking at the above loops it is easy to find startingGID for the ghostColGIDs
  Array<GO> startingGlobalIndices(numDimensions), dimStride(numDimensions);
  Array<GO> startingColIndices(numDimensions), finishingColIndices(numDimensions);
  GO colMinGID = 0;
  Array<LO> colRange(numDimensions);
  dimStride[0] = 1;
  for (int dim = 1; dim < numDimensions; ++dim) {
    dimStride[dim] = dimStride[dim - 1] * gNodesPerDim[dim - 1];
  }
  {
    GO tmp = startingGID;
    for (int dim = numDimensions; dim > 0; --dim) {
      startingGlobalIndices[dim - 1] = tmp / dimStride[dim - 1];
      tmp                            = tmp % dimStride[dim - 1];

      if (startingGlobalIndices[dim - 1] > 0) {
        startingColIndices[dim - 1] = startingGlobalIndices[dim - 1] - 1;
      }
      if (startingGlobalIndices[dim - 1] + range[dim - 1] < gNodesPerDim[dim - 1]) {
        finishingColIndices[dim - 1] = startingGlobalIndices[dim - 1] + range[dim - 1];
      } else {
        finishingColIndices[dim - 1] = startingGlobalIndices[dim - 1] + range[dim - 1] - 1;
      }
      colRange[dim - 1] = finishingColIndices[dim - 1] - startingColIndices[dim - 1] + 1;
      colMinGID += startingColIndices[dim - 1] * dimStride[dim - 1];
    }
  }
  ghostColGIDs.resize(colRange[0] * colRange[1] * colRange[2] * BlkSize);
  for (LO k = 0; k < colRange[2]; ++k) {
    for (LO j = 0; j < colRange[1]; ++j) {
      for (LO i = 0; i < colRange[0]; ++i) {
        for (LO l = 0; l < BlkSize; ++l) {
          ghostColGIDs[(k * colRange[1] * colRange[0] + j * colRange[0] + i) * BlkSize + l] = colMinGID + (k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + i) * BlkSize + l;
        }
      }
    }
  }

  RCP<const Map> ghostedRowMap    = Xpetra::MapFactory<LO, GO, NO>::Build(A->getRowMap()->lib(),
                                                                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                          ghostRowGIDs(),
                                                                          A->getRowMap()->getIndexBase(),
                                                                          A->getRowMap()->getComm());
  RCP<const Map> ghostedColMap    = Xpetra::MapFactory<LO, GO, NO>::Build(A->getRowMap()->lib(),
                                                                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                                                          ghostColGIDs(),
                                                                          A->getRowMap()->getIndexBase(),
                                                                          A->getRowMap()->getComm());
  RCP<const Import> ghostImporter = Xpetra::ImportFactory<LO, GO, NO>::Build(A->getRowMap(),
                                                                             ghostedRowMap);
  RCP<const Matrix> Aghost        = Xpetra::MatrixFactory<SC, LO, GO, NO>::Build(A, *ghostImporter,
                                                                                 ghostedRowMap,
                                                                                 ghostedColMap);

  if (comm->getRank() == 0) {
    // Pick on which coarse element the algorithm is tested
    Array<LO> elemInds(3), elementNodesPerDir(3);

    elemInds[0] = 1;
    elemInds[1] = 1;
    elemInds[2] = 1;

    Array<int> elementFlags(3);
    for (int dim = 0; dim < 3; ++dim) {
      if (elemInds[dim] == 0 && elemInds[dim] == lCoarseElementsPerDir[dim] - 1) {
        elementFlags[dim] = boundaryFlags[dim];
      } else if (elemInds[dim] == 0 && (boundaryFlags[dim] == 1 || boundaryFlags[dim] == 3)) {
        elementFlags[dim] += 1;
      } else if ((elemInds[dim] == lCoarseElementsPerDir[dim] - 1) && (boundaryFlags[dim] == 2 || boundaryFlags[dim] == 3)) {
        elementFlags[dim] += 2;
      } else {
        elementFlags[dim] = 0;
      }

      // Compute the number of nodes in the current element.
      if (dim < numDimensions) {
        if ((elemInds[dim] == lCoarseElementsPerDir[dim]) && (gIndices[dim] + lNodesPerDim[dim] == gNodesPerDim[dim])) {
          elementNodesPerDir[dim] = endRate[dim] + 1;
        } else {
          elementNodesPerDir[dim] = coarseRate[dim] + 1;
        }
      } else {
        elementNodesPerDir[dim] = 1;
      }
    }
    LO numNodesInElement = elementNodesPerDir[0] * elementNodesPerDir[1] * elementNodesPerDir[2];

    Teuchos::SerialDenseMatrix<LO, SC> Pi, Pf, Pe;
    Array<LO> dofType(numNodesInElement * BlkSize), lDofInd(numNodesInElement * BlkSize);
    factTester.TestComputeLocalEntries(Aghost, coarseRate, endRate, BlkSize, elemInds,
                                       lCoarseElementsPerDir, numDimensions, range, gNodesPerDim,
                                       gIndices, lCoarseNodesPerDir, ghostInterface,
                                       elementFlags, "reduced", "coupled", elementNodesPerDir,
                                       numNodesInElement, colGIDs, Pi, Pf, Pe, dofType, lDofInd);

    elemInds[0] = 0;
    elemInds[1] = 0;
    elemInds[2] = 0;

    elementFlags[0] = 0;
    elementFlags[1] = 0;
    elementFlags[2] = 0;
    for (int dim = 0; dim < 3; ++dim) {
      if (elemInds[dim] == 0 && elemInds[dim] == lCoarseElementsPerDir[dim] - 1) {
        elementFlags[dim] = boundaryFlags[dim];
      } else if (elemInds[dim] == 0 && (boundaryFlags[dim] == 1 || boundaryFlags[dim] == 3)) {
        elementFlags[dim] += 1;
      } else if ((elemInds[dim] == lCoarseElementsPerDir[dim] - 1) && (boundaryFlags[dim] == 2 || boundaryFlags[dim] == 3)) {
        elementFlags[dim] += 2;
      } else {
        elementFlags[dim] = 0;
      }

      // Compute the number of nodes in the current element.
      if (dim < numDimensions) {
        if ((elemInds[dim] == lCoarseElementsPerDir[dim]) && (gIndices[dim] + lNodesPerDim[dim] == gNodesPerDim[dim])) {
          elementNodesPerDir[dim] = endRate[dim] + 1;
        } else {
          elementNodesPerDir[dim] = coarseRate[dim] + 1;
        }
      } else {
        elementNodesPerDir[dim] = 1;
      }
    }
    numNodesInElement = elementNodesPerDir[0] * elementNodesPerDir[1] * elementNodesPerDir[2];

    dofType.resize(numNodesInElement * BlkSize);
    lDofInd.resize(numNodesInElement * BlkSize);
    factTester.TestComputeLocalEntries(Aghost, coarseRate, endRate, BlkSize, elemInds,
                                       lCoarseElementsPerDir, numDimensions, range, gNodesPerDim,
                                       gIndices, lCoarseNodesPerDir, ghostInterface,
                                       elementFlags, "reduced", "coupled", elementNodesPerDir,
                                       numNodesInElement, colGIDs, Pi, Pf, Pe, dofType, lDofInd);
  }

}  // Prolongator

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlackBoxPFactory, BBPoisson, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // used Xpetra lib (for maps and smoothers)
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  const LO numDimensions        = 3;
  const std::string stencilType = "reduced";

  RCP<Matrix> A;
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates;
  RCP<Map> map;
  Array<GO> gNodesPerDim(3);
  Array<LO> lNodesPerDim(3);
  GetProblemData(comm, lib, numDimensions, stencilType, A, coordinates, map, gNodesPerDim,
                 lNodesPerDim);

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
  nullSpace->putScalar((Scalar)1.0);
  Teuchos::Array<magnitude_type> norms(1);
  nullSpace->norm1(norms);
  if (comm->getRank() == 0) {
    out << "||NS|| = " << norms[0] << std::endl;
  }

  // fill hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  // create the factory manager and the factories
  FactoryManager M;

  RCP<Factory> Pfact     = rcp(new BlackBoxPFactory());
  RCP<Factory> Rfact     = rcp(new TransPFactory());
  RCP<Factory> Tfact     = rcp(new CoordinatesTransferFactory());
  RCP<RAPFactory> Acfact = rcp(new RAPFactory());
  RCP<Factory> NSfact    = rcp(new NullspaceFactory());

  // Set paramters needed by the factories
  Pfact->SetParameter("Coarsen", Teuchos::ParameterEntry(std::string("{2,2,2}")));
  Pfact->SetParameter("stencil type", Teuchos::ParameterEntry(std::string("reduced")));
  Pfact->SetParameter("block strategy", Teuchos::ParameterEntry(std::string("coupled")));

  Tfact->SetParameter("Geometric", Teuchos::ParameterEntry(true));

  // Set interfactory dependencies
  NSfact->SetFactory("Nullspace", Pfact);

  Acfact->AddTransferFactory(Tfact);

  // Set default factories in the manager
  M.SetFactory("P", Pfact);
  M.SetFactory("R", Rfact);
  M.SetFactory("A", Acfact);
  M.SetFactory("Nullspace", NSfact);
  M.SetFactory("gNodesPerDim", Tfact);
  M.SetFactory("lNodesPerDim", Tfact);
  M.SetFactory("Coordinates", Tfact);
  M.SetFactory("coarseCoordinates", Pfact);
  M.SetFactory("gCoarseNodesPerDim", Pfact);
  M.SetFactory("lCoarseNodesPerDim", Pfact);

  // setup smoothers
  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
  smootherParamList.set("relaxation: sweeps", (LocalOrdinal)1);
  smootherParamList.set("relaxation: damping factor", (Scalar)1.0);
  RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("RELAXATION", smootherParamList));
  RCP<SmootherFactory> SmooFact    = rcp(new SmootherFactory(smooProto));
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

  // Set smoothers and coarse solver in the manager
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarseSolveFact);

  // Populate data on the finest level of the hierarchy
  RCP<Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A", A);                        // set fine level matrix
  Finest->Set("Nullspace", nullSpace);        // set null space information for finest level
  Finest->Set("Coordinates", coordinates);    // set fine level coordinates
  Finest->Set("gNodesPerDim", gNodesPerDim);  // set BlackBoxPFactory specific info
  Finest->Set("lNodesPerDim", lNodesPerDim);  // set BlackBoxPFactory specific info
  Finest->SetFactoryManager(Teuchos::rcpFromRef(M));

  // Setup the hierarchy
  H->SetMaxCoarseSize(10);
  H->Setup(M, 0, 2);

  // Extract the prolongator operator
  RCP<Level> lvl1                          = H->GetLevel(1);
  RCP<Xpetra::Matrix<SC, LO, GO, Node> > P = lvl1->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
  RCP<CrsMatrix> PCrs                      = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

}  // BBPoisson

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlackBoxPFactory, Constructor, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlackBoxPFactory, BlackBoxGhosts, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlackBoxPFactory, GetNodeInfo, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlackBoxPFactory, Prolongator, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlackBoxPFactory, BBPoisson, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
