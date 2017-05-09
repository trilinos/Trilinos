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

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class BlackBoxPFactoryTester {
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor
    BlackBoxPFactoryTester() { }

    //! Destructor.
    virtual ~BlackBoxPFactoryTester() { }
    //@}

    void TestGetGeometricData(Teuchos::RCP<Xpetra::MultiVector<double, LO, GO, NO> >& coordinates, Array<LO> coarseRate, Array<GO> gNodesPerDim, Array<LO> lNodesPerDim, LO BlkSize,
                              Array<GO>& gIndices, Array<GO>& gCoarseNodesPerDir, Array<GO>& ghostGIDs, Array<GO>& coarseNodesGIDs, Array<GO>& colGIDs, GO& gNumCoarseNodes,
                              Array<LO>& myOffset, Array<LO>& lCoarseNodesPerDir, Array<LO>& endRate, LO& lNumCoarseNodes,
                              Array<bool>& ghostInterface, ArrayRCP<Array<double> > coarseNodes) const{
      // Call the method to be tested.
      MueLu::BlackBoxPFactory<SC,LO,GO,Node> mybbmgPFactory;
      mybbmgPFactory.GetGeometricData(coordinates, coarseRate, gNodesPerDim, lNodesPerDim, BlkSize,
                                      gIndices, myOffset, ghostInterface, endRate, gCoarseNodesPerDir,
                                      lCoarseNodesPerDir, ghostGIDs, coarseNodesGIDs, colGIDs,
                                      gNumCoarseNodes, lNumCoarseNodes, coarseNodes);
    };

  };

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  void GetProblemData(RCP<const Teuchos::Comm<int> >& comm, const Xpetra::UnderlyingLib lib, const LocalOrdinal numDimensions,
                      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Op,
                      RCP<Xpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> >& Coordinates,
                      RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal> >& map,
                      Array<GlobalOrdinal>& gNodesPerDim, Array<LocalOrdinal>& lNodesPerDim) {
#include "MueLu_UseShortNames.hpp"

    GO nx = 5;
    GO ny = 5;
    GO nz = (numDimensions == 2 ? 1 : 5);
    GO gNumPoints = nx*ny*nz;

    gNodesPerDim[0] = nx;
    gNodesPerDim[1] = ny;
    gNodesPerDim[2] = nz;

    GO myOffset;
    if(comm->getSize() == 1) {
      myOffset = 0;
      lNodesPerDim[0] = nx;
      lNodesPerDim[1] = ny;
      lNodesPerDim[2] = nz;
    } else if(comm->getSize() == 4) {
      if(comm->getRank() == 0) {
        if(numDimensions == 2) {
          myOffset = 0;
          lNodesPerDim[0] = 3;
          lNodesPerDim[1] = 3;
          lNodesPerDim[2] = 1;
        } else if(numDimensions == 3) {
          myOffset = 0;
          lNodesPerDim[0] = nx;
          lNodesPerDim[1] = 3;
          lNodesPerDim[2] = 3;
        }
      } else if(comm->getRank() == 1) {
        if(numDimensions == 2) {
          myOffset = 3;
          lNodesPerDim[0] = 2;
          lNodesPerDim[1] = 3;
          lNodesPerDim[2] = 1;
        } else if(numDimensions == 3) {
          myOffset = 15;
          lNodesPerDim[0] = nx;
          lNodesPerDim[1] = 2;
          lNodesPerDim[2] = 3;
        }
      } else if(comm->getRank() == 2) {
        if(numDimensions == 2) {
          myOffset = 15;
          lNodesPerDim[0] = 3;
          lNodesPerDim[1] = 2;
          lNodesPerDim[2] = 1;
        } else if(numDimensions == 3) {
          myOffset = 75;
          lNodesPerDim[0] = nx;
          lNodesPerDim[1] = 3;
          lNodesPerDim[2] = 2;
        }
      } else if(comm->getRank() == 3) {
        if(numDimensions == 2) {
          myOffset = 18;
          lNodesPerDim[0] = 2;
          lNodesPerDim[1] = 2;
          lNodesPerDim[2] = 1;
        } else if(numDimensions == 3) {
          myOffset = 90;
          lNodesPerDim[0] = nx;
          lNodesPerDim[1] = 2;
          lNodesPerDim[2] = 2;
        }
      }
    }

    GO myZoffset, myYoffset, myXoffset;
    if(numDimensions == 2) {
      myZoffset = 0;
      myYoffset = myOffset / gNodesPerDim[0];
      myXoffset = myOffset % gNodesPerDim[0];
    } else if(numDimensions == 3) {
      myZoffset = myOffset / (gNodesPerDim[1]*gNodesPerDim[0]);
      myYoffset = (myOffset - myZoffset*gNodesPerDim[1]*gNodesPerDim[0]) / gNodesPerDim[0];
      myXoffset = (myOffset - myZoffset*gNodesPerDim[1]*gNodesPerDim[0]) % gNodesPerDim[0];
    }
    LO lNumPoints = lNodesPerDim[0]*lNodesPerDim[1]*lNodesPerDim[2];

    // Construct map and local coordinates
    Teuchos::Array<GO>     myGIDs(lNumPoints);
    Teuchos::Array<double> myXCoords(lNumPoints);
    Teuchos::Array<double> myYCoords(lNumPoints);
    Teuchos::Array<double> myZCoords(lNumPoints);
    for(LO k = 0; k < lNodesPerDim[2]; ++k) {
      for(LO j = 0; j < lNodesPerDim[1]; ++j) {
        for(LO i = 0; i < lNodesPerDim[0]; ++i) {
          myGIDs[k*lNodesPerDim[1]*lNodesPerDim[0] + j*lNodesPerDim[0] + i] = myOffset + k*gNodesPerDim[1]*gNodesPerDim[0] + j*gNodesPerDim[0] + i;
          myXCoords[k*lNodesPerDim[1]*lNodesPerDim[0] + j*lNodesPerDim[0] + i] = (i + myXoffset) / Teuchos::as<double>(gNodesPerDim[0] - 1);
          myYCoords[k*lNodesPerDim[1]*lNodesPerDim[0] + j*lNodesPerDim[0] + i] = (j + myYoffset) / Teuchos::as<double>(gNodesPerDim[1] - 1);
          myZCoords[k*lNodesPerDim[1]*lNodesPerDim[0] + j*lNodesPerDim[0] + i] = (k + myZoffset) / Teuchos::as<double>(gNodesPerDim[2] - 1);
        }
      }
    }

    Teuchos::Array<Teuchos::ArrayView<const double> > myCoordinates(numDimensions);
    if(numDimensions == 1) {
      myCoordinates[0] = myXCoords();
    } else if(numDimensions == 2) {
      myCoordinates[0] = myXCoords();
      myCoordinates[1] = myYCoords();
    } else if(numDimensions == 3) {
      myCoordinates[0] = myXCoords();
      myCoordinates[1] = myYCoords();
      myCoordinates[2] = myZCoords();
    }

    // Create the map and store coordinates using the above array views
    map         = MapFactory::Build(lib, gNumPoints, myGIDs(), 0, comm);
    Coordinates = Xpetra::MultiVectorFactory<double,LO,GO,NO>::Build(map, myCoordinates(), numDimensions);

    // small parameter list for Galeri
    Teuchos::ParameterList problemParamList;
    if(numDimensions == 1) {
      problemParamList.set("nx",gNodesPerDim[0]);
    } else if(numDimensions == 2) {
      problemParamList.set("nx",gNodesPerDim[0]);
      problemParamList.set("ny",gNodesPerDim[1]);
    } else if(numDimensions == 3) {
      problemParamList.set("nx",gNodesPerDim[0]);
      problemParamList.set("ny",gNodesPerDim[1]);
      problemParamList.set("nz",gNodesPerDim[2]);
    }
    // problemParamList.set("keepBCs", true);

    // create Poisson problem and matrix
    if(numDimensions == 1) {
      Galeri::Xpetra::Laplace1DProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector> Problem(problemParamList, map);
      Op = Problem.BuildMatrix();
    } else if(numDimensions == 2) {
      Galeri::Xpetra::Laplace2DProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector> Problem(problemParamList, map);
      Op = Problem.BuildMatrix();
    } else if(numDimensions == 3) {
      Galeri::Xpetra::Laplace3DProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector> Problem(problemParamList, map);
      Op = Problem.BuildMatrix();
    }

  }

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlackBoxPFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    out << "version: " << MueLu::Version() << std::endl;

    RCP<BlackBoxPFactory> bbmgPFact = rcp(new BlackBoxPFactory);
    TEST_EQUALITY(bbmgPFact != Teuchos::null, true);

  } //Constructor

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlackBoxPFactory, GeometricGhosts, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    RCP<Teuchos::FancyOStream> fancy = getFancyOStream(rcpFromRef(std::cout));
    fancy->setShowAllFrontMatter(false).setShowProcRank(true);
    Teuchos::FancyOStream& out2 = *fancy;

    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // used Xpetra lib (for maps and smoothers)
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    LO numDimensions = 2;

    RCP<Matrix> Op;
    RCP<Xpetra::MultiVector<double,LO,GO,NO> > coordinates;
    RCP<Map> map;
    Array<GO> gNodesPerDim(3);
    Array<LO> lNodesPerDim(3);
    GetProblemData(comm, lib, numDimensions, Op, coordinates, map, gNodesPerDim, lNodesPerDim);

    Array<LO> coarseRate(3);
    coarseRate[0] = 2;
    coarseRate[1] = 2;
    if(numDimensions == 2) {
      coarseRate[2] = 1;
    } else if(numDimensions == 3) {
      coarseRate[2] = 2;
    }

    MueLu::BlackBoxPFactory<SC,LO,GO,Node> bbmgPFact;
    BlackBoxPFactoryTester<SC,LO,GO,Node> factTester;
    Array<GO> gIndices(3), gCoarseNodesPerDir(3), ghostGIDs, coarseNodesGIDs, colGIDs;
    Array<LO> myOffset(3), lCoarseNodesPerDir(3), endRate(3);
    Array<bool> ghostInterface(6);
    Array<ArrayView<const double> > fineNodes(numDimensions);
    ArrayRCP<Array<double> > coarseNodes(numDimensions);
    GO gNumCoarseNodes = 0;
    LO lNumCoarseNodes = 0;
    for(LO dim = 0; dim < numDimensions; ++dim) {
      fineNodes[dim] = coordinates->getData(dim)();
    }
    factTester.TestGetGeometricData(coordinates, coarseRate, gNodesPerDim, lNodesPerDim, 1,
                                    gIndices, gCoarseNodesPerDir, ghostGIDs, coarseNodesGIDs,
                                    colGIDs, gNumCoarseNodes, myOffset, lCoarseNodesPerDir,
                                    endRate, lNumCoarseNodes, ghostInterface, coarseNodes);

    Array<GO> ghostGIDs_check, coarseNodesGIDs_check, colGIDs_check;
    if(map->getComm()->getSize() == 1) {
      coarseNodesGIDs_check.resize(9);
      for(LO i = 0; i < 9; ++i) {
        coarseNodesGIDs_check[i] = i;
      }
      colGIDs_check = coarseNodesGIDs_check;
      ghostGIDs_check.resize(0);
    } else if(map->getComm()->getSize() == 4) {
      colGIDs_check.resize(4);
      if(map->getComm()->getRank() == 0) {
        coarseNodesGIDs_check.resize(4);
        coarseNodesGIDs_check[0] = 0;
        coarseNodesGIDs_check[1] = 1;
        coarseNodesGIDs_check[2] = 3;
        coarseNodesGIDs_check[3] = 4;

        ghostGIDs_check.resize(0);

        colGIDs_check[0] = 0;
        colGIDs_check[1] = 1;
        colGIDs_check[2] = 3;
        colGIDs_check[3] = 4;
      } else if(map->getComm()->getRank() == 1) {
        coarseNodesGIDs_check.resize(2);
        coarseNodesGIDs_check[0] = 2;
        coarseNodesGIDs_check[1] = 5;

        ghostGIDs_check.resize(2);
        ghostGIDs_check[0] = 2;
        ghostGIDs_check[1] = 12;

        colGIDs_check[0] = 2;
        colGIDs_check[1] = 5;
        colGIDs_check[2] = 1;
        colGIDs_check[3] = 4;
      } else if(map->getComm()->getRank() == 2) {
        coarseNodesGIDs_check.resize(2);
        coarseNodesGIDs_check[0] = 6;
        coarseNodesGIDs_check[1] = 7;

        ghostGIDs_check.resize(2);
        ghostGIDs_check[0] = 10;
        ghostGIDs_check[1] = 12;

        colGIDs_check[0] = 6;
        colGIDs_check[1] = 7;
        colGIDs_check[2] = 3;
        colGIDs_check[3] = 4;
      } else if(map->getComm()->getRank() == 3) {
        coarseNodesGIDs_check.resize(1);
        coarseNodesGIDs_check[0] = 8;

        ghostGIDs_check.resize(3);
        ghostGIDs_check[0] = 12;
        ghostGIDs_check[1] = 14;
        ghostGIDs_check[2] = 22;

        colGIDs_check[0] = 8;
        colGIDs_check[1] = 4;
        colGIDs_check[2] = 5;
        colGIDs_check[3] = 7;
      }
    }
    TEST_EQUALITY(coarseNodesGIDs_check == coarseNodesGIDs, true);
    TEST_EQUALITY(ghostGIDs_check == ghostGIDs, true);
    TEST_EQUALITY(colGIDs_check == colGIDs, true);

  } //GeometricGhosts

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlackBoxPFactory, CoarseNodes, Scalar, LocalOrdinal, GlobalOrdinal, Node)
  {
#   include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,Node);

    RCP<Teuchos::FancyOStream> fancy = getFancyOStream(rcpFromRef(std::cout));
    fancy->setShowAllFrontMatter(false).setShowProcRank(true);
    Teuchos::FancyOStream& out2 = *fancy;

    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // used Xpetra lib (for maps and smoothers)
    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

    LO numDimensions = 2;

    RCP<Matrix> Op;
    RCP<Xpetra::MultiVector<double,LO,GO,NO> > coordinates;
    RCP<Map> map;
    Array<GO> gNodesPerDim(3);
    Array<LO> lNodesPerDim(3);
    GetProblemData(comm, lib, numDimensions, Op, coordinates, map, gNodesPerDim, lNodesPerDim);

    Xpetra::IO<double,LO,GO,NO>::Write("/home/lberge/Desktop/fineMap.m", *map);
    Xpetra::IO<double,LO,GO,NO>::Write("/home/lberge/Desktop/fineCoords.m", *coordinates);

    // build nullspace
    RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
    nullSpace->putScalar( (Scalar) 1.0);

    // Setup a fine and a coarse level
    Level fineLevel, coarseLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    fineLevel.setDefaultVerbLevel(Teuchos::VERB_HIGH);
    fineLevel.Set("A", Op);                       // set fine level matrix
    fineLevel.Set("Nullspace",    nullSpace);     // set null space information for finest level
    fineLevel.Set("Coordinates",  coordinates);   // set fine level coordinates
    fineLevel.Set("gNodesPerDim", gNodesPerDim);  // set GeneralGeometricPFactory specific info
    fineLevel.Set("lNodesPerDim", lNodesPerDim);  // set GeneralGeometricPFactory specific info

    // Black Box factory ParameterList
    Teuchos::ParameterList BBParams;
    BBParams.set("Coarsen", "{2,2,2}");
    BBParams.set("axisPermutation", "{0,1,2}");

    // create the black box factory
    RCP<BlackBoxPFactory> bbmgPFact = rcp(new BlackBoxPFactory);
    bbmgPFact->SetParameterList(BBParams);
    coarseLevel.Request("P",            bbmgPFact.get());  // request Ptent
    // coarseLevel.Request("Nullspace",    bbmgPFact.get());
    // coarseLevel.Request("Coordinates",  bbmgPFact.get());
    coarseLevel.Request(*bbmgPFact);
    bbmgPFact->Build(fineLevel,coarseLevel);

    // TEST_EQUALITY(bbmgPFact != Teuchos::null, true);

  } //CoarseNodes


#  define MUELU_ETI_GROUP(Scalar, LO, GO, Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlackBoxPFactory,Constructor,    Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlackBoxPFactory,GeometricGhosts,Scalar,LO,GO,Node) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlackBoxPFactory,CoarseNodes,    Scalar,LO,GO,Node)

#include <MueLu_ETI_4arg.hpp>


} // namespace MueLuTests
