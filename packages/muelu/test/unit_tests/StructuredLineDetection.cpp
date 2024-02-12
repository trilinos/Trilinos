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

#include "MueLu_StructuredLineDetectionFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_GeneralGeometricPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"

namespace MueLuTests {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GetProblemData(RCP<const Teuchos::Comm<int> >& comm, const Xpetra::UnderlyingLib lib, const LocalOrdinal numDimensions,
                    RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Op,
                    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> >& Coordinates,
                    RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& map,
                    Array<GlobalOrdinal>& gNodesPerDim, Array<LocalOrdinal>& lNodesPerDim) {
#include "MueLu_UseShortNames.hpp"

  GO nx         = 5;
  GO ny         = 5;
  GO nz         = (numDimensions == 2 ? 1 : 5);
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
        lNodesPerDim[0] = 3;
        lNodesPerDim[1] = 3;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 0;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 3;
        lNodesPerDim[2] = 3;
      }
    } else if (comm->getRank() == 1) {
      if (numDimensions == 2) {
        myOffset        = 4;
        lNodesPerDim[0] = 2;
        lNodesPerDim[1] = 3;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 15;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 2;
        lNodesPerDim[2] = 3;
      }
    } else if (comm->getRank() == 2) {
      if (numDimensions == 2) {
        myOffset        = 21;
        lNodesPerDim[0] = 3;
        lNodesPerDim[1] = 2;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 75;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 3;
        lNodesPerDim[2] = 2;
      }
    } else if (comm->getRank() == 3) {
      if (numDimensions == 2) {
        myOffset        = 25;
        lNodesPerDim[0] = 2;
        lNodesPerDim[1] = 2;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 90;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 2;
        lNodesPerDim[2] = 2;
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
        myGIDs[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i]    = myOffset + k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + i;
        myXCoords[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] = (i + myXoffset) / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[0] - 1);
        myYCoords[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] = (j + myYoffset) / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[1] - 1);
        myZCoords[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] = (k + myZoffset) / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[2] - 1);
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
  Coordinates = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(map, myCoordinates(), numDimensions);

  // small parameter list for Galeri
  Teuchos::ParameterList problemParamList;
  if (numDimensions == 1) {
    problemParamList.set("nx", gNodesPerDim[0]);
  } else if (numDimensions == 2) {
    problemParamList.set("nx", gNodesPerDim[0]);
    problemParamList.set("ny", gNodesPerDim[1]);
  } else if (numDimensions == 3) {
    problemParamList.set("nx", gNodesPerDim[0]);
    problemParamList.set("ny", gNodesPerDim[1]);
    problemParamList.set("nz", gNodesPerDim[2]);
  }
  // problemParamList.set("keepBCs", true);

  // create Poisson problem and matrix
  if (numDimensions == 1) {
    Galeri::Xpetra::Laplace1DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> Problem(problemParamList, map);
    Op = Problem.BuildMatrix();
  } else if (numDimensions == 2) {
    Galeri::Xpetra::Laplace2DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> Problem(problemParamList, map);
    Op = Problem.BuildMatrix();
  } else if (numDimensions == 3) {
    Galeri::Xpetra::Laplace3DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> Problem(problemParamList, map);
    Op = Problem.BuildMatrix();
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredLineDetectionFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<StructuredLineDetectionFactory> lineDetectionFact = rcp(new StructuredLineDetectionFactory);
  TEST_EQUALITY(lineDetectionFact != Teuchos::null, true);

}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(StructuredLineDetectionFactory, LabelLines, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // used Xpetra lib (for maps and smoothers)
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  LO numDimensions = 3;

  LO maxLevels = 2, maxIter = 10;
  RCP<Matrix> Op;
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates;
  RCP<Xpetra::Map<LO, GO, NO> > map;
  Array<GO> gNodesPerDim(3);
  Array<LO> lNodesPerDim(3);
  GetProblemData<SC, LO, GO, NO>(comm, lib, numDimensions, Op, coordinates, map, gNodesPerDim, lNodesPerDim);

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

  RCP<Factory> Pfact     = rcp(new GeneralGeometricPFactory());
  RCP<Factory> Rfact     = rcp(new TransPFactory());
  RCP<Factory> LDfact    = rcp(new StructuredLineDetectionFactory());
  RCP<Factory> Tfact     = rcp(new CoordinatesTransferFactory());
  RCP<RAPFactory> Acfact = rcp(new RAPFactory());
  RCP<Factory> NSfact    = rcp(new NullspaceFactory());

  // Set paramters needed by the factories
  Pfact->SetParameter("Coarsen", Teuchos::ParameterEntry(std::string("{2,2,2}")));
  Pfact->SetParameter("order", Teuchos::ParameterEntry(1));

  LDfact->SetParameter("orientation", Teuchos::ParameterEntry(std::string("X")));

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
  M.SetFactory("CoarseNumZLayers", LDfact);
  M.SetFactory("LineDetection_VertLineIds", LDfact);

  // setup smoothers
  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "Jacobi");
  smootherParamList.set("relaxation: sweeps", (LocalOrdinal)1);
  smootherParamList.set("relaxation: damping factor", Teuchos::ScalarTraits<Scalar>::one());
  RCP<SmootherPrototype> smooProto = rcp(new TrilinosSmoother("LINESMOOTHING_BANDEDRELAXATION", smootherParamList));
  RCP<SmootherFactory> SmooFact    = rcp(new SmootherFactory(smooProto));
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  RCP<SmootherFactory> coarseSolveFact = rcp(new SmootherFactory(smooProto, Teuchos::null));

  // Set smoothers and coarse solver in the manager
  M.SetFactory("Smoother", SmooFact);
  M.SetFactory("CoarseSolver", coarseSolveFact);

  // Populate data on the finest level of the hierarchy
  RCP<Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A", Op);                       // set fine level matrix
  Finest->Set("Nullspace", nullSpace);        // set null space information for finest level
  Finest->Set("Coordinates", coordinates);    // set fine level coordinates
  Finest->Set("gNodesPerDim", gNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->Set("lNodesPerDim", lNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->SetFactoryManager(Teuchos::rcpFromRef(M));

  // Setup the hierarchy
  H->SetMaxCoarseSize(10);
  H->Setup(M, 0, maxLevels);

  out << "MueLu setup done!" << std::endl;

  // Define RHS and solution vector
  RCP<MultiVector> X   = MultiVectorFactory::Build(map, 1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

  X->putScalar(1.0);
  X->norm2(norms);
  Op->apply(*X, *RHS, Teuchos::NO_TRANS, (Scalar)1.0, (Scalar)0.0);
  X->putScalar((Scalar)0.0);
  out << "Starting solver iterations" << std::endl;
  for (LO i = 0; i < 10; ++i) {
    H->Iterate(*RHS, *X, maxIter);
    X->norm2(norms);
    if (comm->getRank() == 0) {
      out << "||RHS[" << i << "]|| = " << norms << std::endl;
    }
    Op->apply(*X, *RHS, Teuchos::NO_TRANS, (Scalar)-1.0, (Scalar)1.0);
    if (norms[0] < 1.0e-7) {
      break;
    }
  }
  TEST_EQUALITY((norms[0] < 1.0e-7), true);
}  // LabelLines

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredLineDetectionFactory, Constructor, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(StructuredLineDetectionFactory, LabelLines, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
