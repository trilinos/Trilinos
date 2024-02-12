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

#include <complex>

#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_GeneralGeometricPFactory.hpp"
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

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class GeneralGeometricPFactoryTester {
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  GeneralGeometricPFactoryTester() {}

  //! Destructor.
  virtual ~GeneralGeometricPFactoryTester() {}
  //@}

  void TestComputeLinearInterpolationStencil(const GeneralGeometricPFactory& fac,
                                             const LO numDimension,
                                             const Array<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coord,
                                             std::vector<double>& stencil)
      const {
    // Call the method to be tested.
    MueLu::GeneralGeometricPFactory<SC, LO, GO, Node> myGGPFactory;
    myGGPFactory.ComputeLinearInterpolationStencil(numDimension, coord, stencil);
  };
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void GGGetProblemData(RCP<const Teuchos::Comm<int> >& comm, const Xpetra::UnderlyingLib lib,
                      const LocalOrdinal numDimensions, const std::string mode,
                      RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& Op,
                      RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LocalOrdinal, GlobalOrdinal, Node> >& Coordinates,
                      RCP<Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& map,
                      Array<GlobalOrdinal>& gNodesPerDim, Array<LocalOrdinal>& lNodesPerDim) {
#include "MueLu_UseShortNames.hpp"

  GO nx         = 9;
  GO ny         = 9;
  GO nz         = (numDimensions < 3 ? 1 : 9);
  GO gNumPoints = nx * ny * nz;

  gNodesPerDim[0] = nx;
  gNodesPerDim[1] = ny;
  gNodesPerDim[2] = nz;

  GO myOffset = 0, myGIDOffset = 0;
  if (comm->getSize() == 1) {
    myOffset        = 0;
    lNodesPerDim[0] = nx;
    lNodesPerDim[1] = ny;
    lNodesPerDim[2] = nz;
  } else if (comm->getSize() == 4) {
    if (comm->getRank() == 0) {
      if (numDimensions == 2) {
        myOffset        = 0;
        myGIDOffset     = 0;
        lNodesPerDim[0] = 5;
        lNodesPerDim[1] = 5;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 0;
        myGIDOffset     = 0;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 5;
        lNodesPerDim[2] = 5;
      }
    } else if (comm->getRank() == 1) {
      if (numDimensions == 2) {
        myOffset        = 4;
        myGIDOffset     = 16;
        lNodesPerDim[0] = 4;
        lNodesPerDim[1] = 5;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 45;
        myGIDOffset     = 225;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 5;
      }
    } else if (comm->getRank() == 2) {
      if (numDimensions == 2) {
        myOffset        = 21;
        myGIDOffset     = 28;
        lNodesPerDim[0] = 5;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 405;
        myGIDOffset     = 405;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 5;
        lNodesPerDim[2] = 4;
      }
    } else if (comm->getRank() == 3) {
      if (numDimensions == 2) {
        myOffset        = 25;
        myGIDOffset     = 44;
        lNodesPerDim[0] = 4;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 1;
      } else if (numDimensions == 3) {
        myOffset        = 450;
        myGIDOffset     = 585;
        lNodesPerDim[0] = nx;
        lNodesPerDim[1] = 4;
        lNodesPerDim[2] = 4;
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
        if (mode == "Global Lexicographic") {
          myGIDs[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] = myOffset + k * gNodesPerDim[1] * gNodesPerDim[0] + j * gNodesPerDim[0] + i;
        } else if (mode == "Local Lexicographic") {
          myGIDs[k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i] = myGIDOffset + k * lNodesPerDim[1] * lNodesPerDim[0] + j * lNodesPerDim[0] + i;
        }
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
  // problemList.set("keepBCs", true);

  // create Poisson problem and matrix
  if (numDimensions == 1) {
    Galeri::Xpetra::Laplace1DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> Problem(problemList,
                                                                                          map);
    Op = Problem.BuildMatrix();
  } else if (numDimensions == 2) {
    Galeri::Xpetra::Laplace2DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> Problem(problemList,
                                                                                          map);
    Op = Problem.BuildMatrix();
  } else if (numDimensions == 3) {
    Galeri::Xpetra::Laplace3DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> Problem(problemList,
                                                                                          map);
    Op = Problem.BuildMatrix();
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(GeneralGeometricPFactory, Constructor, Scalar, LocalOrdinal,
                                  GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<GeneralGeometricPFactory> ggPFact = rcp(new GeneralGeometricPFactory);
  TEST_EQUALITY(ggPFact != Teuchos::null, true);

}  // Constructor

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(GeneralGeometricPFactory, LinearInterpolation, Scalar,
                                  LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  MueLu::GeneralGeometricPFactory<SC, LO, GO, Node> ggPFact;
  GeneralGeometricPFactoryTester<SC, LO, GO, Node> factTester;

  LO numDimension = 3;
  Array<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coord(9);
  for (int node = 0; node < 9; ++node) {
    coord[node].resize(3);
  }
  coord[0][0] = 0.3;
  coord[0][1] = 0.9;
  coord[0][2] = 0.1;
  coord[1][0] = 0.0;
  coord[1][1] = 0.0;
  coord[1][2] = 0.0;
  coord[2][0] = 1.0;
  coord[2][1] = 0.0;
  coord[2][2] = 0.0;
  coord[3][0] = 0.0;
  coord[3][1] = 1.0;
  coord[3][2] = 0.0;
  coord[4][0] = 1.0;
  coord[4][1] = 1.0;
  coord[4][2] = 0.0;
  coord[5][0] = 0.0;
  coord[5][1] = 0.0;
  coord[5][2] = 1.0;
  coord[6][0] = 1.0;
  coord[6][1] = 0.0;
  coord[6][2] = 1.0;
  coord[7][0] = 0.0;
  coord[7][1] = 1.0;
  coord[7][2] = 1.0;
  coord[8][0] = 1.0;
  coord[8][1] = 1.0;
  coord[8][2] = 1.0;
  std::vector<double> stencil(8);
  factTester.TestComputeLinearInterpolationStencil(ggPFact, numDimension, coord, stencil);

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType x = 0.0, y = 0.0, z = 0.0;
  for (LO i = 0; i < 8; ++i) {
    x += stencil[i] * coord[i + 1][0];
    y += stencil[i] * coord[i + 1][1];
    z += stencil[i] * coord[i + 1][2];
  }

  TEST_EQUALITY((std::abs(x - coord[0][0]) < 1e-5) && (std::abs(y - coord[0][1]) < 1e-5) && (std::abs(z - coord[0][2]) < 1e-5), true);

}  // LinearInterpolation

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(GeneralGeometricPFactory, LinearExtrapolation, Scalar,
                                  LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  MueLu::GeneralGeometricPFactory<SC, LO, GO, Node> ggPFact;
  GeneralGeometricPFactoryTester<SC, LO, GO, Node> factTester;

  LO numDimension = 3;
  Array<Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> > coord(9);
  for (int node = 0; node < 9; ++node) {
    coord[node].resize(3);
  }
  coord[0][0] = 1.1;
  coord[0][1] = 0.3;
  coord[0][2] = 0.8;
  coord[1][0] = 0.0;
  coord[1][1] = 0.0;
  coord[1][2] = 0.0;
  coord[2][0] = 1.0;
  coord[2][1] = 0.0;
  coord[2][2] = 0.0;
  coord[3][0] = 0.0;
  coord[3][1] = 1.0;
  coord[3][2] = 0.0;
  coord[4][0] = 1.0;
  coord[4][1] = 1.0;
  coord[4][2] = 0.0;
  coord[5][0] = 0.0;
  coord[5][1] = 0.0;
  coord[5][2] = 1.0;
  coord[6][0] = 1.0;
  coord[6][1] = 0.0;
  coord[6][2] = 1.0;
  coord[7][0] = 0.0;
  coord[7][1] = 1.0;
  coord[7][2] = 1.0;
  coord[8][0] = 1.0;
  coord[8][1] = 1.0;
  coord[8][2] = 1.0;
  std::vector<double> stencil(8);
  factTester.TestComputeLinearInterpolationStencil(ggPFact, numDimension, coord, stencil);

  typename Teuchos::ScalarTraits<Scalar>::magnitudeType x = 0.0, y = 0.0, z = 0.0;
  for (LO i = 0; i < 8; ++i) {
    x += stencil[i] * coord[i + 1][0];
    y += stencil[i] * coord[i + 1][1];
    z += stencil[i] * coord[i + 1][2];
  }

  TEST_EQUALITY((std::abs(x - coord[0][0]) < 1e-5) && (std::abs(y - coord[0][1]) < 1e-5) && (std::abs(z - coord[0][2]) < 1e-5), true);

}  // LinearExtrapolation

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(GeneralGeometricPFactory, PoissonOnCubeLinear, Scalar,
                                  LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // used Xpetra lib (for maps and smoothers)
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // generate problem
  LocalOrdinal maxLevels = 2;
  // LocalOrdinal its=10;
  GO nx        = 4;
  GO ny        = 4;
  GO nz        = 4;
  GO numPoints = nx * ny * nz;
  Array<GO> gNodesPerDim(3);
  gNodesPerDim[0] = 4;
  gNodesPerDim[1] = 4;
  gNodesPerDim[2] = 4;
  Array<LO> lNodesPerDim(3);
  if (comm->getSize() == 1) {
    lNodesPerDim[0] = nx;
    lNodesPerDim[1] = ny;
    lNodesPerDim[2] = nz;
  } else if (comm->getSize() == 4) {
    lNodesPerDim[0] = nx;
    lNodesPerDim[1] = ny;
    lNodesPerDim[2] = nz / 4;
  }

  const RCP<const Map> map = MapFactory::Build(lib,
                                               gNodesPerDim[0] * gNodesPerDim[1] * gNodesPerDim[2],
                                               lNodesPerDim[0] * lNodesPerDim[1] * lNodesPerDim[2],
                                               0, comm);

  Teuchos::ParameterList problemParamList;
  problemParamList.set("nx", gNodesPerDim[0]);
  problemParamList.set("ny", gNodesPerDim[1]);
  problemParamList.set("nz", gNodesPerDim[2]);
  // problemParamList.set("keepBCs", true);

  // create Poisson problem and matrix
  Galeri::Xpetra::Laplace3DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> PoissonOnCube(problemParamList,
                                                                                              map);
  RCP<Matrix> Op = PoissonOnCube.BuildMatrix();

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
  nullSpace->putScalar((Scalar)1.0);
  Teuchos::Array<magnitude_type> norms(1);
  nullSpace->norm1(norms);
  if (comm->getRank() == 0) {
    out << "||NS|| = " << norms[0] << std::endl;
  }

  // build coordinates on rank 0
  global_size_t numGlobalElements                                                                                 = numPoints;
  size_t numLocalElements                                                                                         = (comm->getRank() == 0 ? numPoints : 0);
  RCP<const Map> exportMap                                                                                        = MapFactory::Build(lib, numGlobalElements, numLocalElements, 0, comm);
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > source_coordinates = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(exportMap, 3);
  if (comm->getRank() == 0) {
    for (LO k = 0; k < gNodesPerDim[2]; ++k) {
      for (LO j = 0; j < gNodesPerDim[1]; ++j) {
        for (LO i = 0; i < gNodesPerDim[0]; ++i) {
          source_coordinates->getDataNonConst(0)[k * ny * nx + j * nx + i] = i / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[0] - 1);
          source_coordinates->getDataNonConst(1)[k * ny * nx + j * nx + i] = j / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[1] - 1);
          source_coordinates->getDataNonConst(2)[k * ny * nx + j * nx + i] = k / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[2] - 1);
        }
      }
    }
  }
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(map, 3);
  RCP<Xpetra::Export<LO, GO, Node> > coords_exporter                                                       = Xpetra::ExportFactory<LO, GO, Node>::Build(exportMap, map);
  coordinates->doExport(*source_coordinates, *coords_exporter, Xpetra::INSERT);

  // fill hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  // create the factory manager and the factories
  FactoryManager M;

  RCP<Factory> Pfact     = rcp(new GeneralGeometricPFactory());
  RCP<Factory> Rfact     = rcp(new TransPFactory());
  RCP<Factory> Tfact     = rcp(new CoordinatesTransferFactory());
  RCP<RAPFactory> Acfact = rcp(new RAPFactory());
  RCP<Factory> NSfact    = rcp(new NullspaceFactory());

  // Set paramters needed by the factories
  Pfact->SetParameter("Coarsen", Teuchos::ParameterEntry(std::string("{2,2,2}")));
  Pfact->SetParameter("order", Teuchos::ParameterEntry(1));

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
  Finest->Set("A", Op);                       // set fine level matrix
  Finest->Set("Nullspace", nullSpace);        // set null space information for finest level
  Finest->Set("Coordinates", coordinates);    // set fine level coordinates
  Finest->Set("gNodesPerDim", gNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->Set("lNodesPerDim", lNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->SetFactoryManager(Teuchos::rcpFromRef(M));

  // Setup the hierarchy
  H->SetMaxCoarseSize(10);
  H->Setup(M, 0, maxLevels);

  // Extract the prolongator operator
  RCP<Level> lvl1                          = H->GetLevel(1);
  RCP<Xpetra::Matrix<SC, LO, GO, Node> > P = lvl1->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
  RCP<CrsMatrix> PCrs                      = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  // Construct vectors to check that a linear vector remains linear after projection
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector0 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(PCrs->getRangeMap(), 1);
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector1 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(PCrs->getDomainMap(), 1);
  {
    ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
    if (comm->getRank() == 0) {
      for (LO i = 0; i < 3; ++i) {
        for (LO j = 0; j < 3; ++j) {
          coarse_data[3 * i + j] = 2.0 * i + j;
        }
      }
    }
  }
  PCrs->apply(*vector1, *vector0, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(),
              Teuchos::ScalarTraits<SC>::zero());

  ArrayRCP<const SC> fine_data = vector0->getData(0);
  bool is_linear = true, is_injected = true;
  Array<LO> fine_inds(9);
  Array<LO> coarse_inds(9);
  ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
  if (comm->getRank() == 0) {
    fine_inds[0] = 0;
    fine_inds[1] = 2;
    fine_inds[2] = 3;
    fine_inds[3] = 8;
    fine_inds[4] = 10;
    fine_inds[5] = 11;
    fine_inds[6] = 12;
    fine_inds[7] = 14;
    fine_inds[8] = 15;
    for (LO i = 0; i < 9; ++i) {
      if (std::abs(fine_data[fine_inds[i]] - coarse_data[i]) > 1.0e-10) {
        is_injected = false;
      }
    }
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType linear_avg = std::abs(coarse_data[0] + coarse_data[1] + coarse_data[3] + coarse_data[4]) / 4.0;
    if (std::abs(fine_data[5] - linear_avg) > 1.0e-10) {
      is_linear = false;
    }
  }

  TEST_EQUALITY(Op != Teuchos::null, true);
  TEST_EQUALITY(H != Teuchos::null, true);
  TEST_EQUALITY(nullSpace != Teuchos::null, true);
  TEST_EQUALITY(is_injected, true);
  TEST_EQUALITY(is_linear, true);

}  // PoissonOnCubeLinear

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(GeneralGeometricPFactory, PoissonOnCubeConstant, Scalar,
                                  LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // used Xpetra lib (for maps and smoothers)
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // generate problem
  LocalOrdinal maxLevels = 2;
  // LocalOrdinal its=10;
  GO nx        = 4;
  GO ny        = 4;
  GO nz        = 4;
  GO numPoints = nx * ny * nz;
  Array<GO> gNodesPerDim(3);
  gNodesPerDim[0] = 4;
  gNodesPerDim[1] = 4;
  gNodesPerDim[2] = 4;
  Array<LO> lNodesPerDim(3);
  if (comm->getSize() == 1) {
    lNodesPerDim[0] = nx;
    lNodesPerDim[1] = ny;
    lNodesPerDim[2] = nz;
  } else if (comm->getSize() == 4) {
    lNodesPerDim[0] = nx;
    lNodesPerDim[1] = ny;
    lNodesPerDim[2] = nz / 4;
  }

  const RCP<const Map> map = MapFactory::Build(lib,
                                               gNodesPerDim[0] * gNodesPerDim[1] * gNodesPerDim[2],
                                               lNodesPerDim[0] * lNodesPerDim[1] * lNodesPerDim[2],
                                               0, comm);

  Teuchos::ParameterList problemParamList;
  problemParamList.set("nx", gNodesPerDim[0]);
  problemParamList.set("ny", gNodesPerDim[1]);
  problemParamList.set("nz", gNodesPerDim[2]);
  // problemParamList.set("keepBCs", true);

  // create Poisson problem and matrix
  Galeri::Xpetra::Laplace3DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> PoissonOnCube(problemParamList,
                                                                                              map);
  RCP<Matrix> Op = PoissonOnCube.BuildMatrix();

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
  nullSpace->putScalar((Scalar)1.0);
  Teuchos::Array<magnitude_type> norms(1);
  nullSpace->norm1(norms);
  if (comm->getRank() == 0) {
    out << "||NS|| = " << norms[0] << std::endl;
  }

  // build coordinates on rank 0
  global_size_t numGlobalElements                                                                                 = numPoints;
  size_t numLocalElements                                                                                         = (comm->getRank() == 0 ? numPoints : 0);
  RCP<const Map> exportMap                                                                                        = MapFactory::Build(lib, numGlobalElements, numLocalElements, 0, comm);
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > source_coordinates = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(exportMap, 3);
  if (comm->getRank() == 0) {
    for (LO k = 0; k < gNodesPerDim[2]; ++k) {
      for (LO j = 0; j < gNodesPerDim[1]; ++j) {
        for (LO i = 0; i < gNodesPerDim[0]; ++i) {
          source_coordinates->getDataNonConst(0)[k * ny * nx + j * nx + i] = i / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[0] - 1);
          source_coordinates->getDataNonConst(1)[k * ny * nx + j * nx + i] = j / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[1] - 1);
          source_coordinates->getDataNonConst(2)[k * ny * nx + j * nx + i] = k / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[2] - 1);
        }
      }
    }
  }
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(map, 3);
  RCP<Xpetra::Export<LO, GO, Node> > coords_exporter                                                       = Xpetra::ExportFactory<LO, GO, Node>::Build(exportMap, map);
  coordinates->doExport(*source_coordinates, *coords_exporter, Xpetra::INSERT);

  // fill hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  // create the factory manager and the factories
  FactoryManager M;

  RCP<Factory> Pfact     = rcp(new GeneralGeometricPFactory());
  RCP<Factory> Rfact     = rcp(new TransPFactory());
  RCP<Factory> Tfact     = rcp(new CoordinatesTransferFactory());
  RCP<RAPFactory> Acfact = rcp(new RAPFactory());
  RCP<Factory> NSfact    = rcp(new NullspaceFactory());

  // Set paramters needed by the factories
  Pfact->SetParameter("Coarsen", Teuchos::ParameterEntry(std::string("{2,2,2}")));
  Pfact->SetParameter("order", Teuchos::ParameterEntry(0));

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
  Finest->Set("A", Op);                       // set fine level matrix
  Finest->Set("Nullspace", nullSpace);        // set null space information for finest level
  Finest->Set("Coordinates", coordinates);    // set fine level coordinates
  Finest->Set("gNodesPerDim", gNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->Set("lNodesPerDim", lNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->SetFactoryManager(Teuchos::rcpFromRef(M));

  // Setup the hierarchy
  H->SetMaxCoarseSize(10);
  H->Setup(M, 0, maxLevels);

  // Extract the prolongator operator
  RCP<Level> lvl1                          = H->GetLevel(1);
  RCP<Xpetra::Matrix<SC, LO, GO, Node> > P = lvl1->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
  RCP<CrsMatrix> PCrs                      = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  // Construct vectors to check that a linear vector remains linear after projection
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector0 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(PCrs->getRangeMap(), 1);
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector1 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(PCrs->getDomainMap(), 1);
  {
    ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
    if (comm->getRank() == 0) {
      for (LO i = 0; i < 3; ++i) {
        for (LO j = 0; j < 3; ++j) {
          coarse_data[3 * i + j] = 2.0 * i + j;
        }
      }
    }
  }

  PCrs->apply(*vector1, *vector0, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(),
              Teuchos::ScalarTraits<SC>::zero());

  ArrayRCP<const SC> fine_data = vector0->getData(0);
  bool is_constant = true, is_injected = true;
  Array<LO> fine_inds(9);
  Array<LO> coarse_inds(9);
  ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
  if (comm->getRank() == 0) {
    fine_inds[0] = 0;
    fine_inds[1] = 2;
    fine_inds[2] = 3;
    fine_inds[3] = 8;
    fine_inds[4] = 10;
    fine_inds[5] = 11;
    fine_inds[6] = 12;
    fine_inds[7] = 14;
    fine_inds[8] = 15;
    for (LO i = 0; i < 9; ++i) {
      if (std::abs(fine_data[fine_inds[i]] - coarse_data[i]) > 1.0e-10) {
        is_injected = false;
      }
    }
    if ((std::abs(fine_data[1] - coarse_data[0]) > 1.0e-10) ||
        (std::abs(fine_data[4] - coarse_data[0]) > 1.0e-10) ||
        (std::abs(fine_data[5] - coarse_data[0]) > 1.0e-10)) {
      is_constant = false;
    }
  }

  TEST_EQUALITY(Op != Teuchos::null, true);
  TEST_EQUALITY(H != Teuchos::null, true);
  TEST_EQUALITY(nullSpace != Teuchos::null, true);
  TEST_EQUALITY(is_injected, true);
  TEST_EQUALITY(is_constant, true);

}  // PoissonOnCubeConstant

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(GeneralGeometricPFactory, LocalLexicographicLinear, Scalar,
                                  LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // used Xpetra lib (for maps and smoothers)
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  LO numDimensions = 3;

  RCP<Matrix> Op;
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates;
  RCP<Map> map;
  Array<GO> gNodesPerDim(3);
  Array<LO> lNodesPerDim(3);
  GGGetProblemData<SC, LO, GO, NO>(comm, lib, numDimensions, "Local Lexicographic", Op,
                                   coordinates, map, gNodesPerDim, lNodesPerDim);

  TEST_EQUALITY(Op != Teuchos::null, true);

  Array<Array<GO> > indsBounds(4);
  Array<GO> minGID(4);
  indsBounds[0].resize(4);
  if (comm->getSize() == 1) {
    indsBounds[0][0] = 0;
    indsBounds[0][1] = 8;
    indsBounds[0][2] = 0;
    indsBounds[0][3] = 8;
  } else {
    indsBounds[0][0] = 0;
    indsBounds[0][1] = 4;
    indsBounds[0][2] = 0;
    indsBounds[0][3] = 4;
  }
  minGID[0] = 0;
  indsBounds[1].resize(4);
  indsBounds[1][0] = 5;
  indsBounds[1][1] = 8;
  indsBounds[1][2] = 0;
  indsBounds[1][3] = 4;
  minGID[1]        = 225;
  indsBounds[2].resize(4);
  indsBounds[2][0] = 0;
  indsBounds[2][1] = 4;
  indsBounds[2][2] = 5;
  indsBounds[2][3] = 8;
  minGID[2]        = 405;
  indsBounds[3].resize(4);
  indsBounds[3][0] = 5;
  indsBounds[3][1] = 8;
  indsBounds[3][2] = 5;
  indsBounds[3][3] = 8;
  minGID[3]        = 585;

  Array<GO> meshData(comm->getSize() * 10);
  for (int rank = 0; rank < comm->getSize(); ++rank) {
    meshData[rank * 10 + 0] = rank;
    meshData[rank * 10 + 1] = rank;
    meshData[rank * 10 + 2] = 0;
    meshData[rank * 10 + 3] = 0;
    meshData[rank * 10 + 4] = gNodesPerDim[0] - 1;
    meshData[rank * 10 + 5] = indsBounds[rank][0];
    meshData[rank * 10 + 6] = indsBounds[rank][1];
    meshData[rank * 10 + 7] = indsBounds[rank][2];
    meshData[rank * 10 + 8] = indsBounds[rank][3];
    meshData[rank * 10 + 9] = minGID[rank];
  }

  // generate problem
  LocalOrdinal maxLevels = 3;

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
  RCP<Factory> Tfact     = rcp(new CoordinatesTransferFactory());
  RCP<RAPFactory> Acfact = rcp(new RAPFactory());
  RCP<Factory> NSfact    = rcp(new NullspaceFactory());

  // Set paramters needed by the factories
  Pfact->SetParameter("Coarsen", Teuchos::ParameterEntry(std::string("{2,2,2}")));
  Pfact->SetParameter("order", Teuchos::ParameterEntry(1));
  Pfact->SetParameter("meshLayout", Teuchos::ParameterEntry(std::string("Local Lexicographic")));

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
  Finest->Set("A", Op);                       // set fine level matrix
  Finest->Set("Nullspace", nullSpace);        // set null space information for finest level
  Finest->Set("Coordinates", coordinates);    // set fine level coordinates
  Finest->Set("gNodesPerDim", gNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->Set("lNodesPerDim", lNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->Set("meshData", meshData);          // set GeneralGeometricPFactory specific info
  Finest->SetFactoryManager(Teuchos::rcpFromRef(M));

  // Setup the hierarchy
  H->SetMaxCoarseSize(10);
  H->Setup(M, 0, maxLevels);

  // Extract the prolongator operator
  RCP<Level> lvl1                           = H->GetLevel(1);
  RCP<Xpetra::Matrix<SC, LO, GO, Node> > P1 = lvl1->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
  RCP<CrsMatrix> P1Crs                      = rcp_dynamic_cast<CrsMatrixWrap>(P1)->getCrsMatrix();

  // Construct vectors to check that a linear vector remains linear after projection
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector0 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(P1Crs->getRangeMap(), 1);
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector1 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(P1Crs->getDomainMap(), 1);
  std::vector<LO> coarse_inds(8);
  {
    ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
    if (comm->getSize() == 1) {
      coarse_inds[0] = 0;
      coarse_inds[1] = 1;
      coarse_inds[2] = 5;
      coarse_inds[3] = 6;
      coarse_inds[4] = 25;
      coarse_inds[5] = 26;
      coarse_inds[6] = 30;
      coarse_inds[7] = 31;
    } else if (comm->getSize() == 4 && comm->getRank() == 0) {
      coarse_inds[0] = 0;
      coarse_inds[1] = 1;
      coarse_inds[2] = 5;
      coarse_inds[3] = 6;
      coarse_inds[4] = 15;
      coarse_inds[5] = 16;
      coarse_inds[6] = 20;
      coarse_inds[7] = 21;
    }
    coarse_data[coarse_inds[0]] = 5;
    coarse_data[coarse_inds[1]] = 1;
    coarse_data[coarse_inds[2]] = 7;
    coarse_data[coarse_inds[3]] = 8;
    coarse_data[coarse_inds[4]] = 0;
    coarse_data[coarse_inds[5]] = 4;
    coarse_data[coarse_inds[6]] = 0;
    coarse_data[coarse_inds[7]] = 9;
  }

  P1Crs->apply(*vector1, *vector0, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(),
               Teuchos::ScalarTraits<SC>::zero());

  ArrayRCP<const SC> fine_data = vector0->getData(0);
  Array<LO> fine_inds(8);
  bool is_linear_lvl1 = true, is_injected_lvl1 = true;
  if (comm->getRank() == 0) {
    if (comm->getSize() == 1) {
      fine_inds[0] = 0;
      fine_inds[1] = 2;
      fine_inds[2] = 18;
      fine_inds[3] = 20;
      fine_inds[4] = 162;
      fine_inds[5] = 164;
      fine_inds[6] = 180;
      fine_inds[7] = 182;
    } else if (comm->getSize() == 4) {
      fine_inds[0] = 0;
      fine_inds[1] = 2;
      fine_inds[2] = 18;
      fine_inds[3] = 20;
      fine_inds[4] = 90;
      fine_inds[5] = 92;
      fine_inds[6] = 108;
      fine_inds[7] = 110;
    }
    ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
    for (LO i = 0; i < 8; ++i) {
      if (std::abs(fine_data[fine_inds[i]] - coarse_data[coarse_inds[i]]) > 1.0e-10) {
        is_injected_lvl1 = false;
      }
    }
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType linear_avg = 0.0;
    for (auto ind : coarse_inds) {
      linear_avg += std::abs(coarse_data[ind]) / 8.0;
    }
    LO fine_iref = 0;
    if (comm->getSize() == 1) {
      fine_iref = 91;
    } else if (comm->getSize() == 4) {
      fine_iref = 55;
    }
    if (std::abs(fine_data[fine_iref] - linear_avg) > 1.0e-10) {
      is_linear_lvl1 = false;
    }
  }

  // Extract the prolongator operator
  RCP<Level> lvl2                           = H->GetLevel(2);
  RCP<Xpetra::Matrix<SC, LO, GO, Node> > P2 = lvl2->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
  RCP<CrsMatrix> P2Crs                      = rcp_dynamic_cast<CrsMatrixWrap>(P2)->getCrsMatrix();

  // Construct vectors to check that a linear vector remains linear after projection
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector2 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(P2Crs->getRangeMap(), 1);
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector3 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(P2Crs->getDomainMap(), 1);
  if (comm->getSize() == 1) {
    coarse_inds[0] = 0;
    coarse_inds[1] = 1;
    coarse_inds[2] = 3;
    coarse_inds[3] = 4;
    coarse_inds[4] = 9;
    coarse_inds[5] = 10;
    coarse_inds[6] = 12;
    coarse_inds[7] = 13;
  } else if (comm->getSize() == 4 && comm->getRank() == 0) {
    coarse_inds[0] = 0;
    coarse_inds[1] = 1;
    coarse_inds[2] = 3;
    coarse_inds[3] = 4;
    coarse_inds[4] = 6;
    coarse_inds[5] = 7;
    coarse_inds[6] = 9;
    coarse_inds[7] = 10;
  }
  {
    ArrayRCP<SC> coarse_data    = vector3->getDataNonConst(0);
    coarse_data[coarse_inds[0]] = 5;
    coarse_data[coarse_inds[1]] = 1;
    coarse_data[coarse_inds[2]] = 7;
    coarse_data[coarse_inds[3]] = 8;
    coarse_data[coarse_inds[4]] = 0;
    coarse_data[coarse_inds[5]] = 4;
    coarse_data[coarse_inds[6]] = 0;
    coarse_data[coarse_inds[7]] = 9;
  }

  P2Crs->apply(*vector3, *vector2, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(),
               Teuchos::ScalarTraits<SC>::zero());

  fine_data           = vector2->getData(0);
  bool is_linear_lvl2 = true, is_injected_lvl2 = true;
  if (comm->getRank() == 0) {
    if (comm->getSize() == 1) {
      fine_inds[0] = 0;
      fine_inds[1] = 2;
      fine_inds[2] = 10;
      fine_inds[3] = 12;
      fine_inds[4] = 50;
      fine_inds[5] = 52;
      fine_inds[6] = 60;
      fine_inds[7] = 62;
    } else if (comm->getSize() == 4) {
      fine_inds[0] = 0;
      fine_inds[1] = 2;
      fine_inds[2] = 10;
      fine_inds[3] = 12;
      fine_inds[4] = 30;
      fine_inds[5] = 32;
      fine_inds[6] = 40;
      fine_inds[7] = 42;
    }
    ArrayRCP<SC> coarse_data = vector3->getDataNonConst(0);
    for (LO i = 0; i < 8; ++i) {
      if (std::abs(fine_data[fine_inds[i]] - coarse_data[coarse_inds[i]]) > 1.0e-10) {
        is_injected_lvl2 = false;
      }
    }
    typename Teuchos::ScalarTraits<Scalar>::magnitudeType linear_avg = 0.0;
    for (auto ind : coarse_inds) {
      linear_avg += std::abs(coarse_data[ind]) / 8.0;
    }
    LO fine_iref = 0;
    if (comm->getSize() == 1) {
      fine_iref = 31;
    } else if (comm->getSize() == 4) {
      fine_iref = 21;
    }
    if (std::abs(fine_data[fine_iref] - linear_avg) > 1.0e-10) {
      is_linear_lvl2 = false;
    }
  }

  TEST_EQUALITY(Op != Teuchos::null, true);
  TEST_EQUALITY(H != Teuchos::null, true);
  TEST_EQUALITY(nullSpace != Teuchos::null, true);
  TEST_EQUALITY(is_injected_lvl1, true);
  TEST_EQUALITY(is_linear_lvl1, true);
  TEST_EQUALITY(is_injected_lvl2, true);
  TEST_EQUALITY(is_linear_lvl2, true);

}  // End LocalLexicographicLinear

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(GeneralGeometricPFactory, LocalLexicographicConstant, Scalar,
                                  LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // used Xpetra lib (for maps and smoothers)
  Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();

  // generate problem
  LocalOrdinal maxLevels = 2;
  // LocalOrdinal its=10;
  GO nx        = 4;
  GO ny        = 4;
  GO nz        = 4;
  GO numPoints = nx * ny * nz;
  Array<GO> gNodesPerDim(3);
  gNodesPerDim[0] = 4;
  gNodesPerDim[1] = 4;
  gNodesPerDim[2] = 4;
  Array<GO> meshData(comm->getSize() * 10);
  Array<LO> lNodesPerDim(3);
  if (comm->getSize() == 1) {
    lNodesPerDim[0] = nx;
    lNodesPerDim[1] = ny;
    lNodesPerDim[2] = nz;

    meshData[0] = comm->getRank();  // my Parent rank
    meshData[1] = comm->getRank();  // my rank
    meshData[2] = 0;                // block ID
    meshData[3] = 0;
    meshData[4] = nx - 1;
    meshData[5] = 0;
    meshData[6] = ny - 1;
    meshData[7] = 0;
    meshData[8] = nz - 1;
    meshData[9] = 0;
  } else if (comm->getSize() == 4) {
    lNodesPerDim[0] = nx;
    lNodesPerDim[1] = ny;
    lNodesPerDim[2] = nz / 4;

    for (int rank = 0; rank < comm->getSize(); ++rank) {
      meshData[rank * 10 + 0] = rank;            // my Parent rank
      meshData[rank * 10 + 1] = rank;            // my rank
      meshData[rank * 10 + 2] = 0;               // block ID
      meshData[rank * 10 + 3] = 0;               // imin
      meshData[rank * 10 + 4] = nx - 1;          // imax
      meshData[rank * 10 + 5] = 0;               // jmin
      meshData[rank * 10 + 6] = ny - 1;          // jmax
      meshData[rank * 10 + 7] = rank;            // kmin
      meshData[rank * 10 + 8] = rank;            // kmax
      meshData[rank * 10 + 9] = rank * ny * nx;  // minGID on rank
    }
  }

  const RCP<const Map> map = MapFactory::Build(lib,
                                               gNodesPerDim[0] * gNodesPerDim[1] * gNodesPerDim[2],
                                               lNodesPerDim[0] * lNodesPerDim[1] * lNodesPerDim[2],
                                               0,
                                               comm);

  Teuchos::ParameterList problemParamList;
  problemParamList.set("nx", gNodesPerDim[0]);
  problemParamList.set("ny", gNodesPerDim[1]);
  problemParamList.set("nz", gNodesPerDim[2]);
  // problemParamList.set("keepBCs", true);

  // create Poisson problem and matrix
  Galeri::Xpetra::Laplace3DProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector> PoissonOnCube(problemParamList,
                                                                                              map);
  RCP<Matrix> Op = PoissonOnCube.BuildMatrix();

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
  nullSpace->putScalar((Scalar)1.0);
  Teuchos::Array<magnitude_type> norms(1);
  nullSpace->norm1(norms);
  if (comm->getRank() == 0) {
    out << "||NS|| = " << norms[0] << std::endl;
  }

  // build coordinates on rank 0
  global_size_t numGlobalElements                                                                                 = numPoints;
  size_t numLocalElements                                                                                         = (comm->getRank() == 0 ? numPoints : 0);
  RCP<const Map> exportMap                                                                                        = MapFactory::Build(lib, numGlobalElements, numLocalElements, 0, comm);
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > source_coordinates = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(exportMap, 3);
  if (comm->getRank() == 0) {
    for (LO k = 0; k < gNodesPerDim[2]; ++k) {
      for (LO j = 0; j < gNodesPerDim[1]; ++j) {
        for (LO i = 0; i < gNodesPerDim[0]; ++i) {
          source_coordinates->getDataNonConst(0)[k * ny * nx + j * nx + i] = i / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[0] - 1);
          source_coordinates->getDataNonConst(1)[k * ny * nx + j * nx + i] = j / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[1] - 1);
          source_coordinates->getDataNonConst(2)[k * ny * nx + j * nx + i] = k / Teuchos::as<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(gNodesPerDim[2] - 1);
        }
      }
    }
  }
  RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> > coordinates = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO>::Build(map, 3);
  RCP<Xpetra::Export<LO, GO, Node> > coords_exporter                                                       = Xpetra::ExportFactory<LO, GO, Node>::Build(exportMap, map);
  coordinates->doExport(*source_coordinates, *coords_exporter, Xpetra::INSERT);

  // fill hierarchy
  RCP<Hierarchy> H = rcp(new Hierarchy());
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  // create the factory manager and the factories
  FactoryManager M;

  RCP<Factory> Pfact     = rcp(new GeneralGeometricPFactory());
  RCP<Factory> Rfact     = rcp(new TransPFactory());
  RCP<Factory> Tfact     = rcp(new CoordinatesTransferFactory());
  RCP<RAPFactory> Acfact = rcp(new RAPFactory());
  RCP<Factory> NSfact    = rcp(new NullspaceFactory());

  // Set paramters needed by the factories
  Pfact->SetParameter("Coarsen", Teuchos::ParameterEntry(std::string("{2,2,2}")));
  Pfact->SetParameter("order", Teuchos::ParameterEntry(0));
  Pfact->SetParameter("meshLayout", Teuchos::ParameterEntry(std::string("Local Lexicographic")));

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
  Finest->Set("A", Op);                       // set fine level matrix
  Finest->Set("Nullspace", nullSpace);        // set null space information for finest level
  Finest->Set("Coordinates", coordinates);    // set fine level coordinates
  Finest->Set("gNodesPerDim", gNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->Set("lNodesPerDim", lNodesPerDim);  // set GeneralGeometricPFactory specific info
  Finest->Set("meshData", meshData);          // set GeneralGeometricPFactory specific info
  Finest->SetFactoryManager(Teuchos::rcpFromRef(M));

  // Setup the hierarchy
  H->SetMaxCoarseSize(10);
  H->Setup(M, 0, maxLevels);

  // Extract the prolongator operator
  RCP<Level> lvl1                          = H->GetLevel(1);
  RCP<Xpetra::Matrix<SC, LO, GO, Node> > P = lvl1->Get<RCP<Matrix> >("P", MueLu::NoFactory::get());
  RCP<CrsMatrix> PCrs                      = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();

  // Construct vectors to check that a linear vector remains linear after projection
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector0 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(PCrs->getRangeMap(), 1);
  RCP<Xpetra::MultiVector<SC, LO, GO, NO> > vector1 = Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(PCrs->getDomainMap(), 1);
  {
    ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
    if (comm->getRank() == 0) {
      for (LO i = 0; i < 3; ++i) {
        for (LO j = 0; j < 3; ++j) {
          coarse_data[3 * i + j] = 2.0 * i + j;
        }
      }
    }
  }

  PCrs->apply(*vector1, *vector0, Teuchos::NO_TRANS, Teuchos::ScalarTraits<SC>::one(),
              Teuchos::ScalarTraits<SC>::zero());

  ArrayRCP<const SC> fine_data = vector0->getData(0);
  bool is_constant = true, is_injected = true;
  Array<LO> fine_inds(9);
  Array<LO> coarse_inds(9);
  ArrayRCP<SC> coarse_data = vector1->getDataNonConst(0);
  if (comm->getRank() == 0) {
    fine_inds[0] = 0;
    fine_inds[1] = 2;
    fine_inds[2] = 3;
    fine_inds[3] = 8;
    fine_inds[4] = 10;
    fine_inds[5] = 11;
    fine_inds[6] = 12;
    fine_inds[7] = 14;
    fine_inds[8] = 15;
    for (LO i = 0; i < 9; ++i) {
      if (std::abs(fine_data[fine_inds[i]] - coarse_data[i]) > 1.0e-10) {
        is_injected = false;
      }
    }
    if ((std::abs(fine_data[1] - coarse_data[0]) > 1.0e-10) ||
        (std::abs(fine_data[4] - coarse_data[0]) > 1.0e-10) ||
        (std::abs(fine_data[5] - coarse_data[0]) > 1.0e-10)) {
      is_constant = false;
    }
  }

  TEST_EQUALITY(Op != Teuchos::null, true);
  TEST_EQUALITY(H != Teuchos::null, true);
  TEST_EQUALITY(nullSpace != Teuchos::null, true);
  TEST_EQUALITY(is_injected, true);
  TEST_EQUALITY(is_constant, true);

}  // LocalLexicographicConstant

#define MUELU_ETI_GROUP(Scalar, LO, GO, Node)                                                                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(GeneralGeometricPFactory, Constructor, Scalar, LO, GO, Node)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(GeneralGeometricPFactory, LinearInterpolation, Scalar, LO, GO, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(GeneralGeometricPFactory, LinearExtrapolation, Scalar, LO, GO, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(GeneralGeometricPFactory, PoissonOnCubeLinear, Scalar, LO, GO, Node)      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(GeneralGeometricPFactory, PoissonOnCubeConstant, Scalar, LO, GO, Node)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(GeneralGeometricPFactory, LocalLexicographicLinear, Scalar, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(GeneralGeometricPFactory, LocalLexicographicConstant, Scalar, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
