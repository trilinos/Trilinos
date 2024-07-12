// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_LineDetectionFactory.hpp"
#include "MueLu_SemiCoarsenPFactory.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(SemiCoarsenPFactory, TestSemiCoarsenP, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);

  out << "version: " << MueLu::Version() << std::endl;

  typedef TestHelpers::TestFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> test_factory;
  Level fineLevel, coarseLevel;
  test_factory::createTwoLevelHierarchy(fineLevel, coarseLevel);
  fineLevel.SetFactoryManager(Teuchos::null);  // factory manager is not used on this test
  coarseLevel.SetFactoryManager(Teuchos::null);

  // Set global geometric data
  const std::string meshLayout = "Global Lexicographic";
  const std::string coupling   = "uncoupled";
  const LO numDimensions       = 3;
  const LO numSize             = 5;
  Array<GO> meshData;
  Array<LO> lNodesPerDir(numDimensions);
  Array<GO> gNodesPerDir(numDimensions);
  for (int dim = 0; dim < numDimensions; ++dim)
    gNodesPerDir[dim] = numSize;

  typedef typename Teuchos::ScalarTraits<Scalar>::coordinateType coord_type;
  typedef Xpetra::MultiVector<coord_type, LO, GO, NO> CoordMV;
  RCP<CoordMV> fineCoords =
      TestHelpers::TestFactory<SC, LO, GO, NO>::BuildGeoCoordinates(numDimensions, gNodesPerDir,
                                                                    lNodesPerDir, meshData,
                                                                    meshLayout);

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", gNodesPerDir[0]);
  matrixList.set("ny", gNodesPerDir[1]);
  matrixList.set("nz", gNodesPerDir[2]);
  matrixList.set("matrixType", "Laplace3D");
  matrixList.set("left boundary", "Neumann");
  matrixList.set("right boundary", "Neumann");
  matrixList.set("front boundary", "Neumann");
  matrixList.set("back boundary", "Neumann");
  matrixList.set("bottom boundary", "Dirichlet");
  matrixList.set("top boundary", "Dirichlet");
  RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr = Galeri::Xpetra::
      BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>("Laplace3D", fineCoords->getMap(),
                                                                matrixList);
  RCP<Matrix> A = Pr->BuildMatrix();

  LocalOrdinal NSdim         = 1;
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(), NSdim);
  nullSpace->putScalar(1.0);

  fineLevel.Set("A", A);
  fineLevel.Set("Nullspace", nullSpace);
  fineLevel.Set("Coordinates", fineCoords);
  fineLevel.Set("CoarseNumZLayers", numSize);

  RCP<LineDetectionFactory> LineDetectionFact = rcp(new LineDetectionFactory());
  LineDetectionFact->SetParameter("linedetection: orientation",
                                  Teuchos::ParameterEntry(std::string("coordinates")));
  LineDetectionFact->SetParameter("linedetection: num layers",
                                  Teuchos::ParameterEntry(numSize));

  RCP<SemiCoarsenPFactory> SemiCoarsenPFact = rcp(new SemiCoarsenPFactory());
  SemiCoarsenPFact->SetParameter("semicoarsen: coarsen rate", Teuchos::ParameterEntry(2));
  SemiCoarsenPFact->SetFactory("LineDetection_VertLineIds", LineDetectionFact);
  SemiCoarsenPFact->SetFactory("LineDetection_Layers", LineDetectionFact);
  SemiCoarsenPFact->SetFactory("CoarseNumZLayers", LineDetectionFact);
  coarseLevel.Request("P", SemiCoarsenPFact.get());
  coarseLevel.Request("Coordinates", SemiCoarsenPFact.get());
  SemiCoarsenPFact->Build(fineLevel, coarseLevel);

  RCP<Matrix> P;
  coarseLevel.Get("P", P, SemiCoarsenPFact.get());
  RCP<CoordMV> coarseCoords = coarseLevel.Get<RCP<CoordMV>>("Coordinates", SemiCoarsenPFact.get());

  coarseLevel.Release("P", SemiCoarsenPFact.get());
  coarseLevel.Release("Coordinates", SemiCoarsenPFact.get());

  // Prolongate coarse coords and compute difference
  using STS                       = Teuchos::ScalarTraits<SC>;
  const auto one                  = STS::one();
  RCP<MultiVector> coarseCoordsSC = Utilities::RealValuedToScalarMultiVector(coarseCoords);
  RCP<MultiVector> fineCoordsDiff = Utilities::RealValuedToScalarMultiVector(fineCoords);
  P->apply(*coarseCoordsSC, *fineCoordsDiff, Teuchos::NO_TRANS, one, -one);

  // check prolongation of coarse coordinates
  // in this special case, the third layer will have the correct coordinates
  const int numNodes    = fineCoordsDiff->getLocalLength();
  const int numVectors  = fineCoordsDiff->getNumVectors();
  const auto fineMap    = fineCoordsDiff->getMap();
  const GO gStartLayer3 = 2 * numSize * numSize;
  const GO gEndLayer3   = 3 * numSize * numSize;
  LO numBadCoords       = 0;
  for (LO node = 0; node < numNodes; ++node) {
    const GO gnode = fineMap->getGlobalElement(node);
    if (gnode >= gStartLayer3 && gnode < gEndLayer3)
      for (int k = 0; k < numVectors; ++k) {
        const auto fineCoordsDiffArray = fineCoordsDiff->getData(k);
        if (STS::magnitude(fineCoordsDiffArray[node]) > 100 * STS::eps())
          numBadCoords += 1;
      }
  }
  const auto comm = Teuchos::DefaultComm<int>::getComm();
  LO gNumBadCoords;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &numBadCoords, &gNumBadCoords);
  TEST_EQUALITY(gNumBadCoords, 0);
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(SemiCoarsenPFactory, TestSemiCoarsenP, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
