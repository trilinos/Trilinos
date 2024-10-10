// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <ostream>

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include <MueLu_InitialBlockNumberFactory.hpp>
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_CoalesceDropFactory_kokkos.hpp"
#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_LWGraph_kokkos.hpp"

#include <Galeri_XpetraParameters.hpp>

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<CoalesceDropFactory_kokkos> coalesceDropFact = rcp(new CoalesceDropFactory_kokkos());
  TEST_EQUALITY(coalesceDropFact != Teuchos::null, true);

  out << *coalesceDropFact << std::endl;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, Build, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);
}  // Build

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 40);

}  // DistanceLaplacian

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (comm->getRank() == 0) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // DistanceLaplacianScaledCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianUnscaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (!comm->getRank()) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // DistanceLaplacianUnscaleCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianCutSym, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);

  // Now we doctor the coordinates so that the off-diagonal pair row 0 will want to keep (0,1) and row 1 will want to drop (1,0)
  if (!comm->getRank()) {
    auto vals = coordinates->getDataNonConst(0);
    vals[0]   = vals[0] - 2000 * 36;
  }

  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.5));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: distance laplacian algo", Teuchos::ParameterEntry(std::string("scaled cut symmetric")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 106);

}  // DistanceLaplacianCutScaled

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicalScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  // Change entry (1,0)
  auto crsA = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true)->getCrsMatrix();
  crsA->resumeFill();
  if (comm->getRank() == 0) {
    Teuchos::Array<GlobalOrdinal> cols(3);
    Teuchos::Array<Scalar> vals(3);
    size_t numEntries;
    crsA->getGlobalRowCopy(1, cols, vals, numEntries);
    vals[0] = 0.5;
    crsA->replaceGlobalValues(1, cols, vals);
  }
  crsA->fillComplete();
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // ClassicalScaledCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicalUnScaledCut, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  // Change entry (1,0)
  auto crsA = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true)->getCrsMatrix();
  crsA->resumeFill();
  if (comm->getRank() == 0) {
    Teuchos::Array<GlobalOrdinal> cols(3);
    Teuchos::Array<Scalar> vals(3);
    size_t numEntries;
    crsA->getGlobalRowCopy(1, cols, vals, numEntries);
    vals[0] = 0.5;
    crsA->replaceGlobalValues(1, cols, vals);
  }
  crsA->fillComplete();
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // ClassicalUnScaledCut

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicalCutSym, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  // Change entry (1,0)
  auto crsA = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true)->getCrsMatrix();
  crsA->resumeFill();
  if (comm->getRank() == 0) {
    Teuchos::Array<GlobalOrdinal> cols(3);
    Teuchos::Array<Scalar> vals(3);
    size_t numEntries;
    crsA->getGlobalRowCopy(1, cols, vals, numEntries);
    vals[0] = 0.5;
    crsA->replaceGlobalValues(1, cols, vals);
  }
  crsA->fillComplete();
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("scaled cut symmetric")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 106);

}  // ClassicalCutSym

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, SignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  // Change entry (1,0)
  auto crsA = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true)->getCrsMatrix();
  crsA->resumeFill();
  if (comm->getRank() == 0) {
    Teuchos::Array<GlobalOrdinal> cols(3);
    Teuchos::Array<Scalar> vals(3);
    size_t numEntries;
    crsA->getGlobalRowCopy(1, cols, vals, numEntries);
    vals[0] *= 2;
    crsA->replaceGlobalValues(1, cols, vals);
  }
  crsA->fillComplete();
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // A_10 = -2
  // A_ij = -1
  // A_ii = 2
  // criterion for dropping is
  // -Re(L_ij) <= tol * max_{k\neq i} Re(-L_ik)
  // -> We drop entry (1,2).
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 105);

}  // SignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, SignedScaledCutClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("scaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  TEST_THROW(coalesceDropFact.Build(fineLevel), MueLu::Exceptions::RuntimeError);

  //    RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &coalesceDropFact);
  //    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  //    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  //    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  //    const RCP<const Map> myDomainMap = graph->GetDomainMap();

  //    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  //    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  //    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
  //    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(36 + (comm->getSize()-1)*2));

  //    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  //    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  //    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
  //    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),36);

  //    TEST_EQUALITY(graph->GetGlobalNumEdges(),36);

}  // SignedScaledCutClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, SignedUnscaledCutClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical")));
  coalesceDropFact.SetParameter("aggregation: classical algo", Teuchos::ParameterEntry(std::string("unscaled cut")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  TEST_THROW(coalesceDropFact.Build(fineLevel), MueLu::Exceptions::RuntimeError);

  //    RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &coalesceDropFact);
  //    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  //    TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  //    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  //    const RCP<const Map> myDomainMap = graph->GetDomainMap();

  //    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  //    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  //    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
  //    TEST_EQUALITY(myImportMap->getGlobalNumElements(),Teuchos::as<size_t>(36 + (comm->getSize()-1)*2));

  //    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  //    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  //    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
  //    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),36);

  //    TEST_EQUALITY(graph->GetGlobalNumEdges(),36);

}  // SignedUnScaledCutClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalColoredSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal colored signed classical")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalColoredSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalNoColoredSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  // this test is only compatible with rank higher than 1
  if (comm->getSize() == 1) {
    return;
  }

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561 = 3^8)
  // Nice size for 1D and perfect aggregation on small numbers of processors. (8748 = 4*3^7)
  Teuchos::CommandLineProcessor clp(false);
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);

  RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  //    getCrsGraph()->getImporter()
  RCP<const Import> importer = ImportFactory::Build(A->getRowMap(), map);
  fineLevel.Set("Importer", importer);
  auto importerTest = A->getCrsGraph()->getImporter();  // NULL
  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal colored signed classical")));
  coalesceDropFact.SetParameter("aggregation: coloring: localize color graph", Teuchos::ParameterEntry(false));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  // Need an importer
  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalNoColoredSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalSignedClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  Teuchos::ParameterList galeriList;
  galeriList.set("nx", Teuchos::as<GlobalOrdinal>(36));
  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), galeriList);
  fineLevel.Set("Coordinates", coordinates);

  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  // We're dropping all the interior off-diagonal entries.
  // dx = 1/36
  // L_ij = -36
  // L_ii = 72
  // criterion for dropping is |L_ij|^2 <= tol^2 * |L_ii*L_jj|
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 0.51));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal signed classical")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);

  const RCP<const Map> myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  const RCP<const Map> myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), Teuchos::as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);

  TEST_EQUALITY(graph->GetGlobalNumEdges(), 36);

}  // BlockDiagonalSignedClassical

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonal, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);

  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalClassical, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0 / 8.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal classical")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalDistanceLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalDistanceDifferentCoordinatesLaplacian, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  GO bnx = 15 * comm->getSize();
  Teuchos::ParameterList bMatrixList;
  matrixList.set("bnx", bnx);
  RCP<Matrix> B = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(bMatrixList, lib);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("1D", B->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, BlockDiagonalDistanceLaplacianWeighted, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildBlockMatrixAsPoint(matrixList, lib);
  A->SetFixedBlockSize(1);  // So we can block diagonalize
  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<InitialBlockNumberFactory> ibFact = rcp(new InitialBlockNumberFactory());
  Teuchos::ParameterList ibList;
  ibList.set("aggregation: block diagonal: interleaved blocksize", 3);
  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetFactory("BlockNumber", ibFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("block diagonal distance laplacian")));
  coalesceDropFact.SetParameter("aggregation: block diagonal: interleaved blocksize", Teuchos::ParameterEntry(3));
  std::vector<double> weights_v{100.0, 1.0, 1.0, 1.0, 100, 1.0, 1.0, 1.0, 100.0};
  Teuchos::Array<double> weights(weights_v);
  coalesceDropFact.SetParameter("aggregation: distance laplacian directional weights", Teuchos::ParameterEntry(weights));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, DistanceLaplacianWeighted, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixList, lib);

  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", A->getRowMap(), matrixList);
  fineLevel.Set("Coordinates", coordinates);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetDefaultVerbLevel(MueLu::Extreme);
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.025));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));
  std::vector<double> weights_v{100.0, 1.0, 1.0};
  Teuchos::Array<double> weights(weights_v);
  coalesceDropFact.SetParameter("aggregation: distance laplacian directional weights", Teuchos::ParameterEntry(weights));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, SignedClassicalSA, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType real_type;

  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  GO nx = 10 * comm->getSize();
  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("ny", (GO)10);
  matrixList.set("nz", (GO)10);
  matrixList.set("matrixType", "Laplace3D");
  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixList, lib);

  Level fineLevel;
  fineLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  CoalesceDropFactory_kokkos coalesceDropFact;
  coalesceDropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  coalesceDropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.0));
  coalesceDropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("signed classical sa")));
  fineLevel.Request("Graph", &coalesceDropFact);
  fineLevel.Request("DofsPerNode", &coalesceDropFact);

  coalesceDropFact.Build(fineLevel);

  RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &coalesceDropFact);
  LO myDofsPerNode          = fineLevel.Get<LO>("DofsPerNode", &coalesceDropFact);
  TEST_EQUALITY(Teuchos::as<int>(myDofsPerNode) == 1, true);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicScalarWithoutFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(36);
  fineLevel.Set("A", A);

  RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact;
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = *graph;
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1) {
          auto v0 = lclLWGraph.getNeighborVertices(0);
          auto v1 = lclLWGraph.getNeighborVertices(1);
          auto v2 = lclLWGraph.getNeighborVertices(2);
          if (v0.length == 2 && ((v0(0) == 0 && v0(1) == 1) || (v0(0) == 1 && v0(1) == 0)) &&
              v1.length == 3 && v2.length == 3)
            correct = true;
        } else {
          if (comm_rank == 0) {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;

          } else {
            if (lclLWGraph.getNeighborVertices(0).length == 3)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(36 + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 35);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), 36);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicScalarWithFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  auto dofMap = MapFactory::Build(lib, 3 * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 1.0, -1.0, -0.0001);

  mtx->SetFixedBlockSize(1, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(0.5));

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == 3 * comm->getSize(), true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = *graph;
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1) {
          auto v0 = lclLWGraph.getNeighborVertices(0);
          auto v1 = lclLWGraph.getNeighborVertices(1);
          auto v2 = lclLWGraph.getNeighborVertices(2);
          if (v0.length == 1 && v0(0) == 0 &&
              v1.length == 2 && ((v1(0) == 0 && v1(1) == 1) || (v1(0) == 1 && v1(1) == 0)) &&
              v2.length == 2 && ((v2(0) == 1 && v2(1) == 2) || (v2(0) == 2 && v2(1) == 1)))
            correct = true;
        } else {
          if (comm_rank == 0) {
            if (lclLWGraph.getNeighborVertices(0).length == 1)
              correct = true;

          } else {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), 3 * comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(3 * comm->getSize() + (comm->getSize() - 1) * 2));

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), 3 * comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(3 * comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 3);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicBlockWithoutFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  int blockSize = 3;

  auto dofMap = MapFactory::Build(lib, blockSize * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);
  mtx->SetFixedBlockSize(blockSize, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);

  bool bCorrectGraph = false;
  int reduction_val  = 0;
  int comm_size = comm->getSize(), comm_rank = comm->getRank();
  auto lclLWGraph = *graph;
  Kokkos::parallel_reduce(
      "MueLu:TentativePF:Build:compute_agg_sizes", Kokkos::RangePolicy<typename NO::execution_space, size_t>(0, 1),
      KOKKOS_LAMBDA(const LO i, int &correct) {
        if (comm_size == 1 && lclLWGraph.getNeighborVertices(0).length == 1) {
          correct = true;
        } else {
          if (comm_rank == 0 || comm_rank == comm_size - 1) {
            if (lclLWGraph.getNeighborVertices(0).length == 2)
              correct = true;

          } else {
            if (static_cast<int>(lclLWGraph.getNeighborVertices(0).length) == blockSize)
              correct = true;
          }
        }
      },
      reduction_val);
  bCorrectGraph = reduction_val;
  TEST_EQUALITY(bCorrectGraph, true);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));

  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 1);
    } else {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 2);
    }
  }
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t maxLocalIndex      = myImportMap->getMaxLocalIndex();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(maxLocalIndex, numLocalRowMapElts * blockSize - 2);
    } else {
      TEST_EQUALITY(maxLocalIndex, numLocalRowMapElts * blockSize - 1);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 1);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, ClassicBlockWithFiltering, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  auto dofMap = MapFactory::Build(lib, 3 * comm->getSize(), 0, comm);
  auto mtx    = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, 0.00001);

  mtx->SetFixedBlockSize(3, 0);
  fineLevel.Set("A", mtx);

  RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory);
  CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
  dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(1.0));

  fineLevel.Request("Graph", &dropFact);
  fineLevel.Request("DofsPerNode", &dropFact);

  dropFact.Build(fineLevel);

  auto graph         = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);
  auto myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);

  TEST_EQUALITY(as<int>(myDofsPerNode) == 3, true);
  TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);

  TEST_EQUALITY(graph->getNeighborVertices(0).size(), 1);

  auto myImportMap = graph->GetImportMap();  // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
  auto myDomainMap = graph->GetDomainMap();

  TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myImportMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myImportMap->getGlobalNumElements(), as<size_t>(comm->getSize() + 2 * (comm->getSize() - 1)));
  if (comm->getSize() > 1) {
    size_t numLocalRowMapElts = graph->GetNodeNumVertices();
    size_t numLocalImportElts = myImportMap->getLocalNumElements();
    if (comm->getRank() == 0 || comm->getRank() == comm->getSize() - 1) {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 1);
    } else {
      TEST_EQUALITY(numLocalImportElts, numLocalRowMapElts + 2);
    }
  }

  TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize() - 1);
  TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMaxLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getMinLocalIndex(), 0);
  TEST_EQUALITY(myDomainMap->getGlobalNumElements(), as<size_t>(comm->getSize()));
  TEST_EQUALITY(myDomainMap->getLocalNumElements(), 1);
}

#if 0
  TEUCHOS_UNIT_TEST(CoalesceDropFactory_kokkos, LaplacianScalarWithoutFiltering)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::Build1DPoisson(36);
    fineLevel.Set("A", A);

    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory_kokkos dropFact;
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop scheme", Teuchos::ParameterEntry(std::string("distance laplacian")));

    fineLevel.Request("Graph",       &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    auto graph         = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph",       &dropFact);
    auto myDofsPerNode = fineLevel.Get<LO>                  ("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(myDofsPerNode) == 1, true);

    bool bCorrectGraph = false;
    if (comm->getSize() == 1) {
      auto v0 = graph->getNeighborVertices(0);
      auto v1 = graph->getNeighborVertices(1);
      auto v2 = graph->getNeighborVertices(2);
      if (v0.size() == 2 && ((v0(0) == 0 && v0(1) == 1) || (v0(0) == 1 && v0(1) == 0)) &&
          v1.size() == 3 && v2.size() == 3)
        bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 ) {
        if (graph->getNeighborVertices(0).size() == 2)
          bCorrectGraph = true;

      } else {
        if (graph->getNeighborVertices(0).size() == 3)
          bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    auto myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    auto myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(),  35);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(),  0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),      0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),  as<size_t>(36 + (comm->getSize()-1)*2));

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(),  35);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(),  0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),      0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),  36);
  }
#endif

#if 0
  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedLW)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    // unit test for block size 3 using a strided map

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    int blockSize=3;

    GO nx = blockSize*comm->getSize();
    RCP<Matrix> A = TestHelpers::TestFactory<SC,LO,GO,NO>::Build1DPoisson(nx);

    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(blockSize));
    LocalOrdinal stridedBlockId = -1;

    RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  A->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  0 /*offset*/
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            A->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            0 /*offset*/
                                            );

    if(A->IsView("stridedMaps") == true) A->RemoveView("stridedMaps");
    A->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", A);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);
    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == blockSize) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t maxLocalIndex = myImportMap->getMaxLocalIndex();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-2), true);
      } else {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-1), true);
      }
    }

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStridedLW

  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStrided2LW)
  {
#include "MueLu_UseShortNames.hpp"
    MUELU_TESTING_SET_OSTREAM;
    MUELU_TESTING_LIMIT_SCOPE(Scalar,GlobalOrdinal,NO);
    out << "version: " << MueLu::Version() << std::endl;

    // unit test for block size 3 = (2,1). wrap block 0

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(2));
    stridingInfo.push_back(as<size_t>(1));
    LocalOrdinal stridedBlockId = 0;

    int blockSize=3;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, blockSize*comm->getSize(), 0,
                                  stridingInfo, comm,
                                  stridedBlockId /*blockId*/, 0 /*offset*/);

    /////////////////////////////////////////////////////

    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  mtx->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  0 /*offset*/
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            mtx->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            0 /*offset*/
                                            );
    if(mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
    mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", mtx);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == blockSize, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == blockSize) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t maxLocalIndex = myImportMap->getMaxLocalIndex();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-2), true);
      } else {
        TEST_EQUALITY(as<bool>(maxLocalIndex==numLocalRowMapElts*blockSize-1), true);
      }
    }

    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getMaxLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStrided2LW


  TEUCHOS_UNIT_TEST(CoalesceDropFactory, AmalgamationStridedOffsetDropping2LW)
  {
    // unit test for block size 9 = (2,3,4). wrap block 1.
    // drop small entries
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = Parameters::getDefaultComm();
    Xpetra::UnderlyingLib lib = TestHelpers::Parameters::getLib();

    // create strided map information
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(as<size_t>(2));
    stridingInfo.push_back(as<size_t>(3));
    stridingInfo.push_back(as<size_t>(4));
    LocalOrdinal stridedBlockId = 1;
    GlobalOrdinal offset = 19;

    RCP<const StridedMap> dofMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 9*comm->getSize(), 0,
                                  stridingInfo, comm,
                                  stridedBlockId, offset);

    /////////////////////////////////////////////////////

    Teuchos::RCP<Matrix> mtx = TestHelpers::TestFactory<SC,LO,GO,NO>::BuildTridiag(dofMap, 2.0, 1.0, 0.0001);

    Level fineLevel;
    TestHelpers::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

    RCP<const Map> stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                                  mtx->getRangeMap(),
                                                  stridingInfo,
                                                  stridedBlockId,
                                                  offset
                                                  );
    RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(
                                            mtx->getDomainMap(),
                                            stridingInfo,
                                            stridedBlockId,
                                            offset
                                            );

    if(mtx->IsView("stridedMaps") == true) mtx->RemoveView("stridedMaps");
    mtx->CreateView("stridedMaps", stridedRangeMap, stridedDomainMap);

    fineLevel.Set("A", mtx);
    RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory);
    CoalesceDropFactory dropFact = CoalesceDropFactory();
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(0.3));

    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    fineLevel.print(out);
    RCP<GraphBase> graph = fineLevel.Get<RCP<GraphBase> >("Graph", &dropFact);

    LO myDofsPerNode = fineLevel.Get<LO>("DofsPerNode", &dropFact);
    TEST_EQUALITY(as<int>(graph->GetDomainMap()->getGlobalNumElements()) == comm->getSize(), true);
    TEST_EQUALITY(as<int>(myDofsPerNode) == 9, true);
    bool bCorrectGraph = false;
    if (comm->getSize() == 1 && graph->getNeighborVertices(0).size() == 1) {
      bCorrectGraph = true;
    } else {
      if (comm->getRank() == 0) {
        if (graph->getNeighborVertices(0).size() == 1) bCorrectGraph = true;
      }
      else {
        if (graph->getNeighborVertices(0).size() == 2) bCorrectGraph = true;
      }
    }
    TEST_EQUALITY(bCorrectGraph, true);

    const RCP<const Map> myImportMap = graph->GetImportMap(); // < note that the ImportMap is built from the column map of the matrix A WITHOUT dropping!
    const RCP<const Map> myDomainMap = graph->GetDomainMap();

    TEST_EQUALITY(myImportMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myImportMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myImportMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myImportMap->getGlobalNumElements(),as<size_t>(comm->getSize()+2*(comm->getSize()-1)));
    if (comm->getSize()>1) {
      size_t numLocalRowMapElts = graph->GetNodeNumVertices();
      size_t numLocalImportElts = myImportMap->getLocalNumElements();
      if (comm->getRank() == 0 || comm->getRank() == comm->getSize()-1) {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+1), true);
      } else {
        TEST_EQUALITY(as<bool>(numLocalImportElts==numLocalRowMapElts+2), true);
      }
    }
    TEST_EQUALITY(myDomainMap->getMaxAllGlobalIndex(), comm->getSize()-1);
    TEST_EQUALITY(myDomainMap->getMinAllGlobalIndex(), 0);
    TEST_EQUALITY(myDomainMap->getMinLocalIndex(),0);
    TEST_EQUALITY(myDomainMap->getGlobalNumElements(),as<size_t>(comm->getSize()));
    TEST_EQUALITY(as<bool>(myDomainMap->getLocalNumElements()==1), true);
  } // AmalgamationStridedOffsetDropping2LW
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, AggresiveDroppingIsMarkedAsBoundary, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
  // Test that when everything but the diagonal is dropped, the node is marked as boundary
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  RCP<const Map> dofMap    = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(lib, 12 * comm->getSize(), 0, comm);
  Teuchos::RCP<Matrix> mtx = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::BuildTridiag(dofMap, 2.0, -1.0, -1.0);

  {
    Level fineLevel;
    TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

    mtx->SetFixedBlockSize(1);
    fineLevel.Set("A", mtx);

    CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
    RCP<AmalgamationFactory> amalgFact  = rcp(new AmalgamationFactory());
    dropFact.SetFactory("UnAmalgamationInfo", amalgFact);
    dropFact.SetParameter("aggregation: drop tol", Teuchos::ParameterEntry(4.1));
    fineLevel.Request("Graph", &dropFact);
    fineLevel.Request("DofsPerNode", &dropFact);

    dropFact.Build(fineLevel);

    RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", &dropFact);

    auto boundaryNodes     = graph->GetBoundaryNodeMap();
    auto boundaryNodesHost = Kokkos::create_mirror_view(boundaryNodes);
    Kokkos::deep_copy(boundaryNodesHost, boundaryNodes);
    bool allNodesAreOnBoundary = true;
    for (LO i = 0; i < Teuchos::as<LO>(boundaryNodesHost.size()); i++)
      allNodesAreOnBoundary &= boundaryNodesHost(i);
    TEST_EQUALITY(allNodesAreOnBoundary, true);
  }

  // {
  //   Level fineLevel;
  //   TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

  //   mtx->SetFixedBlockSize(2);
  //   fineLevel.Set("A", mtx);

  //   CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  //   RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  //   dropFact.SetFactory("UnAmalgamationInfo",amalgFact);
  //   dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(4.1));
  //   fineLevel.Request("Graph", &dropFact);
  //   fineLevel.Request("DofsPerNode", &dropFact);

  //   dropFact.Build(fineLevel);

  //   RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);

  //   auto boundaryNodes = graph->GetBoundaryNodeMap();
  //   bool allNodesAreOnBoundary = true;
  //   for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
  //     allNodesAreOnBoundary &= boundaryNodes[i];
  //   TEST_EQUALITY(allNodesAreOnBoundary, true);
  // }

  // {
  //   Level fineLevel;
  //   TestHelpers_kokkos::TestFactory<SC,LO,GO,NO>::createSingleLevelHierarchy(fineLevel);

  //   mtx->SetFixedBlockSize(3);
  //   fineLevel.Set("A", mtx);

  //   CoalesceDropFactory_kokkos dropFact = CoalesceDropFactory_kokkos();
  //   RCP<AmalgamationFactory> amalgFact = rcp(new AmalgamationFactory());
  //   dropFact.SetFactory("UnAmalgamationInfo",amalgFact);
  //   dropFact.SetParameter("aggregation: drop tol",Teuchos::ParameterEntry(4.1));
  //   fineLevel.Request("Graph", &dropFact);
  //   fineLevel.Request("DofsPerNode", &dropFact);

  //   dropFact.Build(fineLevel);

  //   RCP<LWGraph_kokkos> graph = fineLevel.Get<RCP<LWGraph_kokkos> >("Graph", &dropFact);

  //   auto boundaryNodes = graph->GetBoundaryNodeMap();
  //   bool allNodesAreOnBoundary = true;
  //   for (LO i = 0; i < Teuchos::as<LO>(boundaryNodes.size()); i++)
  //     allNodesAreOnBoundary &= boundaryNodes[i];
  //   TEST_EQUALITY(allNodesAreOnBoundary, true);
  // }

}  // AggresiveDroppingIsMarkedAsBoundary

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoalesceDropFactory_kokkos, 2x2, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include <MueLu_UseShortNames.hpp>
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int>> comm = Parameters::getDefaultComm();
  Xpetra::UnderlyingLib lib          = TestHelpers_kokkos::Parameters::getLib();

  // Set up a 2x2 matrix.
  Teuchos::RCP<Matrix> mtx = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::build2x2(lib,
                                                                                       2.0, -1.0,
                                                                                       -1.5, 2.0);
  mtx->SetFixedBlockSize(1);

  using magnitudeType = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
  Teuchos::RCP<Xpetra::MultiVector<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>> coords;
  coords = Xpetra::MultiVectorFactory<magnitudeType, LocalOrdinal, GlobalOrdinal, Node>::Build(mtx->getRowMap(), 1);
  {
    auto lclCoords  = coords->getHostLocalView(Xpetra::Access::OverwriteAll);
    auto rank       = comm->getRank();
    lclCoords(0, 0) = 2 * rank;
    lclCoords(1, 0) = 2 * rank + 1;
  }

  using TF                = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>;
  using local_matrix_type = typename Matrix::local_matrix_type::HostMirror;

  std::vector<Teuchos::ParameterList> params;
  std::vector<local_matrix_type> expectedFilteredMatrices;
  std::vector<std::vector<bool>> expectedBoundaryNodesVector;

  for (bool reuseGraph : {false, true}) {
    // test case 0
    Teuchos::ParameterList params0 = Teuchos::ParameterList();
    params0.set("aggregation: Dirichlet threshold", 0.);
    // dropFact.SetParameter("aggregation: row sum drop tol", Teuchos::ParameterEntry(0.2));

    params0.set("aggregation: drop scheme", "classical");
    params0.set("aggregation: use ml scaling of drop tol", false);
    params0.set("aggregation: drop tol", 0.);
    params0.set("aggregation: dropping may create Dirichlet", false);
    params0.set("filtered matrix: reuse graph", reuseGraph);
    params0.set("filtered matrix: use lumping", false);
    params.push_back(params0);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});

    // test case 1
    Teuchos::ParameterList params1 = Teuchos::ParameterList(params0);
    params1.set("aggregation: Dirichlet threshold", 0.2);
    params.push_back(params1);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});

    // test case 2
    Teuchos::ParameterList params2 = Teuchos::ParameterList(params0);
    params2.set("aggregation: Dirichlet threshold", 1.1);
    params.push_back(params2);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, false});

    // test case 3
    Teuchos::ParameterList params3 = Teuchos::ParameterList(params0);
    params3.set("aggregation: drop tol", 0.51);
    params.push_back(params3);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, 0.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});

    // test case 4
    Teuchos::ParameterList params4 = Teuchos::ParameterList(params0);
    params4.set("aggregation: Dirichlet threshold", 1.1);
    params4.set("aggregation: drop tol", 0.51);
    params.push_back(params4);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, 0.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, false});

    // test case 5
    Teuchos::ParameterList params5 = Teuchos::ParameterList(params0);
    params5.set("aggregation: Dirichlet threshold", 1.1);
    params5.set("aggregation: dropping may create Dirichlet", true);
    params5.set("aggregation: drop tol", 0.51);
    params.push_back(params5);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, 0.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, false});

    // test case 6
    Teuchos::ParameterList params6 = Teuchos::ParameterList(params4);
    params6.set("aggregation: Dirichlet threshold", 1.1);
    params6.set("filtered matrix: use lumping", true);
    params.push_back(params6);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(1.0, 0.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, false});

    // test case 7
    Teuchos::ParameterList params7 = Teuchos::ParameterList(params0);
    params7.set("aggregation: drop scheme", "distance laplacian");
    params.push_back(params7);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});

    // test case 8
    Teuchos::ParameterList params8 = Teuchos::ParameterList(params0);
    params8.set("aggregation: drop scheme", "distance laplacian");
    params8.set("aggregation: drop tol", 1.01);
    params.push_back(params8);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, 0.0,
                                                             0.0, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({true, true});

    // test case 9
    Teuchos::ParameterList params9 = Teuchos::ParameterList(params0);
    params9.set("aggregation: drop scheme", "classical");
    params9.set("aggregation: classical algo", "unscaled cut");
    params9.set("aggregation: drop tol", 1.0 / 3.6);
    params.push_back(params9);
    expectedFilteredMatrices.push_back(TF::buildLocal2x2Host(2.0, -1.0,
                                                             -1.5, 2.0, reuseGraph));
    expectedBoundaryNodesVector.push_back({false, false});
  }

  for (size_t testNo = 0; testNo < params.size(); ++testNo) {
    out << "\n\nRunning test number " << testNo << std::endl
        << std::endl;

    auto param                  = params[testNo];
    auto expectedFilteredMatrix = expectedFilteredMatrices[testNo];
    auto expectedBoundaryNodes  = expectedBoundaryNodesVector[testNo];

    for (bool useKokkos : {false, true}) {
      RCP<Matrix> filteredA;
      Kokkos::View<bool *, Kokkos::HostSpace> boundaryNodes;
      {
        Level fineLevel;
        TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);
        fineLevel.Set("A", mtx);
        fineLevel.Set("Coordinates", coords);

        RCP<MueLu::SingleLevelFactoryBase> dropFact;
        RCP<MueLu::SingleLevelFactoryBase> filteredAFact;
        if (useKokkos) {
          out << "\nKokkos code path\nparams:\n"
              << param << "\n";
          dropFact = rcp(new CoalesceDropFactory_kokkos());
          dropFact->SetParameterList(param);
          filteredAFact = dropFact;
        } else {
          out << "\nNon-Kokkos code path\nparams:\n"
              << param << "\n";
          dropFact            = rcp(new CoalesceDropFactory());
          filteredAFact       = rcp(new FilteredAFactory());
          auto paramCopy      = Teuchos::ParameterList(param);
          auto paramFilteredA = Teuchos::ParameterList();
          paramFilteredA.set("filtered matrix: use lumping", paramCopy.get<bool>("filtered matrix: use lumping"));
          paramCopy.remove("filtered matrix: use lumping");
          paramFilteredA.set("filtered matrix: reuse graph", paramCopy.get<bool>("filtered matrix: reuse graph"));
          paramCopy.remove("filtered matrix: reuse graph");
          paramFilteredA.set("filtered matrix: Dirichlet threshold", paramCopy.get<double>("aggregation: Dirichlet threshold"));
          dropFact->SetParameterList(paramCopy);
          filteredAFact->SetParameterList(paramFilteredA);
          filteredAFact->SetFactory("Graph", dropFact);
          filteredAFact->SetFactory("Filtering", dropFact);
        }
        fineLevel.Request("A", filteredAFact.get());
        fineLevel.Request("Graph", dropFact.get());
        fineLevel.Request("DofsPerNode", dropFact.get());
        dropFact->Build(fineLevel);

        filteredA = fineLevel.Get<RCP<Matrix>>("A", filteredAFact.get());

        if (useKokkos) {
          auto graph           = fineLevel.Get<RCP<LWGraph_kokkos>>("Graph", dropFact.get());
          auto boundaryNodes_d = graph->GetBoundaryNodeMap();
          boundaryNodes        = Kokkos::View<bool *, Kokkos::HostSpace>("boundaryNodes_host", boundaryNodes_d.extent(0));
          Kokkos::deep_copy(boundaryNodes, boundaryNodes_d);
        } else {
          auto graph    = fineLevel.Get<RCP<LWGraph>>("Graph", dropFact.get());
          boundaryNodes = graph->GetBoundaryNodeMap();
        }
      }

      auto lclFilteredA         = filteredA->getLocalMatrixHost();
      auto lclExpectedFilteredA = expectedFilteredMatrix;

      out << "Filtered A:\n"
          << TF::localMatToString(lclFilteredA);

      TEST_EQUALITY(lclFilteredA.graph.row_map.extent(0), lclExpectedFilteredA.graph.row_map.extent(0));
      for (size_t row = 0; row < std::min(lclFilteredA.graph.row_map.extent(0), lclExpectedFilteredA.graph.row_map.extent(0)); ++row)
        TEST_EQUALITY(lclFilteredA.graph.row_map(row), lclExpectedFilteredA.graph.row_map(row));

      TEST_EQUALITY(lclFilteredA.graph.entries.extent(0), lclExpectedFilteredA.graph.entries.extent(0));
      for (size_t entry = 0; entry < std::min(lclFilteredA.graph.entries.extent(0), lclExpectedFilteredA.graph.entries.extent(0)); ++entry) {
        TEST_EQUALITY(lclFilteredA.graph.entries(entry), lclExpectedFilteredA.graph.entries(entry));
        TEST_EQUALITY(lclFilteredA.values(entry), lclExpectedFilteredA.values(entry));
      }

      auto boundaryNodes_h = Kokkos::create_mirror_view(boundaryNodes);
      Kokkos::deep_copy(boundaryNodes_h, boundaryNodes);
      TEST_EQUALITY(boundaryNodes_h.extent(0), expectedBoundaryNodes.size());
      for (size_t row = 0; row < std::min(boundaryNodes_h.extent(0), expectedBoundaryNodes.size()); ++row)
        TEST_EQUALITY(boundaryNodes_h(row), expectedBoundaryNodes[row]);
    }
    out << "Done with test number " << testNo << std::endl;
  }
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO)                                                                                    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, Constructor, SC, LO, GO, NO)                            \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, Build, SC, LO, GO, NO)                                  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacian, SC, LO, GO, NO)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianScaledCut, SC, LO, GO, NO)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianUnscaledCut, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianCutSym, SC, LO, GO, NO)                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicalScaledCut, SC, LO, GO, NO)                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicalUnScaledCut, SC, LO, GO, NO)                   \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicalCutSym, SC, LO, GO, NO)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, SignedClassical, SC, LO, GO, NO)                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, SignedScaledCutClassical, SC, LO, GO, NO)               \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, SignedUnscaledCutClassical, SC, LO, GO, NO)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalColoredSignedClassical, SC, LO, GO, NO)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalNoColoredSignedClassical, SC, LO, GO, NO)  \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalSignedClassical, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonal, SC, LO, GO, NO)                          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalDistanceLaplacian, SC, LO, GO, NO)         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, BlockDiagonalDistanceLaplacianWeighted, SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, DistanceLaplacianWeighted, SC, LO, GO, NO)              \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, SignedClassicalSA, SC, LO, GO, NO)                      \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicScalarWithoutFiltering, SC, LO, GO, NO)          \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicScalarWithFiltering, SC, LO, GO, NO)             \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicBlockWithoutFiltering, SC, LO, GO, NO)           \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, AggresiveDroppingIsMarkedAsBoundary, SC, LO, GO, NO)    \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, 2x2, SC, LO, GO, NO)

// TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoalesceDropFactory_kokkos, ClassicBlockWithFiltering,     SC, LO, GO, NO) // not implemented yet

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
