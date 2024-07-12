// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Aggregates.hpp"
#include "MueLu_CoarseMapFactory.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CoarseMap_kokkos, StandardCase, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);
  out << "version: " << MueLu::Version() << std::endl;

  Level fineLevel;
  TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::createSingleLevelHierarchy(fineLevel);

  RCP<Matrix> A = TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(15);
  fineLevel.Set("A", A);

  // build dummy aggregate structure
  RCP<Aggregates> aggs = Teuchos::rcp(new Aggregates(A->getRowMap()));
  aggs->SetNumAggregates(10);  // set (local!) number of aggregates
  fineLevel.Set("Aggregates", aggs);

  // build dummy nullspace vector
  RCP<MultiVector> nsp = MultiVectorFactory::Build(A->getRowMap(), 1);
  nsp->putScalar(1.0);
  fineLevel.Set("Nullspace", nsp);

  RCP<CoarseMapFactory> coarseMapFactory = Teuchos::rcp(new CoarseMapFactory());
  coarseMapFactory->SetFactory("Aggregates", MueLu::NoFactory::getRCP());
  coarseMapFactory->SetFactory("Nullspace", MueLu::NoFactory::getRCP());

  fineLevel.Request("CoarseMap", coarseMapFactory.get());
  coarseMapFactory->Build(fineLevel);

  auto myCoarseMap = fineLevel.Get<Teuchos::RCP<const Map> >("CoarseMap", coarseMapFactory.get());

  TEST_EQUALITY(myCoarseMap->getMinAllGlobalIndex() == 0, true);
  TEST_EQUALITY(myCoarseMap->getMaxLocalIndex() == 9, true);
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CoarseMap_kokkos, StandardCase, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
