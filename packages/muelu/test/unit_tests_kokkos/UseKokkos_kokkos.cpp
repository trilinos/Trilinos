// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_UseDefaultTypes.hpp"

#include "MueLu_TestHelpers_kokkos.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_HierarchyManager.hpp"
#include "MueLu_Hierarchy.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(UseKokkos_kokkos, FactoryManager, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, Node);
  out << "version: " << MueLu::Version() << std::endl;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using real_type             = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using RealValuedMultiVector = Xpetra::MultiVector<real_type, LO, GO, NO>;

  const bool useKokkos = false;

  // Build the problem
  RCP<Matrix> A                          = MueLuTests::TestHelpers_kokkos::TestFactory<SC, LO, GO, NO>::Build1DPoisson(2001);
  RCP<const Map> map                     = A->getMap();
  RCP<RealValuedMultiVector> coordinates = Xpetra::MultiVectorFactory<real_type, LO, GO, NO>::Build(map, 1);
  RCP<MultiVector> nullspace             = MultiVectorFactory::Build(map, 1);
  nullspace->putScalar(Teuchos::ScalarTraits<SC>::one());

  // Setting paramters for the hierarchy
  Teuchos::ParameterList paramList;
  paramList.set("verbosity", "test");
  paramList.set("use kokkos refactor", useKokkos);
  RCP<HierarchyManager> mueluFactory = rcp(new ParameterListInterpreter(paramList));
  mueluFactory->CheckConfig();
  const RCP<const FactoryManagerBase> manager = mueluFactory->GetFactoryManager(0);
  RCP<const FactoryManager> factoryManager    = Teuchos::rcp_dynamic_cast<const FactoryManager>(manager, true);

  TEST_EQUALITY(factoryManager->GetKokkosRefactor() == useKokkos, true);
}

#define MUELU_ETI_GROUP(SC, LO, GO, NO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(UseKokkos_kokkos, FactoryManager, SC, LO, GO, NO)

#include <MueLu_ETI_4arg.hpp>

}  // namespace MueLuTests
