// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA) && defined(HAVE_MUELU_THYRATPETRAADAPTERS) && defined(HAVE_MUELU_TEKO)

#include "MueLu_TestHelpers.hpp"
#include "MueLu_TekoMultiPhysBlockPreconditionerFactory.hpp"
#include "Stratimikos_MueLuHelpers.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"

#include "Xpetra_MultiVectorFactory.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST(TekoMultiPhysBlockPreconditionerFactory, Registration) {
  Stratimikos::registerMueLuTekoMultiPhysBlockAMG("Unit Test Block AMG");
  Teuchos::ParameterList settings;
  auto factory = Teko::PreconditionerFactory::buildPreconditionerFactory("Unit Test Block AMG", settings);
  TEST_ASSERT(!Teuchos::is_null(factory));
}

TEUCHOS_UNIT_TEST(TekoMultiPhysBlockPreconditionerFactory, RejectNonContiguousSubblocks) {
  MueLu::TekoMultiPhysBlockPreconditionerFactory factory;
  Teuchos::ParameterList params;
  auto& amg = params.sublist("AMG Settings");
  amg.sublist("subblockList0");
  amg.sublist("subblockList2");

  TEST_THROW(factory.initializeFromParameterList(params), MueLu::Exceptions::RuntimeError);
}

TEUCHOS_UNIT_TEST(TekoMultiPhysBlockPreconditionerFactory, AcceptTopLevelUserData) {
  using ST          = Teko::ST;
  using LO          = Teko::LO;
  using GO          = Teko::GO;
  using NT          = Teko::NT;
  using Map         = Xpetra::Map<LO, GO, NT>;
  using MultiVector = Xpetra::MultiVector<ST, LO, GO, NT>;

  MueLu::TekoMultiPhysBlockPreconditionerFactory factory;
  Teuchos::ParameterList params;
  auto& amg = params.sublist("AMG Settings");
  amg.sublist("subblockList0");
  amg.sublist("subblockList1");

  auto comm   = MueLuTests::TestHelpers::Parameters::getDefaultComm();
  auto map    = Xpetra::MapFactory<LO, GO, NT>::Build(Xpetra::UseTpetra, 4, 0, comm);
  auto coords = Xpetra::MultiVectorFactory<ST, LO, GO, NT>::Build(map, 1);
  params.sublist("user data").set<Teuchos::RCP<MultiVector> >("Coordinates0", coords);

  TEST_NOTHROW(factory.initializeFromParameterList(params));
}

}  // namespace MueLuTests

#endif