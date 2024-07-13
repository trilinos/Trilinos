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

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Utilities.hpp"
#include "MueLu_BlockedRAPFactory.hpp"

#include "MueLu_Exceptions.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(BlockedRAPFactory, Constructor, Scalar, LocalOrdinal, GlobalOrdinal, Node) {
#include "MueLu_UseShortNames.hpp"
  MUELU_TESTING_SET_OSTREAM;
  MUELU_TESTING_LIMIT_SCOPE(Scalar, GlobalOrdinal, NO);

  out << "version: " << MueLu::Version() << std::endl;

  RCP<BlockedRAPFactory> blockedRAPFactory = rcp(new BlockedRAPFactory());
  TEST_EQUALITY(blockedRAPFactory != Teuchos::null, true);
}  // Constructor

#define MUELU_ETI_GROUP(SC, LO, GO, Node) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(BlockedRAPFactory, Constructor, SC, LO, GO, Node)

#include <MueLu_ETI_4arg.hpp>
}  // namespace MueLuTests
