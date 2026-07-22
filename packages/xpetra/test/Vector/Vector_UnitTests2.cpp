// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_UnitTestHelpers.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Tpetra_Vector.hpp>

namespace {

bool testMpi = true;

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
}

//
// UNIT TESTS
//

// Test if the reference returned by toTpetra is valid.
// cf. bug #5603
TEUCHOS_UNIT_TEST(Vector, toXpetra) {
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::RCP<Xpetra::Map<int> > map            = Xpetra::MapFactory<int>::Build(Xpetra::UseTpetra, 1000, 0, comm);
  Teuchos::RCP<Xpetra::Vector<double, int> > vec = Xpetra::VectorFactory<double, int>::Build(map);

  Tpetra::Vector<double, int> *tVec = &Xpetra::toTpetra<double, int>(*vec);
  tVec->putScalar(Teuchos::ScalarTraits<double>::one());

  std::cout << tVec->isDistributed() << std::endl;
  TEST_EQUALITY_CONST(tVec->meanValue(), 1);
}

}  // namespace
