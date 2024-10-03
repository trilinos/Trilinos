// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tParallelInverse_tpetra.hpp"

#include <string>

// Teko-Package includes
#include "Teko_Config.h"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedTpetraOperator.hpp"

#include "Thyra_TpetraLinearOp.hpp"

// Tpetra includes
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraLinearOp.hpp"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Thyra::tpetraLinearOp;

void tParallelInverse_tpetra::initializeTest() { tolerance_ = 1.0e-7; }

void tParallelInverse_tpetra::loadMatrix() {
  // Read in the matrix, store pointer as an RCP
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrA =
      Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<ST, LO, GO, NT> >::readSparseFile(
          "./data/lsc_F_2.mm", GetComm_tpetra());
  F_ = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrA->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrA->getDomainMap()), ptrA);
}

void tParallelInverse_tpetra::loadStridedMatrix() {
  // Read in the matrix, store pointer as an RCP
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrA =
      Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<ST, LO, GO, NT> >::readSparseFile(
          "./data/nsjac.mm", GetComm_tpetra());
  std::vector<int> vec(2);
  vec[0] = 1;
  vec[1] = 2;

  // Block the linear system using a strided epetra operator
  Teko::TpetraHelpers::StridedTpetraOperator sA(vec, ptrA);

  // get 0,0 block
  F_ = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(sA.GetBlock(0, 0)->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(sA.GetBlock(0, 0)->getDomainMap()),
      sA.GetBlock(0, 0));
}

int tParallelInverse_tpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                                     int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tParallelInverse_tpetra";

  status = test_inverse(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"inverse\" ... PASSED", "   \"inverse\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

#ifdef Teko_ENABLE_DEV_MODE  // so the file nsjac.mm isn't required for release mode
  status = test_stridedInverse(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"stridedInverse\" ... PASSED",
                       "   \"stridedInverse\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;
#endif

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tParallelInverse_tpetra...PASSED",
                         "tParallelInverse_tpetra...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tParallelInverse_tpetra...FAILED");
  }

  return failcount;
}

bool tParallelInverse_tpetra::test_inverse(int verbosity, std::ostream& os) {
  // bool status = false;
  bool allPassed = true;

  loadMatrix();

  // build an InverseLibrary
  RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

  // build the inverse factory needed by the example preconditioner
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("Belos");

  Teko::LinearOp inv = invFact->buildInverse(F_);

  return allPassed;
}

bool tParallelInverse_tpetra::test_stridedInverse(int verbosity, std::ostream& os) {
  // bool status = false;
  bool allPassed = true;

  loadStridedMatrix();

  // build an InverseLibrary
  TEST_MSG("\n   Building inverse library");
  RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

  // build the inverse factory needed by the example preconditioner
  TEST_MSG("\n   Building inverse factory");
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("Belos");

  TEST_MSG("\n   Building inverse");
  Teko::LinearOp inv = invFact->buildInverse(F_);

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
