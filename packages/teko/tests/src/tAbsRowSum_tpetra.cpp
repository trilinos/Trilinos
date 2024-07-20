// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tAbsRowSum_tpetra.hpp"

#include <string>

// Teko-Package includes
#include "Teko_Utilities.hpp"

// Thyra includes
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::tpetraLinearOp;

void tAbsRowSum_tpetra::initializeTest() { tolerance_ = 1.0e-15; }

int tAbsRowSum_tpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                               int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tAbsRowSum_tpetra";

  status = test_absRowSum(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"absRowSum\" ... PASSED", "   \"absRowSum\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_invAbsRowSum(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"invAbsRowSum\" ... PASSED",
                       "   \"invAbsRowSum\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tAbsRowSum_tpetra...PASSED", "tAbsRowSum_tpetra...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tAbsRowSum_tpetra...FAILED");
  }

  return failcount;
}

bool tAbsRowSum_tpetra::test_absRowSum(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  Thyra::LinearOpTester<ST> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(tolerance_);

  const RCP<const Tpetra::Map<LO, GO, NT> > map =
      rcp(new const Tpetra::Map<LO, GO, NT>(100, 0, GetComm_tpetra()));

  // A matrix...to be row summed
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > A = Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 5);
  GO indices[5];
  ST values[5] = {1, -2, 3, -4, 5};
  for (size_t i = 0; i < A->getLocalNumRows() - 5; i++) {
    GO index = A->getRowMap()->getGlobalElement(i);
    for (size_t j = 0; j < 5; j++) indices[j] = A->getRowMap()->getGlobalElement(i + j);
    A->insertGlobalValues(index, Teuchos::ArrayView<GO>(indices, 5),
                          Teuchos::ArrayView<ST>(values, 5));
  }
  for (size_t i = A->getLocalNumRows() - 5; i < A->getLocalNumRows(); i++) {
    GO index = A->getRowMap()->getGlobalElement(i);
    for (size_t j = 0; j < 5; j++) indices[j] = A->getRowMap()->getGlobalElement(j);
    A->insertGlobalValues(index, Teuchos::ArrayView<GO>(indices, 5),
                          Teuchos::ArrayView<ST>(values, 5));
  }
  A->fillComplete();

  // B matrix...already row summed
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > B = Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 1);
  ST number[1]                              = {15.0};
  for (size_t i = 0; i < B->getLocalNumRows(); i++) {
    GO index[1] = {B->getRowMap()->getGlobalElement(i)};
    B->insertGlobalValues(index[0], Teuchos::ArrayView<GO>(index, 1),
                          Teuchos::ArrayView<ST>(number, 1));
  }
  B->fillComplete();

  Teko::LinearOp pA = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A->getDomainMap()), A);
  Teko::LinearOp pB = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(B->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(B->getDomainMap()), B);
  Teko::LinearOp absRowA = getAbsRowSumMatrix(pA);

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*pB, *absRowA, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tAbsRowSum_tpetra::test_absRowSum "
                            << ": Testing basic sum functionality");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

bool tAbsRowSum_tpetra::test_invAbsRowSum(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  Thyra::LinearOpTester<ST> tester;
  tester.show_all_tests(true);
  tester.set_all_error_tol(tolerance_);

  Tpetra::Map<LO, GO, NT> map(100, 0, GetComm_tpetra());

  // A matrix...to be row summed
  Tpetra::CrsMatrix<ST, LO, GO, NT> A(rcpFromRef(map), 5);
  GO indices[5];
  ST values[5] = {1, -2, 3, -4, 5};
  for (size_t i = 0; i < A.getLocalNumRows() - 5; i++) {
    GO index = A.getRowMap()->getGlobalElement(i);
    for (size_t j = 0; j < 5; j++) indices[j] = A.getRowMap()->getGlobalElement(i + j);
    A.insertGlobalValues(index, Teuchos::ArrayView<GO>(indices, 5),
                         Teuchos::ArrayView<ST>(values, 5));
  }
  for (size_t i = A.getLocalNumRows() - 5; i < A.getLocalNumRows(); i++) {
    GO index = A.getRowMap()->getGlobalElement(i);
    for (size_t j = 0; j < 5; j++) indices[j] = A.getRowMap()->getGlobalElement(j);
    A.insertGlobalValues(index, Teuchos::ArrayView<GO>(indices, 5),
                         Teuchos::ArrayView<ST>(values, 5));
  }
  A.fillComplete();

  // B matrix...already row summed
  Tpetra::CrsMatrix<ST, LO, GO, NT> B(rcpFromRef(map), 1);
  ST number[1] = {1.0 / 15.0};
  for (size_t i = 0; i < B.getLocalNumRows(); i++) {
    GO index[1] = {B.getRowMap()->getGlobalElement(i)};
    B.insertGlobalValues(index[0], Teuchos::ArrayView<GO>(index, 1),
                         Teuchos::ArrayView<ST>(number, 1));
  }
  B.fillComplete();

  Teko::LinearOp pA = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A.getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(A.getDomainMap()), rcpFromRef(A));
  Teko::LinearOp pB = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(B.getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(B.getDomainMap()), rcpFromRef(B));
  Teko::LinearOp absRowA = getAbsRowSumInvMatrix(pA);

  {
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*pB, *absRowA, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tLummping::test_invAbsRowSum "
                            << ": Testing basic inv abs row sum functionality");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
