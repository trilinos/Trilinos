// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_TestingUtilities.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Reindex_LinearProblem.hpp>
#include "TestCase1_decl.hpp"
#include "TestCase2_decl.hpp"

namespace {

TEUCHOS_STATIC_SETUP() {
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption(
      "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
      "Test MPI (if available) or force test of serial.  In a serial build,"
      " this option is ignored and a serial comm is always used.");
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(ReindexTransform, TestCase_1, LO, GO, Node) {
  using Scalar_t      = Tpetra::CrsMatrix<>::scalar_type;
  using Map_t         = Tpetra::Map<LO, GO, Node>;
  using CrsMatrix_t   = Tpetra::CrsMatrix<Scalar_t, LO, GO, Node>;
  using MultiVector_t = Tpetra::MultiVector<Scalar_t, LO, GO, Node>;
  using Problem_t     = Tpetra::LinearProblem<Scalar_t, LO, GO, Node>;
  auto comm           = Tpetra::getDefaultComm();

  // Create the linear problem related to 'test case 1'
  TestCase1<Scalar_t, LO, GO, Node> testCase(comm);
  Teuchos::RCP<Problem_t> linearProblem = Teuchos::rcp<Problem_t>(testCase.linearProblem(), false);

  // Save the linear problem as 'originalProblem'
  CrsMatrix_t *auxMat = dynamic_cast<CrsMatrix_t *>(testCase.linearProblem()->getMatrix().get());
  CrsMatrix_t originalMat(*auxMat);
  MultiVector_t originalLhs(*(testCase.linearProblem()->getLHS()));
  MultiVector_t originalRhs(*(testCase.linearProblem()->getRHS()));

  Teuchos::RCP<CrsMatrix_t> originalMatRCP(&originalMat, false);
  Teuchos::RCP<MultiVector_t> originalLhsRCP(&originalLhs, false);
  Teuchos::RCP<MultiVector_t> originalRhsRCP(&originalRhs, false);

  Problem_t originalProblem(originalMatRCP, originalLhsRCP, originalRhsRCP);

  // Create the transform
  Teuchos::RCP<Map_t const> newRowMap(Teuchos::null);
  Tpetra::Reindex_LinearProblem<Scalar_t, LO, GO, Node> reindexTransform(newRowMap);

  // Transform the linear problem
  Teuchos::RCP<Problem_t> transformedLinearProblem = reindexTransform(linearProblem);
  reindexTransform.fwd();

  // Check the transformed problem
  {
    bool checkResult = testCase.checkTransformedProblem(transformedLinearProblem.get());
    TEUCHOS_ASSERT(checkResult);
  }

  // Call transform.rvs()
  reindexTransform.rvs();

  // Check 'testCase.linearProblem', after call to rvs(), against 'originalProblem'
  {
    bool checkResult = testCase.checkAfterFwdRvs(&originalProblem);
    TEUCHOS_ASSERT(checkResult);
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(ReindexTransform, TestCase_2, LO, GO, Node) {
  using Scalar_t      = Tpetra::CrsMatrix<>::scalar_type;
  using Map_t         = Tpetra::Map<LO, GO, Node>;
  using CrsMatrix_t   = Tpetra::CrsMatrix<Scalar_t, LO, GO, Node>;
  using MultiVector_t = Tpetra::MultiVector<Scalar_t, LO, GO, Node>;
  using Problem_t     = Tpetra::LinearProblem<Scalar_t, LO, GO, Node>;
  auto comm           = Tpetra::getDefaultComm();

  // Create the linear problem related to 'test case 2'
  TestCase2<Scalar_t, LO, GO, Node> testCase(comm);
  Teuchos::RCP<Problem_t> linearProblem = Teuchos::rcp<Problem_t>(testCase.linearProblem(), false);

  // Save the linear problem as 'originalProblem'
  CrsMatrix_t *auxMat = dynamic_cast<CrsMatrix_t *>(testCase.linearProblem()->getMatrix().get());
  CrsMatrix_t originalMat(*auxMat);
  MultiVector_t originalLhs(*(testCase.linearProblem()->getLHS()));
  MultiVector_t originalRhs(*(testCase.linearProblem()->getRHS()));

  Teuchos::RCP<CrsMatrix_t> originalMatRCP(&originalMat, false);
  Teuchos::RCP<MultiVector_t> originalLhsRCP(&originalLhs, false);
  Teuchos::RCP<MultiVector_t> originalRhsRCP(&originalRhs, false);

  Problem_t originalProblem(originalMatRCP, originalLhsRCP, originalRhsRCP);

  // Create the transform
  Teuchos::RCP<Map_t const> newRowMap(Teuchos::null);
  Tpetra::Reindex_LinearProblem<Scalar_t, LO, GO, Node> reindexTransform(newRowMap);

  // Transform the linear problem
  Teuchos::RCP<Problem_t> transformedLinearProblem = reindexTransform(linearProblem);
  reindexTransform.fwd();

  // Check the transformed problem
  {
    bool checkResult = testCase.checkTransformedProblem(transformedLinearProblem.get());
    TEUCHOS_ASSERT(checkResult);
  }

  // Call transform.rvs()
  reindexTransform.rvs();

  // Check 'testCase.linearProblem', after call to rvs(), against 'originalProblem'
  {
    bool checkResult = testCase.checkAfterFwdRvs(&originalProblem);
    TEUCHOS_ASSERT(checkResult);
  }
}

#define UNIT_TEST_GROUP(LO, GO, NODE)                                              \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(ReindexTransform, TestCase_1, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(ReindexTransform, TestCase_2, LO, GO, NODE)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_LGN(UNIT_TEST_GROUP)

}  // namespace
