// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tBlockLowerTriInverseOp_tpetra.hpp"

#include "Thyra_LinearOpTester.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::RCP;
using Teuchos::rcp;

static const RCP<const Thyra::LinearOpBase<ST> > build2x2op(
    const RCP<const Teuchos::Comm<int> > comm, ST a, ST b, ST c, ST d) {
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));

  GO indices[2];
  ST row0[2];
  ST row1[2];

  indices[0] = 0;
  indices[1] = 1;

  // build a CrsMatrix
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > blk = Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  row0[0]                                     = a;
  row0[1]                                     = c;
  row1[0]                                     = b;
  row1[1]                                     = d;
  blk->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices, 2), Teuchos::ArrayView<ST>(row0, 2));
  blk->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices, 2), Teuchos::ArrayView<ST>(row1, 2));
  blk->fillComplete();

  return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(blk->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(blk->getRangeMap()), blk);
}

void tBlockLowerTriInverseOp_tpetra::initializeTest() {
  std::vector<GO> indices(2);
  std::vector<ST> row0(2), row1(2);

  RCP<const Teuchos::Comm<int> > comm     = GetComm_tpetra();
  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));

  tolerance_ = 1.0e-11;
  RCP<const Thyra::LinearOpBase<ST> > blk;

  // build forward operator
  A_ = rcp(new Thyra::DefaultBlockedLinearOp<ST>());
  A_->beginBlockFill(3, 3);

  // build 0,0 matrix
  blk = build2x2op(comm, 1.0, 3.0, 2.0, -1.0);
  A_->setBlock(0, 0, blk);

  // build 1,0 matrix
  blk = build2x2op(comm, 2.0, 9.0, 8.0, 3.0);
  A_->setBlock(1, 0, blk);

  // build 1,1 matrix
  blk = build2x2op(comm, 7.0, 8.0, -2.0, 4.0);
  A_->setBlock(1, 1, blk);

  // build 2,1 matrix
  blk = build2x2op(comm, -1.0, 6.0, 2.0, 1.0);
  A_->setBlock(2, 1, blk);

  // build 2,2 matrix
  blk = build2x2op(comm, 3.0, 9.0, 7.0, 1.0);
  A_->setBlock(2, 2, blk);

  A_->endBlockFill();

  // build inverse operator
  invA_ = rcp(new Thyra::DefaultBlockedLinearOp<ST>());
  invA_->beginBlockFill(3, 3);

  // build 0,0 matrix
  blk =
      build2x2op(comm, 0.142857142857143, 0.428571428571429, 0.285714285714286, -0.142857142857143);
  invA_->setBlock(0, 0, blk);
  invDiag_.push_back(blk);

  // build 1,0 matrix
  blk = build2x2op(comm, -0.454545454545455, 0.266233766233766, -0.045454545454545,
                   -0.444805194805195);
  invA_->setBlock(1, 0, blk);

  // build 2,0 matrix
  blk =
      build2x2op(comm, 0.303571428571429, -0.271103896103896, 0.069642857142857, 0.090746753246753);
  invA_->setBlock(2, 0, blk);

  // build 1,1 matrix
  blk =
      build2x2op(comm, 0.090909090909091, -0.181818181818182, 0.045454545454545, 0.159090909090909);
  invA_->setBlock(1, 1, blk);
  invDiag_.push_back(blk);

  // build 2,1 matrix
  blk = build2x2op(comm, -0.050000000000000, 0.086363636363636, -0.045833333333333,
                   -0.019318181818182);
  invA_->setBlock(2, 1, blk);

  // build 2,2 matrix
  blk = build2x2op(comm, -0.016666666666667, 0.150000000000000, 0.116666666666667,
                   -0.050000000000000);
  invA_->setBlock(2, 2, blk);
  invDiag_.push_back(blk);

  invA_->endBlockFill();
}

int tBlockLowerTriInverseOp_tpetra::runTest(int verbosity, std::ostream& stdstrm,
                                            std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tLowerTriInverseOp";

  status = test_apply(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"apply\" ... PASSED", "   \"apply\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_alphabeta(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"alphabeta\" ... PASSED", "   \"alphabeta\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tLowerTriInverseOp...PASSED", "tLowerTriInverseOp...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tLowerTriInverseOp...FAILED");
  }

  return failcount;
}

bool tBlockLowerTriInverseOp_tpetra::test_alphabeta(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  ST diff;

  BlockedLinearOp U = getLowerTriBlocks(A_);
  LinearOp invTri   = createBlockLowerTriInverseOp(U, invDiag_);

  RCP<Thyra::VectorBase<ST> > src  = Thyra::createMember(invA_->domain());
  RCP<Thyra::VectorBase<ST> > dste = Thyra::createMember(invA_->range());

  Thyra::randomize<ST>(-10, 10, src.ptr());
  Thyra::randomize<ST>(-10, 10, dste.ptr());

  RCP<Thyra::VectorBase<ST> > dstn = dste->clone_v();

  diff = Teko::Test::Difference(dste, dstn);
  TEST_ASSERT(diff <= 0.0, std::endl
                               << "   tBlockLowerTriInverseOp_tpetra::test_apply "
                               << toString(status) << ": exact copy failed (abserr=" << diff
                               << " <= " << 0.0 << ")");

  MultiVector dste_mv = dste;
  MultiVector dstn_mv = dstn;

  applyOp(invA_, src, dste_mv, 3.2, -1.9);
  applyOp(invTri, src, dstn_mv, 3.2, -1.9);

  diff = Teko::Test::Difference(dste, dstn) / Thyra::norm_2(*dste);
  TEST_ASSERT(diff <= tolerance_, std::endl
                                      << "   tBlockLowerTriInverseOp_tpetra::test_apply "
                                      << toString(status)
                                      << ": alpha/beta apply operation failed (relerr=" << diff
                                      << " <= " << tolerance_ << ")");

  return allPassed;
}

bool tBlockLowerTriInverseOp_tpetra::test_apply(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<Thyra::PhysicallyBlockedLinearOpBase<ST> > U = getLowerTriBlocks(A_);
  RCP<const Thyra::LinearOpBase<ST> > invTri       = createBlockLowerTriInverseOp(U, invDiag_);

  Thyra::LinearOpTester<ST> tester;
  tester.show_all_tests(true);
  std::stringstream ss;
  Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss), "      |||");
  const bool result = tester.compare(*invA_, *invTri, Teuchos::ptrFromRef(fos));
  TEST_ASSERT(result, std::endl
                          << "   tBlockLowerTriInverseOp_tpetra::test_apply "
                          << ": Comparing implicitly generated operator to exact operator");
  if (not result || verbosity >= 10) os << ss.str();

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
