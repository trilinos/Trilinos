// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tJacobi2x2PreconditionerFactory_tpetra.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"

#include "Teko_InverseLibrary.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_SolveInverseFactory.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teko_Utilities.hpp"

#ifdef TEKO_HAVE_EPETRA
#include "Thyra_MLPreconditionerFactory.hpp"
#include "Thyra_IfpackPreconditionerFactory.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#endif

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"

#include <vector>

// This whole test rig is based on inverting the matrix
//
//      [  1  2  1 -1 ]
//  A = [  2  1 -3  1 ]
//      [  1 -3  0  0 ]
//      [ -1  1  0  0 ]
//
// see the matlab file

namespace Teko {
namespace Test {

using Teko::Test::toString;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using namespace Teko;

void tJacobi2x2PreconditionerFactory_tpetra::initializeTest() {
  std::vector<GO> indices(2);
  std::vector<ST> row0(2), row1(2);

  tolerance_ = 1.0e-14;

  comm                                    = GetComm_tpetra();
  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));

  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrF =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrG =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrD =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrC =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));

  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvF =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvC =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));

  indices[0] = 0;
  indices[1] = 1;

  // build F matrix
  row0[0] = 1.0;
  row0[1] = 2.0;
  row1[0] = 2.0;
  row1[1] = 1.0;
  ptrF->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrF->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrF->fillComplete();
  F_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrF->getRangeMap()), ptrF);

  // build D matrix
  row0[0] = 1.0;
  row0[1] = -3.0;
  row1[0] = -1.0;
  row1[1] = 1.0;
  ptrD->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrD->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrD->fillComplete();
  D_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrD->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrD->getRangeMap()), ptrD);

  // build G matrix
  row0[0] = 1.0;
  row0[1] = -1.0;
  row1[0] = -3.0;
  row1[1] = 1.0;
  ptrG->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrG->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrG->fillComplete();
  G_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrG->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrG->getRangeMap()), ptrG);

  // build C matrix
  row0[0] = 9.0;
  row0[1] = 2.0;
  row1[0] = 6.0;
  row1[1] = 5.0;
  ptrC->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrC->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrC->fillComplete();
  C_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrC->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrC->getRangeMap()), ptrC);

  // build inv(F) matrix
  row0[0] = -1.0 / 3.0;
  row0[1] = 2.0 / 3.0;
  row1[0] = 2.0 / 3.0;
  row1[1] = -1.0 / 3.0;
  ptrInvF->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrInvF->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrInvF->fillComplete();
  invF_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvF->getRangeMap()), ptrInvF);

  // build inv(C) matrix
  row0[0] = 0.151515151515151;
  row0[1] = -0.060606060606061;
  row1[0] = -0.181818181818182;
  row1[1] = 0.272727272727273;
  ptrInvC->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrInvC->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrInvC->fillComplete();
  invC_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvC->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvC->getRangeMap()), ptrInvC);

  A_ = Thyra::block2x2<ST>(F_, G_, D_, C_, "A");
}

int tJacobi2x2PreconditionerFactory_tpetra::runTest(int verbosity, std::ostream& stdstrm,
                                                    std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tJacobi2x2PreconditionerFactory_tpetra";

  status = test_createPrec(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"createPrec\" ... PASSED", "   \"createPrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec\" ... PASSED",
                       "   \"initializePrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_uninitializePrec(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"uninitializePrec\" ... PASSED",
                       "   \"uninitializePrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_isCompatable(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"isCompatable\" ... PASSED",
                       "   \"isCompatable\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_identity(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"identity\" ... PASSED", "   \"identity\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_diagonal(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"diagonal\" ... PASSED", "   \"diagonal\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_result(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"result\" ... PASSED", "   \"result\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializeFromParameterList(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializeFromParameterList\" ... PASSED",
                       "   \"initializeFromParameterList\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tJacobi2x2PreconditionedFactory...PASSED",
                         "tJacobi2x2PreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tJacobi2x2PreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tJacobi2x2PreconditionerFactory_tpetra::test_createPrec(int verbosity, std::ostream& os) {
  RCP<JacobiPreconditionerFactory> fact = rcp(new JacobiPreconditionerFactory(invF_, invC_));

  try {
    // preconditioner factory should return a DefaultPreconditionerBase
    rcp_dynamic_cast<Thyra::DefaultPreconditioner<ST> >(fact->createPrec(), true);
  } catch (std::exception& e) {
    // if the dynamic cast fails...so does the test
    os << std::endl
       << "   test_createPrec: dynamic cast to \"DefaultPreconditioner\" FAILED" << std::endl;
    os << "   Descriptive exception \"" << e.what() << "\"" << std::endl;

    return false;
  }

  return true;
}

bool tJacobi2x2PreconditionerFactory_tpetra::test_initializePrec(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  // Build block2x2 preconditioner
  RCP<Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new JacobiPreconditionerFactory(invF_, invC_));
  RCP<Thyra::PreconditionerBase<ST> > prec = precFactory->createPrec();

  // initialize the preconditioner
  precFactory->initializePrec(Thyra::defaultLinearOpSource(A_), &*prec);

  RCP<const Thyra::LinearOpBase<ST> > op;

  op     = prec->getUnspecifiedPrecOp();
  status = (op != Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getUnspecifiedPrecOp\" is null (it should not be!)" << std::endl;
    ;
  }
  allPassed &= status;

  op     = prec->getRightPrecOp();
  status = (op == Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getRightPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  op     = prec->getLeftPrecOp();
  status = (op == Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getLeftPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  return allPassed;
}

bool tJacobi2x2PreconditionerFactory_tpetra::test_uninitializePrec(int verbosity,
                                                                   std::ostream& os) {
  return true;
}

bool tJacobi2x2PreconditionerFactory_tpetra::test_isCompatable(int verbosity, std::ostream& os) {
  return true;
}

bool tJacobi2x2PreconditionerFactory_tpetra::test_identity(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<ST> > LinearOp;

  bool status    = false;
  bool allPassed = true;
  ST diff        = 0.0;

  LinearOp F = Thyra::identity<ST>(invF_->range());
  LinearOp G = Thyra::identity<ST>(invF_->range());
  LinearOp D = Thyra::identity<ST>(invC_->range());
  LinearOp C = Thyra::identity<ST>(invC_->range());

  LinearOp A = Thyra::block2x2(F, G, D, C);
  RCP<Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new JacobiPreconditionerFactory(F, C));
  RCP<Thyra::PreconditionerBase<ST> > prec = Thyra::prec<ST>(*precFactory, A);

  // build linear operator
  RCP<const Thyra::LinearOpBase<ST> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));
  // construct a couple of vectors
  Tpetra::Vector<ST, LO, GO, NT> ea(map), eb(map);
  const RCP<const Thyra::MultiVectorBase<ST> > x = BlockVector(ea, eb, A->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(A->range(), 1);

  // test vector [0 1 1 3]
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 1.0);
  eb.replaceGlobalValue(0, 1.0);
  eb.replaceGlobalValue(1, 3.0);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, x) / Thyra::norm_2(*x->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_identity " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  // test vector [-2 4 7 9]
  ea.replaceGlobalValue(0, -2.0);
  ea.replaceGlobalValue(1, 4.0);
  eb.replaceGlobalValue(0, 7.0);
  eb.replaceGlobalValue(1, 9.0);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, x) / Thyra::norm_2(*x->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_identity " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  // test vector [1 0 0 -5]
  ea.replaceGlobalValue(0, 1.0);
  ea.replaceGlobalValue(1, 0.0);
  eb.replaceGlobalValue(0, 0.0);
  eb.replaceGlobalValue(1, -5.0);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, x) / Thyra::norm_2(*x->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_identity " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  // test vector [4 -4 6 12]
  ea.replaceGlobalValue(0, 4.0);
  ea.replaceGlobalValue(1, -4.0);
  eb.replaceGlobalValue(0, 6.0);
  eb.replaceGlobalValue(1, 12.0);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, x) / Thyra::norm_2(*x->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_identity " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  return allPassed;
}

bool tJacobi2x2PreconditionerFactory_tpetra::test_diagonal(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<ST> > LinearOp;

  bool status    = false;
  bool allPassed = true;
  ST vec[2];
  ST diff = 0.0;

  // build 4x4 matrix with block 2x2 diagonal subblocks
  //
  //            [ 1 0 7 0 ]
  // [ F G ] =  [ 0 2 0 8 ]
  // [ D C ]    [ 5 0 3 0 ]
  //            [ 0 6 0 4 ]
  //

  vec[0]     = 1.0;
  vec[1]     = 2.0;
  LinearOp F = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 7.0;
  vec[1]     = 8.0;
  LinearOp G = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 5.0;
  vec[1]     = 6.0;
  LinearOp D = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]     = 0.0;
  vec[1]     = 0.0;
  LinearOp C = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]      = 1.0;
  vec[1]      = 0.5;
  LinearOp iF = Teko::Test::DiagMatrix_tpetra(2, vec);

  vec[0]      = 1.0 / 3.0;
  vec[1]      = 0.25;
  LinearOp iC = Teko::Test::DiagMatrix_tpetra(2, vec);

  LinearOp A = Thyra::block2x2(F, G, D, C);
  RCP<Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new JacobiPreconditionerFactory(iF, iC));
  RCP<Thyra::PreconditionerBase<ST> > prec = Thyra::prec<ST>(*precFactory, A);

  // build linear operator
  RCP<const Thyra::LinearOpBase<ST> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));
  // construct a couple of vectors
  Tpetra::Vector<ST, LO, GO, NT> ea(map), eb(map);
  Tpetra::Vector<ST, LO, GO, NT> ef(map), eg(map);
  const RCP<const Thyra::MultiVectorBase<ST> > x = BlockVector(ea, eb, A->domain());
  const RCP<const Thyra::MultiVectorBase<ST> > z = BlockVector(ef, eg, A->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(A->range(), 1);

  // now checks of the preconditioner (should be exact!)
  /////////////////////////////////////////////////////////////////////////

  // test vector [0 1 1 3]
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 1.0);
  eb.replaceGlobalValue(0, 1.0);
  eb.replaceGlobalValue(1, 3.0);
  ef.replaceGlobalValue(0, 0.0);
  ef.replaceGlobalValue(1, 0.5);
  eg.replaceGlobalValue(0, 1.0 / 3.0);
  eg.replaceGlobalValue(1, 0.75);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [-2 4 7 9]
  ea.replaceGlobalValue(0, -2.0);
  ea.replaceGlobalValue(1, 4.0);
  eb.replaceGlobalValue(0, 7.0);
  eb.replaceGlobalValue(1, 9.0);
  ef.replaceGlobalValue(0, -2.000000000000000);
  ef.replaceGlobalValue(1, 2.000000000000000);
  eg.replaceGlobalValue(0, 2.333333333333333);
  eg.replaceGlobalValue(1, 2.250000000000000);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [1 0 0 -5]
  ea.replaceGlobalValue(0, 1.0);
  ea.replaceGlobalValue(1, 0.0);
  eb.replaceGlobalValue(0, 0.0);
  eb.replaceGlobalValue(1, -5.0);
  ef.replaceGlobalValue(0, 1.000000000000000);
  ef.replaceGlobalValue(1, 0.000000000000000);
  eg.replaceGlobalValue(0, 0.000000000000000);
  eg.replaceGlobalValue(1, -1.250000000000000);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [4 -4 6 12]
  ea.replaceGlobalValue(0, 4.0);
  ea.replaceGlobalValue(1, -4.0);
  eb.replaceGlobalValue(0, 6.0);
  eb.replaceGlobalValue(1, 12.0);
  ef.replaceGlobalValue(0, 4.000000000000000);
  ef.replaceGlobalValue(1, -2.000000000000000);
  eg.replaceGlobalValue(0, 2.000000000000000);
  eg.replaceGlobalValue(1, 3.000000000000000);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  return allPassed;
}

bool tJacobi2x2PreconditionerFactory_tpetra::test_result(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  ST diff;

  // Build block2x2 preconditioner
  RCP<Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new JacobiPreconditionerFactory(invF_, invC_));
  RCP<Thyra::PreconditionerBase<ST> > prec = Thyra::prec<ST>(*precFactory, A_);

  // build linear operator
  RCP<const Thyra::LinearOpBase<ST> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));
  // construct a couple of vectors
  Tpetra::Vector<ST, LO, GO, NT> ea(map), eb(map);
  Tpetra::Vector<ST, LO, GO, NT> ef(map), eg(map);

  const RCP<const Thyra::MultiVectorBase<ST> > x = BlockVector(ea, eb, A_->domain());
  const RCP<const Thyra::MultiVectorBase<ST> > z = BlockVector(ef, eg, A_->domain());
  const RCP<Thyra::MultiVectorBase<ST> > y       = Thyra::createMembers(A_->range(), 1);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());

  // now checks of the preconditioner (should be exact!)
  /////////////////////////////////////////////////////////////////////////

  // test vector [0 1 1 3]
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 1.0);
  eb.replaceGlobalValue(0, 1.0);
  eb.replaceGlobalValue(1, 3.0);
  ef.replaceGlobalValue(0, 6.6666666666666663e-01);
  ef.replaceGlobalValue(1, -3.3333333333333331e-01);
  eg.replaceGlobalValue(0, -3.0303030303030304e-02);
  eg.replaceGlobalValue(1, 6.3636363636363635e-01);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [-2 4 7 9]
  ea.replaceGlobalValue(0, -2.0);
  ea.replaceGlobalValue(1, 4.0);
  eb.replaceGlobalValue(0, 7.0);
  eb.replaceGlobalValue(1, 9.0);
  ef.replaceGlobalValue(0, 3.3333333333333330e+00);
  ef.replaceGlobalValue(1, -2.6666666666666665e+00);
  eg.replaceGlobalValue(0, 5.1515151515151514e-01);
  eg.replaceGlobalValue(1, 1.1818181818181817e+00);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [1 0 0 -5]
  ea.replaceGlobalValue(0, 1.0);
  ea.replaceGlobalValue(1, 0.0);
  eb.replaceGlobalValue(0, 0.0);
  eb.replaceGlobalValue(1, -5.0);
  ef.replaceGlobalValue(0, -3.3333333333333331e-01);
  ef.replaceGlobalValue(1, 6.6666666666666663e-01);
  eg.replaceGlobalValue(0, 3.0303030303030298e-01);
  eg.replaceGlobalValue(1, -1.3636363636363635e+00);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [4 -4 6 12]
  ea.replaceGlobalValue(0, 4.0);
  ea.replaceGlobalValue(1, -4.0);
  eb.replaceGlobalValue(0, 6.0);
  eb.replaceGlobalValue(1, 12.0);
  ef.replaceGlobalValue(0, -4.0000000000000000e+00);
  ef.replaceGlobalValue(1, 4.0000000000000000e+00);
  eg.replaceGlobalValue(0, 1.8181818181818177e-01);
  eg.replaceGlobalValue(1, 2.1818181818181817e+00);
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory_tpetra::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  return allPassed;
}

template <typename T>
Teuchos::RCP<const T> getLowsFactory(const RCP<const InverseFactory>& invFact) {
  return rcp_dynamic_cast<const T>(
      rcp_dynamic_cast<const SolveInverseFactory>(invFact)->getLowsFactory());
}

template <typename T>
Teuchos::RCP<const T> getPrecFactory(const RCP<const InverseFactory>& invFact) {
  return rcp_dynamic_cast<const T>(
      rcp_dynamic_cast<const PreconditionerInverseFactory>(invFact)->getPrecFactory());
}

bool tJacobi2x2PreconditionerFactory_tpetra::test_initializeFromParameterList(int verbosity,
                                                                              std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
  {
    Teuchos::ParameterList p;
    p.set("Inverse Type", "Ifpack2");

    RCP<PreconditionerFactory> fact =
        PreconditionerFactory::buildPreconditionerFactory("Block Jacobi", p, invLib);
    RCP<JacobiPreconditionerFactory> jFact = rcp_dynamic_cast<JacobiPreconditionerFactory>(fact);

    // check we have the right factory
    status = (jFact != Teuchos::null);
    if (status) os << "Dynnamic cast failed" << std::endl;
    allPassed &= status;

    RCP<const InvFactoryDiagStrategy> diagStrat =
        rcp_dynamic_cast<const InvFactoryDiagStrategy>(jFact->getInvDiagStrategy());

    // check we have the right factory
    status = (diagStrat != Teuchos::null);
    if (!status) os << "Dynnamic to InvFactoryDiagStrategy cast failed" << std::endl;
    allPassed &= status;

    const std::vector<Teuchos::RCP<InverseFactory> >& facts = diagStrat->getFactories();
    status                                                  = (facts.size() == 1);
    allPassed &= status;
  }

  {
#ifdef TEKO_HAVE_EPETRA
    Teuchos::ParameterList p;
    p.set("Inverse Type", "ML");
    p.set("Inverse Type 1", "Amesos");
    p.set("Inverse Type 3", "Ifpack");

    RCP<PreconditionerFactory> fact =
        PreconditionerFactory::buildPreconditionerFactory("Block Jacobi", p, invLib);
    RCP<JacobiPreconditionerFactory> jFact = rcp_dynamic_cast<JacobiPreconditionerFactory>(fact);
    RCP<const InvFactoryDiagStrategy> diagStrat =
        rcp_dynamic_cast<const InvFactoryDiagStrategy>(jFact->getInvDiagStrategy());

    // check we have the right factory
    const std::vector<Teuchos::RCP<InverseFactory> >& facts = diagStrat->getFactories();
    status                                                  = (facts.size() == 3);
    if (!status) os << "Incorrect number of factories constructed" << std::endl;
    allPassed &= status;

    status = (getLowsFactory<Thyra::AmesosLinearOpWithSolveFactory>(facts[0]) != Teuchos::null);
    if (!status) os << "Checking if Amesos inverse factory was consctructed" << std::endl;
    allPassed &= status;

    status = (getPrecFactory<Thyra::MLPreconditionerFactory>(facts[1]) != Teuchos::null);
    if (!status) os << "Checking if ML inverse factory was consctructed" << std::endl;
    allPassed &= status;

    status = (getPrecFactory<Thyra::IfpackPreconditionerFactory>(facts[2]) != Teuchos::null);
    if (!status) os << "Checking if ML inverse factory was consctructed" << std::endl;
    allPassed &= status;
#endif
  }

  return allPassed;
}

}  // end namespace Test
}  // end namespace Teko
