// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tLSCStablePreconditionerFactory_tpetra.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"

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

using namespace Teuchos;
using namespace Thyra;
using namespace Teko::NS;

void tLSCStablePreconditionerFactory_tpetra::initializeTest() {
  std::vector<GO> indices(2);
  std::vector<ST> row0(2), row1(2);

  tolerance_ = 1.0e-13;

  comm                                    = GetComm_tpetra();
  const RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));

  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrF =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrB =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrBt =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));

  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvF =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvS =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 2));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvMass =
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

  // build B matrix
  row0[0] = 1.0;
  row0[1] = -3.0;
  row1[0] = -1.0;
  row1[1] = 1.0;
  ptrB->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrB->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrB->fillComplete();
  B_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrB->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrB->getRangeMap()), ptrB);

  // build Bt matrix
  row0[0] = 1.0;
  row0[1] = -1.0;
  row1[0] = -3.0;
  row1[1] = 1.0;
  ptrBt->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrBt->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrBt->fillComplete();
  Bt_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrBt->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrBt->getRangeMap()), ptrBt);

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

  // build inv(Mass) matrix
  row0[0] = 3.0;
  row0[1] = 0.0;
  row1[0] = 0.0;
  row1[1] = 2.0;
  ptrInvMass->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrInvMass->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrInvMass->fillComplete();
  invMass_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvMass->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvMass->getRangeMap()), ptrInvMass);

  // build inv(Pschur) matrix
  row0[0] = 0.208333333333333;
  row0[1] = 0.375000000000000;
  row1[0] = 0.375000000000000;
  row1[1] = 0.875000000000000;
  ptrInvS->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrInvS->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrInvS->fillComplete();
  invBQBt_ = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvS->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvS->getRangeMap()), ptrInvS);

  A_ = Thyra::block2x2<ST>(
      F_, Bt_, B_,
      Thyra::zero<ST>(Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrBt->getRangeMap()),
                      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrB->getDomainMap())),
      "A");
}

int tLSCStablePreconditionerFactory_tpetra::runTest(int verbosity, std::ostream& stdstrm,
                                                    std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tLSCStablePreconditionerFactory_tpetra";

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

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tLSCStablePreconditionedFactory...PASSED",
                         "tLSCStablePreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tLSCStablePreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tLSCStablePreconditionerFactory_tpetra::test_createPrec(int verbosity, std::ostream& os) {
  // RCP<LSCStablePreconditionerFactory> fact = rcp(new
  // LSCStablePreconditionerFactory(invF_,invBQBt_));
  const RCP<const Thyra::PreconditionerFactoryBase<ST> > fact =
      rcp(new LSCPreconditionerFactory(invF_, invBQBt_, Teuchos::null));

  try {
    // preconditioner factory should return a DefaultPreconditionerBase
    rcp_dynamic_cast<DefaultPreconditioner<ST> >(fact->createPrec(), true);
  } catch (std::exception& e) {
    // if the dynamic cast fails...so does the test
    os << std::endl
       << "   test_createPrec: dynamic cast to \"DefaultPreconditioner\" FAILED" << std::endl;
    os << "   Descriptive exception \"" << e.what() << "\"" << std::endl;

    return false;
  }

  return true;
}

bool tLSCStablePreconditionerFactory_tpetra::test_initializePrec(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  // Build block2x2 preconditioner
  // RCP<Thyra::PreconditionerFactoryBase<double> > precFactory
  //      = rcp(new LSCStablePreconditionerFactory(invF_,invBQBt_));
  const RCP<const Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new LSCPreconditionerFactory(invF_, invBQBt_, Teuchos::null));
  RCP<Thyra::PreconditionerBase<ST> > prec = precFactory->createPrec();

  // initialize the preconditioner
  precFactory->initializePrec(Thyra::defaultLinearOpSource(A_), &*prec);

  RCP<const Thyra::LinearOpBase<ST> > op;

  op     = prec->getUnspecifiedPrecOp();
  status = (op != Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_initializePrec " << toString(status)
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
       << "   tLSCStablePreconditionerFactory_tpetra::test_initializePrec " << toString(status)
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
       << "   tLSCStablePreconditionerFactory_tpetra::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getLeftPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  return allPassed;
}

bool tLSCStablePreconditionerFactory_tpetra::test_uninitializePrec(int verbosity,
                                                                   std::ostream& os) {
  return true;
}

bool tLSCStablePreconditionerFactory_tpetra::test_isCompatable(int verbosity, std::ostream& os) {
  return true;
}

bool tLSCStablePreconditionerFactory_tpetra::test_identity(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<ST> > LinearOp;

  bool status    = false;
  bool allPassed = true;

  LinearOp Iu      = Thyra::identity<ST>(invF_->range());
  LinearOp Ip      = Thyra::identity<ST>(invBQBt_->range());
  LinearOp Zp      = Thyra::zero<ST>(invBQBt_->range(), invF_->domain());
  LinearOp invBQBt = Ip;

  LinearOp A = Thyra::block2x2(Iu, Ip, Iu, Zp);
  const RCP<const Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new LSCPreconditionerFactory(Iu, invBQBt, Teuchos::null));
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

  // test vector [0 1 1 3]
  ea.replaceGlobalValue(0, 0.0);
  ea.replaceGlobalValue(1, 1.0);
  eb.replaceGlobalValue(0, 1.0);
  eb.replaceGlobalValue(1, 3.0);

  ef.replaceGlobalValue(0, 1.0);
  ef.replaceGlobalValue(1, 4.0);
  eg.replaceGlobalValue(0, -1.0);
  eg.replaceGlobalValue(1, -3.0);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(y, z) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_Identity " << toString(status)
       << ": (y=inv(A)*x) != z" << std::endl;
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

  ef.replaceGlobalValue(0, 5.0);
  ef.replaceGlobalValue(1, 13.0);
  eg.replaceGlobalValue(0, -7.0);
  eg.replaceGlobalValue(1, -9.0);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(y, z) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_Identity " << toString(status)
       << ": (y=inv(A)*x) != z" << std::endl;
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

  ef.replaceGlobalValue(0, 1.0);
  ef.replaceGlobalValue(1, -5.0);
  eg.replaceGlobalValue(0, 0.0);
  eg.replaceGlobalValue(1, 5.0);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(y, z) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_Identity " << toString(status)
       << ": (y=inv(A)*x) != z" << std::endl;
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

  ef.replaceGlobalValue(0, 10.0);
  ef.replaceGlobalValue(1, 8.0);
  eg.replaceGlobalValue(0, -6.0);
  eg.replaceGlobalValue(1, -12.0);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(y, z) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_Identity " << toString(status)
       << ": (y=inv(A)*x) != z" << std::endl;
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

bool tLSCStablePreconditionerFactory_tpetra::test_diagonal(int verbosity, std::ostream& os) {
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
  // [ D C ]    [ 5 0 0 0 ]
  //            [ 0 6 0 0 ]
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

  // S = -C+D*iF*G
  vec[0]        = 0.028571428571429;
  vec[1]        = 0.020833333333333;
  LinearOp iBBt = Teko::Test::DiagMatrix_tpetra(2, vec);

  LinearOp A = Thyra::block2x2(F, G, D, C);
  const RCP<const Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new LSCPreconditionerFactory(iF, iBBt, Teuchos::null));
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

  ef.replaceGlobalValue(0, 0.200000000000000);
  ef.replaceGlobalValue(1, 1.000000000000000);
  eg.replaceGlobalValue(0, -0.028571428571429);
  eg.replaceGlobalValue(1, -0.125);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_diagonal " << toString(status)
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

  ef.replaceGlobalValue(0, -0.600000000000000);
  ef.replaceGlobalValue(1, 3.500000000000000);
  eg.replaceGlobalValue(0, -0.200000000000000);
  eg.replaceGlobalValue(1, -0.375000000000000);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_diagonal " << toString(status)
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
  ef.replaceGlobalValue(1, -0.833333333333333);
  eg.replaceGlobalValue(0, 0.000000000000000);
  eg.replaceGlobalValue(1, 0.208333333333333);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_diagonal " << toString(status)
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

  ef.replaceGlobalValue(0, 5.200000000000000);
  ef.replaceGlobalValue(1, 0.000000000000000);
  eg.replaceGlobalValue(0, -0.171428571428571);
  eg.replaceGlobalValue(1, -0.500000000000000);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_diagonal " << toString(status)
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

bool tLSCStablePreconditionerFactory_tpetra::test_result(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  ST diff;

  // Build block2x2 preconditioner
  // RCP<Thyra::PreconditionerFactoryBase<double> > precFactory
  //      = rcp(new LSCStablePreconditionerFactory(invF_,invBQBt_,invMass_));
  const RCP<const Thyra::PreconditionerFactoryBase<ST> > precFactory =
      rcp(new LSCPreconditionerFactory(invF_, invBQBt_, invMass_));
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

  ef.replaceGlobalValue(0, -4.333333333333333);
  ef.replaceGlobalValue(1, -2.333333333333327);
  eg.replaceGlobalValue(0, -10.499999999999991);
  eg.replaceGlobalValue(1, -19.499999999999979);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_result " << toString(status)
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

  ef.replaceGlobalValue(0, -13.666666666666668);
  ef.replaceGlobalValue(1, -10.666666666666645);
  eg.replaceGlobalValue(0, -37.499999999999964);
  eg.replaceGlobalValue(1, -70.499999999999915);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_result " << toString(status)
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

  ef.replaceGlobalValue(0, 7.166666666666667);
  ef.replaceGlobalValue(1, 3.166666666666658);
  eg.replaceGlobalValue(0, 14.999999999999986);
  eg.replaceGlobalValue(1, 27.499999999999968);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_result " << toString(status)
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

  ef.replaceGlobalValue(0, -25.000000000000000);
  ef.replaceGlobalValue(1, -4.999999999999973);
  eg.replaceGlobalValue(0, -44.999999999999964);
  eg.replaceGlobalValue(1, -83.999999999999915);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLSCStablePreconditionerFactory_tpetra::test_result " << toString(status)
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

}  // end namespace Test
}  // end namespace Teko
