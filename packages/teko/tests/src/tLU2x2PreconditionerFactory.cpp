// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tLU2x2PreconditionerFactory.hpp"
#include "Teko_LU2x2PreconditionerFactory.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"

#include <vector>

// This whole test rig is based on inverting the matrix
//
//      [  1  2  1 -1 ]
//  A = [  2  1 -3  1 ]
//      [  1 -3  1  2 ]
//      [ -1  1  2  1 ]
//
// see the matlab file

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;

void tLU2x2PreconditionerFactory::initializeTest() {
  std::vector<int> indicies(2);
  std::vector<double> row0(2), row1(2);

  tolerance_ = 9.0e-15;

  comm                      = rcp(new Epetra_SerialComm());
  const RCP<Epetra_Map> map = rcp(new Epetra_Map(2, 0, *comm));

  const RCP<Epetra_CrsMatrix> ptrF  = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrB  = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrBt = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));

  const RCP<Epetra_CrsMatrix> ptrInvF = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrInvS = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));

  indicies[0] = 0;
  indicies[1] = 1;

  // build F matrix
  row0[0] = 1.0;
  row0[1] = 2.0;
  row1[0] = 2.0;
  row1[1] = 1.0;
  ptrF->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrF->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrF->FillComplete();
  F_ = Thyra::epetraLinearOp(ptrF, "ptrF");

  // build B matrix
  row0[0] = 1.0;
  row0[1] = -3.0;
  row1[0] = -1.0;
  row1[1] = 1.0;
  ptrB->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrB->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrB->FillComplete();
  B_ = Thyra::epetraLinearOp(ptrB, "ptrB");

  // build Bt matrix
  row0[0] = 1.0;
  row0[1] = -1.0;
  row1[0] = -3.0;
  row1[1] = 1.0;
  ptrBt->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrBt->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrBt->FillComplete();
  Bt_ = Thyra::epetraLinearOp(ptrBt, "ptrBt");

  // build inv(F) matrix
  row0[0] = -1.0 / 3.0;
  row0[1] = 2.0 / 3.0;
  row1[0] = 2.0 / 3.0;
  row1[1] = -1.0 / 3.0;
  ptrInvF->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrInvF->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrInvF->FillComplete();
  invF_ = Thyra::epetraLinearOp(ptrInvF, "ptrInvF");

  // build inv(Pschur) matrix
  row0[0] = 0.1428571428571428;
  row0[1] = 0.0952380952380952;
  row1[0] = 0.0952380952380952;
  row1[1] = 0.3968253968253968;
  ptrInvS->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrInvS->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrInvS->FillComplete();
  invS_ = Thyra::scale<double>(-1.0, Thyra::epetraLinearOp(ptrInvS, "ptrInvS"));

  A_ = Thyra::block2x2<double>(F_, Bt_, B_, F_, "A");
}

int tLU2x2PreconditionerFactory::runTest(int verbosity, std::ostream& stdstrm,
                                         std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tLU2x2PreconditionerFactory";

  status = test_createPrec(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"createPrec\" ... PASSED", "   \"createPrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"initializePrec\" ... PASSED", "   \"initializePrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_uninitializePrec(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"uninitializePrec\" ... PASSED",
                "   \"uninitializePrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_isCompatable(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"isCompatable\" ... PASSED", "   \"isCompatable\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_identity(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"identity\" ... PASSED", "   \"identity\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_diagonal(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"diagonal\" ... PASSED", "   \"diagonal\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_result(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"result\" ... PASSED", "   \"result\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_alphabeta(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"alphabeta\" ... PASSED", "   \"alphabeta\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tLU2x2PreconditionedFactory...PASSED",
                  "tLU2x2PreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tLU2x2PreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tLU2x2PreconditionerFactory::test_createPrec(int verbosity, std::ostream& os) {
  RCP<LU2x2PreconditionerFactory> fact = rcp(new LU2x2PreconditionerFactory(invF_, invS_));

  try {
    // preconditioner factory should return a DefaultPreconditionerBase
    rcp_dynamic_cast<DefaultPreconditioner<double> >(fact->createPrec(), true);
  } catch (std::exception& e) {
    // if the dynamic cast fails...so does the test
    os << std::endl
       << "   test_createPrec: dynamic cast to \"DefaultPreconditioner\" FAILED" << std::endl;
    os << "   Descriptive exception \"" << e.what() << "\"" << std::endl;

    return false;
  }

  return true;
}

bool tLU2x2PreconditionerFactory::test_initializePrec(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  // Build block2x2 preconditioner
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new LU2x2PreconditionerFactory(invF_, invS_));
  RCP<Thyra::PreconditionerBase<double> > prec = precFactory->createPrec();

  // initialize the preconditioner
  precFactory->initializePrec(Thyra::defaultLinearOpSource(A_), &*prec);

  RCP<const Thyra::LinearOpBase<double> > op;

  op     = prec->getUnspecifiedPrecOp();
  status = (op != Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
    os << "      "
       << "Preconditioner \"getUnspecifiedPrecOp\" is null (it should not be!)" << std::endl;
    ;
  }
  allPassed &= status;

  op     = prec->getRightPrecOp();
  status = (op == Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
    os << "      "
       << "Preconditioner \"getRightPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  op     = prec->getLeftPrecOp();
  status = (op == Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
    os << "      "
       << "Preconditioner \"getLeftPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  return allPassed;
}

bool tLU2x2PreconditionerFactory::test_uninitializePrec(int verbosity, std::ostream& os) {
  return true;
}

bool tLU2x2PreconditionerFactory::test_isCompatable(int verbosity, std::ostream& os) {
  return true;
}

bool tLU2x2PreconditionerFactory::test_identity(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

  bool status    = false;
  bool allPassed = true;

  LinearOp Iu   = Thyra::identity<double>(invF_->range());
  LinearOp Ip   = Thyra::identity<double>(invS_->range());
  LinearOp Zu   = Thyra::zero<double>(invF_->range(), invS_->domain());
  LinearOp Zp   = Thyra::zero<double>(invS_->range(), invF_->domain());
  LinearOp invS = Thyra::scale(-1.0, Ip);

  LinearOp A = Thyra::block2x2(Iu, Zp, Zu, Ip);
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new LU2x2PreconditionerFactory(Iu, invS));
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A);

  // build linear operator
  RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Epetra_Map> map = rcp(new Epetra_Map(2, 0, *comm));
  // construct a couple of vectors
  Epetra_Vector ea(*map), eb(*map);
  const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea, eb, A->domain());
  const RCP<Thyra::MultiVectorBase<double> > y       = Thyra::createMembers(A->range(), 1);

  // test vector [0 1 1 3]
  ea[0] = 0.0;
  ea[1] = 1.0;
  eb[0] = 1.0;
  eb[1] = 3.0;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(x, y) < tolerance_);
  if (not status) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_Identity " << toString(status) << ": A*x != y"
       << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  // test vector [-2 4 7 9]
  ea[0] = -2.0;
  ea[1] = 4.0;
  eb[0] = 7.0;
  eb[1] = 9.0;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(x, y) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_Identity " << toString(status) << ": A*x != y"
       << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  // test vector [1 0 0 -5]
  ea[0] = 1.0;
  ea[1] = 0.0;
  eb[0] = 0.0;
  eb[1] = -5.0;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(x, y) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_Identity " << toString(status) << ": A*x != y"
       << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  // test vector [4 -4 6 12]
  ea[0] = 4.0;
  ea[1] = -4.0;
  eb[0] = 6.0;
  eb[1] = 12.0;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(x, y) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_Identity " << toString(status) << ": A*x != y"
       << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  return allPassed;
}

bool tLU2x2PreconditionerFactory::test_diagonal(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

  bool status    = false;
  bool allPassed = true;
  double vec[2];

  // build 4x4 matrix with block 2x2 diagonal subblocks
  //
  //            [ 1 0 7 0 ]
  // [ F G ] =  [ 0 2 0 8 ]
  // [ D C ]    [ 5 0 3 0 ]
  //            [ 0 6 0 4 ]
  //

  vec[0]     = 1.0;
  vec[1]     = 2.0;
  LinearOp F = Teko::Test::DiagMatrix(2, vec);

  vec[0]     = 7.0;
  vec[1]     = 8.0;
  LinearOp G = Teko::Test::DiagMatrix(2, vec);

  vec[0]     = 5.0;
  vec[1]     = 6.0;
  LinearOp D = Teko::Test::DiagMatrix(2, vec);

  vec[0]     = 3.0;
  vec[1]     = 4.0;
  LinearOp C = Teko::Test::DiagMatrix(2, vec);

  vec[0]      = 1.0;
  vec[1]      = 0.5;
  LinearOp iF = Teko::Test::DiagMatrix(2, vec);

  // S = -C+D*iF*G
  vec[0]      = 0.03125;
  vec[1]      = 0.05;
  LinearOp iS = Teko::Test::DiagMatrix(2, vec);

  LinearOp A = Thyra::block2x2(F, G, D, C);
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new LU2x2PreconditionerFactory(iF, iS));
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A);

  // build linear operator
  RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Epetra_Map> map = rcp(new Epetra_Map(2, 0, *comm));
  // construct a couple of vectors
  Epetra_Vector ea(*map), eb(*map);
  Epetra_Vector ef(*map), eg(*map);
  const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea, eb, A->domain());
  const RCP<const Thyra::MultiVectorBase<double> > z = BlockVector(ef, eg, A->domain());
  const RCP<Thyra::MultiVectorBase<double> > y       = Thyra::createMembers(A->range(), 1);

  // first some sanity checks of the forward operator
  /////////////////////////////////////////////////////////////////////////

  double diff;

  // test vector [0 1 1 3]
  ea[0] = 0.0;
  ea[1] = 1.0;
  eb[0] = 1.0;
  eb[1] = 3.0;
  ef[0] = 7.0;
  ef[1] = 26.0;
  eg[0] = 3.0;
  eg[1] = 18.0;
  Thyra::apply(*A, Thyra::NOTRANS, *x, y.ptr());
  diff   = Teko::Test::Difference(y, z);
  status = (diff < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_diagonal " << toString(status) << ":  A*y != z"
       << std::endl;
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [-2 4 7 9]
  ea[0] = -2.0;
  ea[1] = 4.0;
  eb[0] = 7.0;
  eb[1] = 9.0;
  ef[0] = 47.0;
  ef[1] = 80.0;
  eg[0] = 11.0;
  eg[1] = 60.0;
  Thyra::apply(*A, Thyra::NOTRANS, *x, y.ptr());
  diff   = Teko::Test::Difference(y, z);
  status = (diff < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_diagonal " << toString(status) << ":  A*y != z"
       << std::endl;
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [1 0 0 -5]
  ea[0] = 1.0;
  ea[1] = 0.0;
  eb[0] = 0.0;
  eb[1] = -5.0;
  ef[0] = 1.0;
  ef[1] = -40.0;
  eg[0] = 5.0;
  eg[1] = -20.0;
  Thyra::apply(*A, Thyra::NOTRANS, *x, y.ptr());
  diff   = Teko::Test::Difference(y, z);
  status = (diff < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_diagonal " << toString(status) << ":  A*y != z"
       << std::endl;
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [4 -4 6 12]
  ea[0] = 4.0;
  ea[1] = -4.0;
  eb[0] = 6.0;
  eb[1] = 12.0;
  ef[0] = 46.0;
  ef[1] = 88.0;
  eg[0] = 38.0;
  eg[1] = 24.0;
  Thyra::apply(*A, Thyra::NOTRANS, *x, y.ptr());
  diff   = Teko::Test::Difference(y, z);
  status = (diff < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_diagonal " << toString(status) << ":  A*y != z"
       << std::endl;
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // now checks of the preconditioner (should be exact!)
  /////////////////////////////////////////////////////////////////////////

  // test vector [0 1 1 3]
  ea[0] = 0.0;
  ea[1] = 1.0;
  eb[0] = 1.0;
  eb[1] = 3.0;
  ef[0] = 0.21875;
  ef[1] = 0.5;
  eg[0] = -0.03125;
  eg[1] = 0.0;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(y, z) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [-2 4 7 9]
  ea[0] = -2.0;
  ea[1] = 4.0;
  eb[0] = 7.0;
  eb[1] = 9.0;
  ef[0] = 1.71875;
  ef[1] = 1.4;
  eg[0] = -0.53125;
  eg[1] = 0.15;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(y, z) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [1 0 0 -5]
  ea[0] = 1.0;
  ea[1] = 0.0;
  eb[0] = 0.0;
  eb[1] = -5.0;
  ef[0] = -0.09375;
  ef[1] = -1.0;
  eg[0] = 0.15625;
  eg[1] = 0.25;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(y, z) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [4 -4 6 12]
  ea[0] = 4.0;
  ea[1] = -4.0;
  eb[0] = 6.0;
  eb[1] = 12.0;
  ef[0] = 0.9375;
  ef[1] = 2.8;
  eg[0] = 0.4375;
  eg[1] = -1.2;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = (Teko::Test::Difference(y, z) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_diagonal " << toString(status)
       << ":  (y=inv(A)*x) != z" << std::endl;
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

bool tLU2x2PreconditionerFactory::test_result(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  double diff;

  // Build block2x2 preconditioner
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new LU2x2PreconditionerFactory(invF_, invS_));
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A_);

  // build linear operator
  RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Epetra_Map> map = rcp(new Epetra_Map(2, 0, *comm));
  // construct a couple of vectors
  Epetra_Vector ea(*map), eb(*map);
  Epetra_Vector ef(*map), eg(*map);

  const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea, eb, A_->domain());
  const RCP<const Thyra::MultiVectorBase<double> > z = BlockVector(ef, eg, A_->domain());
  const RCP<Thyra::MultiVectorBase<double> > y       = Thyra::createMembers(A_->range(), 1);

  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());

  // now checks of the preconditioner (should be exact!)
  /////////////////////////////////////////////////////////////////////////

  // test vector [0 1 1 3]
  ea[0] = 0.0;
  ea[1] = 1.0;
  eb[0] = 1.0;
  eb[1] = 3.0;
  ef[0] = -0.190476190476190;
  ef[1] = 0.714285714285714;
  eg[0] = 0.285714285714286;
  eg[1] = 1.523809523809524;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z)) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [-2 4 7 9]
  ea[0] = -2.0;
  ea[1] = 4.0;
  eb[0] = 7.0;
  eb[1] = 9.0;
  ef[0] = -0.317460317460317;
  ef[1] = 1.523809523809524;
  eg[0] = 0.809523809523810;
  eg[1] = 5.539682539682540;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z)) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [1 0 0 -5]
  ea[0] = 1.0;
  ea[1] = 0.0;
  eb[0] = 0.0;
  eb[1] = -5.0;
  ef[0] = 1.269841269841270;
  ef[1] = -1.095238095238095;
  eg[0] = -0.238095238095238;
  eg[1] = -2.158730158730159;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z)) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [4 -4 6 12]
  ea[0] = 4.0;
  ea[1] = -4.0;
  eb[0] = 6.0;
  eb[1] = 12.0;
  ef[0] = 0.539682539682540;
  ef[1] = 1.809523809523809;
  eg[0] = 3.523809523809524;
  eg[1] = 3.682539682539683;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z)) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_result " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2 = " << diff << ")" << std::endl;
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

bool tLU2x2PreconditionerFactory::test_alphabeta(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  double diff;

  // Build block2x2 preconditioner
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new LU2x2PreconditionerFactory(invF_, invS_));
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A_);

  // build linear operator
  RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

  const RCP<Epetra_Map> map = rcp(new Epetra_Map(2, 0, *comm));
  // construct a couple of vectors
  Epetra_Vector ea(*map), eb(*map);
  Epetra_Vector ef(*map), eg(*map);
  Epetra_Vector ec(*map), ed(*map);

  const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea, eb, A_->domain());
  const RCP<const Thyra::MultiVectorBase<double> > z = BlockVector(ef, eg, A_->domain());
  const RCP<const Thyra::MultiVectorBase<double> > q = BlockVector(ec, ed, A_->domain());
  RCP<Thyra::MultiVectorBase<double> > y             = Thyra::createMembers(q->range(), 1);

  // now checks of the preconditioner (should be exact!)
  /////////////////////////////////////////////////////////////////////////

  double alpha = 0.5;
  double beta  = -1.9;
  ec[0]        = 6.2;
  ec[1]        = -2.4;
  ed[0]        = 9.7;
  ed[1]        = 0.04;  // set q

  // test vector [0 1 1 3]
  Thyra::assign<double>(y.ptr(), *q);
  ea[0] = 0.0;
  ea[1] = 1.0;
  eb[0] = 1.0;
  eb[1] = 3.0;
  ef[0] = -11.875238095238094;
  ef[1] = 4.917142857142856;
  eg[0] = -18.287142857142854;
  eg[1] = 0.685904761904762;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr(), alpha, beta);
  status = ((diff = Teko::Test::Difference(y, z)) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_alphabeta " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [-2 4 7 9]
  Thyra::assign<double>(y.ptr(), *q);
  ea[0] = -2.0;
  ea[1] = 4.0;
  eb[0] = 7.0;
  eb[1] = 9.0;
  ef[0] = -11.938730158730158;
  ef[1] = 5.321904761904761;
  eg[0] = -18.025238095238091;
  eg[1] = 2.693841269841270;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr(), alpha, beta);
  status = ((diff = Teko::Test::Difference(y, z)) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_alphabeta " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [1 0 0 -5]
  Thyra::assign<double>(y.ptr(), *q);
  ea[0] = 1.0;
  ea[1] = 0.0;
  eb[0] = 0.0;
  eb[1] = -5.0;
  ef[0] = -11.145079365079365;
  ef[1] = 4.012380952380952;
  eg[0] = -18.549047619047617;
  eg[1] = -1.155365079365079;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr(), alpha, beta);
  status = ((diff = Teko::Test::Difference(y, z)) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_alphabeta " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
    os << "      ";
    Print(os, "z", z);
  }
  allPassed &= status;

  // test vector [4 -4 6 12]
  Thyra::assign<double>(y.ptr(), *q);
  ea[0] = 4.0;
  ea[1] = -4.0;
  eb[0] = 6.0;
  eb[1] = 12.0;
  ef[0] = -11.510158730158729;
  ef[1] = 5.464761904761904;
  eg[0] = -16.668095238095233;
  eg[1] = 1.765269841269841;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr(), alpha, beta);
  status = ((diff = Teko::Test::Difference(y, z)) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tLU2x2PreconditionerFactory::test_alphabeta " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2 = " << diff << ")" << std::endl;
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
