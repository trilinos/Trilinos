// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tJacobi2x2PreconditionerFactory.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"

#include "Teko_InverseLibrary.hpp"
#include "Teko_PreconditionerInverseFactory.hpp"
#include "Teko_SolveInverseFactory.hpp"

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
#include "Thyra_MLPreconditionerFactory.hpp"
#include "Thyra_IfpackPreconditionerFactory.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"

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

void tJacobi2x2PreconditionerFactory::initializeTest() {
  std::vector<int> indicies(2);
  std::vector<double> row0(2), row1(2);

  tolerance_ = 1.0e-14;

  comm                      = rcp(new Epetra_SerialComm());
  const RCP<Epetra_Map> map = rcp(new Epetra_Map(2, 0, *comm));

  const RCP<Epetra_CrsMatrix> ptrF = rcp(new Epetra_CrsMatrix(Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrG = rcp(new Epetra_CrsMatrix(Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrD = rcp(new Epetra_CrsMatrix(Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrC = rcp(new Epetra_CrsMatrix(Copy, *map, 2));

  const RCP<Epetra_CrsMatrix> ptrInvF = rcp(new Epetra_CrsMatrix(Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrInvC = rcp(new Epetra_CrsMatrix(Copy, *map, 2));

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

  // build D matrix
  row0[0] = 1.0;
  row0[1] = -3.0;
  row1[0] = -1.0;
  row1[1] = 1.0;
  ptrD->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrD->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrD->FillComplete();
  D_ = Thyra::epetraLinearOp(ptrD, "ptrD");

  // build G matrix
  row0[0] = 1.0;
  row0[1] = -1.0;
  row1[0] = -3.0;
  row1[1] = 1.0;
  ptrG->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrG->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrG->FillComplete();
  G_ = Thyra::epetraLinearOp(ptrG, "ptrG");

  // build C matrix
  row0[0] = 9.0;
  row0[1] = 2.0;
  row1[0] = 6.0;
  row1[1] = 5.0;
  ptrC->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrC->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrC->FillComplete();
  C_ = Thyra::epetraLinearOp(ptrC, "ptrC");

  // build inv(F) matrix
  row0[0] = -1.0 / 3.0;
  row0[1] = 2.0 / 3.0;
  row1[0] = 2.0 / 3.0;
  row1[1] = -1.0 / 3.0;
  ptrInvF->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrInvF->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrInvF->FillComplete();
  invF_ = Thyra::epetraLinearOp(ptrInvF, "ptrInvF");

  // build inv(C) matrix
  row0[0] = 0.151515151515151;
  row0[1] = -0.060606060606061;
  row1[0] = -0.181818181818182;
  row1[1] = 0.272727272727273;
  ptrInvC->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrInvC->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrInvC->FillComplete();
  invC_ = Thyra::epetraLinearOp(ptrInvC, "ptrInvC");

  A_ = Thyra::block2x2<double>(F_, G_, D_, C_, "A");
}

int tJacobi2x2PreconditionerFactory::runTest(int verbosity, std::ostream& stdstrm,
                                             std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tJacobi2x2PreconditionerFactory";

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

  status = test_initializeFromParameterList(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"initializeFromParameterList\" ... PASSED",
                "   \"initializeFromParameterList\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tJacobi2x2PreconditionedFactory...PASSED",
                  "tJacobi2x2PreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tJacobi2x2PreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tJacobi2x2PreconditionerFactory::test_createPrec(int verbosity, std::ostream& os) {
  RCP<JacobiPreconditionerFactory> fact = rcp(new JacobiPreconditionerFactory(invF_, invC_));

  try {
    // preconditioner factory should return a DefaultPreconditionerBase
    rcp_dynamic_cast<Thyra::DefaultPreconditioner<double> >(fact->createPrec(), true);
  } catch (std::exception& e) {
    // if the dynamic cast fails...so does the test
    os << std::endl
       << "   test_createPrec: dynamic cast to \"DefaultPreconditioner\" FAILED" << std::endl;
    os << "   Descriptive exception \"" << e.what() << "\"" << std::endl;

    return false;
  }

  return true;
}

bool tJacobi2x2PreconditionerFactory::test_initializePrec(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  // Build block2x2 preconditioner
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new JacobiPreconditionerFactory(invF_, invC_));
  RCP<Thyra::PreconditionerBase<double> > prec = precFactory->createPrec();

  // initialize the preconditioner
  precFactory->initializePrec(Thyra::defaultLinearOpSource(A_), &*prec);

  RCP<const Thyra::LinearOpBase<double> > op;

  op     = prec->getUnspecifiedPrecOp();
  status = (op != Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_initializePrec " << toString(status)
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
       << "   tJacobi2x2PreconditionerFactory::test_initializePrec " << toString(status)
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
       << "   tJacobi2x2PreconditionerFactory::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getLeftPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  return allPassed;
}

bool tJacobi2x2PreconditionerFactory::test_uninitializePrec(int verbosity, std::ostream& os) {
  return true;
}

bool tJacobi2x2PreconditionerFactory::test_isCompatable(int verbosity, std::ostream& os) {
  return true;
}

bool tJacobi2x2PreconditionerFactory::test_identity(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

  bool status    = false;
  bool allPassed = true;
  double diff    = 0.0;

  LinearOp F = Thyra::identity<double>(invF_->range());
  LinearOp G = Thyra::identity<double>(invF_->range());
  LinearOp D = Thyra::identity<double>(invC_->range());
  LinearOp C = Thyra::identity<double>(invC_->range());

  LinearOp A = Thyra::block2x2(F, G, D, C);
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new JacobiPreconditionerFactory(F, C));
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
  status = ((diff = Teko::Test::Difference(y, x) / Thyra::norm_2(*x->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_identity " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
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
  status = ((diff = Teko::Test::Difference(y, x) / Thyra::norm_2(*x->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_identity " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
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
  status = ((diff = Teko::Test::Difference(y, x) / Thyra::norm_2(*x->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_identity " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
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
  status = ((diff = Teko::Test::Difference(y, x) / Thyra::norm_2(*x->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_identity " << toString(status)
       << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = " << diff << ")" << std::endl;
    os << "      ";
    Print(os, "x", x);
    os << "      ";
    Print(os, "y", y);
  }
  allPassed &= status;

  return allPassed;
}

bool tJacobi2x2PreconditionerFactory::test_diagonal(int verbosity, std::ostream& os) {
  // make sure the preconditioner is working by testing against the identity matrix
  typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

  bool status    = false;
  bool allPassed = true;
  double vec[2];
  double diff = 0.0;

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

  vec[0]     = 0.0;
  vec[1]     = 0.0;
  LinearOp C = Teko::Test::DiagMatrix(2, vec);

  vec[0]      = 1.0;
  vec[1]      = 0.5;
  LinearOp iF = Teko::Test::DiagMatrix(2, vec);

  vec[0]      = 1.0 / 3.0;
  vec[1]      = 0.25;
  LinearOp iC = Teko::Test::DiagMatrix(2, vec);

  LinearOp A = Thyra::block2x2(F, G, D, C);
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new JacobiPreconditionerFactory(iF, iC));
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

  // now checks of the preconditioner (should be exact!)
  /////////////////////////////////////////////////////////////////////////

  // test vector [0 1 1 3]
  ea[0] = 0.0;
  ea[1] = 1.0;
  eb[0] = 1.0;
  eb[1] = 3.0;
  ef[0] = 0.0;
  ef[1] = 0.5;
  eg[0] = 1.0 / 3.0;
  eg[1] = 0.75;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_diagonal " << toString(status)
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
  ea[0] = -2.0;
  ea[1] = 4.0;
  eb[0] = 7.0;
  eb[1] = 9.0;
  ef[0] = -2.000000000000000;
  ef[1] = 2.000000000000000;
  eg[0] = 2.333333333333333;
  eg[1] = 2.250000000000000;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_diagonal " << toString(status)
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
  ea[0] = 1.0;
  ea[1] = 0.0;
  eb[0] = 0.0;
  eb[1] = -5.0;
  ef[0] = 1.000000000000000;
  ef[1] = 0.000000000000000;
  eg[0] = 0.000000000000000;
  eg[1] = -1.250000000000000;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_diagonal " << toString(status)
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
  ea[0] = 4.0;
  ea[1] = -4.0;
  eb[0] = 6.0;
  eb[1] = 12.0;
  ef[0] = 4.000000000000000;
  ef[1] = -2.000000000000000;
  eg[0] = 2.000000000000000;
  eg[1] = 3.000000000000000;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_diagonal " << toString(status)
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

bool tJacobi2x2PreconditionerFactory::test_result(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  double diff;

  // Build block2x2 preconditioner
  RCP<Thyra::PreconditionerFactoryBase<double> > precFactory =
      rcp(new JacobiPreconditionerFactory(invF_, invC_));
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
  ef[0] = 6.6666666666666663e-01;
  ef[1] = -3.3333333333333331e-01;
  eg[0] = -3.0303030303030304e-02;
  eg[1] = 6.3636363636363635e-01;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_result " << toString(status)
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
  ea[0] = -2.0;
  ea[1] = 4.0;
  eb[0] = 7.0;
  eb[1] = 9.0;
  ef[0] = 3.3333333333333330e+00;
  ef[1] = -2.6666666666666665e+00;
  eg[0] = 5.1515151515151514e-01;
  eg[1] = 1.1818181818181817e+00;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_result " << toString(status)
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
  ea[0] = 1.0;
  ea[1] = 0.0;
  eb[0] = 0.0;
  eb[1] = -5.0;
  ef[0] = -3.3333333333333331e-01;
  ef[1] = 6.6666666666666663e-01;
  eg[0] = 3.0303030303030298e-01;
  eg[1] = -1.3636363636363635e+00;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_result " << toString(status)
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
  ea[0] = 4.0;
  ea[1] = -4.0;
  eb[0] = 6.0;
  eb[1] = 12.0;
  ef[0] = -4.0000000000000000e+00;
  ef[1] = 4.0000000000000000e+00;
  eg[0] = 1.8181818181818177e-01;
  eg[1] = 2.1818181818181817e+00;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tJacobi2x2PreconditionerFactory::test_result " << toString(status)
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

bool tJacobi2x2PreconditionerFactory::test_initializeFromParameterList(int verbosity,
                                                                       std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  RCP<InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
  {
    Teuchos::ParameterList p;
    p.set("Inverse Type", "ML");

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
  }

  return allPassed;
}

}  // end namespace Test
}  // end namespace Teko
