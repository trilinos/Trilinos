// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_Config.h"
#include "tSIMPLEPreconditionerFactory.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_TimingsSIMPLEPreconditionerFactory.hpp"
#include "Teko_InverseLibrary.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"

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
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include <vector>

// This whole test rig is based on inverting the matrix
//
//      [  1  2  1 -1 ]
//  A = [  2  1 -3  1 ]
//      [  1 -3  1  2 ]
//      [ -1  1  1  2 ]
//
// see the matlab file

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;
using namespace Teko::NS;

void tSIMPLEPreconditionerFactory::initializeTest() {
  std::vector<int> indicies(2);
  std::vector<double> row0(2), row1(2);

  tolerance_ = 1.0e-13;

  comm = GetComm();

  const RCP<Epetra_Map> map = rcp(new Epetra_Map(2, 0, *comm));

  const RCP<Epetra_CrsMatrix> ptrF  = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrB  = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrBt = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrC  = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));

  const RCP<Epetra_CrsMatrix> ptrInvF = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrInvS = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));
  const RCP<Epetra_CrsMatrix> ptrInvMass =
      rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, *map, 2));

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

  // build block info for block diagonal (one block)
  block_starts_     = arcp(new int[2], 0, 2);
  block_gids_       = arcp(new int[2], 0, 2);
  int *block_starts = block_starts_.getRawPtr();
  int *block_gids   = block_gids_.getRawPtr();
  block_starts[0]   = 0;
  block_starts[1]   = 2;
  block_gids[0]     = ptrF->LRID(0);
  block_gids[1]     = ptrF->LRID(1);

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

  // build C matrix
  row0[0] = 1.0;
  row0[1] = 2.0;
  row1[0] = 2.0;
  row1[1] = 1.0;
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
  invF_ = rcp(new StaticOpInverseFactory(Thyra::epetraLinearOp(ptrInvF, "ptrInvF")));

  // build inv(Pschur) matrix
  row0[0] = 0.037037037037037;
  row0[1] = 0.222222222222222;
  row1[0] = 0.222222222222222;
  row1[1] = 0.333333333333333;
  ptrInvS->InsertGlobalValues(0, 2, &row0[0], &indicies[0]);
  ptrInvS->InsertGlobalValues(1, 2, &row1[0], &indicies[0]);
  ptrInvS->FillComplete();
  invS_ = rcp(new StaticOpInverseFactory(Thyra::epetraLinearOp(ptrInvS, "ptrInvS")));

  A_ = Thyra::block2x2<double>(F_, Bt_, B_, C_, "A");
}

int tSIMPLEPreconditionerFactory::runTest(int verbosity, std::ostream &stdstrm,
                                          std::ostream &failstrm, int &totalrun) {
  bool allTests = true;
  bool status   = true;
  int failcount = 0;

  failstrm << "tSIMPLEPreconditionerFactory";

  status = test_createPrec(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"createPrec\" ... PASSED", "   \"createPrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, 0);
  Teko_TEST_MSG(stdstrm, 1, "   \"initializePrec(diag)\" ... PASSED",
                "   \"initializePrec(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

#ifdef Teko_ENABLE_Isorropia
  status = test_initializePrec(verbosity, failstrm, 1);
  Teko_TEST_MSG(stdstrm, 1, "   \"initializePrec(block-1)\" ... PASSED",
                "   \"initializePrec(block-1)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, 2);
  Teko_TEST_MSG(stdstrm, 1, "   \"initializePrec(block-2)\" ... PASSED",
                "   \"initializePrec(block-2)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;
#endif

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

  status = test_diagonal(verbosity, failstrm, 0);
  Teko_TEST_MSG(stdstrm, 1, "   \"diagonal(diag)\" ... PASSED", "   \"diagonal(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

#ifdef Teko_ENABLE_Isorropia
  status = test_diagonal(verbosity, failstrm, 1);
  Teko_TEST_MSG(stdstrm, 1, "   \"diagonal(block-1)\" ... PASSED",
                "   \"diagonal(block-1)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_diagonal(verbosity, failstrm, 2);
  Teko_TEST_MSG(stdstrm, 1, "   \"diagonal(block-2)\" ... PASSED",
                "   \"diagonal(block-2)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;
#endif

  status = test_result(verbosity, failstrm, 0);
  Teko_TEST_MSG(stdstrm, 1, "   \"result(diag)\" ... PASSED", "   \"result(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_result_timings(verbosity, failstrm, 0);
  Teko_TEST_MSG(stdstrm, 1, "   \"result_timings(diag)\" ... PASSED",
                "   \"result_timings(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

#ifdef Teko_ENABLE_Isorropia
  status = test_result(verbosity, failstrm, 1);
  Teko_TEST_MSG(stdstrm, 1, "   \"result(block-1)\" ... PASSED",
                "   \"result(block-1)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_result(verbosity, failstrm, 2);
  Teko_TEST_MSG(stdstrm, 1, "   \"result(block-2)\" ... PASSED",
                "   \"result(block-2)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;
#endif

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tSIMPLEPreconditionedFactory...PASSED",
                  "tSIMPLEPreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tSIMPLEPreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tSIMPLEPreconditionerFactory::test_createPrec(int verbosity, std::ostream &os) {
  const RCP<const Thyra::PreconditionerFactoryBase<double> > fact =
      rcp(new SIMPLEPreconditionerFactory(invF_, 0.9));

  try {
    // preconditioner factory should return a DefaultPreconditionerBase
    rcp_dynamic_cast<DefaultPreconditioner<double> >(fact->createPrec(), true);
  } catch (std::exception &e) {
    // if the dynamic cast fails...so does the test
    os << std::endl
       << "   test_createPrec: dynamic cast to \"DefaultPreconditioner\" FAILED" << std::endl;
    os << "   Descriptive exception \"" << e.what() << "\"" << std::endl;

    return false;
  }

  return true;
}

bool tSIMPLEPreconditionerFactory::test_initializePrec(int verbosity, std::ostream &os,
                                                       int use_blocking) {
  bool status    = false;
  bool allPassed = true;

  // Build block2x2 preconditioner
  RCP<SIMPLEPreconditionerFactory> sFactory = rcp(new SIMPLEPreconditionerFactory(invF_, 0.9));
  const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory = sFactory;
  RCP<Thyra::PreconditionerBase<double> > prec = precFactory->createPrec();

  if (use_blocking == 1) {
    // parameter list for (1,1) block
    Teuchos::ParameterList List, BlkList;
    BlkList.set("number of local blocks", 1);
    BlkList.set("block start index", &*block_starts_);
    BlkList.set("block entry gids", &*block_gids_);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 2) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("contiguous block size", 2);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos");
    sFactory->initializeFromParameterList(List);
  }

  // initialize the preconditioner
  precFactory->initializePrec(Thyra::defaultLinearOpSource(A_), &*prec);

  RCP<const Thyra::LinearOpBase<double> > op;

  op     = prec->getUnspecifiedPrecOp();
  status = (op != Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
    os << "      "
       << "Preconditioner \"getUnspecifiedPrecOp\" is null (it should not be!)" << std::endl;
    ;
  }
  allPassed &= status;

  op     = prec->getRightPrecOp();
  status = (op == Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
    os << "      "
       << "Preconditioner \"getRightPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  op     = prec->getLeftPrecOp();
  status = (op == Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_initializePrec " << toString(status) << std::endl;
    os << "      "
       << "Preconditioner \"getLeftPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  return allPassed;
}

bool tSIMPLEPreconditionerFactory::test_uninitializePrec(int verbosity, std::ostream &os) {
  return true;
}

bool tSIMPLEPreconditionerFactory::test_isCompatable(int verbosity, std::ostream &os) {
  return true;
}

bool tSIMPLEPreconditionerFactory::test_diagonal(int verbosity, std::ostream &os,
                                                 int use_blocking) {
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
  LinearOp F = Teko::Test::DiagMatrix(2, vec, "F");

  vec[0]     = 7.0;
  vec[1]     = 8.0;
  LinearOp G = Teko::Test::DiagMatrix(2, vec, "G");

  vec[0]     = 5.0;
  vec[1]     = 6.0;
  LinearOp D = Teko::Test::DiagMatrix(2, vec, "D");

  vec[0]     = 3.0;
  vec[1]     = 4.0;
  LinearOp C = Teko::Test::DiagMatrix(2, vec, "C");

  vec[0]      = 1.0;
  vec[1]      = 0.5;
  LinearOp iF = Teko::Test::DiagMatrix(2, vec, "inv(F)");

  vec[0]      = -0.03125;
  vec[1]      = -0.05;
  LinearOp iS = Teko::Test::DiagMatrix(2, vec, "inv(S)");

  RCP<Teko::InverseFactory> invF = rcp(new Teko::StaticOpInverseFactory(iF));
  RCP<Teko::InverseFactory> invS = rcp(new Teko::StaticOpInverseFactory(iS));

  LinearOp A                                = Thyra::block2x2(F, G, D, C);
  RCP<SIMPLEPreconditionerFactory> sFactory = rcp(new SIMPLEPreconditionerFactory(invF, invS, 0.9));
  const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory = sFactory;
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A);

  if (use_blocking == 1) {
    // parameter list for (1,1) block
    Teuchos::ParameterList List, BlkList;
    BlkList.set("number of local blocks", 1);
    BlkList.set("block start index", &*block_starts_);
    BlkList.set("block entry gids", &*block_gids_);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 2) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("contiguous block size", 2);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos");
    sFactory->initializeFromParameterList(List);
  }

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
  ef[0] = 0.21875;
  ef[1] = 0.5;
  eg[0] = -0.028125;
  eg[1] = 0.0;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_diagonal " << toString(status)
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
  ef[0] = 1.71875;
  ef[1] = 1.4;
  eg[0] = -0.478125;
  eg[1] = 0.135;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_diagonal " << toString(status)
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
  ef[0] = -0.09375;
  ef[1] = -1.0;
  eg[0] = 0.140625;
  eg[1] = 0.225;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_diagonal " << toString(status)
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
  ef[0] = 0.9375;
  ef[1] = 2.800000000000001;
  eg[0] = 0.39375;
  eg[1] = -1.08;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_diagonal " << toString(status)
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

bool tSIMPLEPreconditionerFactory::test_result(int verbosity, std::ostream &os, int use_blocking) {
  bool status    = false;
  bool allPassed = true;
  double diff    = -1000.0;

  // Build block2x2 preconditioner
  RCP<SIMPLEPreconditionerFactory> sFactory =
      rcp(new SIMPLEPreconditionerFactory(invF_, invS_, 0.9));
  const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory = sFactory;
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A_);

  if (use_blocking) {
    // parameter list for (1,1) block
    Teuchos::ParameterList List, BlkList;
    BlkList.set("number of local blocks", 1);
    BlkList.set("block start index", &*block_starts_);
    BlkList.set("block entry gids", &*block_gids_);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 2) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("contiguous block size", 2);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos");
    sFactory->initializeFromParameterList(List);
  }

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
  ef[0] = 0.987654320987654;
  ef[1] = 1.074074074074074;
  eg[0] = 0.777777777777778;
  eg[1] = 1.066666666666667;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_result " << toString(status)
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
  ef[0] = 4.197530864197531;
  ef[1] = 2.814814814814815;
  eg[0] = 2.855555555555555;
  eg[1] = 3.633333333333334;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_result " << toString(status)
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
  ef[0] = -0.567901234567901;
  ef[1] = -1.592592592592592;
  eg[0] = -1.122222222222222;
  eg[1] = -1.333333333333333;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_result " << toString(status)
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
  ef[0] = 0.518518518518519;
  ef[1] = 2.888888888888889;
  eg[0] = 1.533333333333334;
  eg[1] = 5.600000000000001;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory::test_result " << toString(status)
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

bool tSIMPLEPreconditionerFactory::test_result_timings(int verbosity, std::ostream &os,
                                                       int use_blocking) {
  bool status    = false;
  bool allPassed = true;
  double diff    = -1000.0;

  // Build block2x2 preconditioner
  RCP<TimingsSIMPLEPreconditionerFactory> sFactory =
      rcp(new TimingsSIMPLEPreconditionerFactory(invF_, invS_, 0.9));
  const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory = sFactory;
  RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory, A_);

  if (use_blocking) {
    // parameter list for (1,1) block
    Teuchos::ParameterList List, BlkList;
    BlkList.set("number of local blocks", 1);
    BlkList.set("block start index", &*block_starts_);
    BlkList.set("block entry gids", &*block_gids_);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 2) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("contiguous block size", 2);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Amesos");
    sFactory->initializeFromParameterList(List);
  }

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
  ef[0] = 0.987654320987654;
  ef[1] = 1.074074074074074;
  eg[0] = 0.777777777777778;
  eg[1] = 1.066666666666667;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tTimingsSIMPLEPreconditionerFactory::test_result " << toString(status)
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
  ef[0] = 4.197530864197531;
  ef[1] = 2.814814814814815;
  eg[0] = 2.855555555555555;
  eg[1] = 3.633333333333334;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tTimingsSIMPLEPreconditionerFactory::test_result " << toString(status)
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
  ef[0] = -0.567901234567901;
  ef[1] = -1.592592592592592;
  eg[0] = -1.122222222222222;
  eg[1] = -1.333333333333333;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tTimingsSIMPLEPreconditionerFactory::test_result " << toString(status)
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
  ef[0] = 0.518518518518519;
  ef[1] = 2.888888888888889;
  eg[0] = 1.533333333333334;
  eg[1] = 5.600000000000001;
  Thyra::apply(*precOp, Thyra::NOTRANS, *x, y.ptr());
  status = ((diff = Teko::Test::Difference(y, z) / Thyra::norm_2(*z->col(0))) < tolerance_);
  if (not status || verbosity >= 10) {
    os << std::endl
       << "   tTimingsSIMPLEPreconditionerFactory::test_result " << toString(status)
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
