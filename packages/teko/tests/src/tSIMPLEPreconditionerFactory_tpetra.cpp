// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_Config.h"
#include "tSIMPLEPreconditionerFactory_tpetra.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InverseLibrary.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"

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
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraVectorSpace.hpp"
#include "Teko_Utilities.hpp"

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

void tSIMPLEPreconditionerFactory_tpetra::initializeTest() {
  std::vector<GO> indices(2);
  std::vector<ST> row0(2), row1(2);

  tolerance_ = 1.0e-13;

  comm = GetComm_tpetra();

  const RCP<const Tpetra::Map<LO, GO, NT> > map =
      rcp(new const Tpetra::Map<LO, GO, NT>(2, 0, comm));

  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrF =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrB =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrBt =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrC =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);

  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvF =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvS =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrInvMass =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);

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

  // build block info for block diagonal (one block)
  block_starts_    = arcp(new GO[2], 0, 2);
  block_gids_      = arcp(new GO[2], 0, 2);
  GO *block_starts = block_starts_.getRawPtr();
  GO *block_gids   = block_gids_.getRawPtr();
  block_starts[0]  = 0;
  block_starts[1]  = 2;
  block_gids[0]    = ptrF->getRowMap()->getGlobalElement(0);
  block_gids[1]    = ptrF->getRowMap()->getGlobalElement(1);

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

  // build C matrix
  row0[0] = 1.0;
  row0[1] = 2.0;
  row1[0] = 2.0;
  row1[1] = 1.0;
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
  invF_ = rcp(new StaticOpInverseFactory(Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvF->getRangeMap()), ptrInvF)));

  // build inv(Pschur) matrix
  row0[0] = 0.037037037037037;
  row0[1] = 0.222222222222222;
  row1[0] = 0.222222222222222;
  row1[1] = 0.333333333333333;
  ptrInvS->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row0));
  ptrInvS->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices), Teuchos::ArrayView<ST>(row1));
  ptrInvS->fillComplete();
  invS_ = rcp(new StaticOpInverseFactory(Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvS->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrInvS->getRangeMap()), ptrInvS)));

  A_ = Thyra::block2x2<ST>(F_, Bt_, B_, C_, "A");
}

int tSIMPLEPreconditionerFactory_tpetra::runTest(int verbosity, std::ostream &stdstrm,
                                                 std::ostream &failstrm, int &totalrun) {
  bool allTests = true;
  bool status   = true;
  int failcount = 0;

  failstrm << "tSIMPLEPreconditionerFactory_tpetra";

  status = test_createPrec(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"createPrec\" ... PASSED", "   \"createPrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, -2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec(lumped)\" ... PASSED",
                       "   \"initializePrec(lumped)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, -1);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec(diag)\" ... PASSED",
                       "   \"initializePrec(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm, 0);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"initializePrec(absrowsum)\" ... PASSED",
                       "   \"initializePrec(absrowsum)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  /*
  #ifdef Teko_ENABLE_Isorropia
     status = test_initializePrec(verbosity,failstrm,1);
     Teko_TEST_MSG_tpetra(stdstrm,1,"   \"initializePrec(block-1)\" ... PASSED","
  \"initializePrec(block-1)\" ... FAILED"); allTests &= status; failcount += status ? 0 : 1;
     totalrun++;

     status = test_initializePrec(verbosity,failstrm,2);
     Teko_TEST_MSG_tpetra(stdstrm,1,"   \"initializePrec(block-2)\" ... PASSED","
  \"initializePrec(block-2)\" ... FAILED"); allTests &= status; failcount += status ? 0 : 1;
     totalrun++;
  #endif
  */

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

  status = test_iterativeSolves(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"iterativeSolves\" ... PASSED",
                       "   \"iterativeSolves\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_hierarchicalSolves(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"hierarchicalSolves\" ... PASSED",
                       "   \"hierarchicalSolves\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_diagonal(verbosity, failstrm, 0);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"diagonal(diag)\" ... PASSED",
                       "   \"diagonal(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

#ifdef Teko_ENABLE_Isorropia
  status = test_diagonal(verbosity, failstrm, 1);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"diagonal(block-1)\" ... PASSED",
                       "   \"diagonal(block-1)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_diagonal(verbosity, failstrm, 2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"diagonal(block-2)\" ... PASSED",
                       "   \"diagonal(block-2)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;
#endif

  status = test_result(verbosity, failstrm, 0);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"result(diag)\" ... PASSED",
                       "   \"result(diag)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

#ifdef Teko_ENABLE_Isorropia
  status = test_result(verbosity, failstrm, 1);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"result(block-1)\" ... PASSED",
                       "   \"result(block-1)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_result(verbosity, failstrm, 2);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"result(block-2)\" ... PASSED",
                       "   \"result(block-2)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;
#endif

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tSIMPLEPreconditionedFactory...PASSED",
                         "tSIMPLEPreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tSIMPLEPreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_createPrec(int verbosity, std::ostream &os) {
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

bool tSIMPLEPreconditionerFactory_tpetra::test_initializePrec(int verbosity, std::ostream &os,
                                                              int use_blocking) {
  bool status    = false;
  bool allPassed = true;

  // Build block2x2 preconditioner
  RCP<SIMPLEPreconditionerFactory> sFactory = rcp(new SIMPLEPreconditionerFactory(invF_, 0.9));
  const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory = sFactory;
  RCP<Thyra::PreconditionerBase<double> > prec = precFactory->createPrec();

  if (use_blocking == -2) {
    // parameter list for (1,1) block
    Teuchos::ParameterList List;
    List.set("Explicit Velocity Inverse Type", "Lumped");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == -1) {
    // parameter list for (1,1) block
    Teuchos::ParameterList List;
    List.set("Explicit Velocity Inverse Type", "Diagonal");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 0) {
    // parameter list for (1,1) block
    Teuchos::ParameterList List;
    List.set("Explicit Velocity Inverse Type", "AbsRowSum");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 1) {
    // parameter list for (1,1) block
    Teuchos::ParameterList List, BlkList;
    BlkList.set("number of local blocks", 1);
    BlkList.set("block start index", &*block_starts_);
    BlkList.set("block entry gids", &*block_gids_);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  } else if (use_blocking == 2) {
    Teuchos::ParameterList List, BlkList;
    BlkList.set("contiguous block size", 2);
    List.set("H options", BlkList);
    List.set("Explicit Velocity Inverse Type", "BlkDiag");
    List.set("Inverse Pressure Type", "Ifpack2");
    List.set("Inverse Velocity Type", "Ifpack2");
    sFactory->initializeFromParameterList(List);
  }

  // initialize the preconditioner
  precFactory->initializePrec(Thyra::defaultLinearOpSource(A_), &*prec);

  RCP<const Thyra::LinearOpBase<double> > op;

  op     = prec->getUnspecifiedPrecOp();
  status = (op != Teuchos::null);
  if (not status) {
    os << std::endl
       << "   tSIMPLEPreconditionerFactory_tpetra::test_initializePrec " << toString(status)
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
       << "   tSIMPLEPreconditionerFactory_tpetra::test_initializePrec " << toString(status)
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
       << "   tSIMPLEPreconditionerFactory_tpetra::test_initializePrec " << toString(status)
       << std::endl;
    os << "      "
       << "Preconditioner \"getLeftPrecOp\" is not null (it should be!)" << std::endl;
    ;
  }
  allPassed &= status;

  return allPassed;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_uninitializePrec(int verbosity, std::ostream &os) {
  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_isCompatable(int verbosity, std::ostream &os) {
  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_diagonal(int verbosity, std::ostream &os,
                                                        int use_blocking) {
  /*   // make sure the preconditioner is working by testing against the identity matrix
     typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

     bool status = false;
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

     vec[0] = 1.0; vec[1] = 2.0;
     LinearOp F = Teko::Test::DiagMatrix(2,vec,"F");

     vec[0] = 7.0; vec[1] = 8.0;
     LinearOp G = Teko::Test::DiagMatrix(2,vec,"G");

     vec[0] = 5.0; vec[1] = 6.0;
     LinearOp D = Teko::Test::DiagMatrix(2,vec,"D");

     vec[0] = 3.0; vec[1] = 4.0;
     LinearOp C = Teko::Test::DiagMatrix(2,vec,"C");

     vec[0] = 1.0; vec[1] = 0.5;
     LinearOp iF = Teko::Test::DiagMatrix(2,vec,"inv(F)");

     vec[0] = -0.03125; vec[1] = -0.05;
     LinearOp iS = Teko::Test::DiagMatrix(2,vec,"inv(S)");

     RCP<Teko::InverseFactory> invF = rcp(new Teko::StaticOpInverseFactory(iF));
     RCP<Teko::InverseFactory> invS = rcp(new Teko::StaticOpInverseFactory(iS));

     LinearOp A = Thyra::block2x2(F,G,D,C);
     RCP<SIMPLEPreconditionerFactory> sFactory = rcp(new
     SIMPLEPreconditionerFactory(invF,invS,0.9)); const RCP<const
     Thyra::PreconditionerFactoryBase<double> > precFactory =sFactory;
     RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,A);

     if(use_blocking==1){
       // parameter list for (1,1) block
       Teuchos::ParameterList List,BlkList;
       BlkList.set("number of local blocks",1);
       BlkList.set("block start index",&*block_starts_);
       BlkList.set("block entry gids",&*block_gids_);
       List.set("H options",BlkList);
       List.set("Explicit Velocity Inverse Type","BlkDiag");
       List.set("Inverse Pressure Type","Amesos");
       sFactory->initializeFromParameterList(List);
     }
     else if(use_blocking==2){
       Teuchos::ParameterList List,BlkList;
       BlkList.set("contiguous block size",2);
       List.set("H options",BlkList);
       List.set("Explicit Velocity Inverse Type","BlkDiag");
       List.set("Inverse Pressure Type","Amesos");
       sFactory->initializeFromParameterList(List);
     }

     // build linear operator
     RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

     const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));
     // construct a couple of vectors
     Epetra_Vector ea(*map),eb(*map);
     Epetra_Vector ef(*map),eg(*map);
     const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea,eb,A->domain());
     const RCP<const Thyra::MultiVectorBase<double> > z = BlockVector(ef,eg,A->domain());
     const RCP<Thyra::MultiVectorBase<double> > y = Thyra::createMembers(A->range(),1);

     // now checks of the preconditioner (should be exact!)
     /////////////////////////////////////////////////////////////////////////

     // test vector [0 1 1 3]
     ea[0] = 0.0; ea[1] = 1.0; eb[0] = 1.0; eb[1] = 3.0;
     ef[0] =  0.21875; ef[1]  =  0.5;
     eg[0] = -0.028125; eg[1] =  0.0;
     Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
     status = ((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_);
     if(not status || verbosity>=10 ) {
        os << std::endl << "   tSIMPLEPreconditionerFactory_tpetra::test_diagonal " <<
     toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
                        << diff << ")" << std::endl;
        os << "      "; Print(os,"x",x);
        os << "      "; Print(os,"y",y);
        os << "      "; Print(os,"z",z);
     }
     allPassed &= status;

     // test vector [-2 4 7 9]
     ea[0] =-2.0; ea[1] = 4.0; eb[0] = 7.0; eb[1] = 9.0;
     ef[0] = 1.71875; ef[1] =  1.4;
     eg[0] = -0.478125; eg[1] = 0.135;
     Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
     status = ((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_);
     if(not status || verbosity>=10 ) {
        os << std::endl << "   tSIMPLEPreconditionerFactory_tpetra::test_diagonal " <<
     toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
                        << diff << ")" << std::endl;
        os << "      "; Print(os,"x",x);
        os << "      "; Print(os,"y",y);
        os << "      "; Print(os,"z",z);
     }
     allPassed &= status;

     // test vector [1 0 0 -5]
     ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] =-5.0;
     ef[0] = -0.09375; ef[1] = -1.0;
     eg[0] =  0.140625; eg[1] =  0.225;
     Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
     status = ((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_);
     if(not status || verbosity>=10 ) {
        os << std::endl << "   tSIMPLEPreconditionerFactory_tpetra::test_diagonal " <<
     toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
                        << diff << ")" << std::endl;
        os << "      "; Print(os,"x",x);
        os << "      "; Print(os,"y",y);
        os << "      "; Print(os,"z",z);
     }
     allPassed &= status;

     // test vector [4 -4 6 12]
     ea[0] = 4.0; ea[1] =-4.0; eb[0] = 6.0; eb[1] =12.0;
     ef[0] = 0.9375; ef[1] =  2.800000000000001;
     eg[0] = 0.39375; eg[1] = -1.08;
     Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
     status = ((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_);
     if(not status || verbosity>=10 ) {
        os << std::endl << "   tSIMPLEPreconditionerFactory_tpetra::test_diagonal " <<
     toString(status) << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
                        << diff << ")" << std::endl;
        os << "      "; Print(os,"x",x);
        os << "      "; Print(os,"y",y);
        os << "      "; Print(os,"z",z);
     }
     allPassed &= status;

     return allPassed;*/
  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_result(int verbosity, std::ostream &os,
                                                      int use_blocking) {
  /* bool status = false;
   bool allPassed = true;
   double diff = -1000.0;

   // Build block2x2 preconditioner
   RCP<SIMPLEPreconditionerFactory> sFactory = rcp(new
   SIMPLEPreconditionerFactory(invF_,invS_,0.9)); const RCP<const
   Thyra::PreconditionerFactoryBase<double> > precFactory  = sFactory;
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,A_);

   if(use_blocking){
     // parameter list for (1,1) block
     Teuchos::ParameterList List,BlkList;
     BlkList.set("number of local blocks",1);
     BlkList.set("block start index",&*block_starts_);
     BlkList.set("block entry gids",&*block_gids_);
     List.set("H options",BlkList);
     List.set("Explicit Velocity Inverse Type","BlkDiag");
     List.set("Inverse Pressure Type","Amesos");
     sFactory->initializeFromParameterList(List);
   }
   else if(use_blocking==2){
     Teuchos::ParameterList List,BlkList;
     BlkList.set("contiguous block size",2);
     List.set("H options",BlkList);
     List.set("Explicit Velocity Inverse Type","BlkDiag");
     List.set("Inverse Pressure Type","Amesos");
     sFactory->initializeFromParameterList(List);
   }

   // build linear operator
   RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp();

   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));
   // construct a couple of vectors
   Epetra_Vector ea(*map),eb(*map);
   Epetra_Vector ef(*map),eg(*map);

   const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea,eb,A_->domain());
   const RCP<const Thyra::MultiVectorBase<double> > z = BlockVector(ef,eg,A_->domain());
   const RCP<Thyra::MultiVectorBase<double> > y = Thyra::createMembers(A_->range(),1);

   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());

   // now checks of the preconditioner (should be exact!)
   /////////////////////////////////////////////////////////////////////////

   // test vector [0 1 1 3]
   ea[0] = 0.0; ea[1] = 1.0; eb[0] = 1.0; eb[1] = 3.0;
   ef[0] = 0.987654320987654; ef[1] = 1.074074074074074;
   eg[0] = 0.777777777777778; eg[1] = 1.066666666666667;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   status = ((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_);
   if(not status || verbosity>=10 ) {
      os << std::endl << "   tSIMPLEPreconditionerFactory_tpetra::test_result " << toString(status)
   << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [-2 4 7 9]
   ea[0] =-2.0; ea[1] = 4.0; eb[0] = 7.0; eb[1] = 9.0;
   ef[0] = 4.197530864197531; ef[1] = 2.814814814814815;
   eg[0] = 2.855555555555555; eg[1] = 3.633333333333334;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   status = ((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_);
   if(not status || verbosity>=10 ) {
      os << std::endl << "   tSIMPLEPreconditionerFactory_tpetra::test_result " << toString(status)
   << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [1 0 0 -5]
   ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] =-5.0;
   ef[0] = -0.567901234567901; ef[1] = -1.592592592592592;
   eg[0] = -1.122222222222222; eg[1] = -1.333333333333333;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   status = ((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_);
   if(not status || verbosity>=10 ) {
      os << std::endl << "   tSIMPLEPreconditionerFactory_tpetra::test_result " << toString(status)
   << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   // test vector [4 -4 6 12]
   ea[0] = 4.0; ea[1] =-4.0; eb[0] = 6.0; eb[1] =12.0;
   ef[0] = 0.518518518518519; ef[1] = 2.888888888888889;
   eg[0] = 1.533333333333334; eg[1] = 5.600000000000001;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   status = ((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_);
   if(not status || verbosity>=10 ) {
      os << std::endl << "   tSIMPLEPreconditionerFactory_tpetra::test_result " << toString(status)
   << ":  (y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
                      << diff << ")" << std::endl;
      os << "      "; Print(os,"x",x);
      os << "      "; Print(os,"y",y);
      os << "      "; Print(os,"z",z);
   }
   allPassed &= status;

   return allPassed;*/
  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_iterativeSolves(int verbosity, std::ostream &os) {
  bool status    = false;
  bool allPassed = true;

  RCP<ParameterList> params = Teuchos::rcp(new ParameterList());
  ParameterList &tekoList   = params->sublist("Preconditioner Types").sublist("Teko");
  tekoList.set("Inverse Type", "SIMPLE");
  ParameterList &ifl = tekoList.sublist("Inverse Factory Library");
  ifl.sublist("SIMPLE").set("Type", "NS SIMPLE");
  ifl.sublist("SIMPLE").set("Inverse Velocity Type", "Belos");
  ifl.sublist("SIMPLE").set("Preconditioner Velocity Type", "Ifpack2");
  ifl.sublist("SIMPLE").set("Inverse Pressure Type", "Belos");
  ifl.sublist("SIMPLE").set("Preconditioner Pressure Type", "Ifpack2");

  RCP<Teko::InverseLibrary> invLib        = Teko::InverseLibrary::buildFromParameterList(ifl);
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("SIMPLE");

  Teko::ModifiableLinearOp inv = Teko::buildInverse(*invFact, A_);
  TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");

  if (!inv.is_null()) {
    Teko::rebuildInverse(*invFact, A_, inv);
    TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");
  }

  return true;
}

bool tSIMPLEPreconditionerFactory_tpetra::test_hierarchicalSolves(int verbosity, std::ostream &os) {
  bool status    = false;
  bool allPassed = true;

  const RCP<Thyra::PhysicallyBlockedLinearOpBase<ST> > blkOp = Thyra::defaultBlockedLinearOp<ST>();
  blkOp->beginBlockFill(3, 3);
  blkOp->setBlock(0, 0, F_);
  blkOp->setBlock(0, 1, Bt_);
  blkOp->setBlock(1, 0, B_);
  blkOp->setBlock(1, 1, C_);
  blkOp->setBlock(1, 2, B_);
  blkOp->setBlock(2, 1, Bt_);
  blkOp->setBlock(2, 2, F_);
  blkOp->endBlockFill();

  RCP<ParameterList> params = Teuchos::rcp(new ParameterList());
  ParameterList &tekoList   = params->sublist("Preconditioner Types").sublist("Teko");
  tekoList.set("Inverse Type", "HBGS");
  ParameterList &ifl = tekoList.sublist("Inverse Factory Library");

  ifl.sublist("HBGS").set("Type", "Hierarchical Block Gauss-Seidel");

  ifl.sublist("HBGS").sublist("Hierarchical Block 1").set("Included Subblocks", "1,2");
  ifl.sublist("HBGS").sublist("Hierarchical Block 1").set("Inverse Type", "Belos");
  ifl.sublist("HBGS").sublist("Hierarchical Block 1").set("Preconditioner Type", "SIMPLE");

  ifl.sublist("HBGS").sublist("Hierarchical Block 2").set("Included Subblocks", "3");
  ifl.sublist("HBGS").sublist("Hierarchical Block 2").set("Inverse Type", "Belos");
  ifl.sublist("HBGS").sublist("Hierarchical Block 2").set("Preconditioner Type", "SINGLE_BLOCK");

  ifl.sublist("SIMPLE").set("Type", "NS SIMPLE");
  ifl.sublist("SIMPLE").set("Inverse Velocity Type", "Belos");
  ifl.sublist("SIMPLE").set("Preconditioner Velocity Type", "Ifpack2");
  ifl.sublist("SIMPLE").set("Inverse Pressure Type", "Belos");
  ifl.sublist("SIMPLE").set("Preconditioner Pressure Type", "Ifpack2");

  ifl.sublist("SINGLE_BLOCK").set("Type", "Ifpack2");

  RCP<Teko::InverseLibrary> invLib        = Teko::InverseLibrary::buildFromParameterList(ifl);
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("HBGS");

  Teko::ModifiableLinearOp inv = Teko::buildInverse(*invFact, blkOp);
  TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");

  if (!inv.is_null()) {
    Teko::rebuildInverse(*invFact, blkOp, inv);
    TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");
  }

  return true;
}

}  // end namespace Test
}  // end namespace Teko
