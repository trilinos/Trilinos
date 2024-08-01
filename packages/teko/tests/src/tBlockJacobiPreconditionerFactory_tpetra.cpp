// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tBlockJacobiPreconditionerFactory_tpetra.hpp"
#include "Teko_JacobiPreconditionerFactory.hpp"

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
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_LinearOpTester.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#include <vector>
#include "Teko_StratimikosFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Teuchos_AbstractFactoryStd.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;

void tBlockJacobiPreconditionerFactory_tpetra::initializeTest() {
  const Epetra_Comm& comm_epetra                   = *GetComm();
  const RCP<const Teuchos::Comm<int> > comm_tpetra = GetComm_tpetra();

  tolerance_ = 1.0e-14;

  int nx = 39;  // essentially random values
  int ny = 53;

  // create some big blocks to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm_epetra, false);
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  Epetra_CrsMatrix& epetraF = FGallery.GetMatrixRef();
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > tpetraF =
      Teko::TpetraHelpers::epetraCrsMatrixToTpetra(rcpFromRef(epetraF), comm_tpetra);
  F_ = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getRangeMap()), tpetraF);

  Trilinos_Util::CrsMatrixGallery CGallery("laplace_2d", comm_epetra, false);
  CGallery.Set("nx", nx);
  CGallery.Set("ny", ny);
  Epetra_CrsMatrix& epetraC = CGallery.GetMatrixRef();
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > tpetraC =
      Teko::TpetraHelpers::epetraCrsMatrixToTpetra(rcpFromRef(epetraC), comm_tpetra);
  C_ = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraC->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraC->getRangeMap()), tpetraC);

  Trilinos_Util::CrsMatrixGallery BGallery("diag", comm_epetra, false);
  BGallery.Set("nx", nx * ny);
  BGallery.Set("a", 5.0);
  Epetra_CrsMatrix& epetraB = BGallery.GetMatrixRef();
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > tpetraB =
      Teko::TpetraHelpers::epetraCrsMatrixToTpetra(rcpFromRef(epetraB), comm_tpetra);
  B_ = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraB->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraB->getRangeMap()), tpetraB);

  Trilinos_Util::CrsMatrixGallery BtGallery("diag", comm_epetra, false);
  BtGallery.Set("nx", nx * ny);
  BtGallery.Set("a", 3.0);
  Epetra_CrsMatrix& epetraBt = BtGallery.GetMatrixRef();
  RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> > tpetraBt =
      Teko::TpetraHelpers::epetraCrsMatrixToTpetra(rcpFromRef(epetraBt), comm_tpetra);
  Bt_ = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraBt->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraBt->getRangeMap()), tpetraBt);

  // build some inverse operators
  RCP<Tpetra::Vector<ST, LO, GO, NT> > dF =
      rcp(new Tpetra::Vector<ST, LO, GO, NT>(tpetraF->getRangeMap()));
  RCP<Tpetra::Vector<ST, LO, GO, NT> > dC =
      rcp(new Tpetra::Vector<ST, LO, GO, NT>(tpetraC->getRangeMap()));

  tpetraF->getLocalDiagCopy(*dF);
  dF->reciprocal(*dF);

  tpetraC->getLocalDiagCopy(*dC);
  dC->reciprocal(*dC);

  invF_ = Teko::TpetraHelpers::thyraDiagOp(dF, *tpetraF->getRowMap(), "inv(diag(F))");
  invC_ = Teko::TpetraHelpers::thyraDiagOp(dC, *tpetraC->getRowMap(), "inv(diag(C))");
}

int tBlockJacobiPreconditionerFactory_tpetra::runTest(int verbosity, std::ostream& stdstrm,
                                                      std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tBlockJacobiPreconditionerFactory_tpetra";

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

  status = test_isCompatible(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"isCompatible\" ... PASSED", "   \"isCompatible\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_iterativeSolves(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"iterativeSolves\" ... PASSED",
                "   \"iterativeSolves\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tBlockJacobiPreconditionedFactory...PASSED",
                  "tBlockJacobiPreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tBlockJacobiPreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tBlockJacobiPreconditionerFactory_tpetra::test_createPrec(int verbosity, std::ostream& os) {
  RCP<JacobiPreconditionerFactory> fact = rcp(new JacobiPreconditionerFactory(invF_, invC_));

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

  fact = rcp(new JacobiPreconditionerFactory(rcp(new StaticInvDiagStrategy(invF_, invC_))));

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

bool tBlockJacobiPreconditionerFactory_tpetra::test_initializePrec(int verbosity,
                                                                   std::ostream& os) {
  using Thyra::zero;

  bool status    = false;
  bool allPassed = true;

  std::string constrType[3] = {std::string("Static"), std::string("2x2 Static Strategy"),
                               std::string("3x3 Static Strategy")};

  // three by three bloock diagonal
  std::vector<RCP<const Thyra::LinearOpBase<ST> > > invD;
  invD.push_back(invF_);
  invD.push_back(invC_);
  invD.push_back(invF_);

  // allocate new linear operator
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

  const RCP<Thyra::PhysicallyBlockedLinearOpBase<ST> > invBlkOp =
      Thyra::defaultBlockedLinearOp<ST>();
  invBlkOp->beginBlockFill(3, 3);
  invBlkOp->setBlock(0, 0, invF_);
  invBlkOp->setBlock(1, 1, invC_);
  invBlkOp->setBlock(2, 2, invF_);
  invBlkOp->endBlockFill();

  // build factory array
  RCP<JacobiPreconditionerFactory> fact_array[3] = {
      rcp(new JacobiPreconditionerFactory(invF_, invC_)),
      rcp(new JacobiPreconditionerFactory(rcp(new StaticInvDiagStrategy(invF_, invC_)))),
      rcp(new JacobiPreconditionerFactory(rcp(new StaticInvDiagStrategy(invD))))};

  RCP<const Thyra::LinearOpBase<double> > A[3] = {block2x2(F_, Bt_, B_, C_),
                                                  block2x2(F_, Bt_, B_, C_), blkOp};

  // this is what the factory should build
  RCP<const Thyra::LinearOpBase<double> > invA[3] = {
      block2x2(invF_, zero(Bt_->range(), Bt_->domain()), zero(B_->range(), B_->domain()), invC_),
      block2x2(invF_, zero(Bt_->range(), Bt_->domain()), zero(B_->range(), B_->domain()), invC_),
      invBlkOp};

  // test both constructors
  for (int i = 0; i < 3; i++) {
    RCP<const Thyra::LinearOpBase<double> > op;

    RCP<Thyra::PreconditionerFactoryBase<double> > fact = fact_array[i];
    RCP<Thyra::PreconditionerBase<double> > prec        = fact->createPrec();

    // initialize the preconditioner
    fact->initializePrec(Thyra::defaultLinearOpSource(A[i]), &*prec);

    op = prec->getRightPrecOp();
    TEST_EQUALITY(op, Teuchos::null,
                  std::endl
                      << "   tBlockJacobiPreconditionerFactory_tpetra::test_initializePrec "
                      << "using \"" << constrType[i] << "\" constructor " << toString(status)
                      << ": Preconditioner \"getRightPrecOp\" is not null (it should be!)");

    op = prec->getLeftPrecOp();
    TEST_EQUALITY(op, Teuchos::null,
                  std::endl
                      << "   tBlockJacobiPreconditionerFactory_tpetra::test_initializePrec "
                      << "using \"" << constrType[i] << "\" constructor " << toString(status)
                      << ": Preconditioner \"getLeftPrecOp\" is not null (it should be!)");

    op = prec->getUnspecifiedPrecOp();
    TEST_NOT_EQUAL(op, Teuchos::null,
                   std::endl
                       << "   tBlockJacobiPreconditionerFactory_tpetra::test_initializePrec "
                       << "using \"" << constrType[i] << "\" constructor " << toString(status)
                       << ": Preconditioner \"getUnspecifiedPrecOp\" is null (it should not be!)");

    LinearOpTester<double> tester;
    tester.show_all_tests(true);
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*invA[i], *op, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tBlockJacobiPreconditionerFactory_tpetra::test_initializePrec "
                            << ": Comparing factory generated operator to correct operator");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

bool tBlockJacobiPreconditionerFactory_tpetra::test_uninitializePrec(int verbosity,
                                                                     std::ostream& os) {
  return true;
}

bool tBlockJacobiPreconditionerFactory_tpetra::test_isCompatible(int verbosity, std::ostream& os) {
  // bool status = false;
  bool allPassed = true;

  // with the "new" PreconditionerFactory this test is now meaningless.
#if 0
   {
      // three by three bloock diagonal 
      std::vector<RCP<const Thyra::LinearOpBase<double> > > invD;
      invD.push_back(invF_); invD.push_back(invC_); invD.push_back(invF_);
   
      // allocate new linear operator
      const RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blkOp
           = Thyra::defaultBlockedLinearOp<double>();
      blkOp->beginBlockFill(3,3);
      blkOp->setBlock(0,0,F_);  blkOp->setBlock(0,1,Bt_);
      blkOp->setBlock(1,0,B_);  blkOp->setBlock(1,1,C_);  blkOp->setBlock(1,2,B_);
      blkOp->setBlock(2,1,Bt_); blkOp->setBlock(2,2,F_);
      blkOp->endBlockFill();
   
      // build factory array
      RCP<const Thyra::PreconditionerFactoryBase<double> > fact
            = rcp(new JacobiPreconditionerFactory(rcp(new StaticInvDiagStrategy(invD))));
   
      const RCP<const Thyra::LinearOpBase<double> > A = blkOp;
      bool result = fact->isCompatible(*Thyra::defaultLinearOpSource<double>(A));
      TEST_ASSERT(result,
               std::endl << "   tBlockJacobiPreconditionerFactory_tpetra::test_isCompatible " << toString(status) 
               << ": On failure 3x3 based block operator is incompatible with 3x3 preconditioner!");
   }

   {
      // three by three bloock diagonal 
      std::vector<RCP<const Thyra::LinearOpBase<double> > > invD;
      invD.push_back(invF_); invD.push_back(invC_); invD.push_back(invF_);
   
      // allocate new linear operator
      const RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blkOp
           = Thyra::defaultBlockedLinearOp<double>();
      blkOp->beginBlockFill(2,3);
      blkOp->setBlock(0,0,F_);  blkOp->setBlock(0,1,Bt_);
      blkOp->setBlock(1,0,B_);  blkOp->setBlock(1,1,C_);  blkOp->setBlock(1,2,B_);
      blkOp->endBlockFill();
   
      // build factory array
      RCP<const Thyra::PreconditionerFactoryBase<double> > fact
            = rcp(new JacobiPreconditionerFactory(rcp(new StaticInvDiagStrategy(invD))));
   
      const RCP<const Thyra::LinearOpBase<double> > A = blkOp;
      bool result = fact->isCompatible(*Thyra::defaultLinearOpSource<double>(A));
      TEST_ASSERT(not result,
               std::endl << "   tBlockJacobiPreconditionerFactory_tpetra::test_isCompatible " << toString(status) 
               << ": On failure 2x3 based block operator is compatible with 3x3 preconditioner!");
   }

   {
      // three by three bloock diagonal 
      std::vector<RCP<const Thyra::LinearOpBase<double> > > invD;
      invD.push_back(invF_); invD.push_back(invC_); invD.push_back(invF_);
   
      // allocate new linear operator
      const RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blkOp
           = Thyra::defaultBlockedLinearOp<double>();
      blkOp->beginBlockFill(3,2);
      blkOp->setBlock(0,0,F_);  blkOp->setBlock(0,1,Bt_);
      blkOp->setBlock(1,0,B_);  blkOp->setBlock(1,1,C_);
      blkOp->setBlock(2,1,Bt_); 
      blkOp->endBlockFill();
   
      // build factory array
      RCP<const Thyra::PreconditionerFactoryBase<double> > fact
            = rcp(new JacobiPreconditionerFactory(rcp(new StaticInvDiagStrategy(invD))));
   
      const RCP<const Thyra::LinearOpBase<double> > A = blkOp;
      bool result = fact->isCompatible(*Thyra::defaultLinearOpSource<double>(A));
      TEST_ASSERT(not result,
               std::endl << "   tBlockJacobiPreconditionerFactory_tpetra::test_isCompatible " << toString(status) 
               << ": On failure 3x2 based block operator is compatible with 3x3 preconditioner!");
   }
#endif

  return allPassed;
}

bool tBlockJacobiPreconditionerFactory_tpetra::test_iterativeSolves(int verbosity,
                                                                    std::ostream& os) {
  using Thyra::zero;

  bool status    = false;
  bool allPassed = true;

  // allocate new linear operator
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

  // build stratimikos factory, adding Teko's version
  Stratimikos::DefaultLinearSolverBuilder stratFactory;
  stratFactory.setPreconditioningStrategyFactory(
      Teuchos::abstractFactoryStd<Thyra::PreconditionerFactoryBase<double>,
                                  Teko::StratimikosFactory>(),
      "Teko");
  RCP<ParameterList> params = Teuchos::rcp(new ParameterList(*stratFactory.getValidParameters()));
  ParameterList& tekoList   = params->sublist("Preconditioner Types").sublist("Teko");
  tekoList.set("Write Block Operator", false);
  tekoList.set("Test Block Operator", false);
  tekoList.set("Strided Blocking", "1 1");
  tekoList.set("Inverse Type", "BGS");
  ParameterList& ifl = tekoList.sublist("Inverse Factory Library");
  ifl.sublist("BGS").set("Type", "Block Gauss-Seidel");
  ifl.sublist("BGS").set("Inverse Type", "Belos");
  ifl.sublist("BGS").set("Preconditioner Type", "Ifpack2");

  RCP<Teko::InverseLibrary> invLib        = Teko::InverseLibrary::buildFromParameterList(ifl);
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("BGS");

  Teko::ModifiableLinearOp inv = Teko::buildInverse(*invFact, blkOp);
  TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");

  if (!inv.is_null()) {
    Teko::rebuildInverse(*invFact, blkOp, inv);
    TEST_ASSERT(!inv.is_null(), "Constructed preconditioner is null");
  }

  return allPassed;
}

}  // end namespace Test
}  // end namespace Teko
