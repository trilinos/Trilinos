// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tBlockJacobiPreconditionerFactory.hpp"
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

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;

void tBlockJacobiPreconditionerFactory::initializeTest() {
  const Epetra_Comm& comm = *GetComm();

  tolerance_ = 1.0e-14;

  int nx = 39;  // essentially random values
  int ny = 53;

  // create some big blocks to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm,
                                           false);  // CJ TODO FIXME: change for Epetra64
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  Epetra_CrsMatrix& epetraF = FGallery.GetMatrixRef();
  F_                        = Thyra::epetraLinearOp(rcpFromRef(epetraF));

  Trilinos_Util::CrsMatrixGallery CGallery("laplace_2d", comm,
                                           false);  // CJ TODO FIXME: change for Epetra64
  CGallery.Set("nx", nx);
  CGallery.Set("ny", ny);
  Epetra_CrsMatrix& epetraC = CGallery.GetMatrixRef();
  C_                        = Thyra::epetraLinearOp(rcpFromRef(epetraC));

  Trilinos_Util::CrsMatrixGallery BGallery("diag", comm,
                                           false);  // CJ TODO FIXME: change for Epetra64
  BGallery.Set("nx", nx * ny);
  BGallery.Set("a", 5.0);
  B_ = Thyra::epetraLinearOp(rcpFromRef(BGallery.GetMatrixRef()));

  Trilinos_Util::CrsMatrixGallery BtGallery("diag", comm,
                                            false);  // CJ TODO FIXME: change for Epetra64
  BtGallery.Set("nx", nx * ny);
  BtGallery.Set("a", 3.0);
  Bt_ = Thyra::epetraLinearOp(rcpFromRef(BtGallery.GetMatrixRef()));

  // build some inverse operators
  RCP<Epetra_Vector> dF = rcp(new Epetra_Vector(epetraF.OperatorRangeMap()));
  RCP<Epetra_Vector> dC = rcp(new Epetra_Vector(epetraC.OperatorRangeMap()));

  epetraF.ExtractDiagonalCopy(*dF);
  dF->Reciprocal(*dF);

  epetraC.ExtractDiagonalCopy(*dC);
  dC->Reciprocal(*dC);

  RCP<const Thyra::VectorSpaceBase<double> > vsF = F_->range();
  RCP<const Thyra::VectorSpaceBase<double> > vsC = C_->range();

  invF_ = Thyra::diagonal(Thyra::create_Vector(dF, vsF));
  invC_ = Thyra::diagonal(Thyra::create_Vector(dC, vsC));
}

int tBlockJacobiPreconditionerFactory::runTest(int verbosity, std::ostream& stdstrm,
                                               std::ostream& failstrm, int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tBlockJacobiPreconditionerFactory";

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

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tBlockJacobiPreconditionedFactory...PASSED",
                  "tBlockJacobiPreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tBlockJacobiPreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tBlockJacobiPreconditionerFactory::test_createPrec(int verbosity, std::ostream& os) {
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

bool tBlockJacobiPreconditionerFactory::test_initializePrec(int verbosity, std::ostream& os) {
  using Thyra::zero;

  bool status    = false;
  bool allPassed = true;

  std::string constrType[3] = {std::string("Static"), std::string("2x2 Static Strategy"),
                               std::string("3x3 Static Strategy")};

  // three by three bloock diagonal
  std::vector<RCP<const Thyra::LinearOpBase<double> > > invD;
  invD.push_back(invF_);
  invD.push_back(invC_);
  invD.push_back(invF_);

  // allocate new linear operator
  const RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blkOp =
      Thyra::defaultBlockedLinearOp<double>();
  blkOp->beginBlockFill(3, 3);
  blkOp->setBlock(0, 0, F_);
  blkOp->setBlock(0, 1, Bt_);
  blkOp->setBlock(1, 0, B_);
  blkOp->setBlock(1, 1, C_);
  blkOp->setBlock(1, 2, B_);
  blkOp->setBlock(2, 1, Bt_);
  blkOp->setBlock(2, 2, F_);
  blkOp->endBlockFill();

  const RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > invBlkOp =
      Thyra::defaultBlockedLinearOp<double>();
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
                      << "   tBlockJacobiPreconditionerFactory::test_initializePrec "
                      << "using \"" << constrType[i] << "\" constructor " << toString(status)
                      << ": Preconditioner \"getRightPrecOp\" is not null (it should be!)");

    op = prec->getLeftPrecOp();
    TEST_EQUALITY(op, Teuchos::null,
                  std::endl
                      << "   tBlockJacobiPreconditionerFactory::test_initializePrec "
                      << "using \"" << constrType[i] << "\" constructor " << toString(status)
                      << ": Preconditioner \"getLeftPrecOp\" is not null (it should be!)");

    op = prec->getUnspecifiedPrecOp();
    TEST_NOT_EQUAL(op, Teuchos::null,
                   std::endl
                       << "   tBlockJacobiPreconditionerFactory::test_initializePrec "
                       << "using \"" << constrType[i] << "\" constructor " << toString(status)
                       << ": Preconditioner \"getUnspecifiedPrecOp\" is null (it should not be!)");

    LinearOpTester<double> tester;
    tester.show_all_tests(true);
    std::stringstream ss;
    Teuchos::FancyOStream fos(rcpFromRef(ss), "      |||");
    const bool result = tester.compare(*invA[i], *op, Teuchos::ptrFromRef(fos));
    TEST_ASSERT(result, std::endl
                            << "   tBlockJacobiPreconditionerFactory::test_initializePrec "
                            << ": Comparing factory generated operator to correct operator");
    if (not result || verbosity >= 10) os << ss.str();
  }

  return allPassed;
}

bool tBlockJacobiPreconditionerFactory::test_uninitializePrec(int verbosity, std::ostream& os) {
  return true;
}

bool tBlockJacobiPreconditionerFactory::test_isCompatible(int verbosity, std::ostream& os) {
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
               std::endl << "   tBlockJacobiPreconditionerFactory::test_isCompatible " << toString(status) 
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
               std::endl << "   tBlockJacobiPreconditionerFactory::test_isCompatible " << toString(status) 
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
               std::endl << "   tBlockJacobiPreconditionerFactory::test_isCompatible " << toString(status) 
               << ": On failure 3x2 based block operator is compatible with 3x3 preconditioner!");
   }
#endif

  return allPassed;
}

}  // end namespace Test
}  // end namespace Teko
