// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_single_belos_thyra_solver.hpp"
#include "Stratimikos_InternalConfig.h"
#include "Tpetra_Core.hpp"

#include "Thyra_BelosLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Teuchos_ParameterList.hpp"
#include "MatrixMarket_Tpetra.hpp"

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraParameters.hpp"

#ifdef HAVE_STRATIMIKOS_IFPACK2
#  include "Thyra_Ifpack2PreconditionerFactory.hpp"
#endif

Teuchos::RCP<const Thyra::LinearOpBase<double> > getFwdLinearOp(const std::string &matrixFileName)
{
  using crs_type = Tpetra::CrsMatrix<>;
  using op_type = Tpetra::Operator<>;
  using reader_type = Tpetra::MatrixMarket::Reader<crs_type>;
  auto comm = Tpetra::getDefaultComm();
  auto tpetraCrsMatrix = reader_type::readSparseFile(matrixFileName, comm, true);
  auto fwdLinearOp = Thyra::createLinearOp(Teuchos::rcp_dynamic_cast<op_type>(tpetraCrsMatrix, true));
  return fwdLinearOp;
}

bool Thyra::test_single_belos_thyra_solver(
  const std::string                       matrixFile
  ,const bool                             testTranspose
  ,const bool                             usePreconditioner
  ,const int                              numRhs
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxResid
  ,const double                           maxSolutionError
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::ParameterList                 *belosLOWSFPL
  ,Teuchos::ParameterList                 *precPL
  ,Teuchos::FancyOStream                  *out_arg
  )
{
  using Teuchos::rcp;
  using Teuchos::OSTab;
  bool result, success = true;

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::rcp(out_arg,false);

  try {

    if(out.get()) {
      *out << "\n***"
           << "\n*** Testing Thyra::BelosLinearOpWithSolveFactory (and Thyra::BelosLinearOpWithSolve)"
           << "\n***\n"
           << "\nEchoing input options:"
           << "\n  matrixFile             = " << matrixFile
           << "\n  testTranspose          = " << testTranspose
           << "\n  usePreconditioner      = " << usePreconditioner
           << "\n  numRhs                 = " << numRhs
           << "\n  numRandomVectors       = " << numRandomVectors
           << "\n  maxFwdError            = " << maxFwdError
           << "\n  maxResid               = " << maxResid
           << "\n  showAllTests           = " << showAllTests
           << "\n  dumpAll                = " << dumpAll
           << std::endl;
    }

    auto A = getFwdLinearOp(matrixFile);

    if(out.get() && dumpAll) *out << "\ndescribe(A) =\n" << describe(*A,Teuchos::VERB_EXTREME);

    if(out.get()) *out << "\nB) Creating a BelosLinearOpWithSolveFactory object opFactory ...\n";

    Teuchos::RCP<LinearOpWithSolveFactoryBase<double> >
      lowsFactory;
    {
      Teuchos::RCP<BelosLinearOpWithSolveFactory<double> >
        belosLowsFactory = Teuchos::rcp(new BelosLinearOpWithSolveFactory<double>());
      lowsFactory = belosLowsFactory;
    }

    if(out.get()) {
      *out << "\nlowsFactory.getValidParameters() before setting preconditioner factory:\n";
      lowsFactory->getValidParameters()->print(OSTab(out).o(),0,true,false);
    }

    if(usePreconditioner) {
#ifdef HAVE_STRATIMIKOS_IFPACK2
      if(out.get()) {
        *out << "\nSetting an Ifpack2 preconditioner factory ...\n";
      }
      using crs_type = Tpetra::CrsMatrix<>;
      RCP<PreconditionerFactoryBase<double> >
        precFactory = Teuchos::rcp(new Ifpack2PreconditionerFactory<crs_type>());
      if (precPL)
        precFactory->setParameterList(rcp(precPL,false));
      lowsFactory->setPreconditionerFactory(precFactory,"Ifpack2");
#else
      TEUCHOS_TEST_FOR_EXCEPT(usePreconditioner);
#endif
    }
    
    if(out.get()) {
      *out << "\nlowsFactory.getValidParameters() after setting preconditioner factory:\n";
      lowsFactory->getValidParameters()->print(OSTab(out).o(),0,true,false);
      *out << "\nbelosLOWSFPL before setting parameters:\n";
      belosLOWSFPL->print(OSTab(out).o(),0,true);
    }

    lowsFactory->setParameterList(Teuchos::rcp(belosLOWSFPL,false));

    if(out.get()) {
      *out << "\nbelosLOWSFPL after setting parameters:\n";
      belosLOWSFPL->print(OSTab(out).o(),0,true);
    }

    if(out.get()) *out << "\nC) Creating a BelosLinearOpWithSolve object nsA from A ...\n";

    Teuchos::RCP<LinearOpWithSolveBase<double> > nsA = lowsFactory->createOp();
    Thyra::initializeOp<double>(*lowsFactory,  A, nsA.ptr());

    if(out.get()) *out << "\nD) Testing the LinearOpBase interface of nsA ...\n";

    LinearOpTester<double> linearOpTester;
    linearOpTester.check_adjoint(testTranspose);
    linearOpTester.num_rhs(numRhs);
    linearOpTester.num_random_vectors(numRandomVectors);
    linearOpTester.set_all_error_tol(maxFwdError);
    linearOpTester.set_all_warning_tol(1e-2*maxFwdError);
    linearOpTester.show_all_tests(showAllTests);
    linearOpTester.dump_all(dumpAll);
    Thyra::seed_randomize<double>(0);
    result = linearOpTester.check(*nsA,Teuchos::Ptr<Teuchos::FancyOStream>(out.get()));
    if(!result) success = false;

    if(out.get()) *out << "\nE) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
    LinearOpWithSolveTester<double> linearOpWithSolveTester;
    linearOpWithSolveTester.num_rhs(numRhs);
    linearOpWithSolveTester.turn_off_all_tests();
    linearOpWithSolveTester.check_forward_default(true);
    linearOpWithSolveTester.check_forward_residual(true);
    if(testTranspose) {
      linearOpWithSolveTester.check_adjoint_default(true);
      linearOpWithSolveTester.check_adjoint_residual(true);
    }
    else {
      linearOpWithSolveTester.check_adjoint_default(false);
      linearOpWithSolveTester.check_adjoint_residual(false);
    }
    linearOpWithSolveTester.set_all_solve_tol(maxResid);
    linearOpWithSolveTester.set_all_slack_error_tol(maxResid);
    linearOpWithSolveTester.set_all_slack_warning_tol(1e+1*maxResid);
    linearOpWithSolveTester.forward_default_residual_error_tol(2*maxResid);
    linearOpWithSolveTester.forward_default_solution_error_error_tol(maxSolutionError);
    linearOpWithSolveTester.adjoint_default_residual_error_tol(2*maxResid);
    linearOpWithSolveTester.adjoint_default_solution_error_error_tol(maxSolutionError);
    linearOpWithSolveTester.show_all_tests(showAllTests);
    linearOpWithSolveTester.dump_all(dumpAll);
    Thyra::seed_randomize<double>(0);
    result = linearOpWithSolveTester.check(*nsA,out.get());
    if(!result) success = false;

    if(out.get()) {
      *out << "\nbelosLOWSFPL after solving:\n";
      belosLOWSFPL->print(OSTab(out).o(),0,true);
    }
  }
  catch( const std::exception &excpt ) {
    if(out.get()) *out << std::flush;
    std::cerr << "*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }
   
  return success;
    
}
