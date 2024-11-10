// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "test_single_aztecoo_thyra_solver.hpp"

#ifndef SUN_CXX

#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveFactoryExamples.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_SerialComm.h"
#include "Teuchos_ParameterList.hpp"
#ifdef HAVE_AZTECOO_IFPACK
#  include "Thyra_IfpackPreconditionerFactory.hpp"
#endif
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif

#endif // SUN_CXX

bool Thyra::test_single_aztecoo_thyra_solver(
  const std::string                       matrixFile
  ,const bool                             testTranspose
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxResid
  ,const double                           maxSolutionError
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::ParameterList                 *aztecooLOWSFPL
  ,Teuchos::FancyOStream                  *out_arg
  )
{
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::OSTab;
  typedef Teuchos::ParameterList::PrintOptions PLPrintOptions;
  bool result, success = true;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::rcp(out_arg,false);

  try {

#ifndef SUN_CXX

    if(out.get()) {
      *out
        << "\n***"
        << "\n*** Testing Thyra::AztecOOLinearOpWithSolveFactory (and Thyra::AztecOOLinearOpWithSolve)"
        << "\n***\n"
        << "\nEchoing input options:"
        << "\n  matrixFile             = " << matrixFile
        << "\n  testTranspose          = " << testTranspose
        << "\n  numRandomVectors       = " << numRandomVectors
        << "\n  maxFwdError            = " << maxFwdError
        << "\n  maxResid               = " << maxResid
        << "\n  showAllTests           = " << showAllTests
        << "\n  dumpAll                = " << dumpAll
        << std::endl;
    }

    const bool useAztecPrec = (
      aztecooLOWSFPL
      &&
      aztecooLOWSFPL->sublist("Forward Solve")
      .sublist("AztecOO Settings")
      .get("Aztec Preconditioner","none")!="none"
      );

    if(out.get()) {
      if(useAztecPrec)
        *out << "\nUsing aztec preconditioning so we will not test adjoint solves using internal preconditioning ...\n";
    }

    if(out.get()) *out << "\nA) Reading in an epetra matrix A from the file \'"<<matrixFile<<"\' ...\n";

#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    RCP<Epetra_CrsMatrix> epetra_A;
    EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );

    RCP<const LinearOpBase<double> > A = Thyra::epetraLinearOp(epetra_A);

    if(out.get() && dumpAll) *out << "\ndescribe(A) =\n" << describe(*A,Teuchos::VERB_EXTREME);

    if(out.get()) *out << "\nB) Creating a AztecOOLinearOpWithSolveFactory object opFactory ...\n";

    RCP<LinearOpWithSolveFactoryBase<double> >
      lowsFactory = Teuchos::rcp(new AztecOOLinearOpWithSolveFactory());
    if(out.get()) {
      *out << "\nlowsFactory.getValidParameters() initially:\n";
      lowsFactory->getValidParameters()->print(OSTab(out).o(),PLPrintOptions().showTypes(true).showDoc(true));
    }
    TEUCHOS_ASSERT(aztecooLOWSFPL != NULL);
    aztecooLOWSFPL->sublist("Forward Solve").set("Tolerance",maxResid);
    aztecooLOWSFPL->sublist("Adjoint Solve").set("Tolerance",maxResid);
    if(showAllTests) {
      aztecooLOWSFPL->set("Output Every RHS",bool(true));
    }
    if(out.get()) {
      *out << "\naztecooLOWSFPL before setting parameters:\n";
      aztecooLOWSFPL->print(OSTab(out).o(),0,true);
    }
    if(aztecooLOWSFPL) lowsFactory->setParameterList(Teuchos::rcp(aztecooLOWSFPL,false));

    if(out.get()) *out << "\nC) Creating a AztecOOLinearOpWithSolve object nsA from A ...\n";

    RCP<LinearOpWithSolveBase<double> >
      nsA = lowsFactory->createOp();

    Thyra::initializeOp<double>(*lowsFactory, A, nsA.ptr());

    if(out.get()) *out << "\nD) Testing the LinearOpBase interface of nsA ...\n";

    LinearOpTester<double> linearOpTester;
    linearOpTester.check_adjoint(testTranspose);
    linearOpTester.num_random_vectors(numRandomVectors);
    linearOpTester.set_all_error_tol(maxFwdError);
    linearOpTester.set_all_warning_tol(1e-2*maxFwdError);
    linearOpTester.show_all_tests(showAllTests);
    linearOpTester.dump_all(dumpAll);
    Thyra::seed_randomize<double>(0);
    result = linearOpTester.check(*nsA,out());
    if(!result) success = false;

    if(out.get()) *out << "\nE) Testing the LinearOpWithSolveBase interface of nsA ...\n";

    LinearOpWithSolveTester<double> linearOpWithSolveTester;
    linearOpWithSolveTester.turn_off_all_tests();
    linearOpWithSolveTester.check_forward_default(true);
    linearOpWithSolveTester.check_forward_residual(true);
    if(testTranspose && useAztecPrec) {
      linearOpWithSolveTester.check_adjoint_default(true);
      linearOpWithSolveTester.check_adjoint_residual(true);
    }
    else {
      linearOpWithSolveTester.check_adjoint_default(false);
      linearOpWithSolveTester.check_adjoint_residual(false);
    }
    linearOpWithSolveTester.set_all_solve_tol(maxResid);
    linearOpWithSolveTester.set_all_slack_error_tol(maxResid);
    linearOpWithSolveTester.set_all_slack_warning_tol(10.0*maxResid);
    linearOpWithSolveTester.forward_default_residual_error_tol(2.5*maxResid);
    linearOpWithSolveTester.forward_default_solution_error_error_tol(maxSolutionError);
    linearOpWithSolveTester.adjoint_default_residual_error_tol(2.5*maxResid);
    linearOpWithSolveTester.adjoint_default_solution_error_error_tol(maxSolutionError);
    linearOpWithSolveTester.show_all_tests(showAllTests);
    linearOpWithSolveTester.dump_all(dumpAll);
    Thyra::seed_randomize<double>(0);
    result = linearOpWithSolveTester.check(*nsA,out.get());
    if(!result) success = false;

    if(out.get()) *out << "\nF) Uninitialize nsA, create preconditioner for diagonal scaled by 0.99 and then reinitialize nsA reusing the old preconditioner ...\n";

    // Scale the diagonal of the matrix and then create the preconditioner for it
    Thyra::uninitializeOp<double>(*lowsFactory, nsA.ptr());
    // Above is not required but a good idea since we are changing the matrix
    {
      Epetra_Vector diag(epetra_A->RowMap());
      epetra_A->ExtractDiagonalCopy(diag);
      diag.Scale(0.5);
      epetra_A->ReplaceDiagonalValues(diag);
    }
    Thyra::initializeOp<double>(*lowsFactory, A, nsA.ptr());

    // Scale the matrix back again and then reuse the preconditioner
    Thyra::uninitializeOp<double>(*lowsFactory, nsA.ptr());
    // Above is not required but a good idea since we are changing the matrix
    {
      Epetra_Vector diag(epetra_A->RowMap());
      epetra_A->ExtractDiagonalCopy(diag);
      diag.Scale(1.0/0.5);
      epetra_A->ReplaceDiagonalValues(diag);
    }
    initializeAndReuseOp<double>(*lowsFactory, A, nsA.ptr());

    if(out.get()) *out << "\nG) Testing the LinearOpWithSolveBase interface of nsA ...\n";

    Thyra::seed_randomize<double>(0);
    result = linearOpWithSolveTester.check(*nsA,out.get());
    if(!result) success = false;

    if(useAztecPrec) {

      if(out.get()) *out << "\nH) Reinitialize (A,A,PRECONDITIONER_INPUT_TYPE_AS_MATRIX) => nsA ...\n";

      initializeApproxPreconditionedOp<double>(*lowsFactory, A, A, nsA.ptr());

      if(out.get()) *out << "\nI) Testing the LinearOpWithSolveBase interface of nsA ...\n";

      Thyra::seed_randomize<double>(0);
      result = linearOpWithSolveTester.check(*nsA,out.get());
      if(!result) success = false;

      if(testTranspose && useAztecPrec) {
        linearOpWithSolveTester.check_adjoint_default(true);
        linearOpWithSolveTester.check_adjoint_residual(true);
      }
      else {
        linearOpWithSolveTester.check_adjoint_default(false);
        linearOpWithSolveTester.check_adjoint_residual(false);
      }

    }
    else {

      if(out.get()) *out << "\nSkipping testing steps H and I since we are not using aztec preconditioning and therefore will not test with an external preconditioner matrix!\n";

    }


    RCP<PreconditionerFactoryBase<double> >
      precFactory;

#ifdef HAVE_AZTECOO_IFPACK

    if(useAztecPrec) {

      if(testTranspose) {
        linearOpWithSolveTester.check_adjoint_default(true);
        linearOpWithSolveTester.check_adjoint_residual(true);
      }

      if(out.get()) *out << "\nJ) Create an ifpack preconditioner precA for A ...\n";

      precFactory = Teuchos::rcp(new IfpackPreconditionerFactory());

      if(out.get()) {
        *out << "\nprecFactory.description() = " << precFactory->description() << std::endl;
        *out << "\nprecFactory.getValidParameters() =\n";
        precFactory->getValidParameters()->print(OSTab(out).o(),0,true,false);
      }

      RCP<Teuchos::ParameterList>
        ifpackPFPL = Teuchos::rcp(new Teuchos::ParameterList("IfpackPreconditionerFactory"));
      ifpackPFPL->set("Prec Type","ILUT");
      ifpackPFPL->set("Overlap",int(1));
      if(out.get()) {
        *out << "\nifpackPFPL before setting parameters =\n";
        ifpackPFPL->print(OSTab(out).o(),0,true);
      }

      precFactory->setParameterList(ifpackPFPL);

      RCP<PreconditionerBase<double> >
        precA = precFactory->createPrec();
      Thyra::initializePrec<double>(*precFactory,A,&*precA);

      if(out.get()) {
        *out << "\nifpackPFPL after setting parameters =\n";
        ifpackPFPL->print(OSTab(out).o(),0,true);
        *out << "\nprecFactory.description() = " << precFactory->description() << std::endl;
      }

      if(out.get()) *out << "\nprecA.description() = " << precA->description() << std::endl;
      if(out.get() && dumpAll) *out << "\ndescribe(precA) =\n" << describe(*precA,Teuchos::VERB_EXTREME);

      if(out.get()) *out << "\nK) Reinitialize (A,precA->getUnspecifiedPrecOp(),PRECONDITIONER_INPUT_TYPE_AS_OPERATOR) => nsA ...\n";

      Thyra::initializePreconditionedOp<double>(*lowsFactory,A,precA,&*nsA);

      if(out.get()) *out << "\nL) Testing the LinearOpWithSolveBase interface of nsA ...\n";

      Thyra::seed_randomize<double>(0);
      result = linearOpWithSolveTester.check(*nsA,out.get());
      if(!result) success = false;

      if(testTranspose && useAztecPrec) {
        linearOpWithSolveTester.check_adjoint_default(true);
        linearOpWithSolveTester.check_adjoint_residual(true);
      }
      else {
        linearOpWithSolveTester.check_adjoint_default(false);
        linearOpWithSolveTester.check_adjoint_residual(false);
      }

    }
    else {

      if(out.get()) *out << "\nSkipping testing steps J, K, and L since we are not using aztec preconditioning and therefore will not test with an ifpack preconditioner!\n";

    }

#else // HAVE_AZTECOO_IFPACK

    if(out.get()) *out << "\nSkipping testing steps J, K, and L since they require ifpack and ifpack has not been enabled!\n";

#endif // HAVE_AZTECOO_IFPACK


    if(out.get()) *out << "\nM) Scale the epetra_A object by 2.5, and then reinitialize nsA with epetra_A ...\n";

    Thyra::uninitializeOp<double>(*lowsFactory, nsA.ptr());
    // Not required but a good idea since we are changing the matrix
    epetra_A->Scale(2.5);
    initializeOp<double>(*lowsFactory, A, nsA.ptr());

    if(out.get()) *out << "\nN) Testing the LinearOpWithSolveBase interface of nsA ...\n";

    Thyra::seed_randomize<double>(0);
    result = linearOpWithSolveTester.check(*nsA,out.get());
    if(!result) success = false;

    if(out.get()) *out << "\nO) Create a scaled (by 2.5) copy epetra_A2 of epetra_A, and then reinitialize nsA with epetra_A2 ...\n";

    RCP<Epetra_CrsMatrix>
      epetra_A2 = Teuchos::rcp(new Epetra_CrsMatrix(*epetra_A));
    epetra_A2->Scale(2.5);
    RCP<const LinearOpBase<double> >
      A2 = Thyra::epetraLinearOp(epetra_A2);
    initializeOp<double>(*lowsFactory, A2, nsA.ptr());
    // Note that it was okay not to uninitialize nsA first here since A, which
    // was used to initialize nsA last, was not changed and therefore the
    // state of nsA was fine throughout

    if(out.get()) *out << "\nP) Testing the LinearOpWithSolveBase interface of nsA ...\n";

    Thyra::seed_randomize<double>(0);
    result = linearOpWithSolveTester.check(*nsA,out.get());
    if(!result) success = false;

    if(!useAztecPrec) {

      if(out.get()) *out << "\nQ) Create an implicitly scaled (by 2.5) and transposed matrix A3 = scale(2.5,transpose(A)) and initialize nsA2 ...\n";

      RCP<const LinearOpBase<double> >
        A3 = scale<double>(2.5,transpose<double>(A));
      RCP<LinearOpWithSolveBase<double> >
        nsA2 = linearOpWithSolve(*lowsFactory,A3);

      if(out.get()) *out << "\nR) Testing the LinearOpWithSolveBase interface of nsA2 ...\n";

      Thyra::seed_randomize<double>(0);
      result = linearOpWithSolveTester.check(*nsA2,out.get());
      if(!result) success = false;

      if(out.get()) *out << "\nS) Testing that LinearOpBase interfaces of transpose(nsA) == nsA2 ...\n";

      result = linearOpTester.compare(
        *transpose(Teuchos::rcp_implicit_cast<const LinearOpBase<double> >(nsA)),*nsA2
        ,out()
        );
      if(!result) success = false;

    }
    else {

      if(out.get()) *out << "\nSkipping testing steps Q, R, and S because we are using internal AztecOO preconditioners!\n";

    }

  if(out.get()) *out << "\nT) Running example use cases ...\n";

  nonExternallyPreconditionedLinearSolveUseCases(
    *A,*lowsFactory,testTranspose,*out
    );

  if(precFactory.get()) {
    externallyPreconditionedLinearSolveUseCases(
      *A,*lowsFactory,*precFactory,false,true,*out
      );
  }

#else // SUN_CXX

    if(out.get()) *out << "\nTest failed since is was not even compiled since SUN_CXX was defined!\n";
    success = false;

#endif // SUN_CXX

  }
  catch( const std::exception &excpt ) {
    std::cerr << "\n*** Caught standard exception : " << excpt.what() << std::endl;
    success = false;
  }

  return success;

}
