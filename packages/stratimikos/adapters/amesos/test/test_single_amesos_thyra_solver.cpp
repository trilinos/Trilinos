// @HEADER
// ***********************************************************************
// 
//         Stratimikos: Thyra-based strategies for linear solvers
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "test_single_amesos_thyra_solver.hpp"

#ifndef SUN_CXX

#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveFactoryExamples.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_SerialComm.h"

#endif // SUN_CXX

bool Thyra::test_single_amesos_thyra_solver(
  const std::string                       matrixFile
  ,Teuchos::ParameterList                 *amesosLOWSFPL
  ,const bool                             testTranspose_in
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxError
  ,const double                           maxResid
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,Teuchos::FancyOStream                  *out_arg
  )
{

  using Teuchos::RCP;
  using Teuchos::OSTab;

  bool result, success = true;

  RCP<Teuchos::FancyOStream>
    out = Teuchos::rcp(out_arg,false);

#ifndef SUN_CXX

  if(out.get()) {
    *out
      << "\n***"
      << "\n*** Testing Thyra::AmesosLinearOpWithSolveFactory (and Thyra::AmesosLinearOpWithSolve)"
      << "\n***\n"
      << "\nEchoing input options:"
      << "\n  matrixFile             = " << matrixFile
      << "\n  testTranspose          = " << testTranspose_in
      << "\n  numRandomVectors       = " << numRandomVectors
      << "\n  maxFwdError            = " << maxFwdError
      << "\n  maxError               = " << maxError
      << "\n  maxResid               = " << maxResid
      << "\n  showAllTests           = " << showAllTests
      << "\n  dumpAll                = " << dumpAll
      << std::endl;
    if(amesosLOWSFPL) {
      OSTab tab(out);
      *out
        << "amesosLOWSFPL:\n";
      amesosLOWSFPL->print(OSTab(out).o(),0,true);
    }
  }
  
  if(out.get()) *out << "\nA) Reading in an epetra matrix A from the file \'"<<matrixFile<<"\' ...\n";
  
  Epetra_SerialComm comm;
  RCP<Epetra_CrsMatrix> epetra_A;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );

  RCP<const LinearOpBase<double> >
    A = Thyra::epetraLinearOp(epetra_A);

  if(out.get() && dumpAll) *out << "\ndescribe(A) =\n" << describe(*A,Teuchos::VERB_EXTREME);

  if(out.get()) *out << "\nB) Creating a AmesosLinearOpWithSolveFactory object lowsFactory ...\n";

  RCP<LinearOpWithSolveFactoryBase<double> >
    lowsFactory = Teuchos::rcp(new AmesosLinearOpWithSolveFactory());

  lowsFactory->setParameterList(Teuchos::rcp(amesosLOWSFPL,false));

  if(out.get()) {
    *out << "\nlowsFactory.description() = " << lowsFactory->description() << std::endl;
    *out << "\nlowsFactory.getValidParameters() =\n";
    lowsFactory->getValidParameters()->print(OSTab(out).o(),0,true,false);
  }

  if(out.get()) *out << "\nC) Creating a AmesosLinearOpWithSolve object nsA ...\n";

  RCP<LinearOpWithSolveBase<double> >
    nsA = lowsFactory->createOp();

  initializeOp<double>(*lowsFactory, A, nsA.ptr());

  bool testTranspose = testTranspose_in;
  if (testTranspose_in && !Thyra::solveSupports(*nsA, Thyra::CONJTRANS)) {
    if(out.get()) *out
      << "\nChanging testTranspose=false since "
      << nsA->description() << " does not support adjoint solves!\n";
    testTranspose = false;
  }
  
  if(out.get()) *out << "\nD) Testing the LinearOpBase interface of nsA ...\n";

  Thyra::seed_randomize<double>(0);

  LinearOpTester<double> linearOpTester;
  linearOpTester.check_adjoint(testTranspose);
  linearOpTester.num_random_vectors(numRandomVectors);
  linearOpTester.set_all_error_tol(maxFwdError);
  linearOpTester.set_all_warning_tol(1e-2*maxFwdError);
  linearOpTester.show_all_tests(showAllTests);
  linearOpTester.dump_all(dumpAll);
  result = linearOpTester.check(*nsA,out());
  if(!result) success = false;

  if(out.get()) *out << "\nE) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
  LinearOpWithSolveTester<double> linearOpWithSolveTester;
  linearOpWithSolveTester.turn_off_all_tests();
  linearOpWithSolveTester.check_forward_default(true);
  linearOpWithSolveTester.forward_default_residual_error_tol(1.1*maxResid);
  linearOpWithSolveTester.forward_default_residual_warning_tol(2.0*maxResid);
  linearOpWithSolveTester.check_forward_residual(true);
  linearOpWithSolveTester.forward_residual_solve_tol(maxResid);
  linearOpWithSolveTester.forward_residual_slack_error_tol(1e-1*maxResid);
  linearOpWithSolveTester.forward_residual_slack_warning_tol(maxResid);
  linearOpWithSolveTester.check_adjoint_default(testTranspose);
  linearOpWithSolveTester.adjoint_default_residual_error_tol(1.1*maxResid);
  linearOpWithSolveTester.adjoint_default_residual_warning_tol(2.0*maxResid);
  linearOpWithSolveTester.check_adjoint_residual(testTranspose);
  linearOpWithSolveTester.adjoint_residual_solve_tol(maxResid);
  linearOpWithSolveTester.adjoint_residual_slack_error_tol(1e-1*maxResid);
  linearOpWithSolveTester.adjoint_residual_slack_warning_tol(maxResid);
  linearOpWithSolveTester.num_random_vectors(numRandomVectors);
  linearOpWithSolveTester.show_all_tests(showAllTests);
  linearOpWithSolveTester.dump_all(dumpAll);

  LinearOpWithSolveTester<double> adjLinearOpWithSolveTester(linearOpWithSolveTester);
  adjLinearOpWithSolveTester.check_forward_default(testTranspose);
  adjLinearOpWithSolveTester.check_forward_residual(testTranspose);

  result = linearOpWithSolveTester.check(*nsA,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nF) Uninitialize the matrix object nsA, scale the epetra_A object by 2.5, and then refactor nsA with epetra_A ...\n";

  uninitializeOp<double>(*lowsFactory, nsA.ptr()); // Optional call but a good idea if changing the operator
  epetra_A->Scale(2.5);
  initializeOp<double>(*lowsFactory, A, nsA.ptr());
  
  if(out.get()) *out << "\nG) Testing the LinearOpBase interface of nsA ...\n";

  Thyra::seed_randomize<double>(0);

  result = linearOpTester.check(*nsA, out());
  if(!result) success = false;

  if(out.get()) *out << "\nH) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
  result = linearOpWithSolveTester.check(*nsA,out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nI) Uninitialize the matrix object nsA, create a scaled (by 2.5) copy  epetra_A2 of epetra_A, and then refactor nsA with epetra_A2 ...\n";

  RCP<Epetra_CrsMatrix>
    epetra_A2 = Teuchos::rcp(new Epetra_CrsMatrix(*epetra_A));
  epetra_A2->Scale(2.5);
  RCP<const LinearOpBase<double> >
    A2 = Thyra::epetraLinearOp(epetra_A2);
  initializeOp<double>(*lowsFactory, A2, nsA.ptr());
  
  if(out.get()) *out << "\nJ) Testing the LinearOpBase interface of nsA ...\n";

  Thyra::seed_randomize<double>(0);

  result = linearOpTester.check(*nsA, out());
  if(!result) success = false;

  if(out.get()) *out << "\nK) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
  result = linearOpWithSolveTester.check(*nsA, out.get());
  if(!result) success = false;

  if(out.get()) *out << "\nL) Create an implicitly scaled (by 2.5) and transposed matrix A3 = scale(2.5,transpose(A)) and initialize nsA2 ...\n";

  RCP<const LinearOpBase<double> >
    A3 = scale<double>(2.5,Thyra::transpose<double>(A));
  RCP<LinearOpWithSolveBase<double> >
    nsA2 = linearOpWithSolve(*lowsFactory,A3);
  
  if(out.get()) *out << "\nM) Testing the LinearOpBase interface of nsA2 ...\n";

  Thyra::seed_randomize<double>(0);

  result = linearOpTester.check(*nsA2,out());
  if(!result) success = false;

  if(out.get()) *out << "\nN) Testing the LinearOpWithSolveBase interface of nsA2 ...\n";
    
  result = adjLinearOpWithSolveTester.check(*nsA2,out.get());
  if(!result) success = false;
  
  if(out.get()) *out << "\nO) Testing that LinearOpBase interfaces of transpose(nsA) == nsA2 ...\n";

  result = linearOpTester.compare(
    *transpose(Teuchos::rcp_implicit_cast<const LinearOpBase<double> >(nsA)),*nsA2
    ,out()
    );
  if(!result) success = false;

  if(out.get()) *out << "\nP) Running example use cases ...\n";

  nonExternallyPreconditionedLinearSolveUseCases(
    *A,*lowsFactory,true,*out
    );

  if(out.get()) *out << "\nQ) Creating a DefaultInverseLinearOp object from nsA and testing the LinearOpBase interface ...\n";

  RCP<const LinearOpBase<double> >
    invA = inverse<double>(nsA.getConst());

  result = linearOpTester.check(*invA,out());
  if(!result) success = false;

#else // SUN_CXX
  
  if(out.get()) *out << "\nTest failed since is was not even compiled since SUN_CXX was defined!\n";
  success = false;

#endif // SUN_CXX

  return success;

}
