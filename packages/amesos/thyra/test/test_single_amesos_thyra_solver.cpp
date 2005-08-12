/*
// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "test_single_amesos_thyra_solver.hpp"
#include "Thyra_AmesosLinearOpWithSolveFactory.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_SerialComm.h"

bool Thyra::test_single_amesos_thyra_solver(
  const std::string                       matrixFile
  ,const Amesos::ESolverType              solverType
  ,const Amesos::ERefactorizationPolicy   refactorizationPolicy
  ,const bool                             testTranspose
  ,const int                              numRandomVectors
  ,const double                           maxFwdError
  ,const double                           maxError
  ,const double                           maxResid
  ,const bool                             showAllTests
  ,const bool                             dumpAll
  ,std::ostream                           *out
  )
{
  bool result, success = true;

  if(out) {
    *out << "\n***"
         << "\n*** Testing Thyra::AmesosLinearOpWithSolveFactory (and Thyra::AmesosLinearOpWithSolve)"
         << "\n***\n"
         << "\nEchoing input options:"
         << "\n  matrixFile             = " << matrixFile
         << "\n  solverType             = " << toString(solverType)
         << "\n  refactorizationPolicy  = " << toString(refactorizationPolicy)
         << "\n  testTranspose          = " << testTranspose
         << "\n  numRandomVectors       = " << numRandomVectors
         << "\n  maxFwdError            = " << maxFwdError
         << "\n  maxError               = " << maxError
         << "\n  maxResid               = " << maxResid
         << "\n  showAllTests           = " << showAllTests
         << "\n  dumpAll                = " << dumpAll
         << std::endl;
  }
  
  if(out) *out << "\nA) Reading in an epetra matrix A from the file \'"<<matrixFile<<"\' ...\n";
  
  const std::string indentSpacer = "  ";
  
  Epetra_SerialComm comm;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> epetra_A;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A );

  Teuchos::RefCountPtr<LinearOpBase<double> > A = Teuchos::rcp(new EpetraLinearOp(epetra_A));

  if(out && dumpAll) *out << "\ndescribe(A) =\n" << describe(*A,Teuchos::VERB_EXTREME,indentSpacer,indentSpacer);

  if(out) *out << "\nB) Creating a AmesosLinearOpWithSolveFactory object opFactory ...\n";

  Teuchos::RefCountPtr<const LinearOpWithSolveFactoryBase<double> >
    opFactory = Teuchos::rcp(new AmesosLinearOpWithSolveFactory(solverType,refactorizationPolicy));

  if(out) *out << "\nC) Creating a AmesosLinearOpWithSolve object nsA ...\n";

  Teuchos::RefCountPtr<LinearOpWithSolveBase<double> >
    nsA = opFactory->createOp();

  opFactory->initializeOp( A, &*nsA );

  if(out) *out << "\nD) Testing the LinearOpBase interface of nsA ...\n";

  LinearOpTester<double> linearOpTester;
  linearOpTester.check_adjoint(testTranspose);
  linearOpTester.num_random_vectors(numRandomVectors);
  linearOpTester.set_all_error_tol(maxFwdError);
  linearOpTester.set_all_warning_tol(1e-2*maxFwdError);
  linearOpTester.show_all_tests(showAllTests);
  linearOpTester.dump_all(dumpAll);
  result = linearOpTester.check(*nsA,out,indentSpacer,indentSpacer);
  if(!result) success = false;

  if(out) *out << "\nE) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
  LinearOpWithSolveTester<double> linearOpWithSolveTester;
  linearOpWithSolveTester.turn_off_all_tests();
  linearOpWithSolveTester.check_forward_residual(true);
  linearOpWithSolveTester.forward_residual_solve_tol(maxResid);
  linearOpWithSolveTester.forward_residual_slack_error_tol(1e-1*maxResid);
  linearOpWithSolveTester.forward_residual_slack_warning_tol(maxResid);
  linearOpWithSolveTester.check_forward_solution_error(true);
  linearOpWithSolveTester.forward_solution_error_solve_tol(maxError);
  linearOpWithSolveTester.forward_solution_error_slack_error_tol(1e-1*maxError);
  linearOpWithSolveTester.forward_solution_error_slack_warning_tol(maxError);
  if(testTranspose) {
    linearOpWithSolveTester.check_adjoint_residual(true);
    linearOpWithSolveTester.adjoint_residual_solve_tol(maxResid);
    linearOpWithSolveTester.adjoint_residual_slack_error_tol(1e-1*maxResid);
    linearOpWithSolveTester.adjoint_residual_slack_warning_tol(maxResid);
    linearOpWithSolveTester.check_adjoint_solution_error(true);
    linearOpWithSolveTester.adjoint_solution_error_solve_tol(maxError);
    linearOpWithSolveTester.adjoint_solution_error_slack_error_tol(1e-1*maxError);
    linearOpWithSolveTester.adjoint_solution_error_slack_warning_tol(maxError);
  }
  linearOpWithSolveTester.num_random_vectors(numRandomVectors);
  linearOpWithSolveTester.show_all_tests(showAllTests);
  linearOpWithSolveTester.dump_all(dumpAll);
  result = linearOpWithSolveTester.check(*nsA,out,indentSpacer,indentSpacer);
  if(!result) success = false;

  if(out) *out << "\nF) Uninitialize the matrix object, scale the matrix by 2.5, and then refactor it ...\n";

  opFactory->uninitializeOp(&*nsA );
  epetra_A->Scale(2.5);
  opFactory->initializeOp(A,&*nsA );

  if(out) *out << "\nG) Testing the LinearOpBase interface of nsA ...\n";

  result = linearOpTester.check(*nsA,out,indentSpacer,indentSpacer);
  if(!result) success = false;

  if(out) *out << "\nH) Testing the LinearOpWithSolveBase interface of nsA ...\n";
    
  result = linearOpWithSolveTester.check(*nsA,out,indentSpacer,indentSpacer);
  if(!result) success = false;

  return success;

}
