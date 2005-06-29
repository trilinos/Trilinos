// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ******
                                                                     
// 1D Finite Element Test Problem
/* Solves continuation problem (Parameter c="Right BC")
 *
 * d2u 
 * --- + a * u**3 = 0
 * dx2
 *
 * subject to @ x=0, u=b
 * subject to @ x=1, u=c
 */

// LOCA Objects
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX_TestCompare.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "Tcubed_FiniteElementProblem.H"              

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0;
  int MyPID;

  try {
  
    // scale factor to test arc-length scaling
    double scale = 1.0;

    // Initialize MPI
#ifdef HAVE_MPI
    MPI_Init(&argc,&argv);
#endif

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
    Epetra_SerialComm Comm;
#endif

    // Get the process ID and the total number of processors
    MyPID = Comm.MyPID();
    int NumProc = Comm.NumProc();

    // Check for verbose output
    bool verbose = false;
    if (argc>1) 
      if (argv[1][0]=='-' && argv[1][1]=='v') 
	verbose = true;

    // Get the number of elements from the command line
    int NumGlobalElements = 0;
    if ((argc > 2) && (verbose))
      NumGlobalElements = atoi(argv[2]) + 1;
    else if ((argc > 1) && (!verbose))
      NumGlobalElements = atoi(argv[1]) + 1;
    else 
      NumGlobalElements = 101;

    // The number of unknowns must be at least equal to the 
    // number of processors.
    if (NumGlobalElements < NumProc) {
      cout << "numGlobalBlocks = " << NumGlobalElements 
	   << " cannot be < number of processors = " << NumProc << endl;
      exit(1);
    }

    // Create the FiniteElementProblem class.  This creates all required
    // Epetra objects for the problem and allows calls to the 
    // function (RHS) and Jacobian evaluation routines.
    Tcubed_FiniteElementProblem Problem(NumGlobalElements, Comm, scale);

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();

    // Initialize Solution
    soln.PutScalar(1.0);
  
    // Begin LOCA Solver ************************************

    // Create parameter list
    Teuchos::RefCountPtr<NOX::Parameter::List> paramList = 
      Teuchos::rcp(new NOX::Parameter::List);
  
    // Create LOCA sublist
    NOX::Parameter::List& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    NOX::Parameter::List& locaStepperList = locaParamsList.sublist("Stepper");
    //locaStepperList.setParameter("Continuation Method", "Natural");
    locaStepperList.setParameter("Continuation Method", "Arc Length");
    locaStepperList.setParameter("Bordered Solver Method", "Householder");
    locaStepperList.setParameter("Continuation Parameter", "Right BC");
    //locaStepperList.setParameter("Continuation Parameter", "Nonlinear Factor");
    locaStepperList.setParameter("Initial Value", 0.1/scale);
    locaStepperList.setParameter("Max Value", 100.0/scale);
    locaStepperList.setParameter("Min Value", 0.05/scale);
    locaStepperList.setParameter("Max Steps", 30);
    locaStepperList.setParameter("Max Nonlinear Iterations", 15);
    locaStepperList.setParameter("Enable Arc Length Scaling", false);
    locaStepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5);
    locaStepperList.setParameter("Max Arc Length Parameter Contribution", 0.7);
    locaStepperList.setParameter("Initial Scale Factor", 1.0);
    locaStepperList.setParameter("Min Scale Factor", 1.0e-8);
    locaStepperList.setParameter("Enable Tangent Factor Step Size Scaling",false);
    locaStepperList.setParameter("Min Tangent Factor", 0.8);
    locaStepperList.setParameter("Tangent Factor Exponent",1.5);

    // Create bifurcation sublist
    NOX::Parameter::List& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.setParameter("Method", "None");

    // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
    locaStepperList.setParameter("Compute Eigenvalues",true);
    NOX::Parameter::List& aList = locaStepperList.sublist("Eigensolver");
    aList.setParameter("Method", "Anasazi");
    aList.setParameter("Block Size", 1);
    aList.setParameter("Arnoldi Size", 10);
    aList.setParameter("NEV", 3);
    aList.setParameter("Tol", 2.0e-7);
    aList.setParameter("Convergence Check", 1);
    aList.setParameter("Restarts",2);
    aList.setParameter("Frequency",1);
    aList.setParameter("Debug Level",0);
  
    // Create predictor sublist
    NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
    //predictorList.setParameter("Method", "Constant");
    predictorList.setParameter("Method", "Tangent");
    //predictorList.setParameter("Method", "Secant");

    // Create step size sublist
    NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.setParameter("Method", "Constant");
    stepSizeList.setParameter("Method", "Adaptive");
    stepSizeList.setParameter("Initial Step Size", 0.1/scale);
    stepSizeList.setParameter("Min Step Size", 1.0e-3/scale);
    stepSizeList.setParameter("Max Step Size", 2000.0/scale);
    stepSizeList.setParameter("Aggressiveness", 0.1);
    stepSizeList.setParameter("Failed Step Reduction Factor", 0.5);
    stepSizeList.setParameter("Successful Step Increase Factor", 1.26); // for constant

    // Set the LOCA Utilities
    NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
    locaUtilsList.setParameter("MyPID", MyPID);
    if (verbose)
      locaUtilsList.setParameter("Output Information", 
				 LOCA::Utils::Error +
				 LOCA::Utils::Warning +
				 LOCA::Utils::StepperIteration +
				 LOCA::Utils::StepperDetails +
				 LOCA::Utils::Solver +
				 LOCA::Utils::SolverDetails +
				 LOCA::Utils::Parameters);
    else
      locaUtilsList.setParameter("Output Information", LOCA::Utils::Error);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    NOX::Parameter::List& nlParams = paramList->sublist("NOX");
    nlParams.setParameter("Nonlinear Solver", "Line Search Based");

    // Create the NOX printing parameter list
    NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.setParameter("MyPID", MyPID); 
    if (verbose)
      nlPrintParams.setParameter("Output Information", 
				 NOX::Utils::OuterIteration + 
				 NOX::Utils::OuterIterationStatusTest + 
				 NOX::Utils::InnerIteration +
				 NOX::Utils::Parameters + 
				 NOX::Utils::Details + 
				 NOX::Utils::Warning +
				 NOX::Utils::TestDetails + 
				 NOX::Utils::Error);
    else
      nlPrintParams.setParameter("Output Information", NOX::Utils::Error);

    // Create the "Line Search" sublist for the "Line Search Based" solver
    NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
    searchParams.setParameter("Method", "Full Step");
    searchParams.setParameter("Max Iters", 7);
    searchParams.setParameter("Default Step", 1.0000);
    searchParams.setParameter("Recovery Step", 0.0001);
    searchParams.setParameter("Minimum Step", 0.0001);

    // Create the "Direction" sublist for the "Line Search Based" solver
    NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
    NOX::Parameter::List& newParams = dirParams.sublist("Newton");
    dirParams.setParameter("Method", "Newton");
    newParams.setParameter("Forcing Term Method", "Constant");

    // Create the "Linear Solver" sublist for the "Direction" sublist
    NOX::Parameter::List& lsParams = newParams.sublist("Linear Solver");
    lsParams.setParameter("Aztec Solver", "GMRES");  
    lsParams.setParameter("Max Iterations", 100);  
    lsParams.setParameter("Tolerance", 1e-4);
    if (verbose)
      lsParams.setParameter("Output Frequency", 1);
    else
      lsParams.setParameter("Output Frequency", 0);
    lsParams.setParameter("Scaling", "None");             
    lsParams.setParameter("Preconditioner", "Ifpack");
    //lsParams.setParameter("Preconditioner", "AztecOO");
    //lsParams.setParameter("Jacobian Operator", "Matrix-Free");
    //lsParams.setParameter("Preconditioner Operator", "Finite Difference");
    lsParams.setParameter("Aztec Preconditioner", "ilut"); 
    //lsParams.setParameter("Overlap", 2);   
    //lsParams.setParameter("Fill Factor", 2.0); 
    //lsParams.setParameter("Drop Tolerance", 1.0e-12);

    // Create and initialize the parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("Nonlinear Factor",1.0);
    pVector.addParameter("Left BC", 0.0);
    pVector.addParameter("Right BC", 0.1);

    // Create the interface between the test problem and the nonlinear solver
    // This is created by the user using inheritance of the abstract base class:
    // NLS_PetraGroupInterface
    Problem_Interface interface(Problem);

    // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner by 
    // uncommenting one or more of the following lines:
    // 1. User supplied (Epetra_RowMatrix)
    Epetra_RowMatrix& A = Problem.getJacobian();
    // 2. Matrix-Free (Epetra_Operator)
    //NOX::Epetra::MatrixFree A(interface, soln);
    // 3. Finite Difference (Epetra_RowMatrix)
    //NOX::Epetra::FiniteDifference A(interface, soln);
    // 4. Jacobi Preconditioner
    //NOX::Epetra::JacobiPreconditioner Prec(soln);

    // Create the linear systems
    NOX::EpetraNew::LinearSystemAztecOO linsys(nlPrintParams, lsParams,
					       interface, interface,
					       A, soln);
    //   NOX::EpetraNew::LinearSystemAztecOO linsys(nlPrintParams, lsParams,
    // 					     interface, soln);

    // Create the loca vector
    NOX::Epetra::Vector locaSoln(soln);

    // Create the Group
    Teuchos::RefCountPtr<LOCA::EpetraNew::Group> grp = 
      Teuchos::rcp(new LOCA::EpetraNew::Group(nlPrintParams, interface, 
					      locaSoln, linsys, pVector));
    grp->computeF();

    // Create the Solver convergence test
    NOX::StatusTest::NormF wrms(1.0e-8);
    NOX::StatusTest::MaxIters maxiters(searchParams.getParameter("Max Iters", 10));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(wrms);
    combo->addStatusTest(maxiters);

    // Create the Epetra Factory
    Teuchos::RefCountPtr<LOCA::Epetra::Factory> factory = 
      Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create the stepper  
    LOCA::NewStepper stepper(grp, combo, paramList, factory);
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      ierr = 1;
      if (LOCA::Utils::doPrint(LOCA::Utils::Error))
	cout << "Stepper failed to converge!" << endl;
    }

    // Get the final solution from the stepper
    Teuchos::RefCountPtr<const LOCA::EpetraNew::Group> finalGroup = 
      Teuchos::rcp_dynamic_cast<const LOCA::EpetraNew::Group>(stepper.getSolutionGroup());
    const NOX::Epetra::Vector& finalSolution = 
      dynamic_cast<const NOX::Epetra::Vector&>(finalGroup->getX());

    // Output the parameter list
    if (LOCA::Utils::doPrint(LOCA::Utils::Parameters)) {
      cout << endl << "Final Parameters" << endl
	   << "****************" << endl;
      stepper.getParameterList()->print(cout);
      cout << endl;
    }

    // Check some statistics on the solution
    NOX::Utils utils(nlPrintParams);
    NOX::TestCompare testCompare(cout, utils);
  
    if (utils.isPrintProcessAndType(NOX::Utils::TestDetails))
      cout << endl << "***** Checking solutions statistics *****" << endl;

    // Check number of steps
    int numSteps = stepper.getStepNumber();
    int numSteps_expected = 20;
    ierr += testCompare.testValue(numSteps, numSteps_expected, 0.0,
				  "number of continuation steps",
				  NOX::TestCompare::Absolute);

    // Check number of failed steps
    int numFailedSteps = stepper.getNumFailedSteps();
    int numFailedSteps_expected = 0;
    ierr += testCompare.testValue(numFailedSteps, numFailedSteps_expected, 0.0,
				  "number of failed continuation steps",
				  NOX::TestCompare::Absolute);

    // Check final value of continuation parameter
    double right_bc_final = finalGroup->getParam("Right BC");
    double right_bc_expected = 0.05;
    ierr += testCompare.testValue(right_bc_final, right_bc_expected, 1.0e-14,
				  "final value of continuation parameter", 
				  NOX::TestCompare::Relative);
 
    // Check norm of solution
    double norm_x = finalSolution.norm();
    double norm_x_expected = 25.00498021;
    ierr += testCompare.testValue(norm_x, norm_x_expected, 1.0e-7,
				  "norm of final solution",
				  NOX::TestCompare::Relative);

  }

  catch (const char *s) {
    cout << s << endl;
    ierr = 1;
  }
  catch (...) {
    cout << "Caught unknown exception!" << endl;
    ierr = 1;
  }

  if (MyPID == 0) {
    if (ierr == 0)
      cout << "All tests passed!" << endl;
    else
      cout << ierr << " test(s) failed!" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
  return ierr ;
}
