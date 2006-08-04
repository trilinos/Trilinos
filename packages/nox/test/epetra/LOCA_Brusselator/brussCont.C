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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

// 1D Finite Element Brusselator Test Problem

/* Solves the nonlinear equation:
 *
 * dT       d2T    
 * --- - D1 --- - alpha + (beta+1)*T - C*T**2 = 0
 * dt       dx2   
 *
 * T(t,0) = T(t,1) = alpha = 0.6
 * T(0,x) = alpha + sinusoidal perturbation
 *
 *
 * dC       d2C    
 * --- - D2 --- - beta*T + C*T**2 = 0
 * dt       dx2   
 *
 * C(t,0) = C(t,1) = beta / alpha = 2.0 / 0.6
 * C(0,x) = beta / alpha + sinusoidal perturbation
 *
 * and
 *      D1 = D2 = 0.025
 *
 * with d representing partial differentiation.
 */

// NOX Objects
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX_Common.H"
#include "NOX_TestCompare.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// Added to allow timings
#include "Epetra_Time.h"

// Headers needed for FD coloring 
#include <vector> 

// User's application specific files 
#include "Problem_Interface.H" // Interface file to NOX
#include "Brusselator.H"              

int main(int argc, char *argv[])
{
  int ierr = 0;
  int MyPID = 0;
  double alpha = 0.25;
  double beta = 0.0;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  int maxNewtonIters = 10;

  try {

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
    int NumGlobalNodes = 0;
    if ((argc > 2) && (verbose))
      NumGlobalNodes = atoi(argv[2]) + 1;
    else if ((argc > 1) && (!verbose))
      NumGlobalNodes = atoi(argv[1]) + 1;
    else 
      NumGlobalNodes = 101;

    // The number of unknowns must be at least equal to the 
    // number of processors.
    if (NumGlobalNodes < NumProc) {
      std::cout << "numGlobalNodes = " << NumGlobalNodes 
		<< " cannot be < number of processors = " << NumProc 
		<< std::endl;
      exit(1);
    }

    // Create the Brusselator problem class.  This creates all required
    // Epetra objects for the problem and allows calls to the 
    // function (F) and Jacobian evaluation routines.
    Brusselator Problem(NumGlobalNodes, Comm);

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();

    // Begin LOCA Solver ************************************

    // Create the top level parameter list

    Teuchos::RefCountPtr<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    //stepperList.set("Continuation Method", "Natural");
    stepperList.set("Continuation Method", "Arc Length");
    stepperList.set("Bordered Solver Method", "Householder");
    stepperList.set("Continuation Parameter", "beta");
    stepperList.set("Initial Value", beta);
    stepperList.set("Max Value", 2.0);
    stepperList.set("Min Value", 0.0);
    stepperList.set("Max Steps", 100);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);
    stepperList.set("Enable Arc Length Scaling", true);
    stepperList.set("Goal Arc Length Parameter Contribution", 0.5);
    stepperList.set("Max Arc Length Parameter Contribution", 0.7);
    stepperList.set("Initial Scale Factor", 1.0);
    stepperList.set("Min Scale Factor", 1.0e-8);
    stepperList.set("Enable Tangent Factor Step Size Scaling",false);
    stepperList.set("Min Tangent Factor", -1.0);
    stepperList.set("Tangent Factor Exponent",1.0);
    stepperList.set("Compute Eigenvalues",false);

     // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList = 
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "None");

    // Create predictor sublist
    Teuchos::ParameterList& predictorList = 
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Constant");
    //predictorList.set("Method", "Tangent");
    //predictorList.set("Method", "Secant");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    //stepSizeList.set("Method", "Constant");
    stepSizeList.set("Method", "Adaptive");
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 10.0);
    stepSizeList.set("Aggressiveness", 0.1);
    stepSizeList.set("Failed Step Reduction Factor", 0.5);
    stepSizeList.set("Successful Step Increase Factor", 1.26); // for constant

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
    // Set the nonlinear solver method
    nlParams.set("Nonlinear Solver", "Line Search Based");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
    printParams.set("MyPID", MyPID); 
    printParams.set("Output Precision", 3);
    printParams.set("Output Processor", 0);
     if (verbose)
      printParams.set("Output Information", 
		      NOX::Utils::OuterIteration + 
		      NOX::Utils::OuterIterationStatusTest + 
		      NOX::Utils::InnerIteration +
		      NOX::Utils::Parameters + 
		      NOX::Utils::Details + 
		      NOX::Utils::Warning +
		      NOX::Utils::TestDetails + 
		      NOX::Utils::Error + 
		      NOX::Utils::StepperIteration +
		      NOX::Utils::StepperDetails);
     else
       printParams.set("Output Information", NOX::Utils::Error);

    // Sublist for line search 
    Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
    searchParams.set("Method", "Full Step");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");
    

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");  
    lsParams.set("Max Iterations", 800);  
    lsParams.set("Tolerance", 1e-6);
    lsParams.set("Output Frequency", 50);    
    lsParams.set("Preconditioner", "Ifpack"); 
    //lsParams.set("Preconditioner", "AztecOO"); 
    //lsParams.set("Aztec Preconditioner", "ilu"); 
    //lsParams.set("Overlap", 2);  
    //lsParams.set("Graph Fill", 2); 
    //lsParams.set("Aztec Preconditioner", "ilut"); 
    //lsParams.set("Overlap", 2);   
    //lsParams.set("Fill Factor", 2);   
    //lsParams.set("Drop Tolerance", 1.0e-12);   
    //lsParams.set("Aztec Preconditioner", "Polynomial"); 
    //lsParams.set("Polynomial Order", 6); 

    // Create the interface between the test problem and the nonlinear solver
    Teuchos::RefCountPtr<Problem_Interface> interface = 
      Teuchos::rcp(new Problem_Interface(Problem));

    // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
    Teuchos::RefCountPtr<Epetra_RowMatrix> A = 
      Teuchos::rcp(&Problem.getJacobian(),false);


    // Use an Epetra Scaling object if desired
    Teuchos::RefCountPtr<Epetra_Vector> scaleVec = 
      Teuchos::rcp(new Epetra_Vector(soln));
    NOX::Epetra::Scaling scaling;
    scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);

    Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> iReq = interface;

    // Create the Linear System
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
							iReq, iJac, A, soln));
                                                        //&scaling);

    // Create initial guess
    NOX::Epetra::Vector initialGuess(Teuchos::rcp(&soln,false), 
				     NOX::Epetra::Vector::CreateView,
				     NOX::DeepCopy);

    // Create and initialize the parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("alpha",alpha);
    pVector.addParameter("beta",beta);
    pVector.addParameter("D1",D1);
    pVector.addParameter("D2",D2);

    // Create Epetra factory
    Teuchos::RefCountPtr<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create global data object
    Teuchos::RefCountPtr<LOCA::GlobalData> globalData = 
      LOCA::createGlobalData(paramList, epetraFactory);

    // Create the Group
    Teuchos::RefCountPtr<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams,
					   iReq, initialGuess, linSys, 
					   pVector));

    grp->computeF();

    // Create the convergence tests
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8, 
					   NOX::StatusTest::NormF::Unscaled));
    //NOX::StatusTest::NormF relresid(*grp.get(), 1.0e-2);
    //NOX::StatusTest::NormUpdate update(1.0e-5);
    //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
    //NOX::StatusTest::Combo converged(NOX::StatusTest::Combo::AND);
    //converged.addStatusTest(absresid);
    //converged.addStatusTest(relresid);
    //converged.addStatusTest(wrms);
    //converged.addStatusTest(update);
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(absresid);
    combo->addStatusTest(maxiters);

    // Create stepper
    LOCA::NewStepper stepper(globalData, grp, combo, paramList);

    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      ierr = 1;
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
	globalData->locaUtils->out() 
	  << "Stepper failed to converge!" << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RefCountPtr<const LOCA::Epetra::Group> finalGroup = 
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(stepper.getSolutionGroup());
    const NOX::Epetra::Vector& finalSolution = 
      dynamic_cast<const NOX::Epetra::Vector&>(finalGroup->getX());

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
      globalData->locaUtils->out() 
	<< std::endl << "Final Parameters" << std::endl
	<< "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    // Check some statistics on the solution
    NOX::TestCompare testCompare(globalData->locaUtils->out(), 
				 *(globalData->locaUtils));
  
    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out() 
	<< std::endl 
	<< "***** Checking solution statistics *****" 
	<< std::endl;

    // Check number of continuation steps
    int numSteps = stepper.getStepNumber();
    int numSteps_expected = 14;
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
    double beta_final = finalGroup->getParam("beta");
    double beta_expected = 2.0;
    ierr += testCompare.testValue(beta_final, beta_expected, 1.0e-14,
				  "final value of continuation parameter", 
				  NOX::TestCompare::Relative);

    // Check final of solution
    NOX::Epetra::Vector final_x_expected(finalSolution);
    int n = final_x_expected.getEpetraVector().MyLength()/2;
    for (int i=0; i<n; i++) {
      final_x_expected.getEpetraVector()[2*i] = alpha;
      final_x_expected.getEpetraVector()[2*i+1] = beta_final/alpha;
    }
    ierr += testCompare.testVector(finalSolution, final_x_expected, 
				   1.0e-10, 1.0e-10,
				   "value of final solution");

    destroyGlobalData(globalData);
  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    ierr = 1;
  }

  if (MyPID == 0) {
    if (ierr == 0)
      std::cout << "All tests passed!" << std::endl;
    else
      std::cout << ierr << " test(s) failed!" << std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

return ierr;
}
