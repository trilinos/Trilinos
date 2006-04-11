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
#include "FiniteElementProblem.H"              

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0;

  double nonlinear_factor = 1.0;
  double left_bc = 0.0;
  double right_bc = 2.07;

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
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Get the number of elements from the command line
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalElements = atoi(argv[1]) + 1;

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
  FiniteElementProblem Problem(NumGlobalElements, Comm);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Initialize Solution
  soln.PutScalar(1.0);

  // Create initial guess for the null vector of jacobian
  Teuchos::RefCountPtr<NOX::Abstract::Vector> nullVec = 
    Teuchos::rcp(new NOX::Epetra::Vector(soln));  
  nullVec->init(1.0);             // initial value 1.0
  
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
  locaStepperList.setParameter("Bordered Solver Method", "Nested");
  //locaStepperList.setParameter("Bordered Solver Method", "Householder");
  //locaStepperList.setParameter("Continuation Parameter", "Right BC");
  locaStepperList.setParameter("Continuation Parameter", "Nonlinear Factor");
  locaStepperList.setParameter("Initial Value", nonlinear_factor);
  locaStepperList.setParameter("Max Value", 2.0);
  locaStepperList.setParameter("Min Value", 0.05);
  locaStepperList.setParameter("Max Steps", 20);
  locaStepperList.setParameter("Max Nonlinear Iterations", 15);
  locaStepperList.setParameter("Enable Arc Length Scaling", true);
  locaStepperList.setParameter("Goal Arc Length Parameter Contribution", 0.5);
  locaStepperList.setParameter("Max Arc Length Parameter Contribution", 0.7);
  locaStepperList.setParameter("Initial Scale Factor", 1.0);
  locaStepperList.setParameter("Min Scale Factor", 1.0e-8);
  locaStepperList.setParameter("Enable Tangent Factor Step Size Scaling",false);
  locaStepperList.setParameter("Min Tangent Factor", 0.8);
  locaStepperList.setParameter("Tangent Factor Exponent",1.5);

  NOX::Parameter::List& nestedList = 
    locaStepperList.sublist("Nested Bordered Solver");
  nestedList.setParameter("Bordered Solver Method", "Householder");
  nestedList.setParameter("Include UV In Preconditioner", true);
  //nestedList.setParameter("Use P For Preconditioner", true);

  // Create bifurcation sublist
  NOX::Parameter::List& bifurcationList = 
    locaParamsList.sublist("Bifurcation");
  bifurcationList.setParameter("Type", "Turning Point");
  //bifurcationList.setParameter("Formulation", "Moore-Spence");
  bifurcationList.setParameter("Formulation", "Minimally Augmented");
  //bifurcationList.setParameter("Constraint Method", "Modified");
  bifurcationList.setParameter("Solver Method", "Phipps Bordering");
  //bifurcationList.setParameter("Solver Method", "Salinger Bordering");
  bifurcationList.setParameter("Bordered Solver Method", "Householder");
  //bifurcationList.setParameter("Bordered Solver Method", "Augmented");
  bifurcationList.setParameter("Bifurcation Parameter", "Right BC");
  //bifurcationList.setParameter("Bifurcation Parameter", "Nonlinear Factor");
  bifurcationList.setParameter("Length Normalization Vector", nullVec);
//   bifurcationList.setParameter("Initial Null Vector Computation",
// 			       "Solve df/dp");
  bifurcationList.setParameter("Initial Null Vector", nullVec);
  bifurcationList.setParameter("Initial A Vector", nullVec);
  bifurcationList.setParameter("Initial B Vector", nullVec);
  bifurcationList.setParameter("Update Null Vectors Every Continuation Step", 
			       true);
  bifurcationList.setParameter("Update Null Vectors Every Nonlinear Iteration",
			       false);
  bifurcationList.setParameter("Symmetric Jacobian", false);
  bifurcationList.setParameter("Include UV In Preconditioner", true);
  //bifurcationList.setParameter("Use P For Preconditioner", true);

  // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
  locaStepperList.setParameter("Compute Eigenvalues",false);
  NOX::Parameter::List& aList = locaStepperList.sublist("Anasazi");
  aList.setParameter("Block Size", 1);
  aList.setParameter("Arnoldi Size", 10);
  aList.setParameter("NEV", 3);
  aList.setParameter("Tol", 2.0e-7);
  aList.setParameter("Convergence Check", 1);
  aList.setParameter("Restarts",2);
  aList.setParameter("Frequency",2);
  aList.setParameter("Debug Level",1);
  
  // Create predictor sublist
  NOX::Parameter::List& predictorList = locaParamsList.sublist("Predictor");
  //predictorList.setParameter("Method", "Constant");
  //predictorList.setParameter("Method", "Tangent");
  predictorList.setParameter("Method", "Secant");

  // Create step size sublist
  NOX::Parameter::List& stepSizeList = locaParamsList.sublist("Step Size");
  //stepSizeList.setParameter("Method", "Constant");
  stepSizeList.setParameter("Method", "Adaptive");
  stepSizeList.setParameter("Initial Step Size", 0.1);
  stepSizeList.setParameter("Min Step Size", 1.0e-3);
  stepSizeList.setParameter("Max Step Size", 2000.0);
  stepSizeList.setParameter("Aggressiveness", 0.1);
  stepSizeList.setParameter("Failed Step Reduction Factor", 0.5);
  stepSizeList.setParameter("Successful Step Increase Factor", 1.26); // for constant

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  NOX::Parameter::List& nlParams = paramList->sublist("NOX");
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");

  // Create the NOX printing parameter list
  NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.setParameter("MyPID", MyPID); 
  nlPrintParams.setParameter("Output Precision", 6); 
  nlPrintParams.setParameter("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::LinearSolverDetails +
			     NOX::Utils::Warning + 
			     NOX::Utils::StepperIteration +
			     NOX::Utils::StepperDetails);

  // Create the "Line Search" sublist for the "Line Search Based" solver
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  searchParams.setParameter("Method", "Full Step");
  searchParams.setParameter("Max Iters", 10);
  searchParams.setParameter("Default Step", 1.0000);
  searchParams.setParameter("Recovery Step", 0.0001);
  searchParams.setParameter("Minimum Step", 0.0001);

  // Create the "Direction" sublist for the "Line Search Based" solver
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  NOX::Parameter::List& newParams = dirParams.sublist("Newton");
  dirParams.setParameter("Method", "Newton");
  newParams.setParameter("Forcing Term Method", "Constant");
  //newParams.setParameter("Forcing Term Method", "Type 1");
  //newParams.setParameter("Forcing Term Method", "Type 2");
  //newParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-4);
  //newParams.setParameter("Forcing Term Maximum Tolerance", 0.1);

  // Create the "Linear Solver" sublist for the "Direction" sublist
  NOX::Parameter::List& lsParams = newParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 200);  
  lsParams.setParameter("Tolerance", 1e-6);
  lsParams.setParameter("Output Frequency", 50);    
  //lsParams.setParameter("Scaling", "None");             
  //lsParams.setParameter("Scaling", "Row Sum");  
  lsParams.setParameter("Compute Scaling Manually", false);
  //lsParams.setParameter("Preconditioning", "None");  
  lsParams.setParameter("Preconditioner", "Ifpack");
  //lsParams.setParameter("Preconditioner", "New Ifpack");
  lsParams.setParameter("Ifpack Preconditioner", "ILU");
  Teuchos::ParameterList ifpackParams;
  ifpackParams.set("fact: level-of-fill", 1);
  lsParams.setParameter("Ifpack Teuchos Parameter List", &ifpackParams);
  //lsParams.setParameter("Preconditioning", "AztecOO: Jacobian Matrix");   
  //lsParams.setParameter("Preconditioning", "AztecOO: User RowMatrix"); 
  //lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");
  //lsParams.setParameter("Aztec Preconditioner", "ilu"); 
  //lsParams.setParameter("Overlap", 2);  
  //lsParams.setParameter("Graph Fill", 2); 
  //lsParams.setParameter("Aztec Preconditioner", "ilut"); 
  //lsParams.setParameter("Overlap", 2);   
  //lsParams.setParameter("Fill Factor", 2.0);   
  //lsParams.setParameter("Drop Tolerance", 1.0e-12);   
  //lsParams.setParameter("Aztec Preconditioner", "Polynomial"); 
  //lsParams.setParameter("Polynomial Order", 6); 

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("Nonlinear Factor",nonlinear_factor);
  pVector.addParameter("Left BC", left_bc);
  pVector.addParameter("Right BC", right_bc);

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  Teuchos::RefCountPtr<Problem_Interface> interface = 
    Teuchos::rcp(new Problem_Interface(Problem));
  Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> iReq = interface;
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = interface;
  
  // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
  Teuchos::RefCountPtr<Epetra_RowMatrix> Amat = 
    Teuchos::rcp(&Problem.getJacobian(),false);

  // Create scaling object
  Teuchos::RefCountPtr<NOX::Epetra::Scaling> scaling = Teuchos::null;
//   scaling = Teuchos::rcp(new NOX::Epetra::Scaling);
//   Teuchos::RefCountPtr<Epetra_Vector> scalingVector = 
//     Teuchos::rcp(new Epetra_Vector(soln.Map()));
//   //scaling->addRowSumScaling(NOX::Epetra::Scaling::Left, scalingVector);
//   scaling->addColSumScaling(NOX::Epetra::Scaling::Right, scalingVector);

  // Create transpose scaling object
  Teuchos::RefCountPtr<NOX::Epetra::Scaling> trans_scaling = Teuchos::null;
  trans_scaling = Teuchos::rcp(new NOX::Epetra::Scaling);
  Teuchos::RefCountPtr<Epetra_Vector> transScalingVector = 
    Teuchos::rcp(new Epetra_Vector(soln.Map()));
//   trans_scaling->addRowSumScaling(NOX::Epetra::Scaling::Right, 
// 				  transScalingVector);
  trans_scaling->addColSumScaling(NOX::Epetra::Scaling::Left, 
				  transScalingVector);
  //bifurcationList.setParameter("Transpose Scaling", trans_scaling);
  //bifurcationList.setParameter("Transpose Solver Method","Transpose Preconditioner");
  bifurcationList.setParameter("Transpose Solver Method","Explicit Transpose");
  //bifurcationList.setParameter("Transpose Solver Method","Left Preconditioning");

  // Create the linear systems
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linsys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams,
						      iReq, iJac, Amat, soln,
						      scaling));

  // Create the loca vector
  NOX::Epetra::Vector locaSoln(soln);

  // Create Epetra factory
  Teuchos::RefCountPtr<LOCA::Abstract::Factory> epetraFactory =
    Teuchos::rcp(new LOCA::Epetra::Factory);

  // Create global data object
  Teuchos::RefCountPtr<LOCA::GlobalData> globalData = 
    LOCA::createGlobalData(paramList, epetraFactory);

  // Create the Group
  Teuchos::RefCountPtr<LOCA::Epetra::Group> grp = 
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, iReq, 
					 locaSoln, linsys,
					 pVector));
  grp->computeF();

  // Create the Solver convergence test
  //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> wrms = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12));
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(searchParams.getParameter("Max Iters", 10)));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);
  
  // Create the stepper  
  LOCA::NewStepper stepper(globalData, grp, combo, paramList);
  LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();
  
  if (status != LOCA::Abstract::Iterator::Finished) {
    if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
      globalData->locaUtils->out() 
	<< "Stepper failed to converge!" << std::endl;
  }

  // Output the parameter list
  if (globalData->locaUtils->isPrintType(NOX::Utils::Parameters)) {
    globalData->locaUtils->out() 
      << std::endl << "Final Parameters" << std::endl
      << "****************" << std::endl;
    stepper.getParameterList()->print(globalData->locaUtils->out());
    globalData->locaUtils->out() << std::endl;
  }

  destroyGlobalData(globalData);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
