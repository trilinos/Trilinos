//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
  int ierr = 0, i;

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
  
  // Begin LOCA Solver ************************************

  // Create parameter list
  NOX::Parameter::List paramList;

  // Create LOCA sublist
  NOX::Parameter::List& locaParamsList = paramList.sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  NOX::Parameter::List& locaStepperList = locaParamsList.sublist("Stepper");
  //locaStepperList.setParameter("Continuation Method", "Natural");
  locaStepperList.setParameter("Continuation Method", "Arc Length");
  locaStepperList.setParameter("Continuation Parameter", "Right BC");
  //locaStepperList.setParameter("Continuation Parameter", "Nonlinear Factor");
  locaStepperList.setParameter("Initial Value", 0.1);
  locaStepperList.setParameter("Max Value", 100.0);
  locaStepperList.setParameter("Min Value", 0.05);
  locaStepperList.setParameter("Max Steps", 30);
  locaStepperList.setParameter("Max Nonlinear Iterations", 15);
  locaStepperList.setParameter("Enable Arc Length Scaling", true);
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
  NOX::Parameter::List& aList = locaStepperList.sublist("Anasazi");
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
  stepSizeList.setParameter("Initial Step Size", 0.1);
  stepSizeList.setParameter("Min Step Size", 1.0e-3);
  stepSizeList.setParameter("Max Step Size", 2000.0);
  stepSizeList.setParameter("Aggressiveness", 0.1);
  stepSizeList.setParameter("Failed Step Reduction Factor", 0.5);
  stepSizeList.setParameter("Successful Step Increase Factor", 1.26); // for constant

  // Set the LOCA Utilities
  NOX::Parameter::List& locaUtilsList = locaParamsList.sublist("Utilities");
  locaUtilsList.setParameter("MyPID", MyPID);
  locaUtilsList.setParameter("Output Information", 
			     LOCA::Utils::Warning +
			     LOCA::Utils::StepperIteration +
			     LOCA::Utils::StepperDetails +
			     LOCA::Utils::Solver +
			     LOCA::Utils::SolverDetails +
			     LOCA::Utils::Parameters);

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  NOX::Parameter::List& nlParams = paramList.sublist("NOX");
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");

  // Create the NOX printing parameter list
  NOX::Parameter::List& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.setParameter("MyPID", MyPID); 
  nlPrintParams.setParameter("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning);

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
  //newParams.setParameter("Forcing Term Method", "Type 1");
  //newParams.setParameter("Forcing Term Method", "Type 2");
  //newParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-4);
  //newParams.setParameter("Forcing Term Maximum Tolerance", 0.1);

  // Create the "Linear Solver" sublist for the "Direction" sublist
  NOX::Parameter::List& lsParams = newParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 100);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Output Frequency", 50);    
  lsParams.setParameter("Scaling", "None");             
  //lsParams.setParameter("Scaling", "Row Sum");          
  //lsParams.setParameter("Preconditioning", "None");   
  lsParams.setParameter("Preconditioning", "AztecOO: Jacobian Matrix");   
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

  // Create the loca vector
  NOX::Epetra::Vector locaSoln(soln);

  // Create the Group
  LOCA::Epetra::Group grp(nlPrintParams, lsParams, interface, pVector, locaSoln, A);
  //NOX::Epetra::Group grp(nlPrintParams, lsParams, interface, soln, A); 
  //NOX::Epetra::Group grp(nlPrintParams, lsParams, interface, soln, A, Prec); 
  grp.computeF();

  // Create the Solver convergence test
  //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  NOX::StatusTest::NormF wrms(1.0e-8);
  NOX::StatusTest::MaxIters maxiters(searchParams.getParameter("Max Iters", 10));
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR);
  combo.addStatusTest(wrms);
  combo.addStatusTest(maxiters);

  // Create the stepper  
  LOCA::Stepper stepper(grp, combo, paramList);
  LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

  if (status != LOCA::Abstract::Iterator::Finished)
    if (MyPID==0) 
      cout << "Stepper failed to converge!" << endl;

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(stepper.getSolutionGroup());
  const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  // End Nonlinear Solver **************************************

  // Output the parameter list
  if (LOCA::Utils::doPrint(LOCA::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    stepper.getParameterList().print(cout);
    cout << endl;
  }

  // Print solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = soln.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",MyPID);
  ifp = fopen(file_name, "w");
  for (i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E\n", soln.Map().MinMyGID()+i, finalSolution[i]);
  fclose(ifp);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
