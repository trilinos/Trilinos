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
                                                                                
//  Simple 2 equation test for quadratic and cubic line searches 
//  from Dennis & Schnabel's book, chp 6.  The test problem is from
//  Example 6.5.1
/*  
 *    U0**2 + U1**2 - 2 = 0
 *    exp(U0-1) + U1**3 -2 = 0
 */

// NOX Library
#include "NOX.H"
#include "NOX_Epetra.H"

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
#include "DennisSchnabel.H"              

using namespace std;

int main(int argc, char *argv[])
{
  int i;

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

  bool verbose = false;
  if (argc > 1)
    if (argv[1][0]=='-' && argv[1][1]=='v')
      verbose = true;

  // Get the process ID and the total number of processors
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  int NumGlobalElements = 2;  // Hardcoded for D&S Example problem

  // A maximum of 2 procesors is allowed since there are only 2 equations
  if (NumProc >= 3) {
    cout << "ERROR: Maximum number of processors is 2!" << endl;
    exit(1);
  }

  // Create the Problem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  DennisSchnabel Problem(NumGlobalElements, Comm);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Initialize Solution
  if (MyPID==0) {
    soln[0]=2.0;
    if (NumProc==1) 
      soln[1]=0.5;
  } 
  else 
    soln[0]=0.5;

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  NOX::Parameter::List nlParams;

  // Set the nonlinear solver method
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", MyPID); 
  printParams.setParameter("Output Precision", 5);
  printParams.setParameter("Output Processor", 0);
  if ( verbose )
    printParams.setParameter("Output Information",
                     NOX::Utils::OuterIteration +
                     NOX::Utils::OuterIterationStatusTest +
                     NOX::Utils::InnerIteration +
                     NOX::Utils::Parameters +
                     NOX::Utils::Details +
                     NOX::Utils::Warning +
                     NOX::Utils::TestDetails);
  else
    printParams.setParameter("Output Information", NOX::Utils::Error +
                     NOX::Utils::TestDetails);

  NOX::Utils printing(printParams);

  // Identify the test problem
  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails))
    cout << "Starting epetra/DS6.5.1/DS_6_5_1.exe" << endl;

  // Identify processor information
#ifdef HAVE_MPI
  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    cout << "Parallel Run" << endl;
    cout << "Number of processors = " << NumProc << endl;
    cout << "Print Process = " << MyPID << endl;
  }
  Comm.Barrier();
  if (printing.isPrintType(NOX::Utils::TestDetails))
    cout << "Process " << MyPID << " is alive!" << endl;
  Comm.Barrier();
#else
  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails))
    cout << "Serial Run" << endl;
#endif
  // Sublist for line search 
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  // This test is designed to exercise the following linesearch
  searchParams.setParameter("Method", "Polynomial");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "Newton");
  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 20);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Preconditioner", "None");   
  if( verbose )
    lsParams.setParameter("Output Frequency", 1);    

  // Create the interface between the test problem and the nonlinear solver
  Problem_Interface interface(Problem);
  
  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)
  Epetra_RowMatrix& A = Problem.getJacobian();

  // Create the callback interfaces for filling the residual and Jacbian
  NOX::EpetraNew::Interface::Required& iReq = interface;
  NOX::EpetraNew::Interface::Jacobian& iJac = interface;

  // Create the Linear System
  NOX::EpetraNew::LinearSystemAztecOO linSys(printParams, lsParams,
                                             iReq, iJac, A, soln);

  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::DeepCopy, true);
  NOX::EpetraNew::Group grp(printParams, iReq, initialGuess, linSys); 

  // Establish initial convergence status
  grp.computeF();

  // Create the convergence tests
  NOX::StatusTest::NormF testNormF(1.0e-6);
  NOX::StatusTest::MaxIters testMaxIters(25);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR, testNormF, testMaxIters);

  // Create the method
  NOX::Solver::Manager solver(grp, combo, nlParams);
  NOX::StatusTest::StatusType status = solver.solve();

  int testStatus = 0;

  if (status != NOX::StatusTest::Converged) {
    if (MyPID==0) cout << "Nonlinear solver failed to converge!" << endl;
    testStatus = -1;
  }

  // Test for correct line-search and linear solver behavior
  NOX::Parameter::List & outParams = nlParams.sublist("Line Search").sublist("Output");
  if( (outParams.getParameter("Total Number of Failed Line Searches", -1) != 0) ||
    (outParams.getParameter("Total Number of Line Search Calls", 0) != 7) ||
    (outParams.getParameter("Total Number of Line Search Inner Iterations", 0) != 4) ||
    (outParams.getParameter("Total Number of Non-trivial Line Searches", 0) != 2) ||
    (lsParams.sublist("Output").getParameter("Total Number of Linear Iterations", 0) != 14) ) {
    testStatus = -1;
  }

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::EpetraNew::Group& finalGroup = 
      dynamic_cast<const NOX::EpetraNew::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  // End Nonlinear Solver **************************************

  // Output the parameter list
  if (printing.isPrintProcessAndType(NOX::Utils::Parameters)) {
    cout << endl << "Final Parameters" << endl
	 << "****************" << endl;
    solver.getParameterList().print(cout);
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

  if (printing.isPrintProcessAndType(NOX::Utils::TestDetails)) {
    if (testStatus == 0)
      cout << "Test Successfull!" << endl;
    else
      cout << "Test Failed!" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
 
  return testStatus;

} // end main
