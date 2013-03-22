//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
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
    std::cout << "ERROR: Maximum number of processors is 2!" << std::endl;
    exit(1);
  }

  // Create the Problem class.  This creates all required
  // Epetra objects for the problem and allows calls to the 
  // function (RHS) and Jacobian evaluation routines.
  DennisSchnabel Problem(NumGlobalElements, Comm);

  // Get the vector from the Problem
  Teuchos::RCP<Epetra_Vector> soln = Problem.getSolution();
  NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

  // Initialize Solution
  if (MyPID==0) {
    (*soln)[0]=2.0;
    if (NumProc==1) 
      (*soln)[1]=0.5;
  } 
  else 
    (*soln)[0]=0.5;

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Inexact Trust Region Based");
  nlParams.sublist("Trust Region").
    set("Inner Iteration Method", "Standard Trust Region");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID); 
  printParams.set("Output Precision", 5);
  printParams.set("Output Processor", 0);
  if ( verbose )
    printParams.set("Output Information",
                     NOX::Utils::OuterIteration +
                     NOX::Utils::OuterIterationStatusTest +
                     NOX::Utils::InnerIteration +
                     NOX::Utils::Parameters +
                     NOX::Utils::Details +
                     NOX::Utils::Warning +
                     NOX::Utils::TestDetails);
  else
    printParams.set("Output Information", NOX::Utils::Error +
                     NOX::Utils::TestDetails);

  NOX::Utils printing(printParams);

  // Identify the test problem
  if (printing.isPrintType(NOX::Utils::TestDetails))
    std::cout << "Starting epetra/DS6.5.1/DS_6_5_1.exe" << std::endl;

  // Identify processor information
#ifdef HAVE_MPI
  if (printing.isPrintType(NOX::Utils::TestDetails)) {
    std::cout << "Parallel Run" << std::endl;
    std::cout << "Number of processors = " << NumProc << std::endl;
    std::cout << "Print Process = " << MyPID << std::endl;
  }
  Comm.Barrier();
  if (printing.isPrintType(NOX::Utils::TestDetails))
    std::cout << "Process " << MyPID << " is alive!" << std::endl;
  Comm.Barrier();
#else
  if (printing.isPrintType(NOX::Utils::TestDetails))
    std::cout << "Serial Run" << std::endl;
#endif

  // Sublist for Newton direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "Newton");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  lsParams.set("Max Iterations", 20);  
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Preconditioner", "None");   
  if( verbose )
    lsParams.set("Output Frequency", 1);    

  // Sublist for Cauchy direction
  Teuchos::ParameterList& cauchyDirParams = nlParams.sublist("Cauchy Direction");
  cauchyDirParams.set("Method", "Steepest Descent");
  Teuchos::ParameterList& sdParams = cauchyDirParams.sublist("Steepest Descent");
    sdParams.set("Scaling Type", "Quadratic Model Min");

  // Create the interface between the test problem and the nonlinear solver
  Teuchos::RCP<Problem_Interface> interface = 
    Teuchos::rcp(new Problem_Interface(Problem));
  
  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)
  Teuchos::RCP<Epetra_RowMatrix> A = Problem.getJacobian();

  // Create the callback interfaces for filling the residual and Jacbian
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;

  // Create the Linear System
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      iReq, iJac, A, noxSoln));

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
					iReq, 
					noxSoln, 
					linSys)); 

  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> testNormF = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-6));
  Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(25));
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					    testNormF, testMaxIters));

  // Create the method
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType status = solver->solve();

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = 
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  // End Nonlinear Solver **************************************

  // Output the parameter list
  if (printing.isPrintType(NOX::Utils::Parameters)) {
    std::cout << std::endl << "Final Parameters" << std::endl
	 << "****************" << std::endl;
    solver->getList().print(cout);
    std::cout << std::endl;
  }

  // Print solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = soln->Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",MyPID);
  ifp = fopen(file_name, "w");
  for (i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E\n", soln->Map().MinMyGID()+i, finalSolution[i]);
  fclose(ifp);

  // Report results

  int testStatus = 0; // Converged

  // 1. Convergence
  if (status != NOX::StatusTest::Converged) {
    if (MyPID==0) std::cout << "Nonlinear solver failed to converge!" << std::endl;
    testStatus = 1;
  }
  // 2. Nonlinear Iterations (5)
  if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Output").get("Nonlinear Iterations",0) != 5) {
    testStatus = 2;
  }

  if (testStatus == 0)
    std::cout << "Test passed!" << std::endl;
  else
    std::cout << "Test failed!" << std::endl;

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
 
  return testStatus;

} // end main
