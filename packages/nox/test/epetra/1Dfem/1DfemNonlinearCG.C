//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
                                                                    
// 1D Finite Element NonlinearCG Test Problem
/* Solves the nonlinear equation:
 *
 * d2u 
 * --- - k * u**2 = 0
 * dx2
 *
 * subject to @ x=0, u=1
*
* using nonlinear CG direction and linesearch methods.

  This problem is not challenging but serves primarily to demonstrate
  use of NonlinearCG as well as ensure it is not broken.

  Preconditioning is performed using 5 gmres iterations with the 
  analytic jacobian matrix for the problem.  Note that solving these
  linear systems completely, omitting orthogonalization
  and using a full step (1.0) is equivalent to using Newton's method.

 */

// NOX Objects
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
#include "1DfemInterface.H" 
#include "1DfemPrePostOperator.H"

#include "Teuchos_ParameterList.hpp"

using namespace std;

int main(int argc, char *argv[])
{
 
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

  // Check verbosity level
  bool verbose = false;
  if (argc > 1)
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

  // The number of unknowns must be at least equal to the number of processors.
  if (NumGlobalElements < NumProc) 
  {
    cout << "numGlobalBlocks = " << NumGlobalElements << " cannot be < number of processors = " 
         << NumProc << endl;
    cout << "Test failed!" << endl;
    throw "NOX Error";
  }

  // Create the interface between NOX and the application
  // This object is derived from NOX::Epetra::Interface
  Teuchos::RefCountPtr<Interface> interface = Teuchos::rcp(new Interface(NumGlobalElements, Comm));

  // Get the vector from the Problem
  Teuchos::RefCountPtr<Epetra_Vector> soln = interface->getSolution();
  Teuchos::RefCountPtr<NOX::Epetra::Vector> noxSoln = 
    Teuchos::rcp(new NOX::Epetra::Vector(soln, NOX::Epetra::Vector::CreateView));

  // Set the PDE factor (for nonlinear forcing term).  This could be specified
  // via user input.
  interface->setPDEfactor(1000.0);

  // Set the initial guess 
  soln->PutScalar(1.0);

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RefCountPtr<Teuchos::ParameterList> nlParamsPtr = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

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
			     NOX::Utils::LinearSolverDetails +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning +
                             NOX::Utils::Debug +
			     NOX::Utils::TestDetails +
			     NOX::Utils::Error);
  else
    printParams.set("Output Information", NOX::Utils::Error +
			     NOX::Utils::TestDetails);

  // Create a print class for controlling output below
  NOX::Utils printing(printParams);

  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method", "NonlinearCG"); // "Full Step" can also work well sometimes

  // Sublist for direction
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "NonlinearCG");
  Teuchos::ParameterList& nonlinearcg = dirParams.sublist("Nonlinear CG");
    nonlinearcg.set("Restart Frequency", 100);
    nonlinearcg.set("Precondition", "On");
    nonlinearcg.set("Orthogonalize", "Fletcher-Reeves"); // or "Polak-Ribiere"

  // Sublist for linear solver for the Newton method
  Teuchos::ParameterList& lsParams = nonlinearcg.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");  
  //lsParams.set("Preconditioner Operator", "Use Jacobian");
  lsParams.set("Preconditioner", "AztecOO");
  lsParams.set("AztecOO Preconditioner Iterations", 5);
  lsParams.set("Preconditioner Reuse Policy", "Recompute");

  // Let's force all status tests to do a full check
  nlParams.sublist("Solver Options").set("Status Test Check Type", NOX::StatusTest::Complete);

  // 1. User supplied (Epetra_RowMatrix)
  Teuchos::RefCountPtr<Epetra_RowMatrix> Analytic = interface->getJacobian();

  // Create the linear system
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
						      interface, 
						      iJac, Analytic, 
						      *soln));
  
  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::Epetra::Vector::CreateView);
  Teuchos::RefCountPtr<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
					iReq, 
					initialGuess, 
					linSys));  

  // Create the convergence tests
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(absresid);
  combo->addStatusTest(maxiters);

  // Create the solver
  NOX::Solver::Manager solver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType solvStatus = solver.solve();

  // End Nonlinear Solver **************************************

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = 
    dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = 
    (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).
    getEpetraVector();

  // Output the parameter list
  if (verbose) 
  {
    if (printing.isPrintType(NOX::Utils::Parameters)) {
      printing.out() << endl << "Final Parameters" << endl
	   << "****************" << endl;
      solver.getList().print(printing.out());
      printing.out() << endl;
    }
  }

  // Print solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = soln->Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",MyPID);
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E\n", soln->Map().MinMyGID()+i, finalSolution[i]);
  fclose(ifp);


  // Tests
  int status = 0; // Converged
  
  // 1. Convergence
  if (solvStatus != NOX::StatusTest::Converged) 
  {
      status = 1;
      if (printing.isPrintType(NOX::Utils::Error))
	printing.out() << "Nonlinear solver failed to converge!" << endl;
  }
  // 2. Nonlinear solve iterations (10)
  if (const_cast<Teuchos::ParameterList&>(solver.getList()).sublist("Output").get("Nonlinear Iterations", 0) != 13)
    status = 2;

  // Summarize test results 
  if (status == 0)
    printing.out() << "Test passed!" << endl;
  else 
    printing.out() << "Test failed!" << endl;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  printing.out() << "Status = " << status << endl;

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
