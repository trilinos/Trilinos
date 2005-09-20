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
/* Solves the nonlinear equation:
 *
 * d2u 
 * --- - k * u**2 = 0
 * dx2
 *
 * subject to @ x=0, u=1
 */

// NOX Objects
#include "NOX.H"
#include "NOX_Epetra.H"

#include "NOX_Belos_Group.H"

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
#include "1DfemInterface.H" // Interface file to NOX

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
  Interface interface(NumGlobalElements, Comm);
  interface.setPDEfactor(1000.0);

  // Get the vector from the Problem
  Epetra_Vector& soln = interface.getSolution();

  // Initialize Solution
  soln.PutScalar(1.0);
  
  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RefCountPtr<NOX::Parameter::List> nlParamsPtr =
    Teuchos::rcp(new NOX::Parameter::List);
  NOX::Parameter::List& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");
  //nlParams.setParameter("Nonlinear Solver", "Trust Region Based");

  // Set the printing parameters in the "Printing" sublist
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", MyPID); 
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  if (verbose) {
    printParams.setParameter("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning);
  }
  else
    printParams.setParameter("Output Information",NOX::Utils::Error);

  // Create a print handler for output control
  NOX::Utils utils(printParams);

  // Sublist for line search 
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  searchParams.setParameter("Method", "Full Step");
  //searchParams.setParameter("Method", "Interval Halving");
  //searchParams.setParameter("Method", "Polynomial");
  //searchParams.setParameter("Method", "NonlinearCG");
  //searchParams.setParameter("Method", "Quadratic");
  //searchParams.setParameter("Method", "More'-Thuente");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
//  dirParams.setParameter("Method", "Modified-Newton");
//  NOX::Parameter::List& newtonParams = dirParams.sublist("Modified-Newton");
//    newtonParams.setParameter("Max Age of Jacobian", 2);
  dirParams.setParameter("Method", "Newton");
  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");
    //newtonParams.setParameter("Forcing Term Method", "Type 1");
    //newtonParams.setParameter("Forcing Term Method", "Type 2");
    //newtonParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-4);
    //newtonParams.setParameter("Forcing Term Maximum Tolerance", 0.1);
  //dirParams.setParameter("Method", "Steepest Descent");
  //NOX::Parameter::List& sdParams = dirParams.sublist("Steepest Descent");
    //sdParams.setParameter("Scaling Type", "None");
    //sdParams.setParameter("Scaling Type", "2-Norm");
    //sdParams.setParameter("Scaling Type", "Quadratic Model Min");
  //dirParams.setParameter("Method", "NonlinearCG");
  //NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
    //nlcgParams.setParameter("Restart Frequency", 2000);
    //nlcgParams.setParameter("Precondition", "On");
    //nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

  // Sublist for linear solver for the Newton method
  NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  //lsParams.setParameter("Aztec Solver", "GMRES"); 
  lsParams.setParameter("Belos Solver", "GMRES"); 
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Preconditioner", "Ifpack");
  lsParams.setParameter("Preconditioner Operator", "Use Jacobian");
  lsParams.setParameter("Output Frequency", 0);
  lsParams.setParameter("Verbosity Level", 1);

  // Create the Epetra_RowMatrix.  Uncomment one or more of the following:
  // 1. User supplied (Epetra_RowMatrix)
  Epetra_RowMatrix& A = interface.getJacobian();

  // Create the linear system
  NOX::Epetra::Interface::Required& iReq = interface;
  NOX::Epetra::Interface::Jacobian& iJac = interface;
  NOX::Epetra::LinearSystemAztecOO linSys(printParams, lsParams,
					  iReq, iJac, A, soln);

  // Create the Group
  NOX::Epetra::Vector nv(soln, NOX::DeepCopy, true);
  Teuchos::RefCountPtr<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
					iReq, 
					nv, 
					linSys));  
  NOX::Epetra::Group& grp = *(grpPtr.get());

  // Use an Epetra Scaling object if desired
  Epetra_Vector scaleVec(soln);
  NOX::Epetra::Scaling scaling;
  scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);
  //grp.setLinearSolveScaling(scaling);

  // Create the convergence tests
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::NormF> relresid = 
    Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-2));
  Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate> update =
    Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
  Teuchos::RefCountPtr<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(relresid);
  converged->addStatusTest(wrms);
  converged->addStatusTest(update);
  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create the belos group
  Teuchos::RefCountPtr<NOX::Belos::Group> belos_grp = 
    Teuchos::rcp(new NOX::Belos::Group(grp, printParams));

  // Create the method
  NOX::Solver::Manager solver(belos_grp, combo, nlParamsPtr);
  NOX::StatusTest::StatusType solverStatus = solver.solve();

  if (verbose) {
    if (solverStatus != NOX::StatusTest::Converged)
      if (MyPID==0) 
	utils.out() << "Nonlinear solver failed to converge!" << endl;
  }

  // Get the Epetra_Vector with the final solution from the solver
  /*
  const NOX::Belos::Group& finalGroup = dynamic_cast<const NOX::Belos::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
  */

  // End Nonlinear Solver **************************************

  // Output the parameter list
  if (verbose) {
    if (utils.isPrintType(NOX::Utils::Parameters)) {
      utils.out() << endl << "Final Parameters" << endl
	   << "****************" << endl;
      solver.getParameterList().print(utils.out());
      utils.out() << endl;
    }
  }

  // Tests
  int status = 0; // Converged
  
  // 1. Convergence
  if (solverStatus != NOX::StatusTest::Converged) {
    status = 1;
    if (utils.isPrintType(NOX::Utils::Error))
      utils.out() << "Nonlinear solver failed to converge!" << endl;
  }
  // 2. Nonlinear solve iterations (10)
  if (solver.getParameterList().sublist("Output").getParameter("Nonlinear Iterations", 0) != 10)
    status = 1;
  
  
  // Summarize test results 
  if (status == 0)
    utils.out() << "Test passed!" << endl;
  else 
    utils.out() << "Test failed!" << endl;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

/* end main
*/
return status;
}
