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

  // The number of unknowns must be at least equal to the 
  // number of processors.
  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << endl;
    cout << "Test failed!" << endl;
    throw "NOX Error";
  }

  // Create the interface between NOX and the application
  // This object is derived from NOX::Epetra::Interface
  Interface interface(NumGlobalElements, Comm);

  // Get the vector from the Problem
  Epetra_Vector& soln = interface.getSolution();

  // Set the PDE factor (for nonlinear forcing term).  This could be specified
  // via user input.
  interface.setPDEfactor(1000.0);

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  NOX::Parameter::List nlParams;

  // Set the nonlinear solver method
  nlParams.setParameter("Nonlinear Solver", "Line Search Based");

  // Set the printing parameters in the "Printing" sublist
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", MyPID); 
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  if (verbose)
    printParams.setParameter("Output Information", 
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
    printParams.setParameter("Output Information", NOX::Utils::Error +
			     NOX::Utils::TestDetails);

  // Sublist for line search 
  NOX::Parameter::List& searchParams = nlParams.sublist("Line Search");
  searchParams.setParameter("Method", "Full Step");

  // Sublist for direction
  NOX::Parameter::List& dirParams = nlParams.sublist("Direction");
  dirParams.setParameter("Method", "Newton");
  NOX::Parameter::List& newtonParams = dirParams.sublist("Newton");
    newtonParams.setParameter("Forcing Term Method", "Constant");

  // Sublist for linear solver for the Newton method
  NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);  
  //lsParams.setParameter("Max Age Of Prec", 3);

  // Various Preconditioner options
  //lsParams.setParameter("Preconditioner", "None");
  //lsParams.setParameter("Preconditioner", "AztecOO");
  lsParams.setParameter("Preconditioner", "Ifpack");
  lsParams.setParameter("Max Age Of Prec", 5);

  // New Ifpack interface
  //lsParams.setParameter("Preconditioner", "New Ifpack");
  Teuchos::ParameterList IfpackList;


#ifdef HAVE_NOX_ML_EPETRA
  // Comment out the previous Preconditioner spec and uncomment this line
  // to turn on ML 
  //lsParams.setParameter("Preconditioner", "ML");
  // create a parameter list for ML options
  Teuchos::ParameterList MLList;
#endif

  //lsParams.setParameter("Jacobian Operator", "Finite Difference");
  //lsParams.setParameter("Jacobian Operator", "Matrix-Free");
  //lsParams.setParameter("Preconditioner Operator", "Finite Difference");

  // Set specific Preconditioning Package options if valid
  if( lsParams.getParameter("Preconditioner", "None") == "New Ifpack" ) {
    lsParams.setParameter("Ifpack Preconditioner", "ILU");
    // Can fill Ifpack options into IfpackList here

    lsParams.setParameter("Ifpack Teuchos Parameter List", &IfpackList);
  }

#ifdef HAVE_NOX_ML_EPETRA
  if( lsParams.getParameter("Preconditioner", "None") == "ML" ) {
    // This Teuchos parameter list is needed for ML
  
    // These specifications come straight from the example in 
    // Trilinos/packages/ml/example/ml_example_epetra_preconditioner.cpp
  
    // set defaults for classic smoothed aggregation
    ML_Epetra::SetDefaults("SA",MLList);
    // maximum number of levels
    MLList.set("max levels",5);
    MLList.set("increasing or decreasing","decreasing");
    // use Uncoupled scheme to create the aggregate,
    // from level 3 use the better but more expensive MIS
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("aggregation: type (level 3)", "MIS");
    // smoother is Gauss-Seidel. Example file
    // ml_example_epetra_preconditioner_2level.cpp shows how to use
    // AZTEC's preconditioners as smoothers
    MLList.set("smoother: type","Gauss-Seidel");
    // use both pre and post smoothing
    MLList.set("smoother: pre or post", "both");
    // solve with serial direct solver KLU
    MLList.set("coarse: type","Jacobi");
  
    lsParams.setParameter("ML Teuchos Parameter List", &MLList);
  }
#endif

  // Use an Epetra Scaling object if desired
  Epetra_Vector scaleVec(soln);
  NOX::Epetra::Scaling scaling;
  scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);
  //grp.setLinearSolveScaling(scaling);

  // Create all possible Epetra_Operators.
  // 1. User supplied (Epetra_RowMatrix)
  //Epetra_RowMatrix& Analytic = interface.getJacobian();
  // 2. Matrix-Free (Epetra_Operator)
  NOX::Epetra::MatrixFree MF(interface, soln);
  // 3. Finite Difference (Epetra_RowMatrix)
  NOX::Epetra::FiniteDifference FD(interface, soln);

  // Four constructors to create the Linear System
  NOX::Epetra::Interface::Required& iReq = interface;

  // **** Ctor #1 - No Jac and No Prec
  //NOX::Epetra::LinearSystemAztecOO linSys(printParams, lsParams, 
  //				      iReq, soln);

  // **** Ctor #2 - Jac but no Prec
  //NOX::Epetra::Interface::Jacobian& iJac = FD;
  //NOX::Epetra::Interface::Jacobian& iJac = interface;
  //NOX::Epetra::LinearSystemAztecOO linSys(printParams, lsParams,
  //				     iReq, iJac, FD, soln);

  // **** Ctor #3 - Prec but no Jac
  //NOX::Epetra::Interface::Preconditioner& iPrec = FD;
  //NOX::Epetra::LinearSystemAztecOO linSys(printParams, lsParams,
  //				      iReq, iPrec, FD, soln);

  // **** Ctor #4 - Prec and Jac
  NOX::Epetra::Interface::Jacobian& iJac = MF;
  NOX::Epetra::Interface::Preconditioner& iPrec = FD;
  NOX::Epetra::LinearSystemAztecOO linSys(printParams, lsParams,
					     iJac, MF, iPrec, FD, soln);
					     //&scaling);

  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::DeepCopy, true);
  NOX::Epetra::Group grp(printParams, iReq, initialGuess, linSys);  

  // uncomment the following for loca supergroups
  MF.setGroupForComputeF(grp);
  FD.setGroupForComputeF(grp);

  // Test group accessor
  //const NOX::Epetra::LinearSystem& lins = grp.getLinearSystem();
  //lins.createPreconditioner(soln, lsParams, false);

  // Create the convergence tests
  NOX::StatusTest::NormF absresid(1.0e-8);
  NOX::StatusTest::NormF relresid(grp, 1.0e-2);
  NOX::StatusTest::NormUpdate update(1.0e-5);
  NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  NOX::StatusTest::Combo converged(NOX::StatusTest::Combo::AND);
  converged.addStatusTest(absresid);
  converged.addStatusTest(relresid);
  converged.addStatusTest(wrms);
  converged.addStatusTest(update);
  NOX::StatusTest::MaxIters maxiters(20);
  NOX::StatusTest::Combo combo(NOX::StatusTest::Combo::OR);
  NOX::StatusTest::FiniteValue fv;
  combo.addStatusTest(fv);
  combo.addStatusTest(converged);
  combo.addStatusTest(maxiters);

  // Create the solver
  NOX::Solver::Manager solver(grp, combo, nlParams);
  NOX::StatusTest::StatusType solvStatus = solver.solve();

  // Create a print class for controlling output below
  NOX::Utils printing(printParams);

  //
  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
  //

  // End Nonlinear Solver **************************************

  // Output the parameter list
  if (verbose) {
    if (printing.isPrintProcessAndType(NOX::Utils::Parameters)) {
      cout << endl << "Final Parameters" << endl
	   << "****************" << endl;
      solver.getParameterList().print(cout);
      cout << endl;
    }
  }

  //
  // Print solution
  char file_name[25];
  FILE *ifp;
  int NumMyElements = soln.Map().NumMyElements();
  (void) sprintf(file_name, "output.%d",MyPID);
  ifp = fopen(file_name, "w");
  for (int i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E\n", soln.Map().MinMyGID()+i, finalSolution[i]);
  fclose(ifp);
  //


  // Tests
  int status = 0; // Converged
  
  // 1. Convergence
  if (solvStatus != NOX::StatusTest::Converged) {
      status = 1;
      if (printing.isPrintProcessAndType(NOX::Utils::Error))
	cout << "Nonlinear solver failed to converge!" << endl;
  }
  // 2. Linear solve iterations (53)
  if (solver.getParameterList().sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").getParameter("Total Number of Linear Iterations",0) != 53) {
    status = 1;
  }
  // 3. Nonlinear solve iterations (10)
  if (solver.getParameterList().sublist("Output").getParameter("Nonlinear Iterations", 0) != 10)
    status = 1;


  // Summarize test results 
  if (printing.isPrintProcess()) { 
    if (status == 0)
      cout << "Test passed!" << endl;
    else 
      cout << "Test failed!" << endl;
  }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
