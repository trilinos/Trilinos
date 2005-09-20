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
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

#include <vector>
#ifdef HAVE_NOX_EPETRAEXT 	// Use epetraext package in Trilinos
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#else  // Otherwise use the local version of the needed coloring files
#include "nox_EDT_CrsGraph_MapColoring.H"
#include "nox_EDT_CrsGraph_MapColoringIndex.H"
#endif

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
    cout << "Error: numGlobalBlocks = " << NumGlobalElements 
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
  if (verbose)
    printParams.setParameter("Output Information", 
			     NOX::Utils::OuterIteration + 
			     NOX::Utils::OuterIterationStatusTest + 
			     NOX::Utils::InnerIteration +
			     NOX::Utils::Parameters + 
			     NOX::Utils::Details + 
			     NOX::Utils::Warning);
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
  lsParams.setParameter("Output Frequency", 0); 
  lsParams.setParameter("Preconditioner", "Ifpack");
  lsParams.setParameter("Max Age Of Prec", 5); 

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  // NLS_PetraGroupInterface
  Problem_Interface interface(Problem);

  // Create the Epetra_RowMatrix using Finite Difference with Coloring
#ifndef HAVE_NOX_EPETRAEXT
  bool verbose_ = false;
  Epetra_MapColoring* colorMap = Epetra_Transform::
        CrsGraph_MapColoring( verbose_ )( Problem.getGraph());
  Epetra_Transform::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
  vector<Epetra_IntVector>* columns = 
                                     colorMapIndex(Problem.getGraph());
#else
  bool verbose_ = false;
  EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType = 
    EpetraExt::CrsGraph_MapColoring::GREEDY;
  EpetraExt::CrsGraph_MapColoring tmpMapColoring( algType, verbose_ );
  Epetra_MapColoring* colorMap = &tmpMapColoring(Problem.getGraph());
  EpetraExt::CrsGraph_MapColoringIndex colorMapIndex(*colorMap);
  vector<Epetra_IntVector>* columns = &colorMapIndex(Problem.getGraph());
#endif

  // Use this constructor to create the graph numerically as a means of timing
  // the old way of looping without colors :
  //  NOX::Epetra::FiniteDifferenceColoring A(interface, soln, 
  //                                          *colorMap, *columns);
  // Or use this as the standard way of using finite differencing with coloring
  // where the application is responsible for creating the matrix graph 
  // beforehand, ie as is done in Problem.
  NOX::Epetra::Interface::Required& iReq = interface;
  NOX::Epetra::FiniteDifferenceColoring A(iReq, soln, Problem.getGraph(),
					  *colorMap, *columns);

  // Create the linear system
  NOX::Epetra::Interface::Jacobian& iJac = interface;
  NOX::Epetra::LinearSystemAztecOO linSys(printParams, lsParams,
					  iReq, iJac, A, soln);

  // Create the Group
  NOX::Epetra::Vector initialGuess(soln, NOX::DeepCopy, true);
  Teuchos::RefCountPtr<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
					iReq, 
					initialGuess, 
					linSys));   
  NOX::Epetra::Group& grp = *(grpPtr.get()); 

  // Use an Epetra Scaling object if desired
  Epetra_Vector scaleVec(soln);
  NOX::Epetra::Scaling scaling;
  scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);
  //grp.setLinearSolveScaling(scaling);

  // ATOL vector if using NOX::StatusTest::WRMS
  NOX::Epetra::Vector weights(soln);
  weights.scale(1.0e-8);

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

  // Create the method
  NOX::Solver::Manager solver(grpPtr, combo, nlParamsPtr);
  NOX::StatusTest::StatusType status = solver.solve();

  if (verbose) {
    if (status != NOX::StatusTest::Converged)
      if (MyPID==0) 
	cout << "Nonlinear solver failed to converge!" << endl;
  }

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

  // End Nonlinear Solver **************************************

  // Output the parameter list
  if (verbose) {
    NOX::Utils utils(printParams);
    if (utils.isPrintType(NOX::Utils::Parameters)) {
      cout << endl << "Final Parameters" << endl
	   << "****************" << endl;
      solver.getParameterList().print(cout);
      cout << endl;
    }
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

  // Summarize test results
  if (Comm.MyPID() == 0) {  
    if (status == NOX::StatusTest::Converged)
      cout << "Test passed!" << endl;
    else 
      cout << "Test failed!" << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  // Final return value (0 = successfull, non-zero = failure)
  return ierr ;
}
