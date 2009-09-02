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
                                                                    
// 1D Finite Element Test Problem
/* Solves the nonlinear equation:
 *
 * d2u 
 * --- - k * u**2 = 0
 * dx2
 *
 * subject to @ x=0, u=1
 */

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

// Azteco Objects
#include "AztecOO.h"
#include "Ifpack.h"

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
  Teuchos::RCP<Interface> interface = 
    Teuchos::rcp(new Interface(NumGlobalElements, Comm));

  // Set the PDE factor (for nonlinear forcing term).  This could be specified
  // via user input.
  interface->setPDEfactor(1000.0);

  // Begin Nonlinear Solver ************************************

  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> IfpackParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList printParams;
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
  NOX::Utils p(printParams);


  // *******************************
  // Setup Test Objects
  // *******************************

  // Create Linear Objects
  // Get the vector from the Problem
  if (verbose)
    p.out() << "Creating Vectors and Matrices" << endl;
  Teuchos::RCP<Epetra_Vector> solution_vec = 
    interface->getSolution();
  Teuchos::RCP<Epetra_Vector> rhs_vec = 
    Teuchos::rcp(new Epetra_Vector(*solution_vec));
  Teuchos::RCP<Epetra_Vector> lhs_vec = 
    Teuchos::rcp(new Epetra_Vector(*solution_vec));
  Teuchos::RCP<Epetra_CrsMatrix> jacobian_matrix = 
    interface->getJacobian();


  if (verbose)
    p.out() << "Evaluating F and J" << endl;
  solution_vec->PutScalar(1.0);
  interface->computeF(*solution_vec, *rhs_vec);
  rhs_vec->Scale(-1.0);
  interface->computeJacobian(*solution_vec, *jacobian_matrix);

  double norm =0.0;
  rhs_vec->Norm2(&norm);
  if (verbose)
    p.out() << "Step 0, ||F|| = " << norm << endl;

  if (verbose)
    p.out() << "Creating Aztec Solver" << endl;

  Teuchos::RCP<AztecOO> aztecSolverPtr = Teuchos::rcp(new AztecOO());
  if (verbose)
    aztecSolverPtr->SetAztecOption(AZ_output, AZ_last);
  else
    aztecSolverPtr->SetAztecOption(AZ_output, AZ_none);

  // *******************************
  // Rebuild Test
  // *******************************

  if (verbose) {
    p.out() << "**********************************************" << endl;
    p.out() << "Testing Newton solve rebuilding prec each step" << endl;
    p.out() << "**********************************************" << endl;
  }

  int step_number = 0;
  int max_steps = 20;
  bool converged = false;
  int total_linear_iterations = 0;
  while (norm > 1.0e-8 && step_number < max_steps) {

    step_number++;

    if (verbose)
      p.out() << "Step " << step_number << ", ||F|| = " << norm << endl;

    if (verbose)
      p.out() << "Creating Ifpack preconditioner" << endl;
    
    Ifpack Factory;
    Teuchos::RCP<Ifpack_Preconditioner> PreconditionerPtr;
    PreconditionerPtr = Teuchos::rcp(Factory.Create("ILU",
						    jacobian_matrix.get(),0));
    Teuchos::ParameterList teuchosParams;
    PreconditionerPtr->SetParameters(teuchosParams);
    PreconditionerPtr->Initialize();
    PreconditionerPtr->Compute();

    aztecSolverPtr->SetUserMatrix(jacobian_matrix.get(), false);
    aztecSolverPtr->SetPrecOperator(PreconditionerPtr.get());
    aztecSolverPtr->SetRHS(rhs_vec.get());
    aztecSolverPtr->SetLHS(lhs_vec.get());

    aztecSolverPtr->Iterate(400, 1.0e-4);

    solution_vec->Update(1.0, *lhs_vec, 1.0);
    
    interface->computeF(*solution_vec, *rhs_vec);
    rhs_vec->Scale(-1.0);
    interface->computeJacobian(*solution_vec, *jacobian_matrix);

    rhs_vec->Norm2(&norm);

    total_linear_iterations += aztecSolverPtr->NumIters();

    if (norm < 1.0e-8) 
      converged = true;
  }

  if (verbose) {
    p.out() << "Final Step " << step_number << ", ||F|| = " << norm << endl;
    if (converged)
      p.out() << "Converged!!" << endl;
    else 
      p.out() << "Failed!!" << endl;
  }

  // Tests
  int status = 0; // Converged

  if (verbose)
    p.out() << "Total Number of Linear Iterations = "
	    << total_linear_iterations << endl;

  if (Comm.NumProc() == 1 && total_linear_iterations != 10)
    status = 1;

  if (!converged)
    status = 2;


  // Summarize test results 
  if (converged && status == 0)
    p.out() << "Test passed!" << endl;
  else 
    p.out() << "Test failed!" << endl;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (verbose)
    p.out() << "Status = " << status << endl;

  // Final return value (0 = successfull, non-zero = failure)
  return status;
}
