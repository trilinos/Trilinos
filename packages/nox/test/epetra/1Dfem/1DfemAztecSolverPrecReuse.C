// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "Teuchos_StandardCatchMacros.hpp"

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

  // Check verbosity level
  bool verbose = false;
  bool success = false;
  try {
    // Get the process ID and the total number of processors
    int MyPID = Comm.MyPID();
    int NumProc = Comm.NumProc();

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
      std::cout << "numGlobalBlocks = " << NumGlobalElements
        << " cannot be < number of processors = " << NumProc << std::endl;
      std::cout << "Test failed!" << std::endl;
      throw std::runtime_error("NOX Error");
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
      p.out() << "Creating Vectors and Matrices" << std::endl;
    Teuchos::RCP<Epetra_Vector> solution_vec =
      interface->getSolution();
    Teuchos::RCP<Epetra_Vector> rhs_vec =
      Teuchos::rcp(new Epetra_Vector(*solution_vec));
    Teuchos::RCP<Epetra_Vector> lhs_vec =
      Teuchos::rcp(new Epetra_Vector(*solution_vec));
    Teuchos::RCP<Epetra_CrsMatrix> jacobian_matrix =
      interface->getJacobian();


    if (verbose)
      p.out() << "Evaluating F and J" << std::endl;
    solution_vec->PutScalar(1.0);
    interface->computeF(*solution_vec, *rhs_vec);
    rhs_vec->Scale(-1.0);
    interface->computeJacobian(*solution_vec, *jacobian_matrix);

    double norm =0.0;
    rhs_vec->Norm2(&norm);
    if (verbose)
      p.out() << "Step 0, ||F|| = " << norm << std::endl;


    if (verbose)
      p.out() << "Creating Ifpack preconditioner" << std::endl;

    Ifpack Factory;
    Teuchos::RCP<Ifpack_Preconditioner> PreconditionerPtr;
    PreconditionerPtr = Teuchos::rcp(Factory.Create("ILU",
          jacobian_matrix.get(),0));
    Teuchos::ParameterList teuchosParams;
    PreconditionerPtr->SetParameters(teuchosParams);
    PreconditionerPtr->Initialize();
    PreconditionerPtr->Compute();


    if (verbose)
      p.out() << "Creating Aztec Solver" << std::endl;

    Teuchos::RCP<AztecOO> aztecSolverPtr = Teuchos::rcp(new AztecOO());
    if (verbose)
      aztecSolverPtr->SetAztecOption(AZ_output, AZ_last);
    else
      aztecSolverPtr->SetAztecOption(AZ_output, AZ_none);

    // *******************************
    // Reuse Test
    // *******************************

    if (verbose) {
      p.out() << "**********************************************" << std::endl;
      p.out() << "Testing Newton solve with prec reuse" << std::endl;
      p.out() << "**********************************************" << std::endl;
    }

    int step_number = 0;
    int max_steps = 20;
    bool converged = false;
    int total_linear_iterations = 0;
    while (norm > 1.0e-8 && step_number < max_steps) {

      step_number++;

      if (verbose)
        p.out() << "Step " << step_number << ", ||F|| = " << norm << std::endl;

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
      p.out() << "Final Step " << step_number << ", ||F|| = " << norm << std::endl;
      if (converged)
        p.out() << "Converged!!" << std::endl;
      else
        p.out() << "Failed!!" << std::endl;
    }

    // Tests
    int status = 0; // Converged

    if (verbose)
      p.out() << "Total Number of Linear Iterations = "
        << total_linear_iterations << std::endl;

    if (Comm.NumProc() == 1 && total_linear_iterations != 157)
      status = 1;

    if (!converged)
      status = 2;

    success = converged && status == 0;

    // Summarize test results
    if (success)
      p.out() << "Test passed!" << std::endl;
    else
      p.out() << "Test failed!" << std::endl;

    if (verbose)
      p.out() << "Status = " << status << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
