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
//#include "AztecOO.h"

#include "NOX_Epetra_LinearSystem_Stratimikos.H"

// User's application specific files
#include "1DfemInterface.H"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  // Initialize MPI
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  bool success = false;
  bool verbose = false;
  try {
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

    // Get the vector from the Problem
    Teuchos::RCP<Epetra_Vector> soln = interface->getSolution();
    Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
      Teuchos::rcp(new NOX::Epetra::Vector(soln,
            NOX::Epetra::Vector::CreateView));

    // Set the PDE factor (for nonlinear forcing term).  This could be specified
    // via user input.
    interface->setPDEfactor(1000.0);

    // Set the initial guess
    soln->PutScalar(1.0);

    // Begin Nonlinear Solver ************************************

    // Create the top level parameter list
    Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
      Teuchos::rcp(new Teuchos::ParameterList);
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
    searchParams.set("Method", "Full Step");

    // Sublist for direction
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");

    // Alternative linear solver list for Stratimikos
    Teuchos::ParameterList& stratLinSolParams = newtonParams.sublist("Stratimikos Linear Solver");
    Teuchos::ParameterList& noxStratParams = stratLinSolParams.sublist("NOX Stratimikos Options");
    Teuchos::ParameterList& stratParams = stratLinSolParams.sublist("Stratimikos");

    noxStratParams.set("Preconditioner Reuse Policy","Rebuild");

    stratParams.set("Linear Solver Type", "AztecOO");
    stratParams.set("Preconditioner Type", "ML");
    if (verbose) stratParams.sublist("Linear Solver Types")
      .sublist("AztecOO").sublist("Forward Solve")
        .sublist("AztecOO Settings").set("Output Frequency", 1);

    // Let's force all status tests to do a full check
    nlParams.sublist("Solver Options").set("Status Test Check Type", "Complete");

    // Create all possible Epetra_Operators.
    // 1. User supplied (Epetra_RowMatrix)
    Teuchos::RCP<Epetra_RowMatrix> Analytic = interface->getJacobian();
    // 2. Matrix-Free (Epetra_Operator)
    Teuchos::RCP<NOX::Epetra::MatrixFree> MF =
      Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface,
            *noxSoln));
    // 3. Finite Difference (Epetra_RowMatrix)
    Teuchos::RCP<NOX::Epetra::FiniteDifference> FD =
      Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface,
            *soln));

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = FD;
    Teuchos::RCP<NOX::Epetra::LinearSystem> linSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(printParams, stratLinSolParams,
            iJac, Analytic,
            iPrec, FD,
            *soln, false));

    // Create the Group
    NOX::Epetra::Vector initialGuess(soln, NOX::Epetra::Vector::CreateView);
    Teuchos::RCP<NOX::Epetra::Group> grpPtr =
      Teuchos::rcp(new NOX::Epetra::Group(printParams,
            iReq,
            initialGuess,
            linSys));
    NOX::Epetra::Group& grp = *grpPtr;

    // uncomment the following for loca supergroups
    //MF->setGroupForComputeF(*grpPtr);
    //FD->setGroupForComputeF(*grpPtr);

    // Create the convergence tests
    Teuchos::RCP<NOX::StatusTest::NormF> absresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RCP<NOX::StatusTest::NormF> relresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(grp, 1.0e-2));
    Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(1.0e-5));
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
      Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    Teuchos::RCP<NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    converged->addStatusTest(absresid);
    converged->addStatusTest(relresid);
    converged->addStatusTest(wrms);
    converged->addStatusTest(update);
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
      Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);
    NOX::StatusTest::StatusType solvStatus = solver->solve();

    // End Nonlinear Solver **************************************

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup =
      dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution =
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).
      getEpetraVector();

    // Output the parameter list
    if (verbose) {
      if (printing.isPrintType(NOX::Utils::Parameters)) {
        printing.out() << std::endl << "Final Parameters" << std::endl
          << "****************" << std::endl;
        solver->getList().print(printing.out());
        printing.out() << std::endl;
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
    if (solvStatus != NOX::StatusTest::Converged) {
      status = 1;
      if (printing.isPrintType(NOX::Utils::Error))
        printing.out() << "Nonlinear solver failed to converge!" << std::endl;
    }

    // 3. Nonlinear solve iterations (10)
    if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Output").get("Nonlinear Iterations", 0) != 10)
      status = 3;

    success = status==0;

    // Summarize test results
    if (success)
      printing.out() << "Test passed!" << std::endl;
    else
      printing.out() << "Test failed!" << std::endl;

    printing.out() << "Status = " << status << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
