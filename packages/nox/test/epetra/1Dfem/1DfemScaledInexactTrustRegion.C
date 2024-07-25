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
#include "NOX_SolverStats.hpp"
#include "Teuchos_TestingHelpers.hpp"

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

  bool success = false;

  try {

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
    interface->setPDEfactor(100000.0);

    // Use a scaled vector space.  The scaling must also be registered
    // with the linear solver so the linear system is consistent!
    Teuchos::RCP<Epetra_Vector> scaleVec =
      Teuchos::rcp(new Epetra_Vector( *(interface->getSolution())));
    scaleVec->PutScalar(2.0);
    Teuchos::RCP<NOX::Epetra::Scaling> scaling =
      Teuchos::rcp(new NOX::Epetra::Scaling);
    scaling->addUserScaling(NOX::Epetra::Scaling::Left, scaleVec);

    // Use a weighted vector space for scaling all norms
    Teuchos::RCP<NOX::Epetra::VectorSpace> weightedVectorSpace =
      Teuchos::rcp(new NOX::Epetra::VectorSpaceScaledL2(scaling));

    // Get the vector from the Problem
    Teuchos::RCP<Epetra_Vector> soln = interface->getSolution();
    Teuchos::RCP<NOX::Epetra::Vector> noxSoln =
      Teuchos::rcp(new NOX::Epetra::Vector(soln,
            NOX::Epetra::Vector::CreateCopy,
            NOX::DeepCopy,
            weightedVectorSpace));

    // Initial Guess
    noxSoln->init(2.0);

    // Begin Nonlinear Solver ************************************

    // Create the top level parameter list
    Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

    // Set the nonlinear solver method
    nlParams.set("Nonlinear Solver", "Inexact Trust Region Based");
    nlParams.sublist("Trust Region").
      set("Inner Iteration Method", "Inexact Trust Region");

    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

    // RPP: Commenting this line out.  There is now a default for MPI
    // specific builds.  We are testing that it works here.
    // //printParams.set("MyPID", MyPID);

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

    // Sublist for direction
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Type 1");

    // Sublist for linear solver for the Newton method
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");
    lsParams.set("Max Iterations", 800);
    lsParams.set("Tolerance", 1e-4);

    // Various Preconditioner options
    //lsParams.set("Preconditioner", "AztecOO");
    lsParams.set("Preconditioner", "Ifpack");
    lsParams.set("Preconditioner Reuse Policy", "Rebuild");

    // Sublist for Cauchy direction
    Teuchos::ParameterList& cauchyDirParams = nlParams.sublist("Cauchy Direction");
    cauchyDirParams.set("Method", "Steepest Descent");
    Teuchos::ParameterList& sdParams = cauchyDirParams.sublist("Steepest Descent");
    sdParams.set("Scaling Type", "Quadratic Model Min");

    // Add a user defined pre/post operator object
    Teuchos::RCP<NOX::Abstract::PrePostOperator> ppo =
      Teuchos::rcp(new UserPrePostOperator(printing));
    nlParams.sublist("Solver Options").set("User Defined Pre/Post Operator",
        ppo);

    // Let's force all status tests to do a full check
    nlParams.sublist("Solver Options").set("Status Test Check Type", "Complete");

    // User supplied (Epetra_RowMatrix)
    Teuchos::RCP<Epetra_RowMatrix> Analytic = interface->getJacobian();

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
            iReq,
            iJac, Analytic,
            *noxSoln,
            scaling));

    // Create the Group
    Teuchos::RCP<NOX::Epetra::Group> grpPtr =
      Teuchos::rcp(new NOX::Epetra::Group(printParams,
            iReq,
            *noxSoln,
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
      Teuchos::rcp(new NOX::StatusTest::MaxIters(200));
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
#ifndef HAVE_MPI
    // 2. Linear solve iterations (14) - SERIAL TEST ONLY!
    //    The number of linear iterations changes with # of procs.
    if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Direction").sublist("Newton").sublist("Linear Solver").sublist("Output").get("Total Number of Linear Iterations",0) != 14) {
      status = 2;
    }
#endif
    // 3. Nonlinear solve iterations (17)
    if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Output").get("Nonlinear Iterations", 0) != 17)
      status = 3;
    // 4. Test the pre/post iterate options
    {
      UserPrePostOperator & ppo2 =
        dynamic_cast<UserPrePostOperator&>(*ppo.get());
      if (ppo2.getNumRunPreIterate() != 17)
        status = 4;
      if (ppo2.getNumRunPostIterate() != 17)
        status = 4;
      if (ppo2.getNumRunPreSolve() != 1)
        status = 4;
      if (ppo2.getNumRunPostSolve() != 1)
        status = 4;
    }
    // 5. Number of Cauchy steps (3)
    if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Trust Region").sublist("Output").get("Number of Cauchy Steps", 0) != 3)
      status = 5;
    // 6. Number of Newton steps (14)
    if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Trust Region").sublist("Output").get("Number of Newton Steps", 0) != 14)
      status = 6;

    success = status==0;

    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->numNonlinearSolves,1,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->numNonlinearIterations,17,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->numTotalNonlinearIterations,17,std::cout,success);

    TEUCHOS_TEST_ASSERT(solver->getSolverStatistics()->linearSolve.lastLinearSolve_Converged,std::cout,success);
    if (NumProc == 1) {
      // linear solve iterations changes with processor count, only check for serial
      TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.lastLinearSolve_NumIterations,1,std::cout,success);
      TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.lastNonlinearSolve_NumLinearIterations,14,std::cout,success);
      TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.allNonlinearSolves_NumLinearIterations,14,std::cout,success);
    }
    TEUCHOS_TEST_ASSERT(solver->getSolverStatistics()->linearSolve.lastLinearSolve_AchievedTolerance < 1.0e-4,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.lastNonlinearSolve_NumLinearSolves,14,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.allNonlinearSolves_NumLinearSolves,14,std::cout,success);

    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->trustRegion.numCauchySteps,3,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->trustRegion.numNewtonSteps,14,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->trustRegion.numDoglegSteps,0,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->trustRegion.numTrustRegionInnerIterations,17,std::cout,success);

    // Summarize test results
    if (success)
      printing.out() << "Test passed!" << std::endl;
    else
      printing.out() << "Test failed!" << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
