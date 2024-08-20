// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// 1D Finite Element Test Problem
/* Solves continuation problem (Parameter c="Right BC")
 *
 * d2u
 * --- + lambda * u - alpha * u**2 + beta * u**3 = 0
 * dx2
 *
 * subject to the boundar conditions u(-1) = u(1) = beta.
 */

// LOCA Objects
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX_TestCompare.H"

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
#include "Pitchfork_FiniteElementProblem.H"

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0;
  int MyPID = 0;

  try {

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
    MyPID = Comm.MyPID();
    int NumProc = Comm.NumProc();

    // Check for verbose output
    bool verbose = false;
    if (argc>1)
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
      exit(1);
    }

    int maxNewtonIters = 20;
    double alpha = 1.0;
    double beta = 0.0;
    double lambda = -2.0;

    // Create the FiniteElementProblem class.  This creates all required
    // Epetra objects for the problem and allows calls to the
    // function (RHS) and Jacobian evaluation routines.
    Pitchfork_FiniteElementProblem Problem(NumGlobalElements, Comm);

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();

    // Initialize Solution
    soln.PutScalar(0.0);

    // Begin LOCA Solver ************************************

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& locaStepperList =
      locaParamsList.sublist("Stepper");
    locaStepperList.set("Continuation Parameter", "lambda");
    locaStepperList.set("Initial Value", lambda);
    locaStepperList.set("Max Value", 4.0);
    locaStepperList.set("Min Value", -4.0);
    locaStepperList.set("Max Steps", 50);
    locaStepperList.set("Max Nonlinear Iterations", maxNewtonIters);

#ifdef HAVE_LOCA_ANASAZI
    // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
    locaStepperList.set("Compute Eigenvalues",true);
    Teuchos::ParameterList& aList = locaStepperList.sublist("Eigensolver");
    aList.set("Method", "Anasazi");
     if (!verbose)
      aList.set("Verbosity", Anasazi::Errors);
#else
    locaStepperList.set("Compute Eigenvalues",false);
#endif

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Initial Step Size", -0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 10.0);
    stepSizeList.set("Aggressiveness", 0.1);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

    // Create the NOX printing parameter list
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("MyPID", MyPID);
    if (verbose)
      nlPrintParams.set("Output Information",
            NOX::Utils::OuterIteration +
            NOX::Utils::OuterIterationStatusTest +
            NOX::Utils::InnerIteration +
            NOX::Utils::Details +
            NOX::Utils::Warning +
            NOX::Utils::TestDetails +
            NOX::Utils::Error +
            NOX::Utils::StepperIteration +
            NOX::Utils::StepperDetails +
            NOX::Utils::StepperParameters);
    else
      nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create the "Linear Solver" sublist
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
    Teuchos::ParameterList& lsParams = newParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");
    lsParams.set("Max Iterations", 100);
    lsParams.set("Tolerance", 1e-4);
    if (verbose)
      lsParams.set("Output Frequency", 50);
    else
      lsParams.set("Output Frequency", 0);
    lsParams.set("Scaling", "None");
    lsParams.set("Preconditioner", "Ifpack");

    // Create and initialize the parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("lambda",lambda);
    pVector.addParameter("alpha", alpha);
    pVector.addParameter("beta", beta);

    // Create the interface between the test problem and the nonlinear solver
    // This is created by the user using inheritance of the abstract base
    // class:
    Teuchos::RCP<Problem_Interface> interface =
      Teuchos::rcp(new Problem_Interface(Problem));
    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;

    // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
    Teuchos::RCP<Epetra_RowMatrix> Amat =
      Teuchos::rcp(&Problem.getJacobian(),false);

    // Create the linear systems
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams,
                            lsParams, iReq, iJac,
                            Amat, soln));

    // Create the loca vector
    NOX::Epetra::Vector locaSoln(soln);

    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, epetraFactory);

    // Create the Group
    Teuchos::RCP<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams,
                       iReq, locaSoln,
                       linsys, pVector));
    grp->computeF();

    // Create the Solver convergence test
    Teuchos::RCP<NOX::StatusTest::NormF> wrms =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(wrms);
    combo->addStatusTest(maxiters);

    // Create the stepper
    LOCA::Stepper stepper(globalData, grp, combo, paramList);
    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      ierr = 1;
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
    globalData->locaUtils->out()
      << "Stepper failed to converge!" << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RCP<const LOCA::Epetra::Group> finalGroup =
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(stepper.getSolutionGroup());
    const NOX::Epetra::Vector& finalSolution =
      dynamic_cast<const NOX::Epetra::Vector&>(finalGroup->getX());

    // Output the parameter list
    if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
      globalData->locaUtils->out()
    << std::endl << "Final Parameters" << std::endl
    << "****************" << std::endl;
      stepper.getList()->print(globalData->locaUtils->out());
      globalData->locaUtils->out() << std::endl;
    }

    // Check some statistics on the solution
    NOX::TestCompare testCompare(globalData->locaUtils->out(),
                 *(globalData->locaUtils));

    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out()
    << std::endl
    << "***** Checking solution statistics *****"
    << std::endl;

    // Check number of steps
    int numSteps = stepper.getStepNumber();
    int numSteps_expected = 13;
    ierr += testCompare.testValue(numSteps, numSteps_expected, 0.0,
                  "number of continuation steps",
                  NOX::TestCompare::Absolute);

    // Check number of failed steps
    int numFailedSteps = stepper.getNumFailedSteps();
    int numFailedSteps_expected = 0;
    ierr += testCompare.testValue(numFailedSteps, numFailedSteps_expected, 0.0,
                  "number of failed continuation steps",
                  NOX::TestCompare::Absolute);

    // Check final value of continuation parameter
    double alpha_final = finalGroup->getParam("lambda");
    double alpha_expected = -4.0;
    ierr += testCompare.testValue(alpha_final, alpha_expected, 1.0e-14,
                  "final value of continuation parameter",
                  NOX::TestCompare::Relative);

    // Check norm of solution
    double norm_x = finalSolution.norm();
    double norm_x_expected = 0.0;
    ierr += testCompare.testValue(norm_x, norm_x_expected, 1.0e-7,
                  "norm of final solution",
                  NOX::TestCompare::Relative);

    LOCA::destroyGlobalData(globalData);
  }

  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
    ierr = 1;
  }
  catch (const char *s) {
    std::cout << s << std::endl;
    ierr = 1;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" << std::endl;
    ierr = 1;
  }

  if (MyPID == 0) {
    if (ierr == 0)
      std::cout << "All tests passed!" << std::endl;
    else
      std::cout << ierr << " test(s) failed!" << std::endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr;
}
