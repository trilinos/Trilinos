// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// 1D Finite Element Brusselator Test Problem

/* Solves the nonlinear equation:
 *
 * dT       d2T
 * --- - D1 --- - alpha + (beta+1)*T - C*T**2 = 0
 * dt       dx2
 *
 * T(t,0) = T(t,1) = alpha = 0.6
 * T(0,x) = alpha + sinusoidal perturbation
 *
 *
 * dC       d2C
 * --- - D2 --- - beta*T + C*T**2 = 0
 * dt       dx2
 *
 * C(t,0) = C(t,1) = beta / alpha = 2.0 / 0.6
 * C(0,x) = beta / alpha + sinusoidal perturbation
 *
 * and
 *      D1 = D2 = 0.025
 *
 * with d representing partial differentiation.
 */

/* This version tests the "Shift-Invert 2 Matrix" eigensolver, which
 * requires a second matrix A2 passed to the shiftedLinSys and the
 * declareSeparateMatricMemory() call to the LOCA::Epetra::Group.
 *
 */

// NOX Objects
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX_Common.H"
#include "NOX_TestCompare.H"

// Trilinos Objects
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"

// Added to allow timings
#include "Epetra_Time.h"

// Headers needed for FD coloring
#include <vector>

// User's application specific files
#include "Problem_Interface.H" // Interface file to NOX
#include "Brusselator.H"

int main(int argc, char *argv[])
{
  int ierr = 0;
  int MyPID = 0;
  double alpha = 0.25;
  double beta = 1.5;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  int maxNewtonIters = 10;

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
    int NumGlobalNodes = 0;
    if ((argc > 2) && (verbose))
      NumGlobalNodes = atoi(argv[2]) + 1;
    else if ((argc > 1) && (!verbose))
      NumGlobalNodes = atoi(argv[1]) + 1;
    else
      NumGlobalNodes = 101;

    // The number of unknowns must be at least equal to the
    // number of processors.
    if (NumGlobalNodes < NumProc) {
      std::cout << "numGlobalNodes = " << NumGlobalNodes
        << " cannot be < number of processors = " << NumProc
        << std::endl;
      exit(1);
    }

    // Create the Brusselator problem class.  This creates all required
    // Epetra objects for the problem and allows calls to the
    // function (F) and Jacobian evaluation routines.
    Brusselator Problem(NumGlobalNodes, Comm);

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();

    // Begin LOCA Solver ************************************

    // Create the top level parameter list

    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Bordered Solver Method", "Householder");
    stepperList.set("Continuation Parameter", "beta");
    stepperList.set("Initial Value", beta);
    stepperList.set("Max Value", 1.6);
    stepperList.set("Min Value", 0.0);
    stepperList.set("Max Steps", 100);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

#ifdef HAVE_LOCA_ANASAZI
    // Create Anasazi Eigensolver sublist (needs --with-loca-anasazi)
    stepperList.set("Compute Eigenvalues",true);
    Teuchos::ParameterList& aList = stepperList.sublist("Eigensolver");
    aList.set("Method", "Anasazi");
    if (!verbose)
      aList.set("Verbosity", Anasazi::Errors);
    aList.set("Block Size", 1);        // Size of blocks
    aList.set("Num Blocks", 50);       // Size of Arnoldi factorization
    aList.set("Num Eigenvalues", 3);   // Number of eigenvalues
    aList.set("Convergence Tolerance", 1.0e-7);          // Tolerance
    aList.set("Step Size", 1);         // How often to check convergence
    aList.set("Maximum Restarts",1);   // Maximum number of restarts
    aList.set("Operator", "Shift-Invert 2 Matrix");
    aList.set("Shift", 0.1);
#else
    stepperList.set("Compute Eigenvalues",false);
#endif

    // Create predictor sublist
    Teuchos::ParameterList& predictorList =
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Constant");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Initial Step Size", 0.01);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 0.01);
    stepSizeList.set("Aggressiveness", 0.1);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

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
              NOX::Utils::Details +
              NOX::Utils::Warning +
              NOX::Utils::TestDetails +
              NOX::Utils::Error +
              NOX::Utils::StepperIteration +
              NOX::Utils::StepperDetails +
              NOX::Utils::StepperParameters);
     else
       printParams.set("Output Information", NOX::Utils::Error);

    // Sublist for "Linear Solver"
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");
    lsParams.set("Max Iterations", 800);
    lsParams.set("Tolerance", 1e-6);
    lsParams.set("Output Frequency", 50);
    lsParams.set("Preconditioner", "Ifpack");

    // Create the interface between the test problem and the nonlinear solver
    Teuchos::RCP<Problem_Interface> interface =
      Teuchos::rcp(new Problem_Interface(Problem));

    // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
    Teuchos::RCP<Epetra_RowMatrix> A =
      Teuchos::rcp(&Problem.getJacobian(),false);
    // Clone matrix for separate memory -- used in 2 Matrix eigensolver
    Teuchos::RCP<Epetra_RowMatrix> A2 =
      Teuchos::rcp(new Epetra_CrsMatrix(Problem.getJacobian()));

    // Use an Epetra Scaling object if desired
    Teuchos::RCP<Epetra_Vector> scaleVec =
      Teuchos::rcp(new Epetra_Vector(soln));
    NOX::Epetra::Scaling scaling;
    scaling.addRowSumScaling(NOX::Epetra::Scaling::Left, scaleVec);

    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;

    // Create the Linear System
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                            iReq, iJac, A, soln));
                                                        //&scaling);
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> shiftedLinSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                            iReq, iJac, A2, soln));

    // Create initial guess
    NOX::Epetra::Vector initialGuess(Teuchos::rcp(&soln,false),
                     NOX::Epetra::Vector::CreateView,
                     NOX::DeepCopy);

    // Create and initialize the parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("alpha",alpha);
    pVector.addParameter("beta",beta);
    pVector.addParameter("D1",D1);
    pVector.addParameter("D2",D2);

    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, epetraFactory);

    // Create the Group
    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime = interface;
    Teuchos::RCP<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams,
                       iTime, initialGuess, linSys,
                       shiftedLinSys, pVector));

    // Set flag in group to indicate that linSys and shiftedLinSys in Group constructor
    // (line above) use separate memory for the matrix. This allows algorithms to
    // compute and store 2 different matrices. Currently used by Shift-Invert 2 Matrix
    // eigensolver.
    grp->declareSeparateMatrixMemory();

    grp->computeF();

    // Create the convergence tests
    Teuchos::RCP<NOX::StatusTest::NormF> absresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,
                       NOX::StatusTest::NormF::Unscaled));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(maxNewtonIters));
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(absresid);
    combo->addStatusTest(maxiters);

    // Create stepper
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

    // Check number of continuation steps
    int numSteps = stepper.getStepNumber();
    int numSteps_expected = 11;
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
    double beta_final = finalGroup->getParam("beta");
    double beta_expected = 1.6;
    ierr += testCompare.testValue(beta_final, beta_expected, 1.0e-14,
                  "final value of continuation parameter",
                  NOX::TestCompare::Relative);

    // Check final of solution
    NOX::Epetra::Vector final_x_expected(finalSolution);
    int n = final_x_expected.getEpetraVector().MyLength()/2;
    for (int i=0; i<n; i++) {
      final_x_expected.getEpetraVector()[2*i] = alpha;
      final_x_expected.getEpetraVector()[2*i+1] = beta_final/alpha;
    }
    ierr += testCompare.testVector(finalSolution, final_x_expected,
                   1.0e-6, 1.0e-6,
                   "value of final solution");

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

return ierr;
}
