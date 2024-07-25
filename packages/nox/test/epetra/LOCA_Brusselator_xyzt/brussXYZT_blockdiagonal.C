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

// NOX Objects
#include "LOCA.H"
#include "LOCA_Epetra.H"
#include "NOX_Common.H"
#include "NOX_TestCompare.H"

// Trilinos Objects
#include "Epetra_Map.h"
#include "Epetra_MapColoring.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Teuchos_GlobalMPISession.hpp"

// Added to allow timings
#include "Epetra_Time.h"

// Headers needed for FD coloring
#include <vector>
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"

// User's application specific files
#include "Problem_Interface.H" // Interface file to NOX
#include "Brusselator.H"

// XYZT-specific headers
#ifdef HAVE_MPI
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif
#include "LOCA_Epetra_Interface_xyzt.H"

using namespace std;

#define MAX_NEWTON_ITERS 200
#define BETA_INIT 1.8

// Forward declarations
void getParamList(Teuchos::ParameterList *p, int MyPID);
void getPrecLSParams(Teuchos::ParameterList *p, int MyPID);
void getPrecPrintParams(Teuchos::ParameterList *p, int MyPID);

int main(int argc, char *argv[])
{
  int ierr = 0;
  int MyPID = 0;
  double alpha = 1.00; // stable steady state
  double beta = BETA_INIT;
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  int NumGlobalNodes = 10;  // default
  int spatialProcs   = 1;   // default
  int numTimeSteps   = 100;  // default

  try {

    // Initialize MPI
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);

    // Check for verbose output
    bool verbose = false;
    if (argc>1)
      if (argv[1][0]=='-' && argv[1][1]=='v')
    verbose = true;

    // Get the number of elements from the command line
    if ((argc > 2) && (verbose))
      NumGlobalNodes = atoi(argv[2]) + 1;
    else if ((argc > 1) && (!verbose))
      NumGlobalNodes = atoi(argv[1]) + 1;

    // Get the number of processors to use for the spatial problem
    if ((argc > 3) && (verbose))
      spatialProcs = atoi(argv[3]);
    else if ((argc > 2) && (!verbose))
      spatialProcs = atoi(argv[2]);

    // Get the number of processors to use for the spatial problem
    if ((argc > 4) && (verbose))
      numTimeSteps = atoi(argv[4]);
    else if ((argc > 3) && (!verbose))
      numTimeSteps = atoi(argv[3]);

    // MPI MANIPULATION FOR XYZT PROBLEMS
#ifdef HAVE_MPI
    Teuchos::RCP<EpetraExt::MultiMpiComm> globalComm =
      Teuchos::rcp(new EpetraExt::MultiMpiComm(MPI_COMM_WORLD,
                                               spatialProcs, numTimeSteps));
#else
    Teuchos::RCP<EpetraExt::MultiSerialComm> globalComm =
      Teuchos::rcp(new EpetraExt::MultiSerialComm(numTimeSteps));
#endif
    Epetra_Comm& Comm = globalComm->SubDomainComm();
    MyPID = globalComm->MyPID();

    // Get the total number of processors
    int NumProc = Comm.NumProc();

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
    Brusselator::OverlapType OType = Brusselator::ELEMENTS;
    Brusselator Problem(NumGlobalNodes, Comm, OType);

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();

    // Begin Nonlinear Solver ************************************

    // Create the top level parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);
    getParamList(&(*paramList), MyPID);

    // Sublist for "Linear Solver"
    Teuchos::ParameterList& lsParams =
      paramList->sublist("NOX").sublist("Direction").sublist("Newton").sublist("Linear Solver");

    // Sublist for "Printing"
    Teuchos::ParameterList& printParams =
      paramList->sublist("NOX").sublist("Printing");

    // Change NOX priting settings if "-verbose" not requested
    if (!verbose)
      printParams.set("Output Information", NOX::Utils::Error);

    // Create the interface between the test problem and the nonlinear solver
    Teuchos::RCP<Problem_Interface> interface =
      Teuchos::rcp(new Problem_Interface(Problem));

    // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
    Teuchos::RCP<Epetra_RowMatrix> A =
      Teuchos::rcp(&Problem.getJacobian(),false);

    // Create initial guess
    Epetra_MultiVector initGuess(soln.Map(),
                 globalComm->NumTimeStepsOnDomain());
    for (int i=0; i<globalComm->NumTimeStepsOnDomain(); i++)
      *(initGuess(i)) = soln;

    // Get XYZT preconditioner linear solver parameters
    Teuchos::RCP<Teuchos::ParameterList> precLSParams =
      Teuchos::rcp(new Teuchos::ParameterList);
    getPrecLSParams(&(*precLSParams), MyPID);

    // Get XYZT preconditioner print parameters
    Teuchos::RCP<Teuchos::ParameterList> precPrintParams =
      Teuchos::rcp(new Teuchos::ParameterList);
    getPrecPrintParams(&(*precPrintParams), MyPID);

    // Change xyztPrec priting settings if "-verbose" not requested
    if (!verbose)
      precPrintParams->set("Output Information", NOX::Utils::Error);

    // Create the XYZT object
    Teuchos::RCP<LOCA::Epetra::Interface::xyzt> ixyzt =
      Teuchos::rcp(new LOCA::Epetra::Interface::xyzt(interface,
                                                     initGuess, A,
                                                     globalComm,
                                                     soln, 0.5,
                                                     precPrintParams.get(),
                                                     precLSParams.get()));

    // Create the XYZT operator, solution, and preconditioner
    Teuchos::RCP<Epetra_RowMatrix> Axyzt =
      Teuchos::rcp(&(ixyzt->getJacobian()),false);
    Epetra_Vector& solnxyzt = ixyzt->getSolution();
    Teuchos::RCP<Epetra_Operator> Mxyzt =
      Teuchos::rcp(&(ixyzt->getPreconditioner()),false);

    // Create the Linear System
    Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = ixyzt;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = ixyzt;
    Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec =
      Teuchos::rcp(&(ixyzt->getPreconditioner()),false);
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                            iJac, Axyzt,
                            iPrec, Mxyzt,
                            solnxyzt));

    // Create full xyzt initial guess
    NOX::Epetra::Vector initialGuess(solnxyzt);

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

    Teuchos::RCP<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams,
                       iReq, initialGuess, linSys,
                       pVector));

    grp->computeF();

    // Create the convergence tests
    Teuchos::RCP<NOX::StatusTest::NormF> absresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8,
                          NOX::StatusTest::NormF::Unscaled));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(MAX_NEWTON_ITERS));
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(absresid);
    combo->addStatusTest(maxiters);

    // Create stepper
    LOCA::Stepper stepper(globalData, grp, combo, paramList);

    LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

    if (status != LOCA::Abstract::Iterator::Finished) {
      ierr =1;
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
    globalData->locaUtils->out()
      << "Stepper failed to converge!" << std::endl;
    }

    // Get the final solution from the stepper
    Teuchos::RCP<const LOCA::Epetra::Group> finalGroup =
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(stepper.getSolutionGroup());

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
    int numSteps_expected = 3;
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
    double beta_expected = 2.0;
    ierr += testCompare.testValue(beta_final, beta_expected, 1.0e-14,
                  "final value of continuation parameter",
                  NOX::TestCompare::Relative);

    // Check final value of ||F||
    double Fnorm_final = (const_cast<Teuchos::ParameterList&>(*stepper.getList()).sublist("NOX").sublist("Output").get("2-Norm of Residual",1.0e+10));
    double Fnorm_expected = 0.0;
    ierr += testCompare.testValue(Fnorm_final, Fnorm_expected, 1.0e-8,
                  "final value of ||F||",
                  NOX::TestCompare::Absolute);

    // Check number of nonlinear iterations on last continuation step
    int nonlin_final = (const_cast<Teuchos::ParameterList&>(*stepper.getList()).sublist("NOX").
            sublist("Output").get("Nonlinear Iterations",MAX_NEWTON_ITERS));
    int nonlin_expected = 4;
    ierr += testCompare.testValue(nonlin_final, nonlin_expected, 0.0,
                  "number of nonlinear iterations on last continuation step",
                  NOX::TestCompare::Absolute);


    // initialize solution comparison norm on all procs
    double solution_diff_norm = 0.0;

    // Compute norm of difference of computed solution - expected solution
    if (globalComm->MyPID() == (globalComm->NumProc()-1)) {
      // Get solution at last time step
      ixyzt->getSolution().ExtractBlockValues(soln,
                          (globalComm->NumTimeStepsOnDomain() +
                           globalComm->FirstTimeStepOnDomain() - 1));

      // Check final solution at final time step
      NOX::Epetra::Vector final_solution(soln);
      NOX::Epetra::Vector final_x_expected(final_solution);
      for (int i=0; i<NumGlobalNodes; i++) {
    final_x_expected.getEpetraVector()[2*i] = alpha;
    final_x_expected.getEpetraVector()[2*i+1] = beta_final/alpha;
      }
      solution_diff_norm = testCompare.computeVectorNorm(final_solution,
                             final_x_expected,
                             1.0e-4,
                             1.0e-4);

    }

    // Check final solution at final time step on all procs
    globalComm->Broadcast(&solution_diff_norm, 1, globalComm->NumProc()-1);
    double solution_diff_norm_expected = 0.0;
    ierr += testCompare.testValue(solution_diff_norm, solution_diff_norm_expected, 1.0e-4,
                  "inf norm of (solution_final - solution_expected) at last time step",
                  NOX::TestCompare::Absolute);

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

  return ierr;
}

void getParamList(Teuchos::ParameterList *p, int MyPID) {
  // Create LOCA sublist
  Teuchos::ParameterList& locaParamsList = p->sublist("LOCA");

  // Create the stepper sublist and set the stepper parameters
  Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
  stepperList.set("Bordered Solver Method", "Householder");
  stepperList.set("Continuation Parameter", "beta");
  stepperList.set("Initial Value", BETA_INIT);
  stepperList.set("Max Value", 2.0);
  stepperList.set("Min Value", 0.0);
  stepperList.set("Max Steps", 100);
  stepperList.set("Max Nonlinear Iterations", MAX_NEWTON_ITERS);
  stepperList.set("Compute Eigenvalues",false);

  // Create predictor sublist
  Teuchos::ParameterList& predictorList =
    locaParamsList.sublist("Predictor");
  predictorList.set("Method", "Constant");

  // Create step size sublist
  Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
  stepSizeList.set("Initial Step Size", 0.1);
  stepSizeList.set("Min Step Size", 1.0e-3);
  stepSizeList.set("Max Step Size", 10.0);
  stepSizeList.set("Aggressiveness", 0.1);

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  Teuchos::ParameterList& nlParams = p->sublist("NOX");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", MyPID);
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
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

  // Sublist for "Linear Solver"
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 800);
  lsParams.set("Tolerance", 1e-6);
  lsParams.set("Output Frequency", AZ_all);
  lsParams.set("Output Solver Details", true);
  lsParams.set("Preconditioner", "User Defined");
}

void getPrecLSParams(Teuchos::ParameterList *p, int MyPID) {
  p->set("Aztec Solver", "GMRES");
  p->set("Preconditioner", "Ifpack");
  p->set("Max Iterations", 800);
  p->set("Tolerance", 1e-6);
  p->set("Output Frequency", 50);
  p->set("XYZTPreconditioner", "BlockDiagonal");
  p->set("Use Preconditioner as Solver", false);
  p->set("Outer Use Preconditioner as Solver", false);
  p->set("Periodic", false);
}

void getPrecPrintParams(Teuchos::ParameterList *p, int MyPID) {
  p->set("MyPID", MyPID);
  p->set("Output Precision", 3);
  p->set("Output Processor", 0);
  p->set("Output Information",
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
}
