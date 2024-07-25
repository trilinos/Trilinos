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
  double pi = 4.0*atan(1.0);
  double D1 = 1.0/40.0;
  double D2 = 1.0/40.0;
  double alpha = 0.25;
  double beta = 2.0*pi*pi*D1 + 1.0 + alpha*alpha;
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
    Problem.setParameters(alpha, beta, D1, D2);
    Problem.initializeSoln();

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();
    NOX::Epetra::Vector initialGuess(Teuchos::rcp(&soln,false),
                     NOX::Epetra::Vector::CreateView,
                     NOX::DeepCopy);

    // Create initial eigenvectors
    Epetra_Vector& x = Problem.getMesh();
    NOX::Epetra::Vector y(initialGuess);
    NOX::Epetra::Vector z(initialGuess);
    double lambda_real = (beta - 1.0 - alpha*alpha)/2.0;
    double lambda_imag = std::sqrt(alpha*alpha - lambda_real*lambda_real);
    double v1_real = -alpha*alpha;
    double v1_imag = 0.0;
    double v2_real = beta - 1.0 - lambda_real;
    double v2_imag = -lambda_imag;
    for (int i=0; i<x.MyLength(); i++) {
      y.getEpetraVector()[2*i]   = v1_real*sin(1.0*pi*x[i]);
      y.getEpetraVector()[2*i+1] = v2_real*sin(1.0*pi*x[i]);
      z.getEpetraVector()[2*i]   = v1_imag*sin(1.0*pi*x[i]);
      z.getEpetraVector()[2*i+1] = v2_imag*sin(1.0*pi*x[i]);
    }

    // Initial guess for frequency (valid for |alpha| > (pi^2)*|D1|)
    double w = lambda_imag;

    // Create length scaling vector (phi)
    NOX::Epetra::Vector phi(initialGuess);
    phi.init(1.0);

    Teuchos::RCP<NOX::Abstract::Vector> y_vec =
      Teuchos::rcp(&y,false);
    Teuchos::RCP<NOX::Abstract::Vector> z_vec =
      Teuchos::rcp(&z,false);
    Teuchos::RCP<NOX::Abstract::Vector> phi_vec =
      Teuchos::rcp(&phi,false);

    // Create initial values for a and b for minimally augmented method
    Teuchos::RCP<NOX::Abstract::Vector> a_vec_real =
      Teuchos::rcp(new NOX::Epetra::Vector(initialGuess));
    Teuchos::RCP<NOX::Abstract::Vector> a_vec_imag =
      Teuchos::rcp(new NOX::Epetra::Vector(initialGuess));
    Teuchos::RCP<NOX::Abstract::Vector> b_vec_real =
      Teuchos::rcp(new NOX::Epetra::Vector(initialGuess));
    Teuchos::RCP<NOX::Abstract::Vector> b_vec_imag =
      Teuchos::rcp(new NOX::Epetra::Vector(initialGuess));
    *a_vec_real = *y_vec;
    *a_vec_imag = *z_vec;
    *b_vec_real = *y_vec;
    *b_vec_imag = *z_vec;

    // Begin LOCA Solver ************************************

    // Create the top level parameter list

    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& stepperList = locaParamsList.sublist("Stepper");
    stepperList.set("Continuation Parameter", "alpha");
    stepperList.set("Initial Value", alpha);
    stepperList.set("Max Value", 1.0);
    stepperList.set("Min Value", 0.24);
    stepperList.set("Max Steps", 100);
    stepperList.set("Max Nonlinear Iterations", maxNewtonIters);

     // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList =
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Hopf");
    bifurcationList.set("Bifurcation Parameter", "beta");
    bifurcationList.set("Initial Frequency", w);

//     // For Moore-Spence Formulation
//     bifurcationList.set("Formulation", "Moore-Spence");
//     bifurcationList.set("Solver Method", "Salinger Bordering");
//     bifurcationList.set("Length Normalization Vector", phi_vec);
//     bifurcationList.set("Initial Real Eigenvector", y_vec);
//     bifurcationList.set("Initial Imaginary Eigenvector", z_vec);

    // For minimally augmented formulation
    bifurcationList.set("Formulation", "Minimally Augmented");
    bifurcationList.set("Initial Real A Vector", a_vec_real);       // Must set
    bifurcationList.set("Initial Imaginary A Vector", a_vec_imag);  // Must set
    bifurcationList.set("Initial Real B Vector", b_vec_real);       // Must set
    bifurcationList.set("Initial Imaginary B Vector", b_vec_imag);  // Must set
    bifurcationList.set("Update Null Vectors Every Continuation Step", true);

    // For minimally augmented method, should set these for good performance
    // Direct solve of bordered equations
    bifurcationList.set("Bordered Solver Method",  "Householder");
    bifurcationList.set("Include UV In Preconditioner", true);

    // Combine arc-length and turning point bordered rows & columns
    stepperList.set("Bordered Solver Method", "Nested");
    Teuchos::ParameterList& nestedList =
      stepperList.sublist("Nested Bordered Solver");
    // Direct solve of combined bordered system
    nestedList.set("Bordered Solver Method", "Householder");
    nestedList.set("Include UV In Preconditioner", true);

    // Create predictor sublist
    Teuchos::ParameterList& predictorList =
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");

//     Teuchos::ParameterList& firstStepPredictor =
//       predictorList.sublist("First Step Predictor");
//     firstStepPredictor.set("Method", "Random");
//     firstStepPredictor.set("Epsilon", 1.0e-3);
//     Teuchos::ParameterList& lastStepPredictor =
//       predictorList.sublist("Last Step Predictor");
//     lastStepPredictor.set("Method", "Random");
//     lastStepPredictor.set("Epsilon", 1.0e-3);

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Initial Step Size", 0.02);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 1.0);
    stepSizeList.set("Aggressiveness", 0.5);

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
              NOX::Utils::LinearSolverDetails +
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

    // Create shifted Linear System
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> shiftedLinSys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
                            iReq, iJac, A, soln));
                                                        //&scaling);

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
    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iTime =
      interface;
    Teuchos::RCP<LOCA::Epetra::Group> grp =
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, printParams,
                       iTime, initialGuess, linSys,
                       shiftedLinSys,
                       pVector));
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
    int numSteps_expected = 12;
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
    double alpha_final = finalGroup->getParam("alpha");
    double alpha_expected = 1.0;
    ierr += testCompare.testValue(alpha_final, alpha_expected, 1.0e-14,
                  "final value of continuation parameter",
                  NOX::TestCompare::Relative);

    // Check final value of bifurcation parameter
    double beta_final = finalGroup->getParam("beta");
    double beta_expected = 5.279;
    ierr += testCompare.testValue(beta_final, beta_expected, 1.0e-3,
                  "final value of bifurcation parameter",
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

return ierr;
}
