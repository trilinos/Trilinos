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
 * --- + a * u**3 = 0
 * dx2
 *
 * subject to @ x=0, u=b
 * subject to @ x=1, u=c
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

#include "NOX_Epetra_LinearSystem_Stratimikos.H"

// User's application specific files
#include "Problem_Interface.H" // Interface file to NOX
#include "Tcubed_FiniteElementProblem.H"

using namespace std;

int tcubed_test(int NumGlobalElements, bool verbose, Epetra_Comm& Comm,
        bool includeUV, bool useP, const std::string& prec,
        const std::string& prec_method)
{
  int ierr = 0;
  int MyPID = Comm.MyPID();

  if (MyPID == 0) {
    std::cout << "********** "
     << "Testing includeUV = " << includeUV << " useP = " << useP
     << " Preconditioner = " << prec
     << " Preconditioner method = " << prec_method
     << " **********" << std::endl;
  }

  try {

    double nonlinear_factor = 1.0;
    double left_bc = 0.0;
    double right_bc = 2.07;

    // Create the FiniteElementProblem class.  This creates all required
    // Epetra objects for the problem and allows calls to the
    // function (RHS) and Jacobian evaluation routines.
    Tcubed_FiniteElementProblem Problem(NumGlobalElements, Comm);

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();

    // Initialize Solution
    soln.PutScalar(1.0);

    // Create initial guess for the null vector of jacobian
    // Only needed for Moore-Spence
    Teuchos::RCP<NOX::Abstract::Vector> nullVec =
      Teuchos::rcp(new NOX::Epetra::Vector(soln));
    nullVec->init(1.0);             // initial value 1.0

    // Begin LOCA Solver ************************************

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& locaStepperList =
      locaParamsList.sublist("Stepper");
    locaStepperList.set("Continuation Method", "Arc Length");
    locaStepperList.set("Continuation Parameter", "Nonlinear Factor");
    locaStepperList.set("Initial Value", nonlinear_factor);
    locaStepperList.set("Max Value", 2.0);
    locaStepperList.set("Min Value", 0.05);
    locaStepperList.set("Max Steps", 20);
    locaStepperList.set("Max Nonlinear Iterations", 15);

    locaStepperList.set("Bordered Solver Method", "Nested");
    Teuchos::ParameterList& nestedList =
      locaStepperList.sublist("Nested Bordered Solver");
    nestedList.set("Bordered Solver Method", "Householder");
    nestedList.set("Include UV In Preconditioner", includeUV);
    nestedList.set("Use P For Preconditioner", useP);

    // Create bifurcation sublist
    Teuchos::ParameterList& bifurcationList =
      locaParamsList.sublist("Bifurcation");
    bifurcationList.set("Type", "Turning Point");
    bifurcationList.set("Bifurcation Parameter", "Right BC");

    bifurcationList.set("Formulation", "Minimally Augmented");
    bifurcationList.set("Symmetric Jacobian", false);
    bifurcationList.set("Update Null Vectors Every Continuation Step", true);
    bifurcationList.set("Update Null Vectors Every Nonlinear Iteration",
            false);
    bifurcationList.set("Transpose Solver Method",prec);
    bifurcationList.set("Multiply Null Vectors by Mass Matrix", true);

    // The others don't seem to work with "Solve df/dp"
    if (prec == "Explicit Transpose")
      bifurcationList.set("Initial Null Vector Computation", "Solve df/dp");
    else {
      bifurcationList.set("Initial A Vector", nullVec);
      bifurcationList.set("Initial B Vector", nullVec);
    }

    bifurcationList.set("Bordered Solver Method", "Householder");
    bifurcationList.set("Include UV In Preconditioner", includeUV);
    bifurcationList.set("Use P For Preconditioner", useP);
    bifurcationList.set("Preconditioner Method", prec_method);

    //bifurcationList.set("Formulation", "Moore-Spence");
    //bifurcationList.set("Solver Method", "Phipps Bordering");
    //bifurcationList.set("Solver Method", "Salinger Bordering");
    //bifurcationList.set("Initial Null Vector", nullVec);
    //bifurcationList.set("Length Normalization Vector", nullVec);

    // Create predictor sublist
    Teuchos::ParameterList& predictorList =
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Secant");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Adaptive");
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 1.0e-3);
    stepSizeList.set("Max Step Size", 2000.0);
    stepSizeList.set("Aggressiveness", 0.1);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

    // Create the NOX printing parameter list
    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("MyPID", MyPID);
    nlPrintParams.set("Output Precision", 4);
    if (verbose)
      nlPrintParams.set("Output Information",
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
      nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create the "Linear Solver" sublist for Newton's method
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
    Teuchos::ParameterList& stratLinSolParams = newParams.sublist("Stratimikos Linear Solver");
    Teuchos::ParameterList& noxStratParams = stratLinSolParams.sublist("NOX Stratimikos Options");
    Teuchos::ParameterList& stratParams = stratLinSolParams.sublist("Stratimikos");

    noxStratParams.set("Preconditioner Reuse Policy","Rebuild");

    stratParams.set("Linear Solver Type", "AztecOO");
    stratParams.set("Preconditioner Type", "Ifpack");
    Teuchos::ParameterList& aztecParams =
      stratParams.sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve");
    Teuchos::ParameterList& aztecSettingsParams =
      aztecParams.sublist("AztecOO Settings");
    if (verbose)
      aztecSettingsParams.set("Output Frequency", 50);
    else
      stratParams.sublist("Verbose Object").set("Verbosity Level", "none");
    aztecParams.set("Max Iterations", 200);
    aztecParams.set("Tolerance", 1e-6);

    // Create and initialize the parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("Nonlinear Factor",nonlinear_factor);
    pVector.addParameter("Left BC", left_bc);
    pVector.addParameter("Right BC", right_bc);

    // Create the interface between the test problem and the nonlinear solver
    // This is created by the user using inheritance of the abstract base class
    Teuchos::RCP<Problem_Interface> interface =
      Teuchos::rcp(new Problem_Interface(Problem));
    Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> iReq = interface;
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;

    // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
    Teuchos::RCP<Epetra_RowMatrix> Amat =
      Teuchos::rcp(&Problem.getJacobian(),false);

    // Create scaling object
    Teuchos::RCP<NOX::Epetra::Scaling> scaling = Teuchos::null;
//       scaling = Teuchos::rcp(new NOX::Epetra::Scaling);
//       Teuchos::RCP<Epetra_Vector> scalingVector =
//     Teuchos::rcp(new Epetra_Vector(soln.Map()));
//       //scaling->addRowSumScaling(NOX::Epetra::Scaling::Left, scalingVector);
//       scaling->addColSumScaling(NOX::Epetra::Scaling::Right, scalingVector);

    // Create transpose scaling object
//     Teuchos::RCP<NOX::Epetra::Scaling> trans_scaling = Teuchos::null;
//     trans_scaling = Teuchos::rcp(new NOX::Epetra::Scaling);
//     Teuchos::RCP<Epetra_Vector> transScalingVector =
//       Teuchos::rcp(new Epetra_Vector(soln.Map()));
//     trans_scaling->addRowSumScaling(NOX::Epetra::Scaling::Right,
//                     transScalingVector);
//     trans_scaling->addColSumScaling(NOX::Epetra::Scaling::Left,
//                     transScalingVector);
//     //bifurcationList.set("Transpose Scaling", trans_scaling);

    // Create the linear systems
    Teuchos::RCP<NOX::Epetra::LinearSystemStratimikos> linsys =
      Teuchos::rcp(new NOX::Epetra::LinearSystemStratimikos(nlPrintParams,
                                stratLinSolParams,
                                iJac, Amat,
                                soln,
                                scaling));

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
      Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams, iReq,
                       locaSoln, linsys, linsys,
                       pVector));
    grp->computeF();

    // Create the Solver convergence test
    //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
    Teuchos::RCP<NOX::StatusTest::NormF> wrms =
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-12));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(locaStepperList.get("Max Nonlinear Iterations", 10)));
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
    int numSteps_expected = 7;
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
    double factor_final = finalGroup->getParam("Nonlinear Factor");
    double factor_expected = 2.0;
    ierr += testCompare.testValue(factor_final, factor_expected, 1.0e-14,
                  "final value of continuation parameter",
                  NOX::TestCompare::Relative);

    // Check final value of bifurcation parameter
    double right_bc_final = finalGroup->getParam("Right BC");
    double right_bc_expected = 1.47241293;
    ierr += testCompare.testValue(right_bc_final, right_bc_expected, 1.0e-7,
                  "final value of bifurcation parameter",
                  NOX::TestCompare::Relative);

    // Check norm of solution
    double norm_x = finalSolution.norm();
    double norm_x_expected = 12.038464;
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

  return ierr ;
}

int main(int argc, char *argv[])
{
  int ierr = 0;

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

  // Test includeUV = true, useP = true, prec = "Transpose Preconditioner"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, true, "Transpose Preconditioner",
              "Jacobian");

  // Test includeUV = true, useP = false, prec = "Transpose Preconditioner"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, false, "Transpose Preconditioner",
              "Jacobian");

  // Test includeUV = false, useP = true, prec = "Transpose Preconditioner"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
                 true, true, "Transpose Preconditioner",
                 "SMW");

  // Test includeUV = false, useP = false, prec = "Transpose Preconditioner"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
                 true, false, "Transpose Preconditioner",
                "SMW");

  if (NumProc > 1) {
    // Test includeUV = false, useP = true, prec = "Transpose Preconditioner"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, true, "Transpose Preconditioner",
            "Jacobian");

    // Test includeUV = false, useP = false, prec = "Transpose Preconditioner"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, false, "Transpose Preconditioner",
            "Jacobian");

    // Test includeUV = false, useP = true, prec = "Transpose Preconditioner"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, true, "Transpose Preconditioner",
            "SMW");

    // Test includeUV = false, useP = false, prec = "Transpose Preconditioner"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, false, "Transpose Preconditioner",
            "SMW");
  }

  // Test includeUV = true, useP = true, prec = "Left Preconditioning"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, true, "Left Preconditioning",
              "Jacobian");

  // Test includeUV = true, useP = false, prec = "Left Preconditioning"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, false, "Left Preconditioning",
              "Jacobian");

  // Test includeUV = true, useP = true, prec = "Left Preconditioning"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, true, "Left Preconditioning",
              "SMW");

  // Test includeUV = true, useP = false, prec = "Left Preconditioning"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, false, "Left Preconditioning",
              "SMW");

  if (NumProc > 1) {
    // Test includeUV = false, useP = true, prec = "Left Preconditioning"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, true, "Left Preconditioning",
            "Jacobian");

    // Test includeUV = false, useP = false, prec = "Left Preconditioning"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, false, "Left Preconditioning",
            "Jacobian");

    // Test includeUV = false, useP = true, prec = "Left Preconditioning"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, true, "Left Preconditioning",
            "SMW");

    // Test includeUV = false, useP = false, prec = "Left Preconditioning"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, false, "Left Preconditioning",
            "SMW");
  }

#ifdef HAVE_NOX_EPETRAEXT
  // Test includeUV = true, useP = true, prec = "Explicit Transpose"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, true, "Explicit Transpose",
              "Jacobian");

  // Test includeUV = true, useP = false, prec = "Explicit Transpose"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, false, "Explicit Transpose",
              "Jacobian");

  if (NumProc > 1) {
    // Test includeUV = false, useP = true, prec = "Explicit Transpose"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, true, "Explicit Transpose",
            "Jacobian");

    // Test includeUV = false, useP = false, prec = "Explicit Transpose"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, false, "Explicit Transpose",
            "Jacobian");

    // Test includeUV = false, useP = true, prec = "Explicit Transpose"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, true, "Explicit Transpose",
            "SMW");

    // Test includeUV = false, useP = false, prec = "Explicit Transpose"
    ierr += tcubed_test(NumGlobalElements, verbose, Comm,
            false, false, "Explicit Transpose",
            "SMW");
  }

  // Test includeUV = false, useP = false, prec = "Explicit Transpose"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, true, "Explicit Transpose",
              "SMW");

  // Test includeUV = false, useP = false, prec = "Explicit Transpose"
  ierr += tcubed_test(NumGlobalElements, verbose, Comm,
              true, false, "Explicit Transpose",
              "SMW");
#endif

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
    return ierr ;
}
