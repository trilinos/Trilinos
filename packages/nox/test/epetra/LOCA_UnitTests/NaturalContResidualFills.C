// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Does a few steps of natural continuation and checks the number of
// residual fills for bug #3013
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
#include "AztecOO.h"

// User's application specific files
#include "Problem_Interface.H" // Interface file to NOX
#include "Tcubed_FiniteElementProblem.H"

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0;
  int MyPID = 0;

  try {

    // scale factor to test arc-length scaling
    double scale = 1.0;

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

    // Create the FiniteElementProblem class.  This creates all required
    // Epetra objects for the problem and allows calls to the
    // function (RHS) and Jacobian evaluation routines.
    Tcubed_FiniteElementProblem Problem(NumGlobalElements, Comm, scale);

    // Get the vector from the Problem
    Epetra_Vector& soln = Problem.getSolution();

    // Initialize Solution
    soln.PutScalar(1.0);

    // Begin LOCA Solver ************************************

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create LOCA sublist
    Teuchos::ParameterList& locaParamsList = paramList->sublist("LOCA");

    // Create the stepper sublist and set the stepper parameters
    Teuchos::ParameterList& locaStepperList =
      locaParamsList.sublist("Stepper");
    locaStepperList.set("Continuation Method", "Natural");
    locaStepperList.set("Skip Parameter Derivative", false);
    locaStepperList.set("Continuation Parameter", "Right BC");
    locaStepperList.set("Initial Value", 0.1);
    locaStepperList.set("Max Value", 0.3);
    locaStepperList.set("Min Value", 0.0);
    locaStepperList.set("Max Steps", 5);
    locaStepperList.set("Max Nonlinear Iterations", 15);
    locaStepperList.set("Skip df/dp", false);

    // Create predictor sublist
    Teuchos::ParameterList& predictorList =
      locaParamsList.sublist("Predictor");
    predictorList.set("Method", "Constant");

    // Create step size sublist
    Teuchos::ParameterList& stepSizeList = locaParamsList.sublist("Step Size");
    stepSizeList.set("Method", "Constant");
    stepSizeList.set("Initial Step Size", 0.1);
    stepSizeList.set("Min Step Size", 0.1);
    stepSizeList.set("Max Step Size", 0.1);

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
            //NOX::Utils::Details +
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
      lsParams.set("Output Frequency", 1);
    else
      lsParams.set("Output Frequency", 0);
    lsParams.set("Scaling", "None");
    lsParams.set("Preconditioner", "Ifpack");

    // Create and initialize the parameter vector
    LOCA::ParameterVector pVector;
    pVector.addParameter("Nonlinear Factor",1.0);
    pVector.addParameter("Left BC", 0.0);
    pVector.addParameter("Right BC", 0.1);

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
      Teuchos::rcp(new NOX::StatusTest::MaxIters(15));
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
    NOX::Epetra::Vector finalSolution(dynamic_cast<const NOX::Epetra::Vector&>(finalGroup->getX()));

    // Check some statistics on the solution
    NOX::TestCompare testCompare(globalData->locaUtils->out(),
                 *(globalData->locaUtils));

    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out()
    << std::endl
    << "***** Checking solution statistics *****"
    << std::endl;

    // Check number of residual fills
    int numFills = interface->numFillsF;
    int numFills_expected = 20;
    ierr += testCompare.testValue(numFills, numFills_expected, 0.0,
                  "number of residual fills",
                  NOX::TestCompare::Absolute);

    // Re-run the problem skipping df/dp computations

    locaStepperList.set("Skip Parameter Derivative", true);
    interface->numFillsF = 0;
    stepper.reset(globalData, grp, combo, paramList);
    status = stepper.run();
    if (status != LOCA::Abstract::Iterator::Finished) {
      ierr = 1;
      if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
    globalData->locaUtils->out()
      << "Stepper failed to converge!" << std::endl;
    }

    // Get the final solution from the stepper
    finalGroup =
      Teuchos::rcp_dynamic_cast<const LOCA::Epetra::Group>(stepper.getSolutionGroup());
    NOX::Epetra::Vector finalSolutionS(dynamic_cast<const NOX::Epetra::Vector&>(finalGroup->getX()));

    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out()
    << std::endl
    << "***** Checking solution statistics *****"
    << std::endl;

    // Check number of residual fills
    int numFillsS = interface->numFillsF;
    int numFillsS_expected = 9;
    ierr += testCompare.testValue(numFillsS, numFillsS_expected, 0.0,
                  "number of residual fills skipping df/dp",
                  NOX::TestCompare::Absolute);

    // Check final solution
    ierr += testCompare.testVector(finalSolutionS, finalSolution,
                   1.0e-13, 1.0e-13,
                   "final solution");

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
  return ierr ;
}
