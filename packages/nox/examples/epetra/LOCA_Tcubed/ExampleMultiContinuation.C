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
#include "LOCA_MultiStepper.H"

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
#include "FiniteElementProblem.H"

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0;

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
  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Get the number of elements from the command line
  int NumGlobalElements = 100 + 1;

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
  FiniteElementProblem Problem(NumGlobalElements, Comm, scale);

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
  Teuchos::ParameterList& locaStepperList = locaParamsList.sublist("Stepper");
  locaStepperList.set("Continuation Method", "Arc Length");
  locaStepperList.set("Bordered Solver Method", "Householder");
  locaStepperList.set("Number of Continuation Parameters", 2);
  locaStepperList.set("Epsilon", 0.1);
  locaStepperList.set("Max Charts", 10000);
  locaStepperList.set("Verbosity", 1);
  locaStepperList.set("Page Charts", 1);
  locaStepperList.set("Dump Polyhedra", true);
  locaStepperList.set("Dump Centers", false);
  locaStepperList.set("Filename", "MFresults");
  locaStepperList.set("Enable Arc Length Scaling", false);
  locaStepperList.set("Max Nonlinear Iterations", 15);
  locaStepperList.set("Aggressiveness", 0.0);
  locaStepperList.set("Max Solution Component", 6.0);

  // Create sublist for each continuation parameter
  Teuchos::ParameterList& paramList1 =
    locaStepperList.sublist("Continuation Parameter 1");
  paramList1.set("Parameter Name", "Right BC");
  paramList1.set("Initial Value", 0.1);
  paramList1.set("Max Value", 4.0);
  paramList1.set("Min Value", 0.0);
  paramList1.set("Initial Step Size", 0.1);
  paramList1.set("Max Step Size", 0.2);
  paramList1.set("Min Step Size", 1.0e-3);

  Teuchos::ParameterList& paramList2 =
    locaStepperList.sublist("Continuation Parameter 2");
  paramList2.set("Parameter Name", "Nonlinear Factor");
  paramList2.set("Initial Value", 1.0);
  paramList2.set("Max Value", 4.0);
  paramList2.set("Min Value", 0.0);
  paramList2.set("Initial Step Size", 0.1);
  paramList2.set("Max Step Size", 0.2);
  paramList2.set("Min Step Size", 1.0e-3);

  // Create predictor sublist
  Teuchos::ParameterList& predictorList = locaParamsList.sublist("Predictor");
  predictorList.set("Method", "Tangent");

  // Create the "Solver" parameters sublist to be used with NOX Solvers
  Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

  // Create the NOX printing parameter list
  Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
  nlPrintParams.set("MyPID", MyPID);
  nlPrintParams.set("Output Information",
            NOX::Utils::OuterIteration +
            NOX::Utils::OuterIterationStatusTest +
            NOX::Utils::InnerIteration +
            NOX::Utils::LinearSolverDetails +
            NOX::Utils::Details +
            NOX::Utils::Warning +
            NOX::Utils::StepperIteration +
            NOX::Utils::StepperDetails +
            NOX::Utils::StepperParameters);

  // Create the "Linear Solver" sublist for Newton's method
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newParams.sublist("Linear Solver");
  lsParams.set("Aztec Solver", "GMRES");
  lsParams.set("Max Iterations", 100);
  lsParams.set("Tolerance", 1e-4);
  lsParams.set("Output Frequency", 50);
  lsParams.set("Scaling", "None");
  lsParams.set("Preconditioner", "Ifpack");

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("Right BC", 0.1);
  pVector.addParameter("Nonlinear Factor",1.0);
  pVector.addParameter("Left BC", 0.0);

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of the abstract base class:
  Teuchos::RCP<Problem_Interface> interface =
    Teuchos::rcp(new Problem_Interface(Problem));
  Teuchos::RCP<LOCA::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;

  // Create the Epetra_RowMatrixfor the Jacobian/Preconditioner
  Teuchos::RCP<Epetra_RowMatrix> Amat =
    Teuchos::rcp(&Problem.getJacobian(),false);

  // Create the linear systems
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams, lsParams,
                              iReq, iJac, Amat, soln));

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
                     locaSoln, linsys,
                     pVector));
  grp->computeF();

  // Create the Solver convergence test
  //NOX::StatusTest::NormWRMS wrms(1.0e-2, 1.0e-8);
  Teuchos::RCP<NOX::StatusTest::NormF> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(locaStepperList.get("Max Nonlinear Iterations", 10)));
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(wrms);
  combo->addStatusTest(maxiters);

  // Create the stepper
  LOCA::MultiStepper stepper(globalData, grp, combo, paramList);
  LOCA::Abstract::Iterator::IteratorStatus status = stepper.run();

  if (status == LOCA::Abstract::Iterator::Finished)
    globalData->locaUtils->out() << "All tests passed" << std::endl;
  else {
    if (globalData->locaUtils->isPrintType(NOX::Utils::Error))
      globalData->locaUtils->out()
    << "Stepper failed to converge!" << std::endl;
    }

  // Output the parameter list
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperParameters)) {
    globalData->locaUtils->out()
      << std::endl << "Final Parameters" << std::endl
      << "****************" << std::endl;
    stepper.getList()->print(globalData->locaUtils->out());
    globalData->locaUtils->out() << std::endl;
  }

  LOCA::destroyGlobalData(globalData);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
