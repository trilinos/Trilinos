// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

// User's application specific files
#include "Problem_Interface.H"
#include "Tcubed_FiniteElementProblem.H"

int testTransposeSolve(
         int NumGlobalElements,
         int nRHS,
         double reltol,
         double abstol,
         Epetra_Comm& Comm,
         const Teuchos::RCP<LOCA::GlobalData>& globalData,
         const Teuchos::RCP<Teuchos::ParameterList>& paramList)
{
  int ierr = 0;

  double left_bc = 0.0;
  double right_bc = 1.0;
  double nonlinear_factor = 1.0;

  // Create the FiniteElementProblem class.  This creates all required
  // Epetra objects for the problem and allows calls to the
  // function (RHS) and Jacobian evaluation routines.
  Tcubed_FiniteElementProblem Problem(NumGlobalElements, Comm);

  // Get the vector from the Problem
  Epetra_Vector& soln = Problem.getSolution();

  // Initialize Solution
  soln.PutScalar(0.0);

  // Create and initialize the parameter vector
  LOCA::ParameterVector pVector;
  pVector.addParameter("Nonlinear Factor",nonlinear_factor);
  pVector.addParameter("Left BC", left_bc);
  pVector.addParameter("Right BC", right_bc);

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

  // Get sublists
  Teuchos::ParameterList& nlParams = paramList->sublist("NOX");
  Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newParams.sublist("Linear Solver");

  // Create the linear systems
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linsys =
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(nlPrintParams,
                              lsParams, iReq, iJac,
                              Amat, soln));

  // Create the loca vector
  NOX::Epetra::Vector locaSoln(soln);

  // Create the Group
  Teuchos::RCP<LOCA::Epetra::Group> grp_tp =
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams,
                     iReq, locaSoln,
                     linsys, pVector));

  // Create the Group
  Teuchos::RCP<LOCA::Epetra::Group> grp_lp =
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams,
                     iReq, locaSoln,
                     linsys, pVector));

  // Create the Group
  Teuchos::RCP<LOCA::Epetra::Group> grp_ep =
    Teuchos::rcp(new LOCA::Epetra::Group(globalData, nlPrintParams,
                     iReq, locaSoln,
                     linsys, pVector));

  // Change initial guess to a random vector
  Teuchos::RCP<NOX::Abstract::Vector> xnew =
    grp_tp->getX().clone();
  xnew->random();
  grp_tp->setX(*xnew);
  grp_lp->setX(*xnew);
  grp_ep->setX(*xnew);

  // Check some statistics on the solution
  Teuchos::RCP<NOX::TestCompare> testCompare =
    Teuchos::rcp(new NOX::TestCompare(globalData->locaUtils->out(),
                      *(globalData->locaUtils)));

  // Evaluate blocks
  grp_tp->computeF();
  grp_lp->computeF();
  grp_ep->computeF();
  grp_tp->computeJacobian();
  grp_lp->computeJacobian();
  grp_ep->computeJacobian();

    // Set up left- and right-hand sides
  Teuchos::RCP<NOX::Abstract::MultiVector> F =
    grp_tp->getX().createMultiVector(nRHS);
  F->random();
  Teuchos::RCP<NOX::Abstract::MultiVector> X_tp =
    F->clone(NOX::ShapeCopy);
  X_tp->init(0.0);
  Teuchos::RCP<NOX::Abstract::MultiVector> X_lp =
    F->clone(NOX::ShapeCopy);
  X_lp->init(0.0);
  Teuchos::RCP<NOX::Abstract::MultiVector> X_ep =
    F->clone(NOX::ShapeCopy);
  X_ep->init(0.0);

  // Set up residuals
  Teuchos::RCP<NOX::Abstract::MultiVector> R_tp =
    F->clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> R_lp =
    F->clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> R_ep =
    F->clone(NOX::ShapeCopy);

  // Set up linear solver lists
  Teuchos::ParameterList lsParams_tp = lsParams;
  Teuchos::ParameterList lsParams_lp = lsParams;
  Teuchos::ParameterList lsParams_ep = lsParams;
  lsParams_tp.set("Transpose Solver Method",
               "Transpose Preconditioner");
  lsParams_lp.set("Transpose Solver Method",
               "Left Preconditioning");
  lsParams_ep.set("Transpose Solver Method",
               "Explicit Transpose");

  NOX::Abstract::Group::ReturnType status;

  // Test transpose solve with transposed preconditioner
  if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
    globalData->locaUtils->out() << std::endl <<
      "\t***** " <<
      "Testing Transposed Preconditioner Residual" <<
      " *****" << std::endl;

  status = grp_tp->applyJacobianTransposeInverseMultiVector(lsParams_tp, *F,
                                *X_tp);
  if (status == NOX::Abstract::Group::Failed)
    ++ierr;
  status = grp_tp->applyJacobianTransposeMultiVector(*X_tp, *R_tp);
  ierr += testCompare->testMultiVector(*R_tp, *F, reltol, abstol,
                       "Residual");

  // Test transpose solve with transposed left-preconditioner
  if (lsParams.get("Preconditioner", "None") != "None") {
    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out() << std::endl <<
    "\t***** " <<
    "Testing Transposed Left Preconditioner Residual" <<
    " *****" << std::endl;

    status = grp_lp->applyJacobianTransposeInverseMultiVector(lsParams_lp, *F,
                                  *X_lp);
    if (status == NOX::Abstract::Group::Failed)
      ++ierr;
    status = grp_lp->applyJacobianTransposeMultiVector(*X_lp, *R_lp);
    ierr += testCompare->testMultiVector(*R_lp, *F, reltol, abstol,
                     "Residual");
  }

#ifdef HAVE_NOX_EPETRAEXT
  // Test transpose solve with explicit preconditioner
  if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
    globalData->locaUtils->out() << std::endl <<
      "\t***** " <<
      "Testing Explicit Preconditioner Residual" <<
      " *****" << std::endl;

  status = grp_ep->applyJacobianTransposeInverseMultiVector(lsParams_ep, *F,
                                *X_ep);
  if (status == NOX::Abstract::Group::Failed)
    ++ierr;
  status = grp_ep->applyJacobianTransposeMultiVector(*X_ep, *R_ep);
  ierr += testCompare->testMultiVector(*R_ep, *F, reltol, abstol,
                       "Residual");
#endif


  // Compare solutions
  if (lsParams.get("Preconditioner", "None") != "None") {
    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out() << std::endl <<
    "\t***** " <<
    "Comparing Transposed Preconditioner and Left Preconditioner Solutions"
                   <<
    " *****" << std::endl;
    ierr += testCompare->testMultiVector(*X_lp, *X_tp,
                     reltol, abstol, "Solution");
  }

#ifdef HAVE_NOX_EPETRAEXT
  if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
    globalData->locaUtils->out() << std::endl <<
      "\t***** " <<
      "Comparing Transposed Preconditioner and Explicit Preconditioner Solutions"
                 <<
      " *****" << std::endl;
  ierr += testCompare->testMultiVector(*X_ep, *X_tp,
                       reltol, abstol, "Solution");
#endif

  return ierr;
}
int main(int argc, char *argv[])
{
  int ierr = 0;
  int MyPID = 0;

  int nRHS = 7;
  double reltol = 1.0e-8;
  double abstol = 1.0e-8;
  double lstol = 1.0e-11;
  int ls_its = 100;

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

    // Get the total number of processors
    MyPID = Comm.MyPID();
    int NumProc = Comm.NumProc();

    bool verbose = false;
    // Check for verbose output
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

    // Create parameter list
    Teuchos::RCP<Teuchos::ParameterList> paramList =
      Teuchos::rcp(new Teuchos::ParameterList);

    // Create the "Solver" parameters sublist to be used with NOX Solvers
    Teuchos::ParameterList& nlParams = paramList->sublist("NOX");

    Teuchos::ParameterList& nlPrintParams = nlParams.sublist("Printing");
    nlPrintParams.set("MyPID", MyPID);
    if (verbose)
       nlPrintParams.set("Output Information",
                  NOX::Utils::Error +
                  NOX::Utils::Details +
                  NOX::Utils::OuterIteration +
                  NOX::Utils::InnerIteration +
                  NOX::Utils::Warning +
                  NOX::Utils::TestDetails +
                  NOX::Utils::StepperIteration +
                  NOX::Utils::StepperDetails);
     else
       nlPrintParams.set("Output Information", NOX::Utils::Error);

    // Create the "Direction" sublist for the "Line Search Based" solver
    Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
    Teuchos::ParameterList& newParams = dirParams.sublist("Newton");
    Teuchos::ParameterList& lsParams = newParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "GMRES");
    lsParams.set("Max Iterations", ls_its);
    lsParams.set("Tolerance", lstol);
    if (verbose)
      lsParams.set("Output Frequency", 1);
    else
      lsParams.set("Output Frequency", 0);
    lsParams.set("Scaling", "None");

    //lsParams.set("Overlap", 2);
    //lsParams.set("Fill Factor", 2.0);
    //lsParams.set("Drop Tolerance", 1.0e-12);

    // Create Epetra factory
    Teuchos::RCP<LOCA::Abstract::Factory> epetraFactory =
      Teuchos::rcp(new LOCA::Epetra::Factory);

    // Create global data object
    Teuchos::RCP<LOCA::GlobalData> globalData =
      LOCA::createGlobalData(paramList, epetraFactory);

    // Test transpose solves with Ifpack preconditioner
    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out() << std::endl <<
    "********** " <<
    "Testing Transpose Solves With Ifpack Preconditioner" <<
    " **********" << std::endl;
    lsParams.set("Preconditioner", "Ifpack");
    ierr += testTransposeSolve(NumGlobalElements, nRHS, reltol, abstol,
                   Comm, globalData, paramList);

    // Test transpose solves with no preconditioner
    if (globalData->locaUtils->isPrintType(NOX::Utils::TestDetails))
      globalData->locaUtils->out() << std::endl <<
    "********** " <<
    "Testing Transpose Solves With No Preconditioner" <<
    " **********" << std::endl;
    lsParams.set("Preconditioner", "None");
    lsParams.set("Max Iterations", NumGlobalElements+1);
    ierr += testTransposeSolve(NumGlobalElements, nRHS, reltol, abstol,
                   Comm, globalData, paramList);

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
