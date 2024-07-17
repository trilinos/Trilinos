// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Unit test objects
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

// NOX Objects
#include "NOX.H"
#include "NOX_Thyra.H"

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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "ModelEvaluator2DSim.hpp"

using namespace std;

TEUCHOS_UNIT_TEST(Thyra_NonlinearSolver, WithResetModel)
{
  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Problem size only supports 1 mpi process
  TEST_EQUALITY(Comm.NumProc(), 1);
  
  // Create the model evaluator object
  double d = 10.0;
  double p0 = 2.0;
  double p1 = 0.0;
  double x00 = 0.0;
  double x01 = 1.0;
  Teuchos::RCP<ModelEvaluator2DSim<double> > thyraModel =
    Teuchos::rcp(new ModelEvaluator2DSim<double>(Teuchos::rcp(&Comm,false),
						 d,p0,p1,x00,x01));
  
  ::Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p = 
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "AztecOO");
  p->set("Preconditioner Type", "Ifpack");
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");

  thyraModel->set_W_factory(lowsFactory);
  
  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");
  
  // Create a Thyra nonlinear solver
  Teuchos::RCP< ::Thyra::NonlinearSolverBase<double> > solver =
    Teuchos::rcp(new ::Thyra::NOXNonlinearSolver);

  solver->setParameterList(nl_params);
  solver->setModel(thyraModel);

  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = thyraModel->getNominalValues().get_x()->clone_v();
  
  ::Thyra::SolveCriteria<double> solve_criteria;
  ::Thyra::SolveStatus<double> solve_status;
  
  solve_status = solver->solve(initial_guess.get(), &solve_criteria);
  
  TEST_ASSERT(solve_status.extraParameters->isType<int>("Number of Iterations"));
  TEST_EQUALITY(solve_status.extraParameters->get<int>("Number of Iterations"), 7);  
  TEST_EQUALITY(solve_status.solveStatus, ::Thyra::SOLVE_STATUS_CONVERGED);

  // Test the reset capability for using a new model evalautor with a
  // different set of parameters
  p1 = 2.0;
  thyraModel =
    Teuchos::rcp(new ModelEvaluator2DSim<double>(Teuchos::rcp(&Comm,false),
						 d,p0,p1,x00,x01));
  thyraModel->set_W_factory(lowsFactory);
  
  solver->setModel(thyraModel);
  initial_guess = thyraModel->getNominalValues().get_x()->clone_v();
  solve_status = solver->solve(initial_guess.get(), &solve_criteria);
  TEST_ASSERT(solve_status.extraParameters->isType<int>("Number of Iterations"));
  TEST_EQUALITY(solve_status.extraParameters->get<int>("Number of Iterations"), 9);  
  TEST_EQUALITY(solve_status.solveStatus, ::Thyra::SOLVE_STATUS_CONVERGED);

  Teuchos::TimeMonitor::summarize();
}
