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
#include "NOX_SolverStats.hpp"

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
#include "ModelEvaluator1DPoisson.hpp"

TEUCHOS_UNIT_TEST(SingleStepSolver, WithResetModel)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Create the model evaluator object
  int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  double a = 1.0;
  double b = -1.0;
  Teuchos::RCP<ModelEvaluator1DPoisson<double> > thyraModel =
    Teuchos::rcp(new ModelEvaluator1DPoisson<double>(Teuchos::rcp(&Comm,false),
          num_elements,a,b));

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
  nl_params->set("Nonlinear Solver", "Single Step");
  nl_params->sublist("Single Step Solver").set("Print Norms", true);

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
  TEST_EQUALITY(solve_status.extraParameters->get<int>("Number of Iterations"), 1);
  TEST_EQUALITY(solve_status.solveStatus, ::Thyra::SOLVE_STATUS_CONVERGED);

  nl_params->print(std::cout);

  Teuchos::TimeMonitor::summarize();
}


TEUCHOS_UNIT_TEST(SingleStepSolver, ignoreLinearSolverFailures)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Create the model evaluator object
  int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  double a = 1.0;
  double b = -1.0;
  Teuchos::RCP<ModelEvaluator1DPoisson<double> > thyraModel =
    Teuchos::rcp(new ModelEvaluator1DPoisson<double>(Teuchos::rcp(&Comm,false),
          num_elements,a,b));

  ::Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "Belos");
  p->set("Preconditioner Type", "None");
  // Pick params to make sure linear solver fails (we are testing to ignore this failure at the nonlinear solver level)
  p->sublist("Linear Solver Types").sublist("Belos").set("Solver Type","Block GMRES");
  p->sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES").set("Convergence Tolerance",1.0e-8);
  p->sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES").set("Maximum Iterations",2);
  p->sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES").set("Verbosity",33);
  p->sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES").set("Output Frequency",1);
  p->sublist("Linear Solver Types").sublist("Belos").sublist("VerboseObject").set("Verbosity Level","medium");
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");

  thyraModel->set_W_factory(lowsFactory);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Single Step");
  nl_params->sublist("Single Step Solver").set("Ignore Linear Solver Failures", true);

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
  TEST_EQUALITY(solve_status.extraParameters->get<int>("Number of Iterations"), 1);
  TEST_EQUALITY(solve_status.solveStatus, ::Thyra::SOLVE_STATUS_CONVERGED);

  Teuchos::TimeMonitor::summarize();
}


TEUCHOS_UNIT_TEST(SingleStepSolver, reuseJacobian)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // Create the model evaluator object
  int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  double a = 1.0;
  double b = -1.0;
  Teuchos::RCP<ModelEvaluator1DPoisson<double> > thyraModel =
    Teuchos::rcp(new ModelEvaluator1DPoisson<double>(Teuchos::rcp(&Comm,false),
          num_elements,a,b));

  ::Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "Belos");
  p->set("Preconditioner Type", "None");
  p->sublist("Linear Solver Types").sublist("Belos").set("Solver Type","Block GMRES");
  p->sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES").set("Convergence Tolerance",1.0e-8);
  p->sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES").set("Maximum Iterations",100);
  p->sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES").set("Verbosity",33);
  p->sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Block GMRES").set("Output Frequency",1);
  p->sublist("Linear Solver Types").sublist("Belos").sublist("VerboseObject").set("Verbosity Level","medium");
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");

  thyraModel->set_W_factory(lowsFactory);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Single Step");
  nl_params->sublist("Single Step Solver").set("Ignore Linear Solver Failures", false);
  nl_params->sublist("Single Step Solver").set("Update Jacobian", false);
  nl_params->sublist("Single Step Solver").set("Print Norms", true);
  nl_params->sublist("Single Step Solver").set("Compute Relative Norm", true);
  nl_params->sublist("Single Step Solver").sublist("Linear Solver").set("Tolerance", 1.0e-7);

  // Create a Thyra nonlinear solver
  Teuchos::RCP< ::Thyra::NOXNonlinearSolver> solver =
    Teuchos::rcp(new ::Thyra::NOXNonlinearSolver);

  solver->setParameterList(nl_params);
  solver->setModel(thyraModel);

  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = thyraModel->getNominalValues().get_x()->clone_v();

  // This test exists to check that we can take an externally
  // evaluated Jacobian and set it on the solver and reuse it in the
  // SingleStep solver. So let's create the Jacobian, evaluate it and
  // set it on the solver.
  auto inArgs = thyraModel->createInArgs();
  auto outArgs = thyraModel->createOutArgs();
  auto J = thyraModel->create_W_op();
  inArgs.set_x(initial_guess);
  outArgs.set_W_op(J);
  thyraModel->evalModel(inArgs,outArgs);
  const bool jacobianIsEvaluated = true;
  Teuchos::RCP<NOX::Thyra::Group> group = 
    Teuchos::rcp(new NOX::Thyra::Group(initial_guess,
                                       thyraModel,
                                       J,
                                       lowsFactory,
                                       Teuchos::null,
                                       Teuchos::null,
                                       Teuchos::null,
                                       Teuchos::null,
                                       Teuchos::null,
                                       false,
                                       false,
                                       jacobianIsEvaluated));
  TEST_ASSERT(group->isJacobian());
  solver->setGroup(group); // Set specially constructed group with Jacobian already evaluated

  ::Thyra::SolveCriteria<double> solve_criteria;
  ::Thyra::SolveStatus<double> solve_status;

  solve_status = solver->solve(initial_guess.get(), &solve_criteria);

  TEST_ASSERT(solve_status.extraParameters->isType<int>("Number of Iterations"));
  TEST_EQUALITY(solve_status.extraParameters->get<int>("Number of Iterations"), 1);
  TEST_EQUALITY(solve_status.solveStatus, ::Thyra::SOLVE_STATUS_CONVERGED);
  // Value is 50 for 4 mpi processes. Pad for different process counts
  TEST_ASSERT(solver->getNOXSolver()->getSolverStatistics()->linearSolve.lastLinearSolve_NumIterations < 55);
  TEST_EQUALITY(solver->getNOXSolver()->getSolverStatistics()->linearSolve.lastLinearSolve_Converged, true);
  TEST_ASSERT(solver->getNOXSolver()->getSolverStatistics()->linearSolve.lastLinearSolve_AchievedTolerance < 1.0e-7);

  if (Comm.MyPID() == 0)
    std::cout << "Final Parameters\n****************\n" << *nl_params << std::endl; 

  Teuchos::TimeMonitor::summarize();
}
