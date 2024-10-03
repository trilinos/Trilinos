// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "ModelEvaluator1DFEM.hpp"

#include "NOX_Thyra_MatrixFreeJacobianOperator.hpp"
#include "NOX_MatrixFree_ModelEvaluatorDecorator.hpp"

using namespace std;

TEUCHOS_UNIT_TEST(NOX_Thyra_1DFEM, AnalyticJacobian_NoPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  const int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  // Create the model evaluator object
  double x00 = 0.0;
  double x01 = 1.0;
  Teuchos::RCP<ModelEvaluator1DFEM<double> > model =
    modelEvaluator1DFEM<double>(Teuchos::rcp(&Comm,false),num_elements,x00,x01);

  ::Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "AztecOO");
  p->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",20);
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = model->getNominalValues().get_x()->clone_v();

  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<double>::one());

  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model, model->create_W_op(), lowsFactory, Teuchos::null, Teuchos::null, Teuchos::null));

  nox_group->computeF();

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  Teuchos::TimeMonitor::summarize();
}

TEUCHOS_UNIT_TEST(NOX_Thyra_1DFEM, AnalyticJacobian_IfpackPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  const int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  // Create the model evaluator object
  double x00 = 0.0;
  double x01 = 1.0;
  Teuchos::RCP<ModelEvaluator1DFEM<double> > model =
    modelEvaluator1DFEM<double>(Teuchos::rcp(&Comm,false),num_elements,x00,x01);

  ::Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "AztecOO");
  p->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",20);
  p->set("Preconditioner Type", "Ifpack");
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = model->getNominalValues().get_x()->clone_v();

  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<double>::one());

  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model, model->create_W_op(), lowsFactory, Teuchos::null, Teuchos::null, Teuchos::null));

  nox_group->computeF();

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  Teuchos::TimeMonitor::summarize();
}

TEUCHOS_UNIT_TEST(NOX_Thyra_1DFEM, JFNK_NoPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  const int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  // Create the model evaluator object
  double x00 = 0.0;
  double x01 = 1.0;
  Teuchos::RCP<ModelEvaluator1DFEM<double> > model =
    modelEvaluator1DFEM<double>(Teuchos::rcp(&Comm,false),num_elements,x00,x01);

  ::Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "AztecOO");
  p->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",20);
  p->set("Preconditioner Type", "None");
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = model->getNominalValues().get_x()->clone_v();

  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<double>::one());

  // Create the JFNK operator
  Teuchos::ParameterList printParams;
  Teuchos::RCP<Teuchos::ParameterList> jfnkParams = Teuchos::parameterList();
  jfnkParams->set("Difference Type","Forward");
  jfnkParams->set("Perturbation Algorithm","KSP NOX 2001");
  jfnkParams->set("lambda",1.0e-4);
  Teuchos::RCP<NOX::Thyra::MatrixFreeJacobianOperator<double> > jfnkOp =
    Teuchos::rcp(new NOX::Thyra::MatrixFreeJacobianOperator<double>(printParams));
  jfnkOp->setParameterList(jfnkParams);
  jfnkParams->print(out);

  // Wrap the model evaluator in a JFNK Model Evaluator
  Teuchos::RCP< ::Thyra::ModelEvaluator<double> > thyraModel =
    Teuchos::rcp(new NOX::MatrixFreeModelEvaluatorDecorator<double>(model));

  // Create the NOX::Thyra::Group
  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, thyraModel, jfnkOp, lowsFactory, Teuchos::null, Teuchos::null, Teuchos::null));

  nox_group->computeF();

  // VERY IMPORTANT!!!  jfnk object needs base evaluation objects.
  // This creates a circular dependency, so use a weak pointer.
  jfnkOp->setBaseEvaluationToNOXGroup(nox_group.create_weak());

  // Create a preconditioner

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  Teuchos::TimeMonitor::summarize();
}


TEUCHOS_UNIT_TEST(NOX_Thyra_1DFEM, JFNK_UserPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  const int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  // Create the model evaluator object
  double x00 = 0.0;
  double x01 = 1.0;
  Teuchos::RCP<ModelEvaluator1DFEM<double> > model =
    modelEvaluator1DFEM<double>(Teuchos::rcp(&Comm,false),num_elements,x00,x01);

  ::Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "AztecOO");
  p->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",20);
  p->set("Preconditioner Type", "Ifpack");
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = model->getNominalValues().get_x()->clone_v();

  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<double>::one());

  // Create the JFNK operator
  Teuchos::ParameterList printParams;
  Teuchos::RCP<Teuchos::ParameterList> jfnkParams = Teuchos::parameterList();
  jfnkParams->set("Difference Type","Forward");
  jfnkParams->set("Perturbation Algorithm","KSP NOX 2001");
  jfnkParams->set("lambda",1.0e-4);
  Teuchos::RCP<NOX::Thyra::MatrixFreeJacobianOperator<double> > jfnkOp =
    Teuchos::rcp(new NOX::Thyra::MatrixFreeJacobianOperator<double>(printParams));
  jfnkOp->setParameterList(jfnkParams);
  jfnkParams->print(out);

  // Wrap the model evaluator in a JFNK Model Evaluator
  Teuchos::RCP< ::Thyra::ModelEvaluator<double> > thyraModel =
    Teuchos::rcp(new NOX::MatrixFreeModelEvaluatorDecorator<double>(model));

  // Create the Preconditioner operator
  Teuchos::RCP< ::Thyra::PreconditionerBase<double> > precOp =
    thyraModel->create_W_prec();

  // Create the NOX::Thyra::Group
  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, thyraModel, jfnkOp, lowsFactory, precOp, Teuchos::null, Teuchos::null));

  nox_group->computeF();

  // VERY IMPORTANT!!!  jfnk object needs base evaluation objects.
  // This creates a circular dependency, so use a weak pointer.
  jfnkOp->setBaseEvaluationToNOXGroup(nox_group.create_weak());

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  Teuchos::TimeMonitor::summarize();
}

TEUCHOS_UNIT_TEST(NOX_Thyra_1DFEM, AndersonAcceleration_UserPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  const int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  // Create the model evaluator object
  double x00 = 0.0;
  double x01 = 1.0;
  Teuchos::RCP<ModelEvaluator1DFEM<double> > model =
    modelEvaluator1DFEM<double>(Teuchos::rcp(&Comm,false),num_elements,x00,x01);

  // Create the initial guess
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = model->getNominalValues().get_x()->clone_v();

  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<double>::one());

  // Create the Preconditioner operator
  Teuchos::RCP< ::Thyra::PreconditionerBase<double> > precOp =
    model->create_W_prec();

  // Create the NOX::Thyra::Group
  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model, Teuchos::null, Teuchos::null, precOp, Teuchos::null, Teuchos::null));

  //nox_group->computeF();

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(2000));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params = Teuchos::parameterList();
  nl_params->set("Nonlinear Solver", "Anderson Accelerated Fixed-Point");
  nl_params->sublist("Anderson Parameters").set("Storage Depth", 100);
  nl_params->sublist("Anderson Parameters").set("Mixing Parameter", -0.9);
  nl_params->sublist("Anderson Parameters").sublist("Preconditioning").set("Precondition", true);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_EQUALITY(solvStatus,NOX::StatusTest::Converged);
  TEST_EQUALITY(solver->getNumIterations(),103);

  Teuchos::TimeMonitor::summarize();
}


TEUCHOS_UNIT_TEST(NOX_Thyra_1DFEM, AndersonAcceleration_IfpackPrec)
{
  Teuchos::TimeMonitor::zeroOutTimers();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  const int num_elements = 100;

  // Make sure there is at least one element on each core
  TEST_ASSERT(Comm.NumProc() < num_elements+1);

  // Create the model evaluator object
  double x00 = 0.0;
  double x01 = 1.0;
  Teuchos::RCP<ModelEvaluator1DFEM<double> > model =
    modelEvaluator1DFEM<double>(Teuchos::rcp(&Comm,false),num_elements,x00,x01);

  ::Stratimikos::DefaultLinearSolverBuilder builder;

  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  p->set("Linear Solver Type", "AztecOO");
  p->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",20);
  p->set("Preconditioner Type", "Ifpack");
  builder.setParameterList(p);

  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");

  model->set_W_factory(lowsFactory);

  // Create the initial guess
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = model->getNominalValues().get_x()->clone_v();

  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<double>::one());

  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model, Teuchos::null));

  nox_group->computeF();

  // Create the NOX status tests and the solver
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
    Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  converged->addStatusTest(absresid);
  converged->addStatusTest(wrms);
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params = Teuchos::parameterList();
  nl_params->set("Nonlinear Solver", "Anderson Accelerated Fixed-Point");
  nl_params->sublist("Anderson Parameters").set("Storage Depth", 100);
  nl_params->sublist("Anderson Parameters").set("Mixing Parameter", -0.9);
  nl_params->sublist("Anderson Parameters").sublist("Preconditioning").set("Precondition", true);
  nl_params->sublist("Printing").sublist("Output Information").set("Details",true);
  nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration",true);
  //nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration StatusTest",true);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  TEST_ASSERT(solvStatus == NOX::StatusTest::Converged);

  // Ifpack preconditioning affects convergence!
  if (Comm.NumProc() == 1) {
    TEST_EQUALITY(solver->getNumIterations(),5);
  }
  else if (Comm.NumProc() == 4) {
    TEST_EQUALITY(solver->getNumIterations(),15);
  }

  Teuchos::TimeMonitor::summarize();
}
