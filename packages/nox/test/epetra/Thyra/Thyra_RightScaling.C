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
#include "AztecOO.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_ScaledModelEvaluator.hpp"
#include "ModelEvaluator2DSim.hpp"

#include "NOX_PrePostOperator_RowSumScaling.H"
#include "Thyra_VectorBase.hpp"

using namespace std;

TEUCHOS_UNIT_TEST(NOX_Thyra_RightScaling, CTOR1)
{
  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  
  // Check we have only one processor since this problem doesn't work
  // for more than one proc
  TEST_ASSERT(Comm.NumProc() == 1);
  
  // Create the model evaluator object
  double d = 1.0;
  double p0 = 1.0e6;
  double p1 = -1.0e-6;
  double x00 = 0.0;
  double x01 = 0.0;
  Teuchos::RCP<ModelEvaluator2DSim<double> > thyraModel =
    Teuchos::rcp(new ModelEvaluator2DSim<double>(Teuchos::rcp(&Comm,false),
                                                 d,p0,p1,x00,x01));
  
  ::Stratimikos::DefaultLinearSolverBuilder builder;
  
  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  {
    p->set("Linear Solver Type", "AztecOO");
    //p->set("Preconditioner Type", "Ifpack");
    p->set("Preconditioner Type", "None");
    Teuchos::ParameterList& az = p->sublist("Linear Solver Types").sublist("AztecOO");
    az.sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency", 1);
    az.sublist("VerboseObject").set("Verbosity Level", "high");
    Teuchos::ParameterList& ip = p->sublist("Preconditioner Types").sublist("Ifpack");
    ip.sublist("VerboseObject").set("Verbosity Level", "high");
  }
  
  builder.setParameterList(p);
  
  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");
  
  thyraModel->set_W_factory(lowsFactory);
  
  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Line Search").set("Method", "Full Step");

  Teuchos::ParameterList& printParams = nl_params->sublist("Printing");
  printParams.set("Output Information",
                  NOX::Utils::OuterIteration +
                  NOX::Utils::OuterIterationStatusTest +
                  NOX::Utils::InnerIteration +
                  NOX::Utils::LinearSolverDetails +
                  NOX::Utils::Parameters +
                  NOX::Utils::Details +
                  NOX::Utils::Warning +
                  NOX::Utils::Debug +
                  NOX::Utils::TestDetails +
                  NOX::Utils::Error);
  
  nl_params->sublist("Solver Options").set("Status Test Check Type", "Complete");
  
  // Enable row sum scaling
  nl_params->sublist("Thyra Group Options").set("Function Scaling", "None");
  
  // Create right scaling vector
  Teuchos::RCP<Thyra::VectorBase<double> > values = Thyra::createMember(thyraModel->get_x_space());
  Thyra::put_scalar(0.0,values.ptr());
  
  // Values chosen such that WRMS convergence is achieved in 1 step
  // Unscaled will not converge in less than 20 steps
  Thyra::set_ele(0,1.0e14,values.ptr());
  Thyra::set_ele(1,1.0e2,values.ptr());
  Teuchos::RCP<NOX::Abstract::Vector > right_scaling = Teuchos::rcp<NOX::Abstract::Vector>(new NOX::Thyra::Vector(values));

  // Insert scaling vector into NOX parameter list
  Teuchos::ParameterEntry entry(right_scaling);
  nl_params->setEntry("Right Scaling Vector", entry);
  
  // Create Status Tests
  {
    Teuchos::ParameterList& st = nl_params->sublist("Status Tests");
    st.set("Test Type", "Combo");
    st.set("Combo Type", "OR");
    st.set("Number of Tests", 3);
    
    {
      Teuchos::ParameterList& normWRMS = st.sublist("Test 0");
      normWRMS.set("Test Type", "NormWRMS");
      normWRMS.set("Absolute Tolerance", 2.0e-8);
      normWRMS.set("Relative Tolerance", 1.0e-5);
      normWRMS.set("Tolerance", 1.0);
      normWRMS.set("BDF Multiplier", 1.0);
      normWRMS.set("Alpha", 1.0);
      normWRMS.set("Beta", 0.5);
      normWRMS.set("Disable Implicit Weighting", true);
    }
    
    {
      Teuchos::ParameterList& fv = st.sublist("Test 1");
      fv.set("Test Type", "FiniteValue");
      fv.set("Vector Type", "F Vector");
      fv.set("Norm Type", "Two Norm");
    }

    {
      Teuchos::ParameterList& maxiters = st.sublist("Test 2");
      maxiters.set("Test Type", "MaxIters");
      maxiters.set("Maximum Iterations", 20);
    }
    
  }
  
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
  
  // Test total num iterations
  {
    // Problem converges in >20 nonlinear iterations with NO scaling
    // Problem converges in 1 nonlinear iterations with scaling
    Teuchos::RCP< ::Thyra::NOXNonlinearSolver> thyra_nox_solver =
      Teuchos::rcp_dynamic_cast< ::Thyra::NOXNonlinearSolver>(solver);
    TEST_EQUALITY(thyra_nox_solver->getNOXSolver()->getNumIterations(), 1);
  }
}


TEUCHOS_UNIT_TEST(NOX_Thyra_RightScaling, CTOR1_withRowSum)
{
  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  
  // Check we have only one processor since this problem doesn't work
  // for more than one proc
  TEST_ASSERT(Comm.NumProc() == 1);
  
  // Create the model evaluator object
  double d = 1.0;
  double p0 = 1.0e6;
  double p1 = -1.0e-6;
  double x00 = 0.0;
  double x01 = 0.0;
  Teuchos::RCP<ModelEvaluator2DSim<double> > thyraModel =
    Teuchos::rcp(new ModelEvaluator2DSim<double>(Teuchos::rcp(&Comm,false),
                                                 d,p0,p1,x00,x01));
  
  ::Stratimikos::DefaultLinearSolverBuilder builder;
  
  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  {
    p->set("Linear Solver Type", "AztecOO");
    //p->set("Preconditioner Type", "Ifpack");
    p->set("Preconditioner Type", "None");
    Teuchos::ParameterList& az = p->sublist("Linear Solver Types").sublist("AztecOO");
    az.sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency", 1);
    az.sublist("VerboseObject").set("Verbosity Level", "high");
    Teuchos::ParameterList& ip = p->sublist("Preconditioner Types").sublist("Ifpack");
    ip.sublist("VerboseObject").set("Verbosity Level", "high");
  }
  
  builder.setParameterList(p);
  
  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");
  
  thyraModel->set_W_factory(lowsFactory);
  
  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Line Search").set("Method", "Full Step");

  Teuchos::ParameterList& printParams = nl_params->sublist("Printing");
  printParams.set("Output Information",
                  NOX::Utils::OuterIteration +
                  NOX::Utils::OuterIterationStatusTest +
                  NOX::Utils::InnerIteration +
                  NOX::Utils::LinearSolverDetails +
                  NOX::Utils::Parameters +
                  NOX::Utils::Details +
                  NOX::Utils::Warning +
                  NOX::Utils::Debug +
                  NOX::Utils::TestDetails +
                  NOX::Utils::Error);
  
  nl_params->sublist("Solver Options").set("Status Test Check Type", "Complete");
  
  // Enable row sum scaling
  nl_params->sublist("Thyra Group Options").set("Function Scaling", "Row Sum");
  
  // Create right scaling vector
  Teuchos::RCP<Thyra::VectorBase<double> > values = Thyra::createMember(thyraModel->get_x_space());
  Thyra::put_scalar(0.0,values.ptr());
  
  // Values chosen such that WRMS convergence is achieved in 1 step
  // Unscaled will not converge in less than 20 steps
  Thyra::set_ele(0,1.0e14,values.ptr());
  Thyra::set_ele(1,1.0e2,values.ptr());
  Teuchos::RCP<NOX::Abstract::Vector > right_scaling = Teuchos::rcp<NOX::Abstract::Vector>(new NOX::Thyra::Vector(values));

  // Insert scaling vector into NOX parameter list
  Teuchos::ParameterEntry entry(right_scaling);
  nl_params->setEntry("Right Scaling Vector", entry);
  nl_params->sublist("Thyra Group Options").set<bool>("Do Right Scaling First", true);
  
  // Create Status Tests
  {
    Teuchos::ParameterList& st = nl_params->sublist("Status Tests");
    st.set("Test Type", "Combo");
    st.set("Combo Type", "OR");
    st.set("Number of Tests", 3);
    
    {
      Teuchos::ParameterList& normWRMS = st.sublist("Test 0");
      normWRMS.set("Test Type", "NormWRMS");
      normWRMS.set("Absolute Tolerance", 2.0e-8);
      normWRMS.set("Relative Tolerance", 1.0e-5);
      normWRMS.set("Tolerance", 1.0);
      normWRMS.set("BDF Multiplier", 1.0);
      normWRMS.set("Alpha", 1.0);
      normWRMS.set("Beta", 0.5);
      normWRMS.set("Disable Implicit Weighting", true);
    }
    
    {
      Teuchos::ParameterList& fv = st.sublist("Test 1");
      fv.set("Test Type", "FiniteValue");
      fv.set("Vector Type", "F Vector");
      fv.set("Norm Type", "Two Norm");
    }

    {
      Teuchos::ParameterList& maxiters = st.sublist("Test 2");
      maxiters.set("Test Type", "MaxIters");
      maxiters.set("Maximum Iterations", 20);
    }
    
  }
  
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
  
  // Test total num iterations
  {
    // Problem converges in >20 nonlinear iterations with NO scaling
    // Problem converges in 1 nonlinear iterations with scaling
    Teuchos::RCP< ::Thyra::NOXNonlinearSolver> thyra_nox_solver =
      Teuchos::rcp_dynamic_cast< ::Thyra::NOXNonlinearSolver>(solver);
    TEST_EQUALITY(thyra_nox_solver->getNOXSolver()->getNumIterations(), 1);
  }
}


TEUCHOS_UNIT_TEST(NOX_Thyra_RightScaling, CTOR2)
{
  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif
  
  // Check we have only one processor since this problem doesn't work
  // for more than one proc
  TEST_ASSERT(Comm.NumProc() == 1);
  
  // Create the model evaluator object
  double d = 1.0;
  double p0 = 1.0e6;
  double p1 = -1.0e-6;
  double x00 = 0.0;
  double x01 = 0.0;
  Teuchos::RCP<ModelEvaluator2DSim<double> > thyraModel =
    Teuchos::rcp(new ModelEvaluator2DSim<double>(Teuchos::rcp(&Comm,false),
                                                 d,p0,p1,x00,x01));
  
  ::Stratimikos::DefaultLinearSolverBuilder builder;
  
  Teuchos::RCP<Teuchos::ParameterList> p =
    Teuchos::rcp(new Teuchos::ParameterList);
  {
    p->set("Linear Solver Type", "AztecOO");
    //p->set("Preconditioner Type", "Ifpack");
    p->set("Preconditioner Type", "None");
    Teuchos::ParameterList& az = p->sublist("Linear Solver Types").sublist("AztecOO");
    az.sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency", 1);
    az.sublist("VerboseObject").set("Verbosity Level", "high");
    Teuchos::ParameterList& ip = p->sublist("Preconditioner Types").sublist("Ifpack");
    ip.sublist("VerboseObject").set("Verbosity Level", "high");
  }
  
  builder.setParameterList(p);
  
  Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
    lowsFactory = builder.createLinearSolveStrategy("");
  
  thyraModel->set_W_factory(lowsFactory);
  
  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->sublist("Line Search").set("Method", "Full Step");

  Teuchos::ParameterList& printParams = nl_params->sublist("Printing");
  printParams.set("Output Information",
                  NOX::Utils::OuterIteration +
                  NOX::Utils::OuterIterationStatusTest +
                  NOX::Utils::InnerIteration +
                  NOX::Utils::LinearSolverDetails +
                  NOX::Utils::Parameters +
                  NOX::Utils::Details +
                  NOX::Utils::Warning +
                  NOX::Utils::Debug +
                  NOX::Utils::TestDetails +
                  NOX::Utils::Error);
  
  nl_params->sublist("Solver Options").set("Status Test Check Type", "Complete");
  
  // Enable row sum scaling
  nl_params->sublist("Thyra Group Options").set("Function Scaling", "None");
  
  // Create right scaling vector
  Teuchos::RCP<Thyra::VectorBase<double> > values = Thyra::createMember(thyraModel->get_x_space());
  Thyra::put_scalar(0.0,values.ptr());
  
  // Values chosen such that WRMS convergence is achieved in 1 step
  // Unscaled will not converge in less than 20 steps
  Thyra::set_ele(0,1.0e14,values.ptr());
  Thyra::set_ele(1,1.0e2,values.ptr());
  Teuchos::RCP<NOX::Abstract::Vector > right_scaling = Teuchos::rcp<NOX::Abstract::Vector>(new NOX::Thyra::Vector(values));
  
  // Create Status Tests
  Teuchos::RCP<NOX::StatusTest::Generic> status_tests;
  {
    Teuchos::RCP<Teuchos::ParameterList> st = 
      Teuchos::parameterList("Status Tests");
    st->set("Test Type", "Combo");
    st->set("Combo Type", "OR");
    st->set("Number of Tests", 3);
    
    {
      Teuchos::ParameterList& normWRMS = st->sublist("Test 0");
      normWRMS.set("Test Type", "NormWRMS");
      normWRMS.set("Absolute Tolerance", 2.0e-8);
      normWRMS.set("Relative Tolerance", 1.0e-5);
      normWRMS.set("Tolerance", 1.0);
      normWRMS.set("BDF Multiplier", 1.0);
      normWRMS.set("Alpha", 1.0);
      normWRMS.set("Beta", 0.5);
      normWRMS.set("Disable Implicit Weighting", true);
    }
    
    {
      Teuchos::ParameterList& fv = st->sublist("Test 1");
      fv.set("Test Type", "FiniteValue");
      fv.set("Vector Type", "F Vector");
      fv.set("Norm Type", "Two Norm");
    }

    {
      Teuchos::ParameterList& maxiters = st->sublist("Test 2");
      maxiters.set("Test Type", "MaxIters");
      maxiters.set("Maximum Iterations", 20);
    }
    // Create a print class for controlling output below
    nl_params->sublist("Printing").set("MyPID", Comm.MyPID());
    nl_params->sublist("Printing").set("Output Precision", 3);
    nl_params->sublist("Printing").set("Output Processor", 0);
    NOX::Utils printing(nl_params->sublist("Printing"));

    NOX::StatusTest::Factory st_factory;
    status_tests = st_factory.buildStatusTests(*st,printing);
  }
  
  Teuchos::RCP< ::Thyra::VectorBase<double> >
    initial_guess = thyraModel->getNominalValues().get_x()->clone_v();

  Teuchos::RCP<NOX::Thyra::Vector> noxThyraRightScaleVec = 
    Teuchos::rcp_dynamic_cast<NOX::Thyra::Vector>(right_scaling);

  Teuchos::RCP<NOX::Thyra::Group> nox_group =
    Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, 
                                       thyraModel, 
                                       thyraModel->create_W_op(), 
                                       lowsFactory, 
                                       thyraModel->create_W_prec(), 
                                       lowsFactory->getPreconditionerFactory(),
                                       Teuchos::null,
                                       noxThyraRightScaleVec->getThyraRCPVector()
                                       ));

  Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, status_tests, nl_params);
  NOX::StatusTest::StatusType solveStatus = solver->solve();
  
  TEST_ASSERT(solveStatus == NOX::StatusTest::Converged);
  
  // Test total num iterations
  {
    // Problem converges in >20 nonlinear iterations with NO scaling
    // Problem converges in 1 nonlinear iterations with scaling
    TEST_EQUALITY(solver->getNumIterations(), 1);
  }
}

