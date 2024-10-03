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
#include <Teuchos_TimeMonitor.hpp>

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
#include "Teuchos_Assert.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_ScaledModelEvaluator.hpp"
#include "ModelEvaluatorRosenbrock.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"

#include "NOX_PrePostOperator_RowSumScaling.H"
#include "Thyra_VectorBase.hpp"

using namespace std;

TEUCHOS_UNIT_TEST(dimension, default)
{
  int status = 0;

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  TEST_ASSERT(Comm.NumProc() == 1);

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

  Teuchos::RCP<RosenbrockModelEvaluator> thyraModel =
    Teuchos::rcp(new RosenbrockModelEvaluator(Teuchos::rcp(&Comm,false)));

  thyraModel->set_W_factory(lowsFactory);

  // Create nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> nl_params =
    Teuchos::rcp(new Teuchos::ParameterList);
  //nl_params->set("Nonlinear Solver", "Line Search Based");
  nl_params->set("Nonlinear Solver", "Pseudo-Transient");
  nl_params->sublist("Pseudo-Transient").set("deltaInit",1.0e-1);
  nl_params->sublist("Pseudo-Transient").set<int>("Maximum Number of Pseudo-Transient Iterations",10);
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

  // Create Status Tests
  {
    Teuchos::ParameterList& st = nl_params->sublist("Status Tests");
    st.set("Test Type", "Combo");
    st.set("Combo Type", "OR");
    st.set("Number of Tests", 3);

    {
      Teuchos::ParameterList& conv = st.sublist("Test 0");
      conv.set("Test Type", "Combo");
      conv.set("Combo Type", "AND");
      conv.set("Number of Tests", 2);

      Teuchos::ParameterList& normF_rel = conv.sublist("Test 0");
      normF_rel.set("Test Type", "NormF");
      normF_rel.set("Tolerance", 1.0e-8);

      Teuchos::ParameterList& normWRMS = conv.sublist("Test 1");
      normWRMS.set("Test Type", "NormWRMS");
      normWRMS.set("Absolute Tolerance", 1.0e-8);
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


  Teuchos::RCP< ::Thyra::NOXNonlinearSolver> thyra_nox_solver =
    Teuchos::rcp_dynamic_cast< ::Thyra::NOXNonlinearSolver>(solver);
  TEST_EQUALITY(thyra_nox_solver->getNOXSolver()->getNumIterations(), 13);

  Teuchos::RCP<const Epetra_Vector> x_analytic = thyraModel->get_analytic_solution();

  Teuchos::RCP<const NOX::Abstract::Vector> x = thyra_nox_solver->getNOXSolver()->getSolutionGroup().getXPtr();

  Teuchos::RCP<const NOX::Thyra::Vector> nox_thyra_x =
    Teuchos::rcp_dynamic_cast<const NOX::Thyra::Vector>(x,true);

  Teuchos::RCP<const Thyra::SpmdVectorBase<double> > spmd_x =
    Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(nox_thyra_x->getThyraRCPVector(),true);

  Teuchos::ArrayRCP<const double> local_values;
  spmd_x->getLocalData(outArg(local_values));

  double tol = 1.0e-8;
  TEST_FLOATING_EQUALITY((*x_analytic)[0],local_values[0],tol);
  TEST_FLOATING_EQUALITY((*x_analytic)[1],local_values[1],tol);

  if (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED)
    std::cout << "Test passed!" << std::endl;

//   std::cout << *p << std::endl;

  Teuchos::TimeMonitor::summarize();

  // Final return value (0 = successfull, non-zero = failure)
  TEST_ASSERT(status == 0);
}
