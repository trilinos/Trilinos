// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_ScaledModelEvaluator.hpp"
#include "ModelEvaluator2DSim.hpp"

#include "NOX_PrePostOperator_RowSumScaling.H"
#include "Thyra_VectorBase.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  bool success = true;
  bool verbose = true;
  try {
    // Parse the command line
    using Teuchos::CommandLineProcessor;
    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    bool verbose = false;
    clp.setOption( "v", "disable-verbosity", &verbose, "Enable verbosity" );

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    if (verbose)
      std::cout << "Verbosity Activated" << std::endl;
    else
      std::cout << "Verbosity Disabled" << std::endl;

    // Create a communicator for Epetra objects
#ifdef HAVE_MPI
    Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
    Epetra_SerialComm Comm;
#endif

    // Check we have only one processor since this problem doesn't work
    // for more than one proc
    if (Comm.NumProc() > 1) {
      std::cerr << "Error!  Problem can only be run with at most 1 processor!"
            << std::endl;
      return -1;
    }

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
    nl_params->sublist("Line Search").set("Method", "Polynomial");

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
        normF_rel.set("Test Type", "RelativeNormF");
        normF_rel.set("Tolerance", 1.0e-4);

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

    // Test total num iterations
    {
      // Problem converges in 7 nonlinear iterations with NO scaling
      // Problem converges in 6 nonlinear iterations with RS scaling (bad test problem - too easy)
      Teuchos::RCP< ::Thyra::NOXNonlinearSolver> thyra_nox_solver =
        Teuchos::rcp_dynamic_cast< ::Thyra::NOXNonlinearSolver>(solver);
      TEUCHOS_ASSERT(thyra_nox_solver->getNOXSolver()->getNumIterations() == 6);
    }

    if (solve_status.solveStatus == ::Thyra::SOLVE_STATUS_CONVERGED) {
      std::cout << "Test passed!" << std::endl;
      success = true;
    }
    else {
      success = false;
    }

    std::cout << *p << std::endl;

    Teuchos::TimeMonitor::summarize();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
