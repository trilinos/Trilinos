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
#include "Teuchos_StandardCatchMacros.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "ModelEvaluatorHeq.hpp"

#include "NOX_Thyra_MatrixFreeJacobianOperator.hpp"
#include "NOX_MatrixFree_ModelEvaluatorDecorator.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  bool success = false;
  bool verbose = false;
  try {
    // Parse the command line
    using Teuchos::CommandLineProcessor;
    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
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

    const int num_elements = 400;

    // Check we have at least one unknown per processor
    if (Comm.NumProc() > num_elements) {
      std::cerr << "Error! Number of elements must be greate than number of processors!"
        << std::endl;
      return EXIT_FAILURE;
    }

    // Create the model evaluator object
    double paramC = 0.99;
    Teuchos::RCP<ModelEvaluatorHeq<double> > model =
      modelEvaluatorHeq<double>(Teuchos::rcp(&Comm,false),num_elements,paramC);

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
      Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model));
    //Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, model, model->create_W_op(), lowsFactory, Teuchos::null, Teuchos::null));

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
    nl_params->set("Nonlinear Solver", "Anderson Accelerated Fixed-Point");
    nl_params->sublist("Anderson Parameters").set("Storage Depth", 5);
    nl_params->sublist("Anderson Parameters").set("Mixing Parameter", 1.0);
    nl_params->sublist("Anderson Parameters").set("Acceleration Start Iteration", 5);
    nl_params->sublist("Anderson Parameters").sublist("Preconditioning").set("Precondition", false);
    nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);
    nl_params->sublist("Printing").sublist("Output Information").set("Details",true);
    nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration",true);
    //nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration StatusTest",true);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(nox_group, combo, nl_params);
    NOX::StatusTest::StatusType solvStatus = solver->solve();

    // 1. Convergence
    int status = 0;
    if (solvStatus != NOX::StatusTest::Converged)
      status = 1;
    // 2. Number of iterations
    if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Output").get("Nonlinear Iterations", 0) != 14)
      status = 2;

    success = status==0;

    if (Comm.MyPID() == 0) {
      if (success)
        std::cout << "Test passed!" << std::endl;
      else
        std::cout << "Test failed!" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

