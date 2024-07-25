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
#include "NOX_SolverStats.hpp"
#include "NOX_Observer_Print.hpp"
#include "NOX_Utils.H"
#include "Teuchos_TestingHelpers.hpp"

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
#include "Teuchos_StandardCatchMacros.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "ModelEvaluator2DSim.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  bool success = false;
  bool verbose = false;
  try {
    // Parse the command line
    using Teuchos::CommandLineProcessor;
    CommandLineProcessor clp;
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

    int status = 0;

    // Check we have only one processor since this problem doesn't work
    // for more than one proc
    if (Comm.NumProc() > 1)
      throw "Error!  Problem can only be run with at most 1 processor!";

    // Create the model evaluator object
    double d = 10.0;
    double p0 = 2.0;
    double p1 = 0.0;
    double x00 = 0.0;
    double x01 = 1.0;
    Teuchos::RCP<ModelEvaluator2DSim<double> > thyraModel =
      Teuchos::rcp(new ModelEvaluator2DSim<double>(Teuchos::rcp(&Comm,false),
            d,p0,p1,x00,x01));

    // Create the linear solver type with Stratimikos
    //Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
    //lowsFactory = rcp(new Thyra::AmesosLinearOpWithSolveFactory());

    ::Stratimikos::DefaultLinearSolverBuilder builder;

    Teuchos::RCP<Teuchos::ParameterList> p =
      Teuchos::rcp(new Teuchos::ParameterList);
    p->set("Linear Solver Type", "AztecOO");
    p->set("Preconditioner Type", "Ifpack");
    //p->set("Enable Delayed Solver Construction", true);
    builder.setParameterList(p);

    Teuchos::RCP< ::Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = builder.createLinearSolveStrategy("");

    thyraModel->set_W_factory(lowsFactory);

    // Create the initial guess
    Teuchos::RCP< ::Thyra::VectorBase<double> >
      initial_guess = thyraModel->getNominalValues().get_x()->clone_v();

    // Create the NOX::Thyra::Group
    Teuchos::RCP<NOX::Thyra::Group> nox_group =
      Teuchos::rcp(new NOX::Thyra::Group(*initial_guess, thyraModel));

    //   nox_group->computeF();
    //   std::cout << "ComputedF!" << std::endl;
    //   const NOX::Thyra::Vector& t_vec =
    //     dynamic_cast<const NOX::Thyra::Vector&>(nox_group->getF());

    //   t_vec.print(std::cout);
    //   exit(0);

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

    // Register Print Observer (and disable normal printing)
    nl_params->sublist("Printing").set("Output Information",NOX::Utils::Error);
    Teuchos::RCP<NOX::Utils> os = Teuchos::rcp(new NOX::Utils(nl_params->sublist("Printing")));
    Teuchos::RCP<NOX::Observer> printObserver = Teuchos::rcp(new NOX::ObserverPrint(os));
    nl_params->sublist("Solver Options").set("Observer",printObserver);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(nox_group, combo, nl_params);
    NOX::StatusTest::StatusType solvStatus = solver->solve();

    success = status==0;

    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->numNonlinearSolves,1,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->numNonlinearIterations,7,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->numTotalNonlinearIterations,7,std::cout,success);

    TEUCHOS_TEST_ASSERT(solver->getSolverStatistics()->linearSolve.lastLinearSolve_Converged,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.lastLinearSolve_NumIterations,1,std::cout,success);
    TEUCHOS_TEST_ASSERT(solver->getSolverStatistics()->linearSolve.lastLinearSolve_AchievedTolerance < 1.0e-4,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.lastNonlinearSolve_NumLinearSolves,7,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.lastNonlinearSolve_NumLinearIterations,7,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.allNonlinearSolves_NumLinearSolves,7,std::cout,success);
    TEUCHOS_TEST_EQUALITY(solver->getSolverStatistics()->linearSolve.allNonlinearSolves_NumLinearIterations,7,std::cout,success);

    if (solvStatus == NOX::StatusTest::Converged)
      std::cout << "Test passed!" << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
