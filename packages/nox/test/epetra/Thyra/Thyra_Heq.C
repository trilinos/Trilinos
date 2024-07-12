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
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_TestingHelpers.hpp"

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
  const Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  bool success = false;
  bool verbose = false;
  int StorageDepth = 2;
  double ParamC = 0.999;
  std::string LineSearch = "Full Step";
  int Reorthogonalize = 0;
  try {
    // Parse the command line
    using Teuchos::CommandLineProcessor;
    CommandLineProcessor  clp;
    clp.throwExceptions(false);
    clp.addOutputSetupOptions(true);
    clp.setOption( "v", "disable-verbosity", &verbose, "Enable verbosity" );
    clp.setOption("storage-depth", &StorageDepth, "Anderson storage depth parameter");
    clp.setOption("c-parameter", &ParamC, "H equation parameter, c");
    clp.setOption("line-search", &LineSearch, "Line search type: Full Step, Backtracking, or Polynomial");
    clp.setOption("ortho-frequency", &Reorthogonalize, "Reorthogonalization frequency");

    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv,&std::cerr);

    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL )
      return parse_return;

    if (verbose)
      *out << "Verbosity Activated" << std::endl;
    else
      *out << "Verbosity Disabled" << std::endl;

    const int num_elements = 400;

    // Check we have at least one unknown per processor
    if (Comm.NumProc() > num_elements)
      throw "Error! Number of elements must be greater than number of processors!";

    // Create the model evaluator object
    Teuchos::RCP<ModelEvaluatorHeq<double> > model =
      modelEvaluatorHeq<double>(Teuchos::rcp(&Comm,false),num_elements,ParamC);

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
      Teuchos::rcp(new NOX::StatusTest::MaxIters(100));
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
    nl_params->sublist("Anderson Parameters").set("Storage Depth", StorageDepth);
    nl_params->sublist("Anderson Parameters").set("Mixing Parameter", 1.0);
    nl_params->sublist("Anderson Parameters").set("Acceleration Start Iteration", 1);
    nl_params->sublist("Anderson Parameters").set("Reorthogonalization Frequency", Reorthogonalize);
    nl_params->sublist("Anderson Parameters").sublist("Preconditioning").set("Precondition", false);
    nl_params->sublist("Direction").sublist("Newton").sublist("Linear Solver").set("Tolerance", 1.0e-4);
    nl_params->sublist("Printing").sublist("Output Information").set("Details",true);
    nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration",true);

    if (verbose){
      nl_params->sublist("Printing").sublist("Output Information").set("Error",true);
      nl_params->sublist("Printing").sublist("Output Information").set("Test Details",true);
      nl_params->sublist("Printing").sublist("Output Information").set("Debug",true);
      nl_params->sublist("Printing").sublist("Output Information").set("Warning",true);
      nl_params->sublist("Printing").sublist("Output Information").set("Parameters",true);
      nl_params->sublist("Printing").sublist("Output Information").set("Linear Solver Details",true);
      nl_params->sublist("Printing").sublist("Output Information").set("Inner Iteration",true);
      nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration StatusTest",true);
    }

    // Line search parameters
    nl_params->sublist("Line Search").set("Method", LineSearch);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(nox_group, combo, nl_params);
    NOX::StatusTest::StatusType solvStatus = solver->solve();
    Teuchos::RCP<NOX::Abstract::Vector> diff = solver->getSolutionGroup().getX().clone();

    // Test the reset function by resolving the same problem
    Teuchos::RCP<NOX::Abstract::Vector> newGuess = diff->clone(NOX::ShapeCopy);
    newGuess->init(1.0);
    solver->reset(*newGuess);
    solver->solve();
    diff->update(1.0, solver->getSolutionGroup().getX(), -1.0);

    // Create a print class for controlling output below
    nl_params->sublist("Printing").set("MyPID", Comm.MyPID());
    nl_params->sublist("Printing").set("Output Precision", 3);
    nl_params->sublist("Printing").set("Output Processor", 0);
    NOX::Utils printing(nl_params->sublist("Printing"));

    // Output the parameter list
    if (verbose) {
      if (printing.isPrintType(NOX::Utils::Parameters)) {
        printing.out() << std::endl << "Final Parameters" << std::endl
          << "****************" << std::endl;
        solver->getList().print(printing.out());
        printing.out() << std::endl;
      }
    }

    *out << "\nCheck for test pass/fail:\n";

    bool loc_success = true;
 
    // 1. Convergence
    TEUCHOS_TEST_EQUALITY_CONST(solvStatus, NOX::StatusTest::Converged, *out, loc_success);
    // 2. Number of iterations
    int numIterations = 0;
    const int *numItersPtr = nullptr;
    if ( nullptr != (numItersPtr = Teuchos::getParameterPtr<int>(
                        solver->getList().sublist("Output"), "Nonlinear Iterations")) ) 
    {
      numIterations = *numItersPtr;
    } 
    TEUCHOS_TEST_EQUALITY_CONST(numIterations, 11, *out, loc_success);
    // 3. Same reset solution
    TEUCHOS_TEST_COMPARE_CONST(diff->norm(), <=, 1.0e-14, *out, loc_success);

    success = loc_success;

    if (success)
      *out << "\nTest passed!" << std::endl;
    else
      *out << "\nTest failed!" << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
