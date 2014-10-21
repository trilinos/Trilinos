//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
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

  // Create a communicator for Epetra objects
#ifdef HAVE_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  bool success = false;
  bool verbose = false;
  int StorageDepth = 10;
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
      std::cout << "Verbosity Activated" << std::endl;
    else
      std::cout << "Verbosity Disabled" << std::endl;

    int status = 0;

    const int num_elements = 400;

    // Check we have only one processor since this problem doesn't work
    // for more than one proc
    if (Comm.NumProc() > num_elements)
      throw "Error! Number of elements must be greate than number of processors!";

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
    //nl_params->sublist("Printing").sublist("Output Information").set("Outer Iteration StatusTest",true);

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

    // 1. Convergence
    if (solvStatus != NOX::StatusTest::Converged)
      status = 1;
    // 2. Number of iterations
    if (const_cast<Teuchos::ParameterList&>(solver->getList()).sublist("Output").get("Nonlinear Iterations", 0) != 19)
      status = 2;

    success = status==0;

    if (success)
      std::cout << "Test passed!" << std::endl;
    else
      std::cout << "Test failed!" << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
