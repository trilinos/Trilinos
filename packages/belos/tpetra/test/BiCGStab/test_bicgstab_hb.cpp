// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a problem from a file, which can only be a Harwell-Boeing (*.hb)
// matrix.  The right-hand side is generated using a random solution vector that has
// the matrix applied to it.  The initial guesses are all set to zero.  
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBiCGStabSolMgr.hpp"
#include "BelosTpetraTestFramework.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include "Teuchos_StandardCatchMacros.hpp"
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

int main (int argc, char *argv[])
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;

  typedef Tpetra::MultiVector<>::scalar_type ST;
  typedef Teuchos::ScalarTraits<ST>       SCT;
  typedef SCT::magnitudeType               MT;
  typedef Tpetra::Operator<ST>             OP;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;

  Teuchos::GlobalMPISession mpisess(&argc,&argv,&cout);

  const ST one  = SCT::one();

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);

  bool verbose = false;
  bool success = true;

  // This "try" relates to TEUCHOS_STANDARD_CATCH_STATEMENTS near the
  // bottom of main().  That macro has the corresponding "catch".
  try {
    bool proc_verbose = false;
    int frequency = -1;        // frequency of status test output.
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxiters = -1;         // maximum number of iterations allowed per linear system
    std::string filename ("cage4.hb");
    MT tol = 1.0e-5;           // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption ("verbose", "quiet", &verbose, "Whether to print messages "
                    "and results.");
    cmdp.setOption ("frequency", &frequency, "Frequency of solver output "
                    "(1 means every iteration; -1 means never).");
    cmdp.setOption ("filename", &filename, "Test matrix filename.  "
                    "Allowed file extensions: *.hb, *.mtx, *.triU, *.triS");
    cmdp.setOption ("tol", &tol, "Relative residual tolerance for solver.");
    cmdp.setOption ("num-rhs", &numrhs, "Number of right-hand sides to solve.");
    cmdp.setOption ("max-iters", &maxiters, "Maximum number of iterations per "
                    "linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse (argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (! verbose) {
      frequency = -1;  // reset frequency if test is not verbose
    }

    proc_verbose = ( verbose && (MyPID==0) );

    //
    // *************Get the problem*********************
    //
    Belos::Tpetra::HarwellBoeingReader<Tpetra::CrsMatrix<ST> > reader( comm );
    RCP<Tpetra::CrsMatrix<ST> > A = reader.readFromFile( filename );
    RCP<const Tpetra::Map<> > map = A->getDomainMap();

    // Create initial vectors
    RCP<MV> B, X;
    X = rcp( new MV(map,numrhs) );
    MVT::MvRandom( *X );
    B = rcp( new MV(map,numrhs) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );

    //
    // Create parameter list for the Belos solver
    //
    const int NumGlobalElements = B->getGlobalLength ();
    if (maxiters == -1) {
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    }
    RCP<ParameterList> belosList (new ParameterList ("Belos"));
    belosList->set ("Maximum Iterations", maxiters);
    belosList->set ("Convergence Tolerance", tol);
    if (verbose) {
      belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
                     Belos::TimingDetails + Belos::StatusTestDetails);
      if (frequency > 0) {
        belosList->set ("Output Frequency", frequency);
      }
    }
    else {
      belosList->set ("Verbosity", Belos::Errors + Belos::Warnings +
                      Belos::FinalSummary);
    }
    //
    // Construct a preconditioned linear problem
    //
    RCP<Belos::LinearProblem<ST,MV,OP> > problem
      = rcp (new Belos::LinearProblem<ST,MV,OP> (A, X, B));
    bool set = problem->setProblem ();
    if (! set) {
      if (proc_verbose) {
        cout << endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << endl;
      }
      return -1;
    }

    // Create a Belos solver.
    RCP<Belos::SolverManager<ST,MV,OP> > solver
      = rcp (new Belos::BiCGStabSolMgr<ST,MV,OP> (problem, belosList));

    if (proc_verbose) {
      cout << endl << endl;
      cout << "Dimension of matrix: " << NumGlobalElements << endl;
      cout << "Number of right-hand sides: " << numrhs << endl;
      cout << "Max number of BiCGStab iterations: " << maxiters << endl;
      cout << "Relative residual tolerance: " << tol << endl;
      cout << endl;
    }

    // Ask Belos to solve the linear system.
    Belos::ReturnType ret = solver->solve();

    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<MT> actual_resids( numrhs );
    std::vector<MT> rhs_norm( numrhs );
    MV resid(map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -one, resid, one, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    }
    for ( int i=0; i<numrhs; i++) {
      MT actRes = actual_resids[i]/rhs_norm[i];
      if (proc_verbose) {
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      }
      if (actRes > tol) badRes = true;
    }

    if ( ret!=Belos::Converged || badRes) {
      if (proc_verbose) {
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
      }
      return -1;
    }
    //
    // Default return value
    //
    if (proc_verbose) {
      std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
