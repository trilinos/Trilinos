//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockStochasticCGSolMgr.hpp"

// I/O for Harwell-Boeing files
#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

int
main (int argc, char *argv[])
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using std::endl;
  using std::cout;

  typedef double                           ST;
  typedef Teuchos::ScalarTraits<ST>       STS;
  typedef STS::magnitudeType               MT;
  typedef Tpetra::Operator<ST>             OP;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  typedef MV::node_type             node_type;

  typedef Teuchos::CommandLineProcessor   CLP;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  bool success = false;
  bool verbose = false;
  try {
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
    RCP<node_type> node; // only for type deduction; null ok

    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool debug = false;
    int frequency = -1;
    int numrhs = 1;
    int blocksize = 1;
    int maxiters = -1;
    std::string filename ("bcsstk14.hb");
    MT tol = 1.0e-5;

    CLP cmdp (false, true);
    cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and "
                    "results.");
    cmdp.setOption ("debug", "nodebug", &debug, "Run debugging checks.");
    cmdp.setOption ("frequency", &frequency, "Solver's frequency for printing "
                    "residuals (#iters).  -1 means \"never\"; 1 means \"every "
                    "iteration.\"");
    cmdp.setOption ("tol", &tol, "Relative residual tolerance used by solver.");
    cmdp.setOption ("filename", &filename, "Filename for Harwell-Boeing test "
                    "matrix.");
    cmdp.setOption ("num-rhs", &numrhs, "Number of right-hand sides for which "
                    "to solve.");
    cmdp.setOption ("max-iters", &maxiters, "Maximum number of iterations per "
                    "linear system (-1 := adapted to problem/block size).");
    cmdp.setOption ("block-size", &blocksize, "Block size to be used by the "
                    "solver.");
    if (cmdp.parse (argc, argv) != CLP::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (debug) {
      verbose = true;
    }
    if (!verbose) {
      frequency = -1;  // reset frequency if test is not verbose
    }

    const int MyPID = comm->getRank ();
    proc_verbose = verbose && MyPID == 0;

    if (proc_verbose) {
      cout << Belos::Belos_Version () << endl << endl;
    }

    //
    // Get the data from the HB file and build the Map,Matrix
    //
    RCP<Tpetra::CrsMatrix<ST> > A;
    Tpetra::Utils::readHBMatrix (filename, comm, node, A);
    RCP<const Tpetra::Map<> > map = A->getDomainMap ();

    // Create initial vectors
    RCP<MV> B, X;
    X = rcp (new MV (map, numrhs));
    MVT::MvRandom (*X);
    B = rcp (new MV (map, numrhs));
    OPT::Apply (*A, *X, *B);
    MVT::MvInit (*X, STS::zero ());

    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1) {
      maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
    }
    //
    ParameterList belosList;
    belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    int verbLevel = Belos::Errors + Belos::Warnings;
    if (debug) {
      verbLevel += Belos::Debug;
    }
    if (verbose) {
      verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
    }
    belosList.set( "Verbosity", verbLevel );
    if (verbose) {
      if (frequency > 0) {
        belosList.set( "Output Frequency", frequency );
      }
    }
    //
    // Construct an unpreconditioned linear problem instance.
    //
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    bool set = problem.setProblem ();
    if (! set) {
      if (proc_verbose) {
        cout << endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << endl;
      }
      return EXIT_FAILURE;
    }
    //
    // *******************************************************************
    // *************Start the block CG iteration***********************
    // *******************************************************************
    //
    typedef Belos::PseudoBlockStochasticCGSolMgr<ST,MV,OP> solver_type;
    solver_type solver (rcpFromRef (problem), rcpFromRef (belosList));

    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      cout << endl << endl;
      cout << "Dimension of matrix: " << NumGlobalElements << endl;
      cout << "Number of right-hand sides: " << numrhs << endl;
      cout << "Block size used by solver: " << blocksize << endl;
      cout << "Max number of CG iterations: " << maxiters << endl;
      cout << "Relative residual tolerance: " << tol << endl;
      cout << endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = solver.solve();
    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<MT> actual_resids (numrhs);
    std::vector<MT> rhs_norm (numrhs);
    MV resid (map, numrhs);
    OPT::Apply (*A, *X, resid);
    MVT::MvAddMv (-STS::one (), resid, STS::one (), *B, resid);
    MVT::MvNorm (resid, actual_resids);
    MVT::MvNorm (*B, rhs_norm);
    if (proc_verbose) {
      cout << "---------- Actual Residuals (normalized) ----------" << endl << endl;
    }
    for (int i=0; i<numrhs; i++) {
      MT actRes = actual_resids[i]/rhs_norm[i];
      if (proc_verbose) {
        cout << "Problem " << i << " : \t" << actRes << endl;
      }
      if (actRes > tol) badRes = true;
    }

    success = ret == Belos::Converged && ! badRes;
    if (success) {
      if (proc_verbose) {
        cout << "\nEnd Result: TEST PASSED" << endl;
      }
    } else {
      if (proc_verbose) {
        cout << "\nEnd Result: TEST FAILED" << endl;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
