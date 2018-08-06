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
#include "BelosPseudoBlockCGSolMgr.hpp"

// I/O for Harwell-Boeing files
#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

int main (int argc, char *argv[])
{
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::tuple;
  using std::cout;
  using std::endl;

  typedef double                           ST;
  typedef Teuchos::ScalarTraits<ST>       SCT;
  typedef SCT::magnitudeType               MT;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Tpetra::Operator<ST>             OP;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;

  Teuchos::GlobalMPISession mpisess (&argc, &argv, &cout);

  bool success = false;
  bool verbose = false;
  try {
    const ST one  = SCT::one();

    int MyPID = 0;

    typedef Tpetra::Map<>::node_type Node;
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    RCP<Node>             node; // only for type deduction; null ok

    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool debug = false;
    int frequency = -1;  // how often residuals are printed by solver
    int numrhs = 1;      // total number of right-hand sides to solve for
    int blocksize = 1;   // blocksize used by solver
    int maxiters = -1;   // maximum number of iterations for solver to use
    bool condest = true; // compute the condition estimate
    std::string filename("bcsstk14.hb");
    MT tol = 1.0e-5;     // relative residual tolerance

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by CG solver.");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
    cmdp.setOption("block-size",&blocksize,"Block size to be used by the CG solver.");
    cmdp.setOption("condest","nocond",&condest,"Compute a condition number estimate");

    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (debug) {
      verbose = true;
    }
    if (!verbose) {
      frequency = -1;  // reset frequency if test is not verbose
    }

    MyPID = comm->getRank ();
    proc_verbose = ( verbose && (MyPID==0) );

    //
    // Get the data from the HB file and build the Map,Matrix
    //
    RCP<Tpetra::CrsMatrix<ST> > A;
    Tpetra::Utils::readHBMatrix(filename,comm,node,A);
    RCP<const Tpetra::Map<> > map = A->getDomainMap();

    // Create initial vectors
    RCP<MV> B, X;
    X = rcp (new MV (map, numrhs));
    MVT::MvRandom (*X);
    B = rcp (new MV (map, numrhs));
    OPT::Apply (*A, *X, *B);
    MVT::MvInit (*X, 0.0);

    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1) {
      maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
    }
    //
    RCP<ParameterList> belosList = parameterList ("Belos");
    belosList->set ("Block Size", blocksize); // Blocksize to be used by iterative solver
    belosList->set ("Maximum Iterations", maxiters); // Maximum number of iterations allowed
    belosList->set ("Convergence Tolerance", tol); // Relative convergence tolerance requested
    int verbLevel = Belos::Errors + Belos::Warnings;
    if (debug) {
      verbLevel += Belos::Debug;
    }
    if (verbose) {
      verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
    }
    belosList->set ("Verbosity", verbLevel);
    if (verbose) {
      if (frequency > 0) {
        belosList->set ("Output Frequency", frequency);
      }
    }
    if (condest) {
      belosList->set("Estimate Condition Number", true);
    }

    //
    // Construct an unpreconditioned linear problem instance.
    //
    typedef Belos::LinearProblem<ST, MV, OP> problem_type;
    RCP<problem_type> problem = rcp (new problem_type (A, X, B));
    bool set = problem->setProblem();
    if (! set) {
      if (proc_verbose)
        cout << endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << endl;
      return -1;
    }
    //
    // *******************************************************************
    // *************Start the block CG iteration***********************
    // *******************************************************************
    //
    Belos::PseudoBlockCGSolMgr<ST,MV,OP> solver;
    solver.setParameters (belosList);
    solver.setProblem (problem);

    if (proc_verbose) {
      cout << "Input parameter list: " << endl;
      belosList->print (cout);
      cout << endl;

      cout << "Belos parameter list: " << endl;
      solver.getCurrentParameters ()->print (cout);
      cout << endl;
    }

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

    {
      std::ostringstream os;
      os << "Proc " << comm->getRank () << ": Finished solve" << endl;
      std::cerr << os.str ();
    }

    if (proc_verbose) {
      cout << "Input parameter list after solve: " << endl;
      belosList->print (cout);
      cout << endl;

      cout << "Belos parameter list after solve: " << endl;
      solver.getCurrentParameters ()->print (cout);
      cout << endl;
    }

    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<MT> actual_resids (numrhs);
    std::vector<MT> rhs_norm (numrhs);
    MV resid (map, numrhs);
    OPT::Apply (*A, *X, resid);
    MVT::MvAddMv (-one, resid, one, *B, resid);
    MVT::MvNorm (resid, actual_resids);
    MVT::MvNorm (*B, rhs_norm);
    if (proc_verbose) {
      cout << "---------- Actual Residuals (normalized) ----------"
           << endl << endl;
    }
    for (int i = 0; i < numrhs; ++i) {
      MT actRes = actual_resids[i]/rhs_norm[i];
      if (proc_verbose) {
        cout << "Problem " << i << " : \t" << actRes << endl;
      }
      if (actRes > tol) {
        badRes = true;
      }
    }

    success = (ret==Belos::Converged && ! badRes);

    // Check condition estimate
    const ST cond = solver.getConditionEstimate ();
    if (condest && cond < SCT::one ()) {
      success = false;
    }
    if (proc_verbose) {
      cout << "Condition Estimate = "<<cond <<endl;
    }

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

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_bl_cg_hb.cpp
