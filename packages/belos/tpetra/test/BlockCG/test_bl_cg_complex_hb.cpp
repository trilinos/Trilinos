// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosTpetraTestFramework.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Harwell-Boeing files
#include <Trilinos_Util_iohb.h>

template <typename ScalarType>
int run(int argc, char *argv[])
{
  using BSC = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using ST  = typename std::complex<BSC>;
  
  using OP = typename Tpetra::Operator<ST>;
  using MV = typename Tpetra::MultiVector<ST>;
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST>;
  
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;

  using SCT = typename Teuchos::ScalarTraits<ST>;
  using MT  = typename SCT::magnitudeType;

  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::MultiVector;

  using Teuchos::GlobalMPISession;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;

  const ST one = SCT::one();

  GlobalMPISession mpisess(&argc,&argv,&std::cout);
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);

  // Get test parameters from command-line processor
  bool verbose = false, proc_verbose = false, debug = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 1;      // total number of right-hand sides to solve for
  int blocksize = 1;   // blocksize used by solver
  int maxiters = -1;   // maximum number of iterations for solver to use
  std::string filename("mhd1280b.cua");
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
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }
  if (!verbose) {
    frequency = -1;  // reset frequency if test is not verbose
  }

  proc_verbose = ( verbose && (MyPID==0) );
  if (proc_verbose) {
    std::cout << Belos::Belos_Version() << std::endl << std::endl;
  }

  Belos::Tpetra::HarwellBoeingReader<tcrsmatrix_t> reader( comm );
  RCP<tcrsmatrix_t> A = reader.readFromFile( filename );
  RCP<const Tpetra::Map<> > map = A->getMap();

  // Create initial vectors
  RCP<MV> B, X;
  X = rcp( new MV(map,numrhs) );
  MVT::MvRandom( *X );
  B = rcp( new MV(map,numrhs) );
  OPT::Apply( *A, *X, *B );
  MVT::MvInit( *X, 0.0 );

  // Other information used by block solver (can be user specified)
  const int NumGlobalElements = B->getGlobalLength();
  if (maxiters == -1) {
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  }

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
  
  // Construct an unpreconditioned linear problem instance.
  Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }

  // Start the block CG iteration
  Belos::BlockCGSolMgr<ST,MV,OP> solver( rcpFromRef(problem), rcpFromRef(belosList) );

  // Print out information about problem
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Max number of CG iterations: " << maxiters << std::endl;
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
  }
  
  // Perform solve
  Belos::ReturnType ret = solver.solve();

  // Compute actual residuals.
  bool badRes = false;
  std::vector<MT> actual_resids( numrhs );
  std::vector<MT> rhs_norm( numrhs );
  MV resid(map, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -one, resid, one, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  if (proc_verbose) {
    std::cout << "---------- Actual Residuals (normalized) ----------" << std::endl<<std::endl;
  }
  for ( int i=0; i<numrhs; i++) {
    MT actRes = actual_resids[i]/rhs_norm[i];
    if (proc_verbose) {
      std::cout << "Problem "<<i<<" : \t"<< actRes << std::endl;
    }
    if (actRes > tol) badRes = true;
  }

  if (ret!=Belos::Converged || badRes) {
    if (proc_verbose) {
      std::cout << "\nEnd Result: TEST FAILED" << std::endl;
    }
    return -1;
  }
  
  // Default return value
  if (proc_verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }

  return 0;
} // end test_bl_cg_complex_hb.cpp

int main(int argc, char *argv[]) {
  return run<double>(argc,argv);

  // wrapped with a check: CMake option Trilinos_ENABLE_FLOAT=ON
  // return run<float>(argc,argv);
}

