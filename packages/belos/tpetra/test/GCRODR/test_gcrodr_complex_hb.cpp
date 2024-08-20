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
//

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosTpetraTestFramework.hpp"

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>

template <typename ScalarType>
int run(int argc, char *argv[])
{ 
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using Teuchos::tuple;  
  
  using BST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using ST = typename std::complex<BST>;

  using SCT = typename Teuchos::ScalarTraits<ST>;
  using MT = typename SCT::magnitudeType;

  using OP = typename Tpetra::Operator<ST>;
  using MV = typename Tpetra::MultiVector<ST>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;

  using tcrsmatrix_t = Tpetra::CrsMatrix<ST>;

  Teuchos::GlobalMPISession mpisess(&argc, &argv, &std::cout);

  const ST one  = SCT::one();

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);

  // Get test parameters from command-line processor
  bool verbose = false, proc_verbose = false, debug = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 1;      // total number of right-hand sides to solve for
  std::string filename("mhd1280b.cua");
  MT tol = 1.0e-5;     // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GCRODR solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
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
  RCP<const Tpetra::Map<> > map = A->getDomainMap();

  // Create initial vectors
  RCP<MV> B, X;
  X = rcp( new MV(map,numrhs) );
  MVT::MvRandom( *X );
  B = rcp( new MV(map,numrhs) );
  OPT::Apply( *A, *X, *B );
  MVT::MvInit( *X, 0.0 );

  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  const int NumGlobalElements = B->getGlobalLength();
  int numIters1, numIters2, numIters3;
  int maxits = NumGlobalElements; // maximum number of iterations to run
  int numBlocks = 100;
  int numRecycledBlocks = 20;
  Teuchos::ParameterList belosList;
  belosList.set( "Maximum Iterations", maxits );         // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  belosList.set( "Num Blocks", numBlocks );
  belosList.set( "Num Recycled Blocks", numRecycledBlocks );

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

  // Start the GCRODR iteration
  Belos::GCRODRSolMgr<ST,MV,OP> solver( rcpFromRef(problem), rcpFromRef(belosList) );

  // Print out information about problem
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Max number of GCRODR iterations: " << maxits << std::endl;
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
  }
  
  // Perform solve
  Belos::ReturnType ret = solver.solve();
  numIters1=solver.getNumIters();
  
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
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  }
  for ( int i=0; i<numrhs; i++) {
    MT actRes = actual_resids[i]/rhs_norm[i];
    if (proc_verbose) {
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
    }
    if (actRes > tol) badRes = true;
  }

  if (proc_verbose) { std::cout << "First solve took " << numIters1 << " iterations." << std::endl; }

  // Resolve linear system with same rhs and recycled space
  MVT::MvInit( *X, SCT::zero() );
  solver.reset(Belos::Problem);
  ret = solver.solve();
  numIters2=solver.getNumIters();

  if (proc_verbose) { std::cout << "Second solve took " << numIters2 << " iterations." << std::endl; }

  // Resolve linear system (again) with same rhs and recycled space
  MVT::MvInit( *X, SCT::zero() );
  solver.reset(Belos::Problem);
  ret = solver.solve();
  numIters3=solver.getNumIters();

  if (proc_verbose) { std::cout << "Third solve took " << numIters3 << " iterations." << std::endl; }

  if ( ret!=Belos::Converged || badRes || numIters1 < numIters2 || numIters2 < numIters3 ) {
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

} // end test_gcrodr_complex_hb.cpp

int main(int argc, char *argv[]) {
  return run<double>(argc, argv);
  // return run<float>(argc, argv);
}

