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

template<typename ScalarType>
int run(int argc, char *argv[])
{
  using Teuchos::Comm;
  using Teuchos::RCP;

  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  
  using SCT = typename Teuchos::ScalarTraits<ST>;
  using MT = typename SCT::magnitudeType;
  
  using OP = typename Tpetra::Operator<ST>;
  using MV = typename Tpetra::MultiVector<ST>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;

  using tcrsmatrix_t = Tpetra::CrsMatrix<ST>;

  Teuchos::GlobalMPISession mpisess(&argc,&argv,&std::cout);

  const ST one  = SCT::one();

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);

  // Get test parameters from command-line processor
  bool verbose = false, debug = false, proc_verbose = false;
  int frequency = -1;        // frequency of status test output.
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  int maxsubspace = 250;      // maximum number of blocks the solver can use for the subspace
  int recycle = 50;      // maximum number of blocks the solver can use for the subspace
  int maxrestarts = 15;      // number of restarts allowed
  std::string filename("orsirr1.hb");
  std::string ortho("IMGS");
  MT tol = 1.0e-10;           // relative residual tolerance

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information from the solver.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
  cmdp.setOption("max-subspace",&maxsubspace,"Maximum number of vectors in search space (including recycle space).");
  cmdp.setOption("recycle",&recycle,"Number of vectors in recycle space.");
  cmdp.setOption("max-cycles",&maxrestarts,"Maximum number of cycles allowed for GCRO-DR solver.");
  cmdp.setOption("ortho-type",&ortho,"Orthogonalization type. Must be one of DGKS, ICGS, IMGS.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose

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

  // Other information used by block solver (can be user specified)
  const int NumGlobalElements = B->getGlobalLength();
  if (maxiters == -1)
    maxiters = NumGlobalElements - 1; // maximum number of iterations to run
  //
  Teuchos::ParameterList belosList;
  belosList.set( "Num Blocks", maxsubspace);             // Maximum number of blocks in Krylov factorization
  belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
  belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  belosList.set( "Num Recycled Blocks", recycle );       // Number of vectors in recycle space
  belosList.set( "Orthogonalization", ortho );           // Orthogonalization type

  int verbosity = Belos::Errors + Belos::Warnings;
  if (verbose) {
    verbosity += Belos::TimingDetails + Belos::StatusTestDetails;
    if (frequency > 0)
      belosList.set( "Output Frequency", frequency );
  }
  if (debug) {
    verbosity += Belos::Debug;
  }
  belosList.set( "Verbosity", verbosity );
  
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
    std::cout << "Max number of GCRODR iterations: " << maxiters << std::endl;
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
  
  // Default return value
  if (proc_verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
  //
} // end test_gcrodr_hb.cpp

int main(int argc, char *argv[]) {
  return run<double>(argc, argv);
  // return run<float>(argc, argv);
}

