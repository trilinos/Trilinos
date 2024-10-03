// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// Update of belos/epetra/example/BlockGmres/BlockGmresEpetraExFile.cpp
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this example.
//

// Tpetra
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraTestFramework.hpp"


template <typename ScalarType>
int run(int argc, char *argv[]) {

  using ST = typename Tpetra::Vector<ScalarType>::scalar_type;
  using LO = typename Tpetra::Vector<>::local_ordinal_type;
  using GO = typename Tpetra::Vector<>::global_ordinal_type;
  using NT = typename Tpetra::Vector<>::node_type;

  using OP  = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV  = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using MT = typename Teuchos::ScalarTraits<ST>::magnitudeType;

  using tmap_t       = Tpetra::Map<LO,GO,NT>;
  using tvector_t    = Tpetra::Vector<ST,LO,GO,NT>;
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST,LO,GO,NT>;

  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);
  const auto comm = Tpetra::getDefaultComm();
  const int myPID = comm->getRank();

  bool verbose = false;
  bool success = true;

  try {
    bool procVerbose = false;
    int frequency = -1;        // frequency of status test output.
    int blockSize = 1;         // blockSize
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxIters = -1;         // maximum number of iterations allowed per linear system
    int maxSubspace = 50;      // maximum number of blocks the solver can use for the subspace
    int maxRestarts = 15;      // number of restarts allowed
    std::string filename("orsirr1.hb");
    MT tol = 1.0e-5;           // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("block-size",&blockSize,"Block size used by GMRES.");
    cmdp.setOption("max-iters",&maxIters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    cmdp.setOption("max-subspace",&maxSubspace,"Maximum number of blocks the solver can use for the subspace.");
    cmdp.setOption("max-restarts",&maxRestarts,"Maximum number of restarts allowed for GMRES solver.");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    procVerbose = ( verbose && (myPID==0) );

    if (procVerbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }

    // Get the problem
    Belos::Tpetra::HarwellBoeingReader<tcrsmatrix_t> reader( comm );
    RCP<tcrsmatrix_t> A = reader.readFromFile( filename );
    RCP<const tmap_t> Map = A->getDomainMap();

    // Create initial vectors
    RCP<MV> B, X;
    X = rcp( new MV(Map,numrhs) );
    MVT::MvRandom( *X );
    B = rcp( new MV(Map,numrhs) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );

    // ********Other information used by block solver***********
    // *****************(can be user specified)******************

    const int numGlobalElements = B->getGlobalLength();
    if (maxIters == -1)
      maxIters = numGlobalElements - 1; // maximum number of iterations to run

    ParameterList belosList;
    belosList.set( "Num Blocks", maxSubspace);             // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", blockSize );              // BlockSize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxIters );       // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxRestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
        Belos::TimingDetails + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

    // Construct an unpreconditioned linear problem instance.
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    bool set = problem.setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // *******************************************************************
    // *************Start the block Gmres iteration*************************
    // *******************************************************************

    // Create an iterative solver manager.
    RCP< Belos::SolverManager<ST,MV,OP> > newSolver
      = rcp( new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));

    // **********Print out information about problem*******************

    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << numGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Block size used by solver: " << blockSize << std::endl;
      std::cout << "Max number of restarts allowed: " << maxRestarts << std::endl;
      std::cout << "Max number of Gmres iterations per restart cycle: " << maxIters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }

    // Perform solve
    Belos::ReturnType ret = newSolver->solve();

    // Compute actual residuals.
    bool badRes = false;
    std::vector<ST> actualResids( numrhs );
    std::vector<ST> rhsNorm( numrhs );
    MV resid(Map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actualResids );
    MVT::MvNorm( *B, rhsNorm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        ST actRes = actualResids[i]/rhsNorm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
    }

    if (ret!=Belos::Converged || badRes) {
      success = false;
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
    } else {
      success = true;
      if (procVerbose)
        std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char *argv[]) {
  // run with different ST
  return run<double>(argc,argv);
  // run<float>(argc,argv); // FAILS
}
