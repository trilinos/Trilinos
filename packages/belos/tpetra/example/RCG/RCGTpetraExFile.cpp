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
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MatrixIO.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosRCGSolMgr.hpp"

// Teuchos
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

template<typename ScalarType>
int run(int argc, char *argv[])
{
  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT = typename Tpetra::MultiVector<>::node_type;

  using OP = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV = typename Tpetra::MultiVector<ST,LO,GO,NT>;

  using tcrsmatrix_t    = Tpetra::CrsMatrix<ST,LO,GO,NT>;
  using tmap_t          = Tpetra::Map<LO,GO,NT>;
  using tmultivector_t  = Tpetra::MultiVector<ST,LO,GO,NT>;

  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // This calls MPI_Init and MPI_Finalize as necessary.
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

  bool verbose = true;
  bool success = true;

  try
  {
    int MyPID = rank(*comm);

    // Get test parameters from command-line processor
    bool proc_verbose = false;
    int frequency = -1;                  // frequency of status test output.
    std::string filename("bcsstk14.hb"); // default input filename
    double tol = 1.0e-6;                 // relative residual tolerance
    int numBlocks = 100;                 // maximum number of blocks the solver can use for the Krylov subspace
    int recycleBlocks = 10;              // maximum number of blocks the solver can use for the recycle space
    int numrhs = 2;                      // number of right-hand sides to solve for
    int maxiters = 4000;                 // maximum number of iterations allowed per linear system

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by the RCG solver.");
    cmdp.setOption("max-subspace",&numBlocks,"Maximum number of vectors in search space (not including recycle space).");
    cmdp.setOption("recycle",&recycleBlocks,"Number of vectors in recycle space.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose
    
    // Get the problem
    RCP<tcrsmatrix_t> A;
    Tpetra::Utils::readHBMatrix(filename, comm, A);
    RCP<const tmap_t> rowMap = A->getDomainMap();

    proc_verbose = ( verbose && (MyPID==0) );

    // Construct initial guess and right-hand sides
    RCP<tmultivector_t> B, X;
    X = rcp(new MV(rowMap, numrhs));
    MVT::MvRandom(*X);
    B = rcp(new MV( rowMap, numrhs));
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );
    
    // Other information used by block solver (can be user specified)
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Num Blocks", numBlocks);               // Maximum number of blocks in Krylov space
    belosList.set( "Num Recycled Blocks", recycleBlocks ); // Number of vectors in recycle space
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    
    // Construct an unpreconditioned linear problem instance.
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    
    // Create an iterative solver manager.
    RCP< Belos::SolverManager<ST,MV,OP> > newSolver
      = rcp( new Belos::RCGSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&belosList,false)) );
    
    // Print out information about problem
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of RCG iterations: " << maxiters << std::endl;
      std::cout << "Max number of vectors in Krylov space: " << numBlocks << std::endl;
      std::cout << "Number of vectors in recycle space: " << recycleBlocks << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    
    // Perform solve
    Belos::ReturnType ret = newSolver->solve();
    
    // Compute actual residuals
    bool badRes = false;
    std::vector<ST> actual_resids( numrhs );
    std::vector<ST> rhs_norm( numrhs );
    MV resid(rowMap, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
    }

    if (ret!=Belos::Converged || badRes) {
      success = false;
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
    }
    else {
      success = true;
      if (proc_verbose)
        std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

int main(int argc, char *argv[]) {
  return run<double>(argc,argv);

  // wrapped with a check: CMake option Trilinos_ENABLE_FLOAT=ON
  // return run<float>(argc,argv);
}
