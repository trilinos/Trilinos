// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// Update of belos/epetra/example/GCRODR/GCRODREpetraExFile.cpp
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this example.

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
#include "BelosGCRODRSolMgr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraTestFramework.hpp"


template <typename ScalarType>
int run(int argc, char *argv[]) {

  using ST = typename Tpetra::Vector<ScalarType>::scalar_type;
  using LO = typename Tpetra::Vector<>::local_ordinal_type;
  using GO = typename Tpetra::Vector<>::global_ordinal_type;
  using NT = typename Tpetra::Vector<>::node_type;

  using OP = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV = typename Tpetra::MultiVector<ST,LO,GO,NT>;
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
  const auto Comm = Tpetra::getDefaultComm();
  const int MyPID = Comm->getRank();

  bool verbose = false;
  bool success = true;

  try {
    bool proc_verbose = false;
    bool debug = false;
    int frequency = -1;        // frequency of status test output.
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxiters = -1;         // maximum number of iterations allowed per linear system
    int maxsubspace = 250;     // maximum number of blocks the solver can use for the subspace
    int recycle = 50;          // maximum number of blocks the solver can use for the subspace
    int maxrestarts = 15;      // number of restarts allowed
    std::string filename("sherman5.hb");
    std::string ortho("IMGS");
    MT tol = 1.0e-8;           // relative residual tolerance

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
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    proc_verbose = ( verbose && (MyPID==0) );

    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }

    // Get the problem
    Belos::Tpetra::HarwellBoeingReader<tcrsmatrix_t> reader( Comm );
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

    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run

    ParameterList belosList;
    belosList.set( "Num Blocks", maxsubspace);        // Maximum number of blocks in Krylov factorization
    belosList.set( "Maximum Iterations", maxiters );  // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxrestarts ); // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );    // Relative convergence tolerance requested
    belosList.set( "Num Recycled Blocks", recycle );  // Number of vectors in recycle space
    belosList.set( "Orthogonalization", ortho );      // Orthogonalization type

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

    // *******************************************************************
    // *************Start the GCRO-DR iteration*************************
    // *******************************************************************

    // Create an iterative solver manager.
    RCP< Belos::SolverManager<ST,MV,OP> > newSolver
      = rcp( new Belos::GCRODRSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));

    // **********Print out information about problem*******************

    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of restarts allowed: " << maxrestarts << std::endl;
      std::cout << "Max number of iterations per restart cycle: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }

    // Perform solve
    Belos::ReturnType ret = newSolver->solve();

    // Get the number of iterations for this solve.
    int numIters = newSolver->getNumIters();
    std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;

    // Compute actual residuals.
    bool badRes = false;
    std::vector<ST> actual_resids( numrhs );
    std::vector<ST> rhs_norm( numrhs );
    MV resid(Map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        ST actRes = actual_resids[i]/rhs_norm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
    }

    if (ret!=Belos::Converged || badRes) {
      success = false;
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
    } else {
      success = true;
      if (proc_verbose)
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
