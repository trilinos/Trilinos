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
// Multiple right-hand-sides are created randomly.
// The initial guesses are all set to zero.
//
// This test is for testing the deflation in the pseudo-block Gmres solver.
// One set of linear systems is solved and then augmented with additional
// linear systems and resolved.  The already solved linear systems should be
// deflated immediately, leaving only the augmented systems to be solved.
//
//

// Tpetra
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraTestFramework.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"


template <typename ScalarType>
int run(int argc, char *argv[]) {

  using ST = typename Tpetra::CrsMatrix<ScalarType>::scalar_type;
  using LO = typename Tpetra::CrsMatrix<>::local_ordinal_type;
  using GO = typename Tpetra::CrsMatrix<>::global_ordinal_type;
  using NT = typename Tpetra::CrsMatrix<>::node_type;

  using OP  = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV  = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using MT = typename Teuchos::ScalarTraits<ST>::magnitudeType;

  using tmap_t       = Tpetra::Map<LO,GO,NT>;
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST,LO,GO,NT>;

  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  Teuchos::GlobalMPISession session(&argc, &argv, nullptr);

  bool verbose = false;
  bool success = true;

  try {

    const auto comm = Tpetra::getDefaultComm();
    const int myPID = comm->getRank();

    bool procVerbose = false;
    int frequency = -1;    // how often residuals are printed by solver
    int initNumRHS = 5;   // how many right-hand sides get solved first
    int augNumRHS = 10;   // how many right-hand sides are augmented to the first group
    int maxRestarts = 15;  // number of restarts allowed
    int length = 100;
    int initBlockSize = 5;// blocksize used for the initial pseudo-block GMRES solve
    int augBlockSize = 3; // blocksize used for the augmented pseudo-block GMRES solve
    int maxIters = -1;     // maximum iterations allowed
    std::string filename("orsirr1.hb");
    MT tol = 1.0e-5;       // relative residual tolerance
    MT aug_tol = 1.0e-5;   // relative residual tolerance for augmented system

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("aug-tol",&aug_tol,"Relative residual tolerance used by GMRES solver for augmented systems.");
    cmdp.setOption("init-num-rhs",&initNumRHS,"Number of right-hand sides to be initially solved for.");
    cmdp.setOption("aug-num-rhs",&augNumRHS,"Number of right-hand sides augmenting the initial solve.");
    cmdp.setOption("max-restarts",&maxRestarts,"Maximum number of restarts allowed for GMRES solver.");
    cmdp.setOption("block-size",&initBlockSize,"Block size used by GMRES for the initial solve.");
    cmdp.setOption("aug-block-size",&augBlockSize,"Block size used by GMRES for the augmented solve.");
    cmdp.setOption("max-iters",&maxIters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    cmdp.setOption("subspace-size",&length,"Dimension of Krylov subspace used by GMRES.");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    procVerbose = verbose && (myPID==0); /* Only print on zero processor */

    // Get the problem
    Belos::Tpetra::HarwellBoeingReader<tcrsmatrix_t> reader( comm );
    RCP<tcrsmatrix_t> A = reader.readFromFile( filename );
    RCP<const tmap_t> Map = A->getDomainMap();

    // ********Other information used by block solver***********
    // *****************(can be user specified)******************

    const int numGlobalElements = Map->getGlobalNumElements();
    if (maxIters == -1)
      maxIters = numGlobalElements - 1; // maximum number of iterations to run

    ParameterList belosList;
    belosList.set( "Num Blocks", length );                 // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", initBlockSize );         // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxIters );       // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxRestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Deflation Quorum", initBlockSize  );  // Number of converged linear systems before deflation
    belosList.set( "Timer Label", "Belos Init" );          // Label used by the timers in this solver
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

    // *****Construct solution std::vector and random right-hand-sides *****

    RCP<MV> initX = rcp( new MV(Map, initNumRHS) );
    RCP<MV> initB = rcp( new MV(Map, initNumRHS) );
    MVT::MvRandom( *initX );
    OPT::Apply( *A, *initX, *initB );
    initX->putScalar( 0.0 );
    Belos::LinearProblem<ST,MV,OP> initProblem( A, initX, initB );
    initProblem.setLabel("Belos Init");

    bool set = initProblem.setProblem();

    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Initial Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // *******************************************************************
    // *********************Perform initial solve*************************
    // *******************************************************************

    RCP< Belos::SolverManager<ST,MV,OP> > initSolver
      = rcp( new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>( rcp(&initProblem,false), rcp(&belosList,false) ) );

    // Perform solve
    Belos::ReturnType ret = initSolver->solve();

    // Compute actual residuals
    bool badRes = false;
    std::vector<ST> actualResids( initNumRHS );
    std::vector<ST> rhsNorm( initNumRHS );
    MV initR( Map, initNumRHS );
    OPT::Apply( *A, *initX, initR );
    MVT::MvAddMv( -1.0, initR, 1.0, *initB, initR );
    MVT::MvNorm( initR, actualResids );
    MVT::MvNorm( *initB, rhsNorm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for (int i=0; i<initNumRHS; i++) {
        ST actRes = actualResids[i]/rhsNorm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
    }

    if (ret != Belos::Converged || badRes==true) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Initial solve did not converge to solution!" << std::endl;
      return -1;
    }

    // ***************Construct augmented linear system****************

    RCP<MV> augX = rcp( new MV(Map, initNumRHS+augNumRHS) );
    RCP<MV> augB = rcp( new MV(Map, initNumRHS+augNumRHS) );
    if (augNumRHS) {
      MVT::MvRandom( *augX );
      OPT::Apply( *A, *augX, *augB );
      augX->putScalar( 0.0 );
    }

    // Copy previous linear system into
    RCP<MV> tmpX = rcp( new MV( *augX ) );
    RCP<MV> tmpB = rcp( new MV( *augB ) );
    tmpX->scale( 1.0, *augX );
    tmpB->scale( 1.0, *augB );

    Belos::LinearProblem<ST,MV,OP> augProblem( A, augX, augB );
    augProblem.setLabel("Belos Aug");

    set = augProblem.setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Augmented Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // *******************************************************************
    // *******************Perform augmented solve*************************
    // *******************************************************************

    belosList.set( "Block Size", augBlockSize );                // Blocksize to be used by iterative solver
    belosList.set( "Convergence Tolerance", aug_tol );           // Relative convergence tolerance requested
    belosList.set( "Deflation Quorum", augBlockSize );          // Number of converged linear systems before deflation
    belosList.set( "Timer Label", "Belos Aug" );          // Label used by the timers in this solver
    belosList.set( "Implicit Residual Scaling", "Norm of RHS" ); // Implicit residual scaling for convergence
    belosList.set( "Explicit Residual Scaling", "Norm of RHS" ); // Explicit residual scaling for convergence
    RCP< Belos::SolverManager<ST,MV,OP> > augSolver
      = rcp( new Belos::PseudoBlockGmresSolMgr<ST,MV,OP>( rcp(&augProblem,false), rcp(&belosList,false) ) );

    // Perform solve
    ret = augSolver->solve();

    if (ret != Belos::Converged) {
      if (procVerbose)
        std::cout << std::endl << "ERROR: Augmented solver did not converge to solution!" << std::endl;
      return -1;
    }

    // **********Print out information about problem*******************

    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << numGlobalElements << std::endl;
      std::cout << "Number of initial right-hand sides: " << initNumRHS << std::endl;
      std::cout << "Number of augmented right-hand sides: " << augNumRHS << std::endl;
      std::cout << "Number of restarts allowed: " << maxRestarts << std::endl;
      std::cout << "Length of block Arnoldi factorization: " << length <<std::endl;
      std::cout << "Max number of Gmres iterations: " << maxIters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      if (aug_tol != tol)
        std::cout << "Relative residual tolerance for augmented systems: " << aug_tol << std::endl;
      std::cout << std::endl;
    }

    // Compute actual residuals.
    badRes = false;
    int total_numrhs = initNumRHS + augNumRHS;
    actualResids.resize( total_numrhs );
    rhsNorm.resize( total_numrhs );
    MV augR( Map, total_numrhs );
    OPT::Apply( *A, *augX, augR );
    MVT::MvAddMv( -1.0, augR, 1.0, *augB, augR );
    MVT::MvNorm( augR, actualResids );
    MVT::MvNorm( *augB, rhsNorm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<total_numrhs; i++) {
        ST actRes = actualResids[i]/rhsNorm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol ) badRes = true;
      }
    }
    if (ret!=Belos::Converged || badRes==true) {
      success = false;
      if (procVerbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
    } else {
      success = true;
      if (procVerbose)
        std::cout << "End Result: TEST PASSED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose,std::cerr,success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
} // run

int main(int argc, char *argv[]) {
  // run with different ST
  return run<double>(argc,argv);
  // run<float>(argc,argv); // FAILS
}
