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
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosTFQMRSolMgr.hpp"
#include "BelosPseudoBlockTFQMRSolMgr.hpp"
#include "BelosTpetraTestFramework.hpp"
#include "BelosKokkosDenseAdapter.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>


int main(int argc, char *argv[]) {
  //
  typedef Tpetra::MultiVector<>::scalar_type ST;
  typedef Tpetra::MultiVector<>::impl_scalar_type IST;
  typedef Kokkos::DualView<IST**,Kokkos::LayoutLeft> DM;
  typedef Teuchos::ScalarTraits<ST>       SCT;
  typedef SCT::magnitudeType               MT;
  typedef Tpetra::Operator<ST>             OP;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV,DM> MVT;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const ST one  = SCT::one();

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);

  //
  bool success = false;
  bool verbose = false;
  try {
    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool explicit_test = false;
    bool comp_recursive = false;
    bool pseudo = false;
    int frequency = -1;  // how often residuals are printed by solver
    int numrhs = 1;  // total number of right-hand sides to solve for
    int maxiters = -1;  // maximum number of iterations for solver to use
    std::string filename("orsirr1.hb");
    double tol = 1.0e-5;  // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");

    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("explicit","implicit-only",&explicit_test,"Compute explicit residuals.");
    cmdp.setOption("recursive","native",&comp_recursive,"Compute recursive residuals.");
    cmdp.setOption("pseudo","not-pseudo",&pseudo,"Use pseudo-block TFQMR solver.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by TFQMR solver.");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem size).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose
    //
    // Get the problem
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

    proc_verbose = ( verbose && (MyPID==0) );
    //
    // Solve using Belos
    //
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (explicit_test)
    {
      belosList.set( "Explicit Residual Test", true );       // Scale by norm of right-hand side vector."
      belosList.set( "Explicit Residual Scaling", "Norm of RHS" ); // Scale by norm of right-hand side vector."
    }
    if (comp_recursive)
    {
      belosList.set( "Compute Recursive Residuals", true );
    }
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    //
    // Construct an unpreconditioned linear problem instance.
    //
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // Create an iterative solver manager.
    //
    RCP< Belos::SolverManager<ST,MV,OP,DM> > solver;

    if (pseudo)
      solver = rcp( new Belos::PseudoBlockTFQMRSolMgr<ST,MV,OP,DM>(rcp(&problem,false), rcp(&belosList,false)) );
    else
      solver = rcp( new Belos::TFQMRSolMgr<ST,MV,OP,DM>(rcp(&problem,false), rcp(&belosList,false)) );
    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of TFQMR iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Perform solve
    //
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
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );

} // end test_tfqmr_hb.cpp
