// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver tests variations of Hybrid GMRES using the Belos GMRES 
// Polynomial solver manager.  One can specify an 'outer solver'
// to which the GMRES Polynomial is applied as a preconditioner.  
// If the 'outer solver' string is empty, the polynomial will be applied
// alone.
//
// Alternately, one can apply the GMRES polynomial without the 'outer
// solver' functionality by using the BelosTpetraOperator interface.
// See the FGMRES test driver for an example of its use.  
//
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero.
//
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosGmresPolySolMgr.hpp"

// I/O for Harwell-Boeing files
#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

using namespace Teuchos;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using std::endl;
using std::cout;
using std::vector;
using Teuchos::tuple;

int main(int argc, char *argv[]) {
  typedef Tpetra::MultiVector<>::scalar_type ST;
  typedef ScalarTraits<ST>                SCT;
  typedef SCT::magnitudeType               MT;
  typedef Tpetra::Operator<ST>             OP;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;

  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  bool success = false;
  bool verbose = false;
  try {
    const ST one  = SCT::one();

    int MyPID = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool userandomrhs = false;                // use linear problem RHS or random RHS to generate poly
    int frequency = -1;                       // frequency of status test output.
    int blocksize = 1;                        // blocksize
    int numrhs = 1;                           // number of right-hand sides to solve for
    int maxiters = -1;                        // maximum number of iterations allowed per linear system
    int maxdegree = 25;                       // maximum degree of polynomial
    int maxsubspace = 50;                     // maximum number of blocks the solver can use for the subspace
    int maxrestarts = 15;                     // number of restarts allowed
    std::string outersolver("Block Gmres");   // name of outer solver
    std::string polytype("Roots");            // polynomial configuration.  
    std::string filename("bcsstk14.hb");      // name of matrix file
    MT tol = 1.0e-5;                          // relative residual tolerance
    MT polytol = tol/10;                      // relative residual tolerance for polynomial construction

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("use-random-rhs","use-rhs",&userandomrhs,"Use linear system RHS or random RHS to generate polynomial.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
    cmdp.setOption("outersolver",&outersolver,"Name of outer solver to be used with GMRES poly");
    cmdp.setOption("poly-type",&polytype,"Name of the polynomial to be generated. Arnoldi, Gmres, or Roots.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("poly-tol",&polytol,"Relative residual tolerance used to construct the GMRES polynomial.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("block-size",&blocksize,"Block size used by GMRES.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    cmdp.setOption("max-degree",&maxdegree,"Maximum degree of the GMRES polynomial.");
    cmdp.setOption("max-subspace",&maxsubspace,"Maximum number of blocks the solver can use for the subspace.");
    cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GMRES solver.");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return EXIT_FAILURE;
    }
    if (!verbose) {
      frequency = -1;  // reset frequency if test is not verbose
    }

    MyPID = rank(*comm);
    proc_verbose = ( verbose && (MyPID==0) );

    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }

    //
    // Get the data from the HB file and build the Map,Matrix
    //
    RCP<CrsMatrix<ST> > A;
    Tpetra::Utils::readHBMatrix(filename,comm,A);
    RCP<const Tpetra::Map<> > map = A->getDomainMap();

    // Create initial vectors
    RCP<MV> B, X;
    X = rcp( new MV(map,numrhs) );
    MVT::MvRandom( *X );
    B = rcp( new MV(map,numrhs) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );

    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1) {
      maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
    }
    //
    //
    //Parameter list used by the outer solver:
    ParameterList belosList;
    belosList.set( "Num Blocks", maxsubspace);             // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    int verbosity = Belos::Errors + Belos::Warnings;
    if (verbose) {
      verbosity += Belos::FinalSummary + Belos::TimingDetails + Belos::StatusTestDetails;
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    belosList.set( "Verbosity", verbosity );

    // Parameter list used by the GMRES Polynomial (inner) solver
    ParameterList polyList;
    polyList.set( "Polynomial Type", polytype );          // Type of polynomial to be generated
    polyList.set( "Maximum Degree", maxdegree );          // Maximum degree of the GMRES polynomial
    polyList.set( "Polynomial Tolerance", polytol );      // Polynomial convergence tolerance requested
    polyList.set( "Verbosity", verbosity );               // Verbosity for polynomial construction
    polyList.set( "Random RHS", userandomrhs );           // Use RHS from linear system or random vector
                                                          // to generate the polynomial.
                                                          
    //
    // Pass the outer solver name and parameters to the polynomial solver manager.
    //
    if ( outersolver != "" ) {
      polyList.set( "Outer Solver", outersolver );
      polyList.set( "Outer Solver Params", belosList );
    }

    //
    // Construct an unpreconditioned linear problem instance. 
    //
    // The polynomial preconditioning will be handled by the solver manager.
    // One could pass another preconditioner such as ILU to the linear problem
    // to be used in combination with the polynomial preconditioner.
    //
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    problem.setInitResVec( B );
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }
    //
    // *******************************************************************
    // *************Start the Hybrid Gmres iteration***********************
    // ********* Using the GMRES Poly Solver Manager *********************
    // *******************************************************************
    //
    Belos::GmresPolySolMgr<ST,MV,OP> solver( rcpFromRef(problem), rcpFromRef(polyList) );

    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Block size used by solver: " << blocksize << std::endl;
      std::cout << "Max number of Gmres iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }

    //
    // Perform solve
    //
    Belos::ReturnType ret = solver.solve();

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

    success = (ret==Belos::Converged && !badRes);

    if (success) {
      if (proc_verbose)
        std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_hybrid_gmres_hb.cpp
