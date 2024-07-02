// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This test is for SIRTR/IRTR solving a standard (Ax=xl) complex Hermitian
// eigenvalue problem where the operator (A) is the 1D finite-differenced Laplacian
// operator.  In this test, the eigenproblem is declared non-Hermitian even though
// the operator is Hermitian.  This illustrated the bug found in bug 2840.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziFactory.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// templated multivector and sparse matrix classes
#include "MyMultiVec.hpp"
#include "MyOperator.hpp"

using namespace Teuchos;

int main(int argc, char *argv[])
{
  using std::cout;
  using std::endl;

  int MyPID = 0;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
#endif

  bool success = false;
  bool verbose = false;
  try {
    bool boolret;
    bool debug = false;
    bool skinny = true;
    bool fakeM = true;
    std::string which("LR");

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("skinny","hefty",&skinny,"Use a skinny (low-mem) or hefty (higher-mem) implementation of IRTR.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
    cmdp.setOption("sort",&which,"Targetted eigenvalues (SR or LR).");
    cmdp.setOption("fakeM","noM",&fakeM,"Use an explicit identity matrix for the mass matrix.");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return -1;
    }

    typedef std::complex<double>                ST;
    typedef ScalarTraits<ST>                   SCT;
    typedef SCT::magnitudeType                  MT;
    typedef Anasazi::MultiVec<ST>               MV;
    typedef Anasazi::Operator<ST>               OP;
    typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
    typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
    ST ONE  = SCT::one();

    if (verbose && MyPID == 0) {
      cout << Anasazi::Anasazi_Version() << endl << endl;
    }

    // -- Set finite difference grid
    int dim = 10;

    // Build the problem matrix
    RCP< const MyOperator<ST> > K
      = rcp( new MyOperator<ST>(dim) );
    // build the identity matrix
    RCP< const MyOperator<ST> > I;
    if (fakeM) {
      std::vector<ST> diag(3);
      diag[0] = Teuchos::as<ST>(0); // lower
      diag[1] = Teuchos::as<ST>(1); // diag
      diag[2] = Teuchos::as<ST>(0); // upper
      I = rcp( new MyOperator<ST>(dim,diag) );
    }

    // Create initial vectors
    const int blockSize = 2;
    RCP<MyMultiVec<ST> > ivec = rcp( new MyMultiVec<ST>(dim,blockSize) );
    ivec->MvRandom();

    // Create eigenproblem
    const int nev = 1;
    RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
      rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(K,I,ivec) );
    //
    // Inform the eigenproblem that the operator K is non-Hermitian (even when it truly is Hermitian)
    problem->setHermitian(true);
    //
    // Set the number of eigenvalues requested
    problem->setNEV( nev );
    //
    // Inform the eigenproblem that you are done passing it information
    boolret = problem->setProblem();
    if (boolret != true) {
      if (verbose && MyPID == 0) {
        cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << endl
          << "End Result: TEST FAILED" << endl;
      }
#ifdef HAVE_MPI
      MPI_Finalize() ;
#endif
      return -1;
    }

    // Set verbosity level
    int verbosity = Anasazi::Errors + Anasazi::Warnings;
    if (verbose) {
      verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
    }
    if (debug) {
      verbosity += Anasazi::Debug;
    }

    // Eigensolver parameters
    int maxIters = 200;
    MT tol = 1.0e-6;
    //
    // Create parameter list to pass into the solver manager
    ParameterList MyPL;
    MyPL.set( "Skinny Solver", skinny);
    MyPL.set( "Verbosity", verbosity );
    MyPL.set( "Which", which );
    MyPL.set( "Block Size", blockSize );
    MyPL.set( "Maximum Iterations", maxIters );
    MyPL.set( "Convergence Tolerance", tol );
    //
    // Create the solver manager
    auto MySolverMgr = Anasazi::Factory::create("RTR", problem, MyPL);

    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMgr->solve();
    success = (returnCode == Anasazi::Converged);

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
    RCP<MV> evecs = sol.Evecs;
    int numev = sol.numVecs;

    if (numev > 0) {

      std::ostringstream os;
      os.setf(std::ios::scientific, std::ios::floatfield);
      os.precision(6);

      // Compute the direct residual
      std::vector<MT> normV( numev );
      SerialDenseMatrix<int,ST> T(numev,numev);
      for (int i=0; i<numev; i++) {
        T(i,i) = sol.Evals[i].realpart;
      }
      RCP<MV> Kvecs = MVT::Clone( *evecs, numev );

      OPT::Apply( *K, *evecs, *Kvecs );

      MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Kvecs );
      MVT::MvNorm( *Kvecs, normV );

      os << "Direct residual norms computed in IRTRComplex_test.exe" << endl
        << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual  " << endl
        << "----------------------------------------" << endl;
      for (int i=0; i<numev; i++) {
        if ( SCT::magnitude(sol.Evals[i].realpart) != SCT::zero() ) {
          normV[i] = SCT::magnitude(normV[i]/sol.Evals[i].realpart);
        }
        os << std::setw(20) << sol.Evals[i].realpart << std::setw(20) << normV[i] << endl;
        success = ( normV[i] < tol );
      }
      if (verbose && MyPID==0) {
        cout << endl << os.str() << endl;
      }

    }

    if (verbose && MyPID==0) {
      if (success)
        cout << "End Result: TEST PASSED" << endl;
      else
        cout << "End Result: TEST FAILED" << endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
