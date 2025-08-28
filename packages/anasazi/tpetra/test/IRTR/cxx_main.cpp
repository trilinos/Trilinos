// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This test is for IRTR solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the RTRSolMgr solver manager.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziRTRSolMgr.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <ModeLaplace1DQ1.hpp>

using namespace Teuchos;

int main(int argc, char *argv[]) 
{

  using std::cout;
  using std::endl;

  bool boolret;

  Tpetra::ScopeGuard tpetraScope (&argc,&argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  int MyPID = comm->getRank();

  bool testFailed = false;
  bool verbose = false;
  bool debug = false;
  bool shortrun = false;
  bool skinny = true;
  bool useSA = false;
  std::string which("SR");
  int numElements = 100;
  int user_verbosity = 0;
  int numicgs = 2;

  bool success = true;
  try {

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
    cmdp.setOption("sort",&which,"Targetted eigenvalues (SR or LR).");
    cmdp.setOption("shortrun","longrun",&shortrun,"Allow only a small number of iterations.");
    cmdp.setOption("skinny","hefty",&skinny,"Use a skinny (low-mem) or hefty (higher-mem) implementation of IRTR.");
    cmdp.setOption("numElements",&numElements,"Number of elements in discretization.");
    cmdp.setOption("useSA","noSA",&useSA,"Use subspace acceleration.");
    cmdp.setOption("verbosity",&user_verbosity,"Additional verbosity falgs.");
    cmdp.setOption("numICGS",&numicgs,"Num ICGS iterations");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (debug) verbose = true;

    typedef double                              ST;
    typedef Teuchos::ScalarTraits<ST>          SCT;
    typedef SCT::magnitudeType                  MT;
    typedef Tpetra::MultiVector<ST>             MV;
    typedef MV::global_ordinal_type             GO;
    typedef MV::local_ordinal_type              LO;
    typedef MV::node_type                     Node;
    typedef Tpetra::Operator<ST>                OP;
    typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
    typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
    const ST ONE  = SCT::one();

    if (verbose && MyPID == 0) {
      cout << Anasazi::Anasazi_Version() << endl << endl;
    }

    //  Problem information
    int space_dim = 1;
    std::vector<double> brick_dim( space_dim );
    brick_dim[0] = 1.0;
    std::vector<int> elements( space_dim );
    elements[0] = numElements;

    // Create problem
    ModeLaplace1DQ1<ST,LO,GO,Node> testCase( comm, brick_dim[0], elements[0]);
    //
    // Get the stiffness and mass matrices
    RCP<const Tpetra::CrsMatrix<ST,LO,GO,Node> > K = testCase.getStiffness();
    RCP<const Tpetra::CrsMatrix<ST,LO,GO,Node> > M = testCase.getMass();
    //
    // Create the initial vectors
    int blockSize = 5;
    RCP<MV> ivec = rcp (new MV (K->getDomainMap(),blockSize));
    ivec->randomize();

    // Create eigenproblem
    const int nev = 5;
    RCP<Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
      rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(K,M,ivec) );
    //
    // Inform the eigenproblem that the operator K is symmetric
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
    verbosity |= user_verbosity;


    // Eigensolver parameters
    int maxIters;
    if (shortrun) {
      maxIters = 50;
    }
    else {
      maxIters = 450;
    }
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
    MyPL.set( "Use SA", useSA );
    MyPL.set( "Num ICGS", numicgs );
    //
    // Create the solver manager
    Anasazi::RTRSolMgr<ST,MV,OP> MySolverMan(problem, MyPL);

    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = MySolverMan.solve();
    testFailed = false;
    if (returnCode != Anasazi::Converged) {
      if (shortrun==false) {
        testFailed = true;
      }
    }

    /*
    // 
    // Check that the parameters were all consumed
    if (MyPL.getEntryPtr("Verbosity")->isUsed() == false ||
        MyPL.getEntryPtr("Which")->isUsed() == false ||
        MyPL.getEntryPtr("Skinny Solver")->isUsed() == false ||
        MyPL.getEntryPtr("Block Size")->isUsed() == false ||
        MyPL.getEntryPtr("Maximum Iterations")->isUsed() == false ||
        MyPL.getEntryPtr("Convergence Tolerance")->isUsed() == false) {
      if (verbose && MyPID==0) {
        cout << "Failure! Unused parameters: " << endl;
        MyPL.unused(cout);
      }
    }
    */

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
    std::vector<Anasazi::Value<ST> > evals = sol.Evals;
    RCP<MV> evecs = sol.Evecs;
    int numev = sol.numVecs;

    if (numev > 0) {

      std::ostringstream os;
      os.setf(std::ios::scientific, std::ios::floatfield);
      os.precision(6);

      // Compute the direct residual
      std::vector<ST> normV( numev );
      SerialDenseMatrix<int,ST> T(numev,numev);
      for (int i=0; i<numev; i++) {
        T(i,i) = evals[i].realpart;
      }
      RCP<MV> Mvecs = MVT::Clone( *evecs, numev ),
        Kvecs = MVT::Clone( *evecs, numev );
      OPT::Apply( *K, *evecs, *Kvecs );
      OPT::Apply( *M, *evecs, *Mvecs );
      MVT::MvTimesMatAddMv( -ONE, *Mvecs, T, ONE, *Kvecs );
      // compute M-norm of residuals
      OPT::Apply( *M, *Kvecs, *Mvecs );
      MVT::MvDot( *Mvecs, *Kvecs, normV );

      os << "Number of iterations performed in RTR_test.exe: " << MySolverMan.getNumIters() << endl
        << "Direct residual norms computed in RTR_test.exe" << endl
        << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual(M)" << endl
        << "----------------------------------------" << endl;
      for (int i=0; i<numev; i++) {
        if ( SCT::magnitude(evals[i].realpart) != SCT::zero() ) {
          normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) / evals[i].realpart );
        }
        else {
          normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) );
        }
        os << std::setw(20) << evals[i].realpart << std::setw(20) << normV[i] << endl;
        if ( normV[i] > tol ) {
          testFailed = true;
        }
      }
      if (verbose && MyPID==0) {
        cout << endl << os.str() << endl;
      }

    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,cout,success);

  if (testFailed || success==false) {
    if (verbose && MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;

}
