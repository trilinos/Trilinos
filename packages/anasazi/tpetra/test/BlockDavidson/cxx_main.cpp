// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This test is for BlockDavidson solving a generalized (Ax=Bxl) real Hermitian
// eigenvalue problem, using the BlockDavidsonSolMgr solver manager.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <ModeLaplace1DQ1.hpp>

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using Tpetra::CrsMatrix;
  using std::cout;
  using std::endl;

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

  Tpetra::ScopeGuard tpetraScope (&argc,&argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  int MyPID = comm->getRank();

  bool testFailed;
  bool verbose = false;
  bool debug = false;
  bool shortrun = false;
  bool insitu = false;
  bool locking = true;
  std::string which("LM");
  int  rblocks = 1;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("insitu","exsitu",&insitu,"Perform in situ restarting.");
  cmdp.setOption("locking","nolocking",&locking,"Perform locking.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("shortrun","longrun",&shortrun,"Allow only a small number of iterations.");
  cmdp.setOption("rblocks",&rblocks,"Number of blocks after restart.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) verbose = true;

  if (verbose && MyPID == 0) {
    std::cout << Anasazi::Anasazi_Version() << std::endl << std::endl;
  }

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;

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
  ivec->randomize ();

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
  bool boolret = problem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      std::cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << std::endl
           << "End Result: TEST FAILED" << std::endl;	
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


  // Eigensolver parameters
  int numBlocks = 8;
  int maxRestarts;
  if (shortrun) {
    maxRestarts = 10;
  }
  else {
    maxRestarts = 100;
  }
  MT tol = 1.0e-6;
  //
  // Create parameter list to pass into the solver manager
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Use Locking", locking );
  MyPL.set( "Locking Tolerance", tol/10 );
  MyPL.set( "In Situ Restarting", insitu );
  MyPL.set( "Num Restart Blocks", rblocks );
  //
  // Create the solver manager
  Anasazi::BlockDavidsonSolMgr<ST,MV,OP> MySolverMan(problem, MyPL);
  // 
  // Check that the parameters were all consumed
  if (MyPL.getEntryPtr("Verbosity")->isUsed() == false ||
      MyPL.getEntryPtr("Which")->isUsed() == false ||
      MyPL.getEntryPtr("Block Size")->isUsed() == false ||
      MyPL.getEntryPtr("Num Blocks")->isUsed() == false ||
      MyPL.getEntryPtr("Maximum Restarts")->isUsed() == false ||
      MyPL.getEntryPtr("Convergence Tolerance")->isUsed() == false ||
      MyPL.getEntryPtr("Use Locking")->isUsed() == false ||
      MyPL.getEntryPtr("In Situ Restarting")->isUsed() == false ||
      MyPL.getEntryPtr("Num Restart Blocks")->isUsed() == false || 
      MyPL.getEntryPtr("Locking Tolerance")->isUsed() == false) {
    if (verbose && MyPID==0) {
      std::cout << "Failure! Unused parameters: " << std::endl;
      MyPL.unused(std::cout);
    }
  }

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMan.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged && shortrun==false) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ST,MV> sol = problem->getSolution();
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<ST> normV( numev );
    Teuchos::SerialDenseMatrix<int,ST> T(numev,numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = sol.Evals[i].realpart;
    }
    RCP<MV> Mvecs = MVT::Clone( *evecs, numev ),
                    Kvecs = MVT::Clone( *evecs, numev );
    OPT::Apply( *K, *evecs, *Kvecs );
    OPT::Apply( *M, *evecs, *Mvecs );
    MVT::MvTimesMatAddMv( -ONE, *Mvecs, T, ONE, *Kvecs );
    // compute M-norm of residuals
    OPT::Apply( *M, *Kvecs, *Mvecs );
    MVT::MvDot( *Mvecs, *Kvecs, normV );

    os << "Number of iterations performed in BlockDavidson_test.exe: " << MySolverMan.getNumIters() << std::endl
       << "Direct residual norms computed in BlockDavidson_test.exe" << std::endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual(M)" << std::endl
       << "----------------------------------------" << std::endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(sol.Evals[i].realpart) != SCT::zero() ) {
        normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) / sol.Evals[i].realpart );
      }
      else {
        normV[i] = SCT::magnitude( SCT::squareroot( normV[i] ) );
      }
      os << std::setw(20) << sol.Evals[i].realpart << std::setw(20) << normV[i] << std::endl;
      if ( normV[i] > tol ) {
        testFailed = true;
      }
    }
    if (verbose && MyPID==0) {
      std::cout << std::endl << os.str() << std::endl;
    }

  }

  if (testFailed) {
    if (verbose && MyPID==0) {
      std::cout << "End Result: TEST FAILED" << std::endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  return 0;

}	
