// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This test is for BlockDavidson solving a standard (Ax=xl) complex Hermitian
// eigenvalue problem, using the BlockDavidsonSolMgr solver manager.
//
// The matrix used is from MatrixMarket:
// Name: MHD1280B: Alfven Spectra in Magnetohydrodynamics
// Source: Source: A. Booten, M.N. Kooper, H.A. van der Vorst, S. Poedts and J.P. Goedbloed University of Utrecht, the Netherlands
// Discipline: Plasma physics
// URL: http://math.nist.gov/MatrixMarket/data/NEP/mhd/mhd1280b.html
// Size: 1280 x 1280
// NNZ: 22778 entries

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziFactory.hpp"

#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// I/O for Harwell-Boeing files
#ifdef HAVE_ANASAZI_TRIUTILS
#include "Trilinos_Util_iohb.h"
#endif

// templated multivector and sparse matrix classes
#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"

using namespace Teuchos;

int main(int argc, char *argv[]) 
{
  using std::cout;
  using std::endl;
  int info = 0;
  bool boolret;
  int MyPID;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &MyPID );
#else
  MyPID = 0;
#endif

  bool testFailed;
  bool verbose = false;
  bool debug = false;
  bool shortrun = false;
  bool insitu = false;
  bool locking = true;
  std::string filename("mhd1280b.cua");
  std::string which("LM");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("insitu","exsitu",&insitu,"Perform in situ restarting.");
  cmdp.setOption("locking","nolocking",&locking,"Perform locking.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("shortrun","longrun",&shortrun,"Allow only a small number of iterations.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (debug) verbose = true;

#ifndef HAVE_ANASAZI_TRIUTILS
  cout << "This test requires Triutils. Please configure with --enable-triutils." << endl;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  if (verbose && MyPID == 0) {
    cout << "End Result: TEST FAILED" << endl;	
  }
  return -1;
#endif

  typedef std::complex<double>                ScalarType;
  typedef ScalarTraits<ScalarType>                   SCT;
  typedef SCT::magnitudeType               MagnitudeType;
  typedef Anasazi::MultiVec<ScalarType>               MV;
  typedef Anasazi::Operator<ScalarType>               OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;
  const ScalarType ONE  = SCT::one();

  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }

  //  Problem information
  int dim,dim2,nnz;
  double *dvals;
  int *colptr,*rowind;
  nnz = -1;
  info = readHB_newmat_double(filename.c_str(),&dim,&dim2,&nnz,
                              &colptr,&rowind,&dvals);
  if (info == 0 || nnz < 0) {
    if (verbose && MyPID == 0) {
      cout << "Error reading '" << filename << "'" << endl
           << "End Result: TEST FAILED" << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  // Convert interleaved doubles to complex values
  std::vector<ScalarType> cvals(nnz);
  for (int ii=0; ii<nnz; ii++) {
    cvals[ii] = ScalarType(dvals[ii*2],dvals[ii*2+1]);
  }
  // Build the problem matrix
  RCP< const MyBetterOperator<ScalarType> > K 
    = rcp( new MyBetterOperator<ScalarType>(dim,colptr,nnz,rowind,&cvals[0]) );

  // Create initial vectors
  int blockSize = 5;
  RCP<MyMultiVec<ScalarType> > ivec = rcp( new MyMultiVec<ScalarType>(dim,blockSize) );
  ivec->MvRandom();

  // Create eigenproblem
  const int nev = 4;
  RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(K,ivec) );
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

  // Get a formatted output stream
  RCP<Anasazi::OutputManager<ScalarType> > MyOM = rcp( new Anasazi::BasicOutputManager<ScalarType>() );
  RCP<Teuchos::FancyOStream> fos = MyOM->getFancyOStream();

  // Eigensolver parameters
  int numBlocks = 8;
  int maxRestarts;
  if (shortrun) {
    maxRestarts = 10;
  }
  else {
    maxRestarts = 100;
  }
  MagnitudeType tol = 1.0e-6;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Output Stream", fos );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Use Locking", locking );
  MyPL.set( "Locking Tolerance", tol/10 );
  MyPL.set( "In Situ Restarting", insitu );
  //
  // Create the solver manager
  RCP<Anasazi::SolverManager<ScalarType,MV,OP>> MySolverMan =
          Anasazi::Factory::create("BLOCK_DAVIDSON", problem, MyPL);
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
      MyPL.getEntryPtr("Locking Tolerance")->isUsed() == false) {
    if (verbose && MyPID==0) {
      cout << "Failure! Unused parameters: " << endl;
      MyPL.unused(cout);
    }
  }

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMan->solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged && shortrun==false) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = problem->getSolution();
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<MagnitudeType> normV( numev );
    SerialDenseMatrix<int,ScalarType> T(numev,numev);
    for (int i=0; i<numev; i++) {
      T(i,i) = sol.Evals[i].realpart;
    }
    RCP<MV> Kvecs = MVT::Clone( *evecs, numev );

    OPT::Apply( *K, *evecs, *Kvecs );

    MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Kvecs );
    MVT::MvNorm( *Kvecs, normV );

    os << "Direct residual norms computed in BlockDavidsonComplex_test.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual(M)" << endl
       << "----------------------------------------" << endl;
    for (int i=0; i<numev; i++) {
      if ( SCT::magnitude(sol.Evals[i].realpart) != SCT::zero() ) {
        normV[i] = SCT::magnitude(normV[i]/sol.Evals[i].realpart);
      }
      os << std::setw(20) << sol.Evals[i].realpart << std::setw(20) << normV[i] << endl;
      if ( normV[i] > tol ) {
        testFailed = true;
      }
    }
    if (verbose && MyPID==0) {
      cout << endl << os.str() << endl;
    }
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  // Clean up.
  std::free( dvals );
  std::free( colptr );
  std::free( rowind );
  
  if (testFailed) {
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
