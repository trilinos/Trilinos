// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  This test is for the internal utilities that are used by the modal analysis solver.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSimpleLOBPCGSolMgr.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include "ModeLaplace1DQ1.h"

using namespace Teuchos;

int main(int argc, char *argv[]) 
{
  bool boolret;
  int MyPID;

#ifdef HAVE_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  MyPID = Comm.MyPID();

  bool testFailed;
  bool verbose = false;
  bool debug = false;
  bool shortrun = false;
  std::string filename("mhd1280b.cua");
  std::string which("LM");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM or LM).");
  cmdp.setOption("shortrun","longrun",&shortrun,"Allow only a small number of iterations.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (debug) verbose = true;

  typedef double ScalarType;
  typedef ScalarTraits<ScalarType>                   SCT;
  typedef SCT::magnitudeType               MagnitudeType;
  typedef Epetra_MultiVector                          MV;
  typedef Epetra_Operator                             OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;
  const ScalarType ONE  = SCT::one();

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
  RCP<ModalProblem> testCase = rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );
  //
  // Get the stiffness and mass matrices
  RCP<const Epetra_CrsMatrix> K = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getStiffness()), false );
  RCP<const Epetra_CrsMatrix> M = rcp( const_cast<Epetra_CrsMatrix *>(testCase->getMass()), false );
  //
  // Create the initial vectors
  int blockSize = 2;
  RCP<Anasazi::EpetraMultiVec> ivec = rcp( new Anasazi::EpetraMultiVec(K->OperatorDomainMap(), blockSize) );
  ivec->Random();
  //
  // Create eigenproblem
  const int nev = 5;
  RCP<Anasazi::BasicEigenproblem<ScalarType,MV,OP> > problem =
    rcp( new Anasazi::BasicEigenproblem<ScalarType,MV,OP>(K,M,ivec) );
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
      std::cout << "Anasazi::BasicEigenproblem::SetProblem() returned with error." << std::endl
           << "End Result: TEST FAILED" << std::endl;	
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
  int maxIters;
  if (shortrun) {
    maxIters = 35;
    blockSize = nev;
  }
  else {
    maxIters = 450;
  }
  MagnitudeType tol = 1.0e-6;
  //
  // Create parameter list to pass into the solver manager
  ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIters );
  MyPL.set( "Convergence Tolerance", tol );
  //
  // Create the solver manager
  Anasazi::SimpleLOBPCGSolMgr<ScalarType,MV,OP> MySolverMan(problem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMan.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged && shortrun==false) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = problem->getSolution();
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  if (numev > 0) {

    // Check the problem against the analytical solutions
    if (verbose) {
      double *revals = new double[numev];
      for (int i=0; i<numev; i++) {
        revals[i] = sol.Evals[i].realpart;
      }
      bool smallest = false;
      if (which == "SM" || which == "SR") {
        smallest = true;
      }
      testCase->eigenCheck( *evecs, revals, 0, smallest );
      delete [] revals;
    }

    // Compute the direct residual
    std::vector<ScalarType> normV( numev );
    SerialDenseMatrix<int,ScalarType> T(numev,numev);
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
  
    for (int i=0; i<numev; i++) {
      normV[i] = SCT::squareroot( normV[i] );
      if ( SCT::magnitude(normV[i]/sol.Evals[i].realpart) > tol ) {
        testFailed = true;
      }
    }

  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

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
