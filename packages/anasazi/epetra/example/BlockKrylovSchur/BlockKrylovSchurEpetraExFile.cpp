// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//  This example compute the eigenvalues of a Harwell-Boeing matrix using the block Arnoldi
//  method.  The matrix is passed to the example routine through the command line, and 
//  converted to an Epetra matrix through some utilty routines.  This matrix is passed to the
//  eigensolver and then used to construct the Krylov decomposition.  The specifics of the 
//  block Arnoldi method can be set by the user.
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_CrsMatrix.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
  using std::cout;

  //
  bool haveM = false;

#ifdef EPETRA_MPI  
  // Initialize MPI  
  MPI_Init(&argc,&argv);   
  Epetra_MpiComm Comm( MPI_COMM_WORLD );  
#else  
  Epetra_SerialComm Comm;  
#endif
  
  int MyPID = Comm.MyPID();

  bool verbose=false;
  bool isHermitian=false;
  std::string k_filename = "";
  std::string m_filename = "";
  std::string which = "LM";
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("hermitian","non-hermitian",&isHermitian,"The eigenproblem being read in is Hermitian.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,or LR).");
  cmdp.setOption("K-filename",&k_filename,"Filename and path of the stiffness matrix.");
  cmdp.setOption("M-filename",&m_filename,"Filename and path of the mass matrix.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  //
  //**********************************************************************
  //******************Set up the problem to be solved*********************
  //**********************************************************************
  //
  // *****Read in matrix from file******
  //
  Teuchos::RCP<Epetra_Map> Map;
  Teuchos::RCP<Epetra_CrsMatrix> K, M;
  EpetraExt::readEpetraLinearSystem( k_filename, Comm, &K, &Map );

  if (haveM) {
    EpetraExt::readEpetraLinearSystem( m_filename, Comm, &M, &Map );
  }
  //
  //************************************
  // Start the block Arnoldi iteration
  //***********************************         
  //
  //  Variables used for the Block Arnoldi Method
  // 
  int nev = 5;
  int blockSize = 5;
  int numBlocks = 10;
  int maxRestarts = 10;
  int step = 5;
  double tol = 1.0e-8;
  
  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Step Size", step );

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;
  //
  // Create the eigenproblem to be solved.
  //
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(*Map, blockSize) );
  ivec->Random();
  
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem;
  if (haveM) {
    MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(K, M, ivec) );
  }
  else { 
    MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(K, ivec) );
  }

  // Inform the eigenproblem that the (K,M) is Hermitian
  MyProblem->setHermitian( isHermitian );

  // Set the number of eigenvalues requested 
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      std::cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << std::endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Initialize the Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);
    
  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && MyPID==0 && verbose) {
    std::cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
  }
  
  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  std::vector<int> index = sol.index;
  int numev = sol.numVecs;

  if (numev > 0) {
    // Compute residuals.
    Teuchos::LAPACK<int,double> lapack;
    std::vector<double> normR(numev);
    
    if (MyProblem->isHermitian()) {
      // Get storage
      Epetra_MultiVector Kevecs(*Map,numev);
      Teuchos::RCP<Epetra_MultiVector> Mevecs;
      Teuchos::SerialDenseMatrix<int,double> B(numev,numev);
      B.putScalar(0.0); 
      for (int i=0; i<numev; i++) {B(i,i) = evals[i].realpart;}
      
      // Compute A*evecs
      OPT::Apply( *K, *evecs, Kevecs );
      if (haveM) {
        Mevecs = Teuchos::rcp( new Epetra_MultiVector(*Map,numev) );
        OPT::Apply( *M, *evecs, *Mevecs );
      } 
      else {
        Mevecs = evecs;
      }
      
      // Compute A*evecs - lambda*evecs and its norm
      MVT::MvTimesMatAddMv( -1.0, *Mevecs, B, 1.0, Kevecs );
      MVT::MvNorm( Kevecs, normR );
      
      // Scale the norms by the eigenvalue
      for (int i=0; i<numev; i++) {
        normR[i] /= Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
      }
    } else {
      // The problem is non-Hermitian.
      int i=0;
      std::vector<int> curind(1);
      std::vector<double> resnorm(1), tempnrm(1);
      Teuchos::RCP<MV> tempKevec, Mevecs;
      Teuchos::RCP<const MV> tempeveci, tempMevec;
      Epetra_MultiVector Kevecs(*Map,numev);
      
      // Compute K*evecs
      OPT::Apply( *K, *evecs, Kevecs );
      if (haveM) {
        Mevecs = Teuchos::rcp( new Epetra_MultiVector(*Map,numev) );
        OPT::Apply( *M, *evecs, *Mevecs );
      }
      else {
        Mevecs = evecs;
      }
      
      Teuchos::SerialDenseMatrix<int,double> Breal(1,1), Bimag(1,1);
      while (i<numev) {
        if (index[i]==0) {
          // Get a view of the M*evecr
          curind[0] = i;
          tempMevec = MVT::CloneView( *Mevecs, curind );

          // Get a copy of A*evecr
          tempKevec = MVT::CloneCopy( Kevecs, curind );

          // Compute K*evecr - lambda*M*evecr
          Breal(0,0) = evals[i].realpart;
          MVT::MvTimesMatAddMv( -1.0, *tempMevec, Breal, 1.0, *tempKevec );

          // Compute the norm of the residual and increment counter
          MVT::MvNorm( *tempKevec, resnorm );
          normR[i] = resnorm[0]/Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
          i++;
        } else {
          // Get a view of the real part of M*evecr
          curind[0] = i;
          tempMevec = MVT::CloneView( *Mevecs, curind );

          // Get a copy of K*evecr
          tempKevec = MVT::CloneCopy( Kevecs, curind );

          // Get a view of the imaginary part of the eigenvector (eveci)
          curind[0] = i+1;
          tempeveci = MVT::CloneView( *Mevecs, curind );

          // Set the eigenvalue into Breal and Bimag
          Breal(0,0) = evals[i].realpart;
          Bimag(0,0) = evals[i].imagpart;

          // Compute K*evecr - M*evecr*lambdar + M*eveci*lambdai
          MVT::MvTimesMatAddMv( -1.0, *tempMevec, Breal, 1.0, *tempKevec );
          MVT::MvTimesMatAddMv( 1.0, *tempeveci, Bimag, 1.0, *tempKevec );
          MVT::MvNorm( *tempKevec, tempnrm );

          // Get a copy of K*eveci
          tempKevec = MVT::CloneCopy( Kevecs, curind );

          // Compute K*eveci - M*eveci*lambdar - M*evecr*lambdai
          MVT::MvTimesMatAddMv( -1.0, *tempMevec, Bimag, 1.0, *tempKevec );
          MVT::MvTimesMatAddMv( -1.0, *tempeveci, Breal, 1.0, *tempKevec );
          MVT::MvNorm( *tempKevec, resnorm );

          // Compute the norms and scale by magnitude of eigenvalue
          normR[i] = lapack.LAPY2( tempnrm[0], resnorm[0] ) /
            lapack.LAPY2( evals[i].realpart, evals[i].imagpart );
          normR[i+1] = normR[i];

          i=i+2;
        }
      }
    }

    // Output computed eigenvalues and their direct residuals
    if (verbose && MyPID==0) {
      std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);	
      std::cout<<std::endl<< "Actual Residuals"<<std::endl;
      if (MyProblem->isHermitian()) {
        std::cout<< std::setw(16) << "Real Part"
          << std::setw(20) << "Direct Residual"<< std::endl;
        std::cout<<"-----------------------------------------------------------"<<std::endl;
        for (int i=0; i<numev; i++) {
          std::cout<< std::setw(16) << evals[i].realpart 
            << std::setw(20) << normR[i] << std::endl;
        }  
        std::cout<<"-----------------------------------------------------------"<<std::endl;
      } 
      else {
        std::cout<< std::setw(16) << "Real Part"
          << std::setw(16) << "Imag Part"
          << std::setw(20) << "Direct Residual"<< std::endl;
        std::cout<<"-----------------------------------------------------------"<<std::endl;
        for (int i=0; i<numev; i++) {
          std::cout<< std::setw(16) << evals[i].realpart 
            << std::setw(16) << evals[i].imagpart 
            << std::setw(20) << normR[i] << std::endl;
        }  
        std::cout<<"-----------------------------------------------------------"<<std::endl;
      }  
    }
  }

#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return 0;

} // end BlockKrylovSchurEpetraExFile.cpp
