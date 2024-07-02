// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \example LOBPCGEpetraFile.cpp
/// \brief Use LOBPCG with Epetra test problem loaded from file.
///
/// This example computes the eigenvalues of largest magnitude of an
/// eigenvalue problem $A x = \lambda x$, using Anasazi's
/// implementation of the LOBPCG method, with Epetra linear algebra.
/// The example loads the matrix from a file whose name is specified
/// at the command line.

#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_InvOperator.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"

// Include header for Ifpack incomplete Cholesky preconditioner
#include "Ifpack.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
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

  //************************************
  // Get the parameters from the command line
  //************************************
  //
  int nev = 5;
  int blockSize = 5;
  int maxIterations = 1000;
  double tol = 1.0e-8;
  bool verbose = true;
  bool locking=false, fullOrtho=true;
  std::string k_filename = "";
  std::string m_filename = "";
  std::string which = "SM";
  bool usePrec = true;
  double prec_dropTol = 1e-4;
  int prec_lofill = 0;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("blocksize",&blockSize,"Block size used in LOBPCG.");
  cmdp.setOption("maxiters",&maxIterations,"Maximum number of iterations used in LOBPCG.");
  cmdp.setOption("tol",&tol,"Convergence tolerance requested for computed eigenvalues.");
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("locking","nolocking",&locking,"Use locking of converged eigenvalues.");
  cmdp.setOption("fullortho","nofullortho",&fullOrtho,"Use full orthogonalization.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,or LR).");
  cmdp.setOption("usePrec","noPrec",&usePrec,"Use Ifpack for preconditioning.");
  cmdp.setOption("prec_dropTol",&prec_dropTol,"Preconditioner: drop tolerance.");
  cmdp.setOption("prec_lofill",&prec_lofill,"Preconditioner: level of fill.");
  cmdp.setOption("K-filename",&k_filename,"Filename and path of the stiffness matrix.");
  cmdp.setOption("M-filename",&m_filename,"Filename and path of the mass matrix.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  if (k_filename=="") {
    std::cout << "The matrix K must be supplied through an input file!!!" << std::endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if (m_filename!="") {
    haveM = true;
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


  //************************************
  // Select the Preconditioner
  //************************************
  //
  Teuchos::RCP<Ifpack_Preconditioner> prec;
  Teuchos::RCP<Epetra_Operator> PrecOp;
  if (usePrec) {
    Ifpack precFactory;
    // additive-Schwartz incomplete Cholesky with thresholding; see IFPACK documentation
    std::string precType = "IC stand-alone";
    int overlapLevel = 0;
    prec = Teuchos::rcp( precFactory.Create(precType,K.get(),overlapLevel) );
    // parameters for preconditioner
    Teuchos::ParameterList precParams;
    precParams.set("fact: drop tolerance",prec_dropTol);
    precParams.set("fact: level-of-fill",prec_lofill);
    IFPACK_CHK_ERR(prec->SetParameters(precParams));
    IFPACK_CHK_ERR(prec->Initialize());
    IFPACK_CHK_ERR(prec->Compute());
    //
    // encapsulate this preconditioner into a IFPACKPrecOp class
    PrecOp = Teuchos::rcp( new Epetra_InvOperator(&*prec) );
  }

  //
  //************************************
  // Start the block Davidson iteration
  //***********************************
  //
  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary;
  }
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Iterations", maxIterations );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Full Ortho", true );
  MyPL.set( "Use Locking", true );

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
    MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>( K, M, ivec ) );
  }
  else {
    MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>( K, ivec ) );
  }

  // Inform the eigenproblem that (K,M) is Hermitian
  MyProblem->setHermitian(true);

  // Pass the preconditioner to the eigenproblem
  if (usePrec) {
    MyProblem->setPrec(PrecOp);
  }

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

  // Initialize the LOBPCG solver
  Anasazi::LOBPCGSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

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

    // Get storage
    Teuchos::RCP<Epetra_MultiVector> Kevecs, Mevecs;
    Teuchos::SerialDenseMatrix<int,double> B(numev,numev);
    B.putScalar(0.0);
    for (int i=0; i<numev; i++) {B(i,i) = evals[i].realpart;}

    // Compute K*evecs
    Kevecs = Teuchos::rcp(new Epetra_MultiVector(*Map,numev) );
    OPT::Apply( *K, *evecs, *Kevecs );

    // Compute M*evecs
    if (haveM) {
      Mevecs = Teuchos::rcp(new Epetra_MultiVector(*Map,numev) );
      OPT::Apply( *M, *evecs, *Mevecs );
    }
    else {
      Mevecs = evecs;
    }

    // Compute K*evecs - lambda*M*evecs and its norm
    MVT::MvTimesMatAddMv( -1.0, *Mevecs, B, 1.0, *Kevecs );
    MVT::MvNorm( *Kevecs, normR );

    // Scale the norms by the eigenvalue
    for (int i=0; i<numev; i++) {
      normR[i] /= Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
    }

    // Output computed eigenvalues and their direct residuals
    if (verbose && MyPID==0) {
      std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      std::cout<<std::endl<< "Actual Residuals"<<std::endl;
      std::cout<< std::setw(16) << "Real Part"
        << std::setw(20) << "Direct Residual"<< std::endl;
      std::cout<<"-----------------------------------------------------------"<<std::endl;
      for (int i=0; i<numev; i++) {
        std::cout<< std::setw(16) << evals[i].realpart
          << std::setw(20) << normR[i] << std::endl;
      }
      std::cout<<"-----------------------------------------------------------"<<std::endl;
    }
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif
  return 0;

} // end BlockDavidsonEpetraExFile.cpp
