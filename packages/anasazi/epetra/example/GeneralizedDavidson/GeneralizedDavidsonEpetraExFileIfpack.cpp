//  This example compute the eigenvalues of a Harwell-Boeing matrix using the block
//  Generalized Davidson method.  The matrix is passed to the example routine through the command line,
//  converted to an Epetra matrix through some utilty routines and passed to the
//  eigensolver.  An Ifpack ILUT factorization of K is used as preconditioner.
#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_CrsMatrix.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"

#include "Ifpack.h"
#include "Ifpack_Preconditioner.h"
#include "Epetra_InvOperator.h"

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
  bool haveM = true;

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
  std::string k_filename = "bfw782a.mtx";
  std::string m_filename = "bfw782b.mtx";
  std::string which = "LR";
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
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
  // Build Preconditioner
  //
  Ifpack factory;
  std::string ifpack_type = "ILUT";
  int overlap = 0;
  Teuchos::RCP<Ifpack_Preconditioner> ifpack_prec = Teuchos::rcp(
          factory.Create( ifpack_type, K.get(), overlap ) );

  //
  // Set parameters and compute preconditioner
  //
  Teuchos::ParameterList ifpack_params;
  double droptol = 1e-2;
  double fill = 2.0;
  ifpack_params.set("fact: drop tolerance",droptol);
  ifpack_params.set("fact: ilut level-of-fill",fill);
  ifpack_prec->SetParameters(ifpack_params);
  ifpack_prec->Initialize();
  ifpack_prec->Compute();

  //
  // GeneralizedDavidson expects preconditioner to be applied with
  //  "Apply" rather than "Apply_Inverse"
  //
  Teuchos::RCP<Epetra_Operator> prec = Teuchos::rcp(
          new Epetra_InvOperator(ifpack_prec.get()) );

  //
  //************************************
  // Start the block Davidson iteration
  //***********************************
  //
  //  Variables used for the Block Arnoldi Method
  //
  int nev = 5;
  int blockSize = 5;
  int maxDim = 40;
  int maxRestarts = 10;
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
  MyPL.set( "Maximum Subspace Dimension", maxDim );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Initial Guess", "User" );

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
    MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>() );
    MyProblem->setA(K);
    MyProblem->setM(M);
    MyProblem->setPrec(prec);
    MyProblem->setInitVec(ivec);
  }
  else {
    MyProblem = Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>() );
    MyProblem->setA(K);
    MyProblem->setPrec(prec);
    MyProblem->setInitVec(ivec);
  }

  // Inform the eigenproblem that the (K,M) is Hermitian
  MyProblem->setHermitian( isHermitian );

  // Set the number of eigenvalues requested
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
#ifdef HAVE_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Initialize the Block Arnoldi solver
  Anasazi::GeneralizedDavidsonSolMgr<double, MV, OP> MySolverMgr(MyProblem, MyPL);

  // Solve the problem to the specified tolerances or length
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && MyPID==0 && verbose) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
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
        normR[i] = resnorm[0];
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
        normR[i] = lapack.LAPY2( tempnrm[i], resnorm[i] );
        normR[i+1] = normR[i];

        i=i+2;
      }
    }

    // Output computed eigenvalues and their direct residuals
    if (verbose && MyPID==0) {
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<endl<< "Actual Residuals"<<endl;
      cout<< std::setw(16) << "Real Part"
        << std::setw(16) << "Imag Part"
        << std::setw(20) << "Direct Residual"<< endl;
      cout<<"-----------------------------------------------------------"<<endl;
      for (int i=0; i<numev; i++) {
        cout<< std::setw(16) << evals[i].realpart
          << std::setw(16) << evals[i].imagpart
          << std::setw(20) << normR[i] << endl;
      }
      cout<<"-----------------------------------------------------------"<<endl;
    }
  }

#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return 0;

} // end BlockKrylovSchurEpetraExFile.cpp
