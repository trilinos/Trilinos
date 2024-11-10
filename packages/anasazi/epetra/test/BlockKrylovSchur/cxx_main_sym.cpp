// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//  This example computes the specified eigenvalues of the discretized 2D Diffusion
//  equation using the block Krylov-Schur method.  This discretized operator is constructed as an
//  Epetra matrix, then passed into the Anasazi::EpetraOp to be used in the construction of the
//  Krylov decomposition.  The specifics of the block Krylov-Schur method can be set by the user.
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicSort.hpp"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Assert.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {

  using std::cout;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool boolret;
  int MyPID = Comm.MyPID();

  bool verbose = true;
  bool debug = false;
  bool dynXtraNev = false;
  std::string which("SM");
  int nx = 10;        // Discretization points in any one direction.
  int nev = 10;
  int blockSize = 3;
  int numBlocks = 10;
  double tol = 1e-5;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,LR,SI,or LI).");
  cmdp.setOption("nx",&nx,"Number of discretization points in each direction; n=nx*nx.");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("blocksize",&blockSize,"Block Size.");
  cmdp.setOption("tol",&tol,"Solver Tolerance.");
  cmdp.setOption("numblocks",&numBlocks,"Number of blocks for the Krylov-Schur form.");
  cmdp.setOption("dynrestart","nodynrestart",&dynXtraNev,"Use dynamic restart boundary to accelerate convergence.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  typedef double ScalarType;
  typedef Teuchos::ScalarTraits<ScalarType>          SCT;
  typedef SCT::magnitudeType               MagnitudeType;
  typedef Epetra_MultiVector                          MV;
  typedef Epetra_Operator                             OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OPT;

  //  Dimension of the matrix
  int NumGlobalElements = nx*nx;  // Size of matrix nx*nx

  // Construct a Map that puts approximately the same number of
  // equations on each processor.

  Epetra_Map Map(NumGlobalElements, 0, Comm);

  // Get update list and number of local equations from newly created Map.

  int NumMyElements = Map.NumMyElements();

  std::vector<int> MyGlobalElements(NumMyElements);
  Map.MyGlobalElements(&MyGlobalElements[0]);

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation
  // on this processor
  std::vector<int> NumNz(NumMyElements);

  /* We are building a matrix of block structure:
  
      | T -I          |
      |-I  T -I       |
      |   -I  T       |
      |        ...  -I|
      |           -I T|

   where each block is dimension nx by nx and the matrix is on the order of
   nx*nx.  The block T is a tridiagonal matrix. 
  */

  for (int i=0; i<NumMyElements; i++) {
    if (MyGlobalElements[i] == 0 || MyGlobalElements[i] == NumGlobalElements-1 || 
        MyGlobalElements[i] == nx-1 || MyGlobalElements[i] == nx*(nx-1) ) {
      NumNz[i] = 3;
    }
    else if (MyGlobalElements[i] < nx || MyGlobalElements[i] > nx*(nx-1) || 
             MyGlobalElements[i]%nx == 0 || (MyGlobalElements[i]+1)%nx == 0) {
      NumNz[i] = 4;
    }
    else {
      NumNz[i] = 5;
    }
  }

  // Create an Epetra_Matrix

  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, Map, &NumNz[0]) );

  // Compute coefficients for discrete convection-diffution operator
  const double one = 1.0;
  std::vector<double> Values(4);
  std::vector<int> Indices(4);
  double h = one /(nx+1);
  double h2 = h*h;
  Values[0] = -one/h2; Values[1] = -one/h2; Values[2] = -one/h2; Values[3]= -one/h2;
  double diag = 4.0 / h2;
  int NumEntries, info;

  for (int i=0; i<NumMyElements; i++)
  {
    if (MyGlobalElements[i]==0)
    {
      Indices[0] = 1;
      Indices[1] = nx;
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] == nx*(nx-1))
    {
      Indices[0] = nx*(nx-1)+1;
      Indices[1] = nx*(nx-2);
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] == nx-1)
    {
      Indices[0] = nx-2;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
      Indices[0] = 2*nx-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
    {
      Indices[0] = NumGlobalElements-2;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
      Indices[0] = nx*(nx-1)-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] < nx)
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i] > nx*(nx-1))
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if (MyGlobalElements[i]%nx == 0)
    {
      Indices[0] = MyGlobalElements[i]+1;
      Indices[1] = MyGlobalElements[i]-nx;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else if ((MyGlobalElements[i]+1)%nx == 0)
    {
      Indices[0] = MyGlobalElements[i]-nx;
      Indices[1] = MyGlobalElements[i]+nx;
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
      Indices[0] = MyGlobalElements[i]-1;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    else
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      Indices[3] = MyGlobalElements[i]+nx;
      NumEntries = 4;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      TEUCHOS_ASSERT( info==0 );
    }
    // Put in the diagonal entry
    info = A->InsertGlobalValues(MyGlobalElements[i], 1, &diag, &MyGlobalElements[i]);
    TEUCHOS_ASSERT( info==0 );
  }

  // Finish up
  info = A->FillComplete();
  TEUCHOS_ASSERT( info==0 );
  A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

  //************************************
  // Start the block Arnoldi iteration
  //***********************************
  //
  //  Variables used for the Block Krylov Schur Method
  //    
  int maxRestarts = 500;

  // Create a sort manager to pass into the block Krylov-Schur solver manager
  // -->  Make sure the reference-counted pointer is of type Anasazi::SortManager<>
  // -->  The block Krylov-Schur solver manager uses Anasazi::BasicSort<> by default,
  //      so you can also pass in the parameter "Which", instead of a sort manager.
  Teuchos::RCP<Anasazi::SortManager<MagnitudeType> > MySort =     
    Teuchos::rcp( new Anasazi::BasicSort<MagnitudeType>( which ) );

  // Set verbosity level
  int verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbose) {
    verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  }
  if (debug) {
    verbosity += Anasazi::Debug;
  }
  //
  // Create parameter list to pass into solver manager
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );
  MyPL.set( "Sort Manager", MySort );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set("Dynamic Extra NEV",dynXtraNev);

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();

  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );
  
  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->setHermitian(true); 
  
  // Set the number of eigenvalues requested
  MyProblem->setNEV( nev );
  
  // Inform the eigenproblem that you are finishing passing it information
  boolret = MyProblem->setProblem();
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
  bool testFailed = false;
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged) {
    testFailed = true;
    if (MyPID==0 && verbose)
      std::cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
  }

  // Get the Ritz values from the eigensolver
  std::vector<Anasazi::Value<double> > ritzValues = MySolverMgr.getRitzValues();
  
  // Output computed eigenvalues and their direct residuals
  if (verbose && MyPID==0) {
    int numritz = (int)ritzValues.size();
    std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
    std::cout<<std::endl<< "Computed Ritz Values"<< std::endl;
    std::cout<< std::setw(16) << "Real Part"
      << std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;
    for (int i=0; i<numritz; i++) {
      std::cout<< std::setw(16) << ritzValues[i].realpart 
        << std::endl;
    }  
    std::cout<<"-----------------------------------------------------------"<<std::endl;
  } 

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  std::vector<int> index = sol.index;
  int numev = sol.numVecs;
  
  if (numev > 0) {
    // Compute residuals.
    std::vector<double> normA(numev);
    
    // Get storage
    Epetra_MultiVector Aevecs(Map,numev);
    Teuchos::SerialDenseMatrix<int,double> B(numev,numev);
    B.putScalar(0.0); 
    for (int i=0; i<numev; i++) {B(i,i) = evals[i].realpart;}
      
    // Compute A*evecs
    OPT::Apply( *A, *evecs, Aevecs );
      
    // Compute A*evecs - lambda*evecs and its norm
    MVT::MvTimesMatAddMv( -1.0, *evecs, B, 1.0, Aevecs );
    MVT::MvNorm( Aevecs, normA );
      
    // Scale the norms by the eigenvalue
    for (int i=0; i<numev; i++) {
      normA[i] /= Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart );
      if ( normA[i] > tol ) {
        testFailed = true;
      }
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
          << std::setw(20) << normA[i] << std::endl;
      }  
      std::cout<<"-----------------------------------------------------------"<<std::endl;
    }
  }  

#ifdef EPETRA_MPI
    MPI_Finalize();
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
