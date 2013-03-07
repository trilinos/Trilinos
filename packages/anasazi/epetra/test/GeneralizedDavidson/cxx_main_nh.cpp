// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
//  This test computes the specified eigenvalues of the discretized 2D Convection-Diffusion
//  equation using the generalized Davidson method.  This discretized operator is constructed as an
//  Epetra matrix.  No preconditioner is applied to the problem.

#include "AnasaziGeneralizedDavidsonSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicSort.hpp"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  bool testFailed;
  bool boolret;
  int MyPID = Comm.MyPID();

  bool verbose = true;
  bool debug = false;
  std::string which("SM");

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,LR,SI,or LI).");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

  typedef double ScalarType;
  typedef Teuchos::ScalarTraits<ScalarType>          ScalarTypeTraits;
  typedef ScalarTypeTraits::magnitudeType            MagnitudeType;
  typedef Epetra_MultiVector                         MV;
  typedef Epetra_Operator                            OP;
  typedef Anasazi::MultiVecTraits<ScalarType,MV>     MVTraits;
  typedef Anasazi::OperatorTraits<ScalarType,MV,OP>  OpTraits;

  //  Dimension of the matrix
  int nx = 10;        // Discretization points in any one direction.
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

  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Copy, Map, &NumNz[0]) );

  // Diffusion coefficient, can be set by user.
  // When rho*h/2 <= 1, the discrete convection-diffusion operator has real eigenvalues.
  // When rho*h/2 > 1, the operator has complex eigenvalues.
  double rho = 2*(nx+1);

  // Compute coefficients for discrete convection-diffution operator
  const double one = 1.0;
  std::vector<double> Values(4);
  std::vector<int> Indices(4);
  double h = one /(nx+1);
  double h2 = h*h;
  double c = 5.0e-01*rho/ h;
  Values[0] = -one/h2 - c; Values[1] = -one/h2 + c; Values[2] = -one/h2; Values[3]= -one/h2;
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
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == nx*(nx-1))
    {
      Indices[0] = nx*(nx-1)+1;
      Indices[1] = nx*(nx-2);
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == nx-1)
    {
      Indices[0] = nx-2;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
      Indices[0] = 2*nx-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] == NumGlobalElements-1)
    {
      Indices[0] = NumGlobalElements-2;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
      Indices[0] = nx*(nx-1)-1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] < nx)
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i] > nx*(nx-1))
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else if (MyGlobalElements[i]%nx == 0)
    {
      Indices[0] = MyGlobalElements[i]+1;
      Indices[1] = MyGlobalElements[i]-nx;
      Indices[2] = MyGlobalElements[i]+nx;
      NumEntries = 3;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[1], &Indices[0]);
      assert( info==0 );
    }
    else if ((MyGlobalElements[i]+1)%nx == 0)
    {
      Indices[0] = MyGlobalElements[i]-nx;
      Indices[1] = MyGlobalElements[i]+nx;
      NumEntries = 2;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[2], &Indices[0]);
      assert( info==0 );
      Indices[0] = MyGlobalElements[i]-1;
      NumEntries = 1;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    else
    {
      Indices[0] = MyGlobalElements[i]-1;
      Indices[1] = MyGlobalElements[i]+1;
      Indices[2] = MyGlobalElements[i]-nx;
      Indices[3] = MyGlobalElements[i]+nx;
      NumEntries = 4;
      info = A->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0], &Indices[0]);
      assert( info==0 );
    }
    // Put in the diagonal entry
    info = A->InsertGlobalValues(MyGlobalElements[i], 1, &diag, &MyGlobalElements[i]);
    assert( info==0 );
  }

  // Finish up
  info = A->FillComplete();
  assert( info==0 );
  A->SetTracebackMode(1); // Shutdown Epetra Warning tracebacks

  //************************************
  // Start the block Davidson iteration
  //***********************************
  //
  //  Variables used for the Generalized Davidson Method
  //
  int nev = 4;
  int blockSize = 1;
  int maxDim = 50;
  int restartDim = 10;
  int maxRestarts = 500;
  double tol = 1e-10;

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
  MyPL.set( "Which", which );
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Maximum Subspace Dimension", maxDim);
  MyPL.set( "Restart Dimension", restartDim);
  MyPL.set( "Maximum Restarts", maxRestarts );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "Relative Convergence Tolerance", true );
  MyPL.set( "Initial Guess", "User" );

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();

  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem = Teuchos::rcp(
    new Anasazi::BasicEigenproblem<double,MV,OP>() );
  MyProblem->setA(A);
  MyProblem->setInitVec(ivec);

  // Inform the eigenproblem that the operator A is non-Hermitian
  MyProblem->setHermitian(false);

  // Set the number of eigenvalues requested
  MyProblem->setNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information
  boolret = MyProblem->setProblem();
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
  testFailed = false;
  if (returnCode != Anasazi::Converged && MyPID==0 && verbose) {
    testFailed = true;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<ScalarType,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<ScalarType> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  std::vector<int> index = sol.index;
  int numev = sol.numVecs;

  // Output computed eigenvalues and their direct residuals
  if (verbose && MyPID==0) {
    int numritz = (int)evals.size();
    cout.setf(std::ios_base::right, std::ios_base::adjustfield);
    cout<<endl<< "Computed Ritz Values"<< endl;
    cout<< std::setw(16) << "Real Part"
        << std::setw(16) << "Imag Part"
        << endl;
    cout<<"-----------------------------------------------------------"<<endl;
    for (int i=0; i<numritz; i++) {
      cout<< std::setw(16) << evals[i].realpart
          << std::setw(16) << evals[i].imagpart
          << endl;
    }
    cout<<"-----------------------------------------------------------"<<endl;
  }

  if (numev > 0) {
    // Compute residuals.
    Teuchos::LAPACK<int,double> lapack;
    std::vector<double> normA(numev);

    // The problem is non-Hermitian.
    int i=0;
    std::vector<int> curind(1);
    std::vector<double> resnorm(1), tempnrm(1);
    Teuchos::RCP<MV> tempAevec;
    Teuchos::RCP<const MV> evecr, eveci;
    Epetra_MultiVector Aevec(Map,numev);

    // Compute A*evecs
    OpTraits::Apply( *A, *evecs, Aevec );

    Teuchos::SerialDenseMatrix<int,double> Breal(1,1), Bimag(1,1);
    while (i<numev) {
      if (index[i]==0) {
        // Get a view of the current eigenvector (evecr)
        curind[0] = i;
        evecr = MVTraits::CloneView( *evecs, curind );

        // Get a copy of A*evecr
        tempAevec = MVTraits::CloneCopy( Aevec, curind );

        // Compute A*evecr - lambda*evecr
        Breal(0,0) = evals[i].realpart;
        MVTraits::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );

        // Compute the norm of the residual and increment counter
        MVTraits::MvNorm( *tempAevec, resnorm );
        normA[i] = resnorm[0] / Teuchos::ScalarTraits<MagnitudeType>::magnitude( evals[i].realpart );
        i++;
      } else {
        // Get a view of the real part of the eigenvector (evecr)
        curind[0] = i;
        evecr = MVTraits::CloneView( *evecs, curind );

        // Get a copy of A*evecr
        tempAevec = MVTraits::CloneCopy( Aevec, curind );

        // Get a view of the imaginary part of the eigenvector (eveci)
        curind[0] = i+1;
        eveci = MVTraits::CloneView( *evecs, curind );

        // Set the eigenvalue into Breal and Bimag
        Breal(0,0) = evals[i].realpart;
        Bimag(0,0) = evals[i].imagpart;

        // Compute A*evecr - evecr*lambdar + eveci*lambdai
        MVTraits::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
        MVTraits::MvTimesMatAddMv( 1.0, *eveci, Bimag, 1.0, *tempAevec );
        MVTraits::MvNorm( *tempAevec, tempnrm );

        // Get a copy of A*eveci
        tempAevec = MVTraits::CloneCopy( Aevec, curind );

        // Compute A*eveci - eveci*lambdar - evecr*lambdai
        MVTraits::MvTimesMatAddMv( -1.0, *evecr, Bimag, 1.0, *tempAevec );
        MVTraits::MvTimesMatAddMv( -1.0, *eveci, Breal, 1.0, *tempAevec );
        MVTraits::MvNorm( *tempAevec, resnorm );

        // Compute the norms and scale by magnitude of eigenvalue
        normA[i] = lapack.LAPY2( tempnrm[0], resnorm[0] ) /
          lapack.LAPY2( evals[i].realpart, evals[i].imagpart );
        normA[i+1] = normA[i];

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
      for (int j=0; j<numev; j++) {
        cout<< std::setw(16) << evals[j].realpart
            << std::setw(16) << evals[j].imagpart
            << std::setw(20) << normA[j] << endl;
        if ( normA[j] > tol ) {
          testFailed = true;
        }
      }
      cout<<"-----------------------------------------------------------"<<endl;
    }
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

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
