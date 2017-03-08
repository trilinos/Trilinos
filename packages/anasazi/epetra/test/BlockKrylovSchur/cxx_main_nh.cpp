// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
//  This test computes the specified eigenvalues of the discretized 2D Convection-Diffusion
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
  bool insitu = false;
  std::string which("SM");

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
  cmdp.setOption("insitu","exsitu",&insitu,"Perform in situ restarting.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,LR,SI,or LI).");
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

  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp( new Epetra_CrsMatrix(Epetra_DataAccess::Copy, Map, &NumNz[0]) );

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
  // Start the block Arnoldi iteration
  //***********************************
  //
  //  Variables used for the Block Krylov Schur Method
  //    
  int nev = 4;
  int blockSize = 1;
  int numBlocks = 20;
  int maxRestarts = 500;
  //int stepSize = 5;
  double tol = 1e-8;

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
  //MyPL.set( "Which", which );  
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Num Blocks", numBlocks );
  MyPL.set( "Maximum Restarts", maxRestarts );
  //MyPL.set( "Step Size", stepSize );
  MyPL.set( "Convergence Tolerance", tol );
  MyPL.set( "In Situ Restarting", insitu );

  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  Teuchos::RCP<Epetra_MultiVector> ivec = Teuchos::rcp( new Epetra_MultiVector(Map, blockSize) );
  ivec->Random();

  // Create the eigenproblem.
  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<double, MV, OP>(A, ivec) );

  // Inform the eigenproblem that the operator A is non-Hermitian
  MyProblem->setHermitian(false); 

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
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  testFailed = false;
  if (returnCode != Anasazi::Converged && MyPID==0 && verbose) {
    testFailed = true;
  }

  // Get the Ritz values from the eigensolver
  std::vector<Anasazi::Value<double> > ritzValues = MySolverMgr.getRitzValues();

  // Output computed eigenvalues and their direct residuals
  if (verbose && MyPID==0) {
    int numritz = (int)ritzValues.size();
    std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);        
    std::cout<<std::endl<< "Computed Ritz Values"<< std::endl;
    std::cout<< std::setw(16) << "Real Part"
        << std::setw(16) << "Imag Part"
        << std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;
    for (int i=0; i<numritz; i++) {
      std::cout<< std::setw(16) << ritzValues[i].realpart 
          << std::setw(16) << ritzValues[i].imagpart 
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
    OPT::Apply( *A, *evecs, Aevec );

    Teuchos::SerialDenseMatrix<int,double> Breal(1,1), Bimag(1,1);
    while (i<numev) {
      if (index[i]==0) {
        // Get a view of the current eigenvector (evecr)
        curind[0] = i;
        evecr = MVT::CloneView( *evecs, curind );

        // Get a copy of A*evecr
        tempAevec = MVT::CloneCopy( Aevec, curind );

        // Compute A*evecr - lambda*evecr
        Breal(0,0) = evals[i].realpart;
        MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );

        // Compute the norm of the residual and increment counter
        MVT::MvNorm( *tempAevec, resnorm );
        normA[i] = resnorm[0] / Teuchos::ScalarTraits<MagnitudeType>::magnitude( evals[i].realpart );
        i++;
      } else {
        // Get a view of the real part of the eigenvector (evecr)
        curind[0] = i;
        evecr = MVT::CloneView( *evecs, curind );

        // Get a copy of A*evecr
        tempAevec = MVT::CloneCopy( Aevec, curind );

        // Get a view of the imaginary part of the eigenvector (eveci)
        curind[0] = i+1;
        eveci = MVT::CloneView( *evecs, curind );

        // Set the eigenvalue into Breal and Bimag
        Breal(0,0) = evals[i].realpart;
        Bimag(0,0) = evals[i].imagpart;

        // Compute A*evecr - evecr*lambdar + eveci*lambdai
        MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
        MVT::MvTimesMatAddMv( 1.0, *eveci, Bimag, 1.0, *tempAevec );
        MVT::MvNorm( *tempAevec, tempnrm );

        // Get a copy of A*eveci
        tempAevec = MVT::CloneCopy( Aevec, curind );

        // Compute A*eveci - eveci*lambdar - evecr*lambdai
        MVT::MvTimesMatAddMv( -1.0, *evecr, Bimag, 1.0, *tempAevec );
        MVT::MvTimesMatAddMv( -1.0, *eveci, Breal, 1.0, *tempAevec );
        MVT::MvNorm( *tempAevec, resnorm );

        // Compute the norms and scale by magnitude of eigenvalue
        normA[i] = lapack.LAPY2( tempnrm[0], resnorm[0] ) /
          lapack.LAPY2( evals[i].realpart, evals[i].imagpart );
        normA[i+1] = normA[i];

        i=i+2;
      }
    }

    // Output computed eigenvalues and their direct residuals
    if (verbose && MyPID==0) {
      std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);        
      std::cout<<std::endl<< "Actual Residuals"<<std::endl;
      std::cout<< std::setw(16) << "Real Part"
          << std::setw(16) << "Imag Part"
          << std::setw(20) << "Direct Residual"<< std::endl;
      std::cout<<"-----------------------------------------------------------"<<std::endl;
      for (int j=0; j<numev; j++) {
        std::cout<< std::setw(16) << evals[j].realpart 
            << std::setw(16) << evals[j].imagpart 
            << std::setw(20) << normA[j] << std::endl;
        if ( normA[j] > tol ) {
          testFailed = true;
        }
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
