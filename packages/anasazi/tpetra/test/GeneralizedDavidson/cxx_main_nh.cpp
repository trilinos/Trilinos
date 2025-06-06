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
//  equation using the generalized Davidson method.  This discretized operator is constructed as a
//  Tpetra matrix.  No preconditioner is applied to the problem.

// Tpetra
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

// Teuchos
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// Anasazi
#include "AnasaziBasicSort.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziGeneralizedDavidsonSolMgr.hpp"


template <typename ScalarType>
int run (int argc, char *argv[]) {

  using ST  = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO  = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO  = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT  = typename Tpetra::MultiVector<>::node_type;

  using MT = typename Teuchos::ScalarTraits<ScalarType>::magnitudeType;

  using OP  = Tpetra::Operator<ST,LO,GO,NT>;
  using MV  = Tpetra::MultiVector<ST,LO,GO,NT>;
  using OPT = Anasazi::OperatorTraits<ST,MV,OP>;
  using MVT = Anasazi::MultiVecTraits<ST,MV>;

  using tmap_t = Tpetra::Map<LO,GO,NT>;
  using tcrsmatrix_t = Tpetra::CrsMatrix<ST,LO,GO,NT>;

  using starray_t = Teuchos::Array<ST>;
  using goarray_t = Teuchos::Array<GO>;

  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &std::cout);

  bool verbose = true;
  bool success = false;

  try {

    const auto comm = Tpetra::getDefaultComm();

    bool boolret;
    const int MyPID = comm->getRank();

    bool debug = false;
    std::string which("SM");

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Print debugging information.");
    cmdp.setOption("sort",&which,"Targetted eigenvalues (SM,LM,SR,LR,SI,or LI).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      throw -1;
    }

    //  Dimension of the matrix
    int nx = 10;        // Discretization points in any one direction.
    int NumGlobalElements = nx*nx;  // Size of matrix nx*nx

    // Construct a Map that puts approximately the same number of
    // equations on each processor.
    if (MyPID == 0) {
      std::cout << "Building Map" << std::endl;
    }

    RCP<tmap_t> Map = rcp(new tmap_t(NumGlobalElements, 0, comm));

    if (MyPID == 0) {
      std::cout << "Setting up info for filling matrix" << std::endl;
    }

    // Get update list and number of local equations from newly created Map.

    int NumMyElements = Map->getLocalNumElements();

    auto MyGlobalElements = Map->getMyGlobalIndices();

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

    if (MyPID == 0) {
      std::cout << "Creating matrix" << std::endl;
    }

    // Create an Tpetra::CrsMatrix
    RCP<tcrsmatrix_t> A = rcp( new tcrsmatrix_t(Map, NumGlobalElements / 2) );

    if (MyPID == 0) {
      std::cout << "Filling matrix" << std::endl;
    }

    // Diffusion coefficient, can be set by user.
    // When rho*h/2 <= 1, the discrete convection-diffusion operator has real eigenvalues.
    // When rho*h/2 > 1, the operator has complex eigenvalues.
    ST rho = 2*(nx+1);

    // Compute coefficients for discrete convection-diffution operator
    const ST one = 1.0;
    goarray_t Indices(4);
    ST h = one /(nx+1);
    ST h2 = h*h;
    ST c = 5.0e-01*rho/ h;
    ST val0 = -one/h2 - c;
    ST val1 = -one/h2 + c;
    ST val2 = -one/h2;
    ST val3 = -one/h2;
    starray_t Values0(4,val0), Values1(4,val1), Values2(4,val2), Values3(4,val3);
    ST diag = 4.0 / h2;
    // starray_t diag(NumMyElements, diag_val);
    LO NumEntries;

    for (GO i=0; i<NumMyElements; i++)
    {
      if (MyGlobalElements[i]==0)
      {
        Indices[0] = 1;
        Indices[1] = nx;
        NumEntries = 2;

        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values3.view(0,NumEntries));
      }
      else if (MyGlobalElements[i] == nx*(nx-1))
      {
        Indices[0] = nx*(nx-1)+1;
        Indices[1] = nx*(nx-2);
        NumEntries = 2;

        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values1.view(0,NumEntries));
      }
      else if (MyGlobalElements[i] == nx-1)
      {
        Indices[0] = nx-2;
        NumEntries = 1;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values0.view(0,NumEntries));
        Indices[0] = 2*nx-1;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values2.view(0,NumEntries));
      }
      else if (MyGlobalElements[i] == NumGlobalElements-1)
      {
        Indices[0] = NumGlobalElements-2;
        NumEntries = 1;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values0.view(0,NumEntries));
        Indices[0] = nx*(nx-1)-1;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values2.view(0,NumEntries));
      }
      else if (MyGlobalElements[i] < nx)
      {
        Indices[0] = MyGlobalElements[i]-1;
        Indices[1] = MyGlobalElements[i]+1;
        Indices[2] = MyGlobalElements[i]+nx;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values0.view(0,NumEntries));
      }
      else if (MyGlobalElements[i] > nx*(nx-1))
      {
        Indices[0] = MyGlobalElements[i]-1;
        Indices[1] = MyGlobalElements[i]+1;
        Indices[2] = MyGlobalElements[i]-nx;
        NumEntries = 3;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values0.view(0,NumEntries));
      }
      else if (MyGlobalElements[i]%nx == 0)
      {
        Indices[0] = MyGlobalElements[i]+1;
        Indices[1] = MyGlobalElements[i]-nx;
        Indices[2] = MyGlobalElements[i]+nx;
        NumEntries = 3;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values1.view(0,NumEntries));
      }
      else if ((MyGlobalElements[i]+1)%nx == 0)
      {
        Indices[0] = MyGlobalElements[i]-nx;
        Indices[1] = MyGlobalElements[i]+nx;
        NumEntries = 2;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values2.view(0,NumEntries));
        Indices[0] = MyGlobalElements[i]-1;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values0.view(0,NumEntries));
      }
      else
      {
        Indices[0] = MyGlobalElements[i]-1;
        Indices[1] = MyGlobalElements[i]+1;
        Indices[2] = MyGlobalElements[i]-nx;
        Indices[3] = MyGlobalElements[i]+nx;
        NumEntries = 4;
        A->insertGlobalValues(MyGlobalElements[i], Indices.view(0,NumEntries), Values0.view(0,NumEntries));
      }
      // Put in the diagonal entry
      NumEntries = 1;
      A->insertGlobalValues(MyGlobalElements[i], 1, &diag, &MyGlobalElements[i]);
    }

    if (MyPID == 0) {
      std::cout << "Calling fillComplete on matrix" << std::endl;
    }

    // Finish up
    A->fillComplete ();

    if (MyPID == 0) {
      std::cout << "Setting Anasazi parameters" << std::endl;
    }

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
    MT tol = 1e-08;

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

    if (MyPID == 0) {
      std::cout << "Creating initial vector for solver" << std::endl;
    }

    // Create a Tpetra::MultiVector for an initial vector to start the solver.
    // Note:  This needs to have the same number of columns as the blocksize.
    RCP<MV> ivec = rcp( new MV(Map, blockSize) );
    MVT::MvRandom(*ivec);

    if (MyPID == 0) {
      std::cout << "Creating eigenproblem" << std::endl;
    }

    // Create the eigenproblem.
    RCP<Anasazi::BasicEigenproblem<ST, MV, OP> > MyProblem = rcp(
        new Anasazi::BasicEigenproblem<ST,MV,OP>() );
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
        std::cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << std::endl;
      }
      return -1;
    }

    if (MyPID == 0) {
      std::cout << "Creating eigensolver (GeneralizedDavidsonSolMgr)" << std::endl;
    }

    // Initialize the Block Arnoldi solver
    Anasazi::GeneralizedDavidsonSolMgr<ST, MV, OP> MySolverMgr(MyProblem, MyPL);

    if (MyPID == 0) {
      std::cout << "Solving eigenproblem" << std::endl;
    }

    // Solve the problem to the specified tolerances or length
    Anasazi::ReturnType returnCode = Anasazi::Unconverged;
    returnCode = MySolverMgr.solve();
    success = (returnCode == Anasazi::Converged);

    std::cout << "Getting eigenvalues and eigenvectors from eigenproblem" << std::endl;

    // Get the eigenvalues and eigenvectors from the eigenproblem
    Anasazi::Eigensolution<ST,MV> sol = MyProblem->getSolution();
    std::vector<Anasazi::Value<ST> > evals = sol.Evals;
    RCP<MV> evecs = sol.Evecs;
    std::vector<int> index = sol.index;
    int numev = sol.numVecs;

    // Output computed eigenvalues and their direct residuals
    if (verbose && MyPID==0) {
      int numritz = (int)evals.size();
      std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      std::cout<<std::endl<< "Computed Ritz Values"<< std::endl;
      std::cout<< std::setw(16) << "Real Part"
        << std::setw(16) << "Imag Part"
        << std::endl;
      std::cout<<"-----------------------------------------------------------"<<std::endl;
      for (int i=0; i<numritz; i++) {
        std::cout<< std::setw(16) << evals[i].realpart
          << std::setw(16) << evals[i].imagpart
          << std::endl;
      }
      std::cout<<"-----------------------------------------------------------"<<std::endl;
    }

    if (numev > 0) {
      // Compute residuals.
      Teuchos::LAPACK<int,ST> lapack;
      std::vector<ST> normA(numev);

      // The problem is non-Hermitian.
      int i=0;
      std::vector<int> curind(1);
      std::vector<ST> resnorm(1), tempnrm(1);
      RCP<MV> tempAevec;
      RCP<const MV> evecr, eveci;
      MV Aevec(Map,numev);

      // Compute A*evecs
      OPT::Apply( *A, *evecs, Aevec );

      Teuchos::SerialDenseMatrix<int,ST> Breal(1,1), Bimag(1,1);
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
          normA[i] = resnorm[0] / Teuchos::ScalarTraits<MT>::magnitude( evals[i].realpart );
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
          success &= ( normA[j] < tol );
        }
        std::cout<<"-----------------------------------------------------------"<<std::endl;
      }
    }

    if (verbose && MyPID==0) {
      if (success)
        std::cout << "End Result: TEST PASSED" << std::endl;
      else
        std::cout << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

int main(int argc, char *argv[]) {
  return run<double>(argc,argv);
  // run<float>(argc,argv);
}
