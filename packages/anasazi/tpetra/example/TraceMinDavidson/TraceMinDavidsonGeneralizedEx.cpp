// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//  This example demonstrates how to use TraceMin-Davidson to solve
//  a generalized eigenvalue problem

// Include autoconfigured header
#include "AnasaziConfigDefs.hpp"

// Include header for TraceMin-Davidson solver
#include "AnasaziTraceMinDavidsonSolMgr.hpp"

// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"

// Include header to provide Anasazi with Tpetra adapters
#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziOperator.hpp"

// Include header for Tpetra compressed-row storage matrix
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"

// Include headers for reading and writing matrix-market files
#include <MatrixMarket_Tpetra.hpp>

// Include header for sparse matrix operations
//#include <TpetraExt_MatrixMatrix_def.hpp>

// Include header for Teuchos serial dense matrix
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Teuchos_ArrayViewDecl.hpp"

#include "Teuchos_ParameterList.hpp"


int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;

  //
  // Specify types used in this example
  //
  typedef double                                       Scalar;
  typedef Tpetra::CrsMatrix<Scalar>                    CrsMatrix;
  typedef Tpetra::MultiVector<Scalar>                  MV;
  typedef Tpetra::Operator<Scalar>                     OP;
  typedef Tpetra::MatrixMarket::Reader<CrsMatrix>      Reader;
  typedef Anasazi::MultiVecTraits<Scalar, MV>          MVT;
  typedef Anasazi::OperatorTraits<Scalar, MV, OP>      OPT;

  //
  // Initialize the MPI session
  //
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  //
  // Get the default communicator
  //
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();
  const int myRank = comm->getRank ();

  //
  // Get parameters from command-line processor
  //
  // FIMME (mfh 12 Feb 2015) The defaults shouldn't point to a
  // specific path.  I don't think the test uses the defaults, though.
  // I don't want to change this because it might break the author's
  // workflow.
  std::string filenameA ("/home/amklinv/matrices/bcsstk06.mtx");
  std::string filenameB ("/home/amklinv/matrices/bcsstm06.mtx");
  Scalar tol = 1e-6;
  int nev = 4;
  int blockSize = 1;
  bool verbose = true;
  std::string whenToShift = "Always";
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("fileA",&filenameA, "Filename for the Matrix-Market stiffness matrix.");
  cmdp.setOption("fileB",&filenameB, "Filename for the Matrix-Market mass matrix.");
  cmdp.setOption("tolerance",&tol, "Relative residual used for solver.");
  cmdp.setOption("nev",&nev, "Number of desired eigenpairs.");
  cmdp.setOption("blocksize",&blockSize, "Number of vectors to add to the subspace at each iteration.");
  cmdp.setOption("verbose","quiet",&verbose, "Whether to print a lot of info or a little bit.");
  cmdp.setOption("whenToShift",&whenToShift, "When to perform Ritz shifts. Options: Never, After Trace Levels, Always.");
  if(cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  //
  // Read the matrices from a file
  //
  RCP<const CrsMatrix> K = Reader::readSparseFile(filenameA, comm);
  RCP<const CrsMatrix> M = Reader::readSparseFile(filenameB, comm);

  //
  // Compute the norm of the matrix
  //
  Scalar mat_norm = std::max(K->getFrobeniusNorm(),M->getFrobeniusNorm());

  //
  // ************************************
  // Start the block Arnoldi iteration
  // ************************************
  //
  //  Variables used for the Block Arnoldi Method
  //
  int verbosity;
  int numRestartBlocks = 2*nev/blockSize;
  int numBlocks = 10*nev/blockSize;
  if(verbose)
    verbosity = Anasazi::TimingDetails + Anasazi::IterationDetails + Anasazi::Debug + Anasazi::FinalSummary;
  else
    verbosity = Anasazi::TimingDetails;
  //
  // Create parameter list to pass into solver
  //
  Teuchos::ParameterList MyPL;
  MyPL.set( "Verbosity", verbosity );                  // How much information should the solver print?
  MyPL.set( "Saddle Solver Type", "Projected Krylov"); // Use projected minres/gmres to solve the saddle point problem
  MyPL.set( "Block Size", blockSize );                 // Add blockSize vectors to the basis per iteration
  MyPL.set( "Convergence Tolerance", tol*mat_norm );   // How small do the residuals have to be
  MyPL.set( "Relative Convergence Tolerance", false);  // Don't scale residuals by eigenvalues (when checking for convergence)
  MyPL.set( "Use Locking", true);                      // Use deflation
  MyPL.set( "Relative Locking Tolerance", false);      // Don't scale residuals by eigenvalues (when checking whether to lock a vector)
  MyPL.set("Num Restart Blocks", numRestartBlocks);    // When we restart, we start back up with 2*nev blocks
  MyPL.set("Num Blocks", numBlocks);                   // Maximum number of blocks in the subspace
  MyPL.set("When To Shift", whenToShift);

  //
  // Create an Epetra_MultiVector for an initial vector to start the solver.
  // Note:  This needs to have the same number of columns as the blocksize.
  //
  RCP<MV> ivec = rcp (new MV (K->getRowMap (), blockSize));
  MVT::MvRandom (*ivec);

  //
  // Create the eigenproblem
  //
  RCP<Anasazi::BasicEigenproblem<Scalar,MV,OP> > MyProblem =
    rcp (new Anasazi::BasicEigenproblem<Scalar,MV,OP> (K, M, ivec));

  //
  // Inform the eigenproblem that the matrix pencil (K,M) is symmetric
  //
  MyProblem->setHermitian(true);

  //
  // Set the number of eigenvalues requested
  //
  MyProblem->setNEV( nev );

  //
  // Inform the eigenproblem that you are finished passing it information
  //
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (myRank == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << std::endl;
    }
    return -1;
  }

  //
  // Initialize the TraceMin-Davidson solver
  //
  Anasazi::Experimental::TraceMinDavidsonSolMgr<Scalar, MV, OP> MySolverMgr(MyProblem, MyPL);

  //
  // Solve the problem to the specified tolerances
  //
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  if (returnCode != Anasazi::Converged && myRank == 0) {
    cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << std::endl;
  }
  else if (myRank == 0)
    cout << "Anasazi::EigensolverMgr::solve() returned converged." << std::endl;

  //
  // Get the eigenvalues and eigenvectors from the eigenproblem
  //
  Anasazi::Eigensolution<Scalar,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<Scalar> > evals = sol.Evals;
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  //
  // Compute the residual, just as a precaution
  //
  if (numev > 0) {

    Teuchos::SerialDenseMatrix<int,Scalar> T(numev,numev);
    for(int i=0; i < numev; i++)
      T(i,i) = evals[i].realpart;
    std::vector<Scalar> normR(sol.numVecs);
    MV Kvec( K->getRowMap(), MVT::GetNumberVecs( *evecs ) );
    MV Mvec( M->getRowMap(), MVT::GetNumberVecs( *evecs ) );

    OPT::Apply( *K, *evecs, Kvec );
    OPT::Apply( *M, *evecs, Mvec );
    MVT::MvTimesMatAddMv( -1.0, Mvec, T, 1.0, Kvec );
    MVT::MvNorm( Kvec, normR );

    if (myRank == 0) {
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<"Actual Eigenvalues: "<<std::endl;
      cout<<"------------------------------------------------------"<<std::endl;
      cout<<std::setw(16)<<"Real Part"
        <<std::setw(16)<<"Error"<<std::endl;
      cout<<"------------------------------------------------------"<<std::endl;
      for (int i=0; i<numev; i++) {
        cout<<std::setw(16)<<evals[i].realpart
          <<std::setw(16)<<normR[i]/mat_norm
          <<std::endl;
      }
      cout<<"------------------------------------------------------"<<std::endl;
    }
  }

  return 0;
}
