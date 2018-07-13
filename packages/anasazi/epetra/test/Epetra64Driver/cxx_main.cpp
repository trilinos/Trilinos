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
//  This program computes the eigenvalues of a laplacian matrix using Anasazi.
//  NOTE:  The laplacian matrix generated has row and column ids larger than INT_MAX
//         to test Epetra64 functionality with Anasazi solvers.
//
//  Karen Devine, April 2012
//  Erik Boman, May 2012
//  K. Devine & H. Thornquist, January 2013
//

#define TINYMATRIX 26
#include "Epetra_ConfigDefs.h"
#ifdef EPETRA_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#define FINALIZE MPI_Finalize()
#else
#include "Epetra_SerialComm.h"
#define FINALIZE
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_InvOperator.h"

#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziRTRSolMgr.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziEpetraAdapter.hpp"

#include "Epetra_Time.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"

#include "build_simple_matrices.hpp"

#include <stdio.h>

////////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[])
{
  using std::cout;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&narg,&arg);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  bool verbose = true;
  int verbosity = 1;

  bool testEpetra64 = true;

  // Matrix properties
  bool isHermitian = true;

  // Multivector properties

  std::string initvec = "random";

  // Eigenvalue properties

  std::string which = "SR";
  std::string method = "LOBPCG";
  std::string precond = "none";
  std::string ortho = "SVQB";
  bool lock = true;
  bool relconvtol = false;
  bool rellocktol = false;
  int nev = 5;

  // Block-Arnoldi properties

  int blockSize = -1;
  int numblocks = -1;
  int maxrestarts = -1;
  int maxiterations = -1;
  int extrablocks = 0;
  int gensize = 25;  // Needs to be long long to test with > INT_MAX rows
  double tol = 1.0e-5;

  // Echo the command line
  if (MyPID == 0)  {
    for (int i = 0; i < narg; i++)
      cout << arg[i] << " ";
    cout << endl;
  }

  // Command-line processing

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("Epetra64", "no-Epetra64", &testEpetra64,
                 "Force code to use Epetra64, even if the problem size does "
                 "not require it. (Epetra64 will be used automatically for "
                 "sufficiently large problems, or not used if Epetra does not have built in support.)");
  cmdp.setOption("gen",&gensize,
                 "Generate a simple Laplacian matrix of size n.");
  cmdp.setOption("verbosity", &verbosity, "0=quiet, 1=low, 2=medium, 3=high.");
  cmdp.setOption("method",&method,
                 "Solver method to use:  LOBPCG, BD, BKS or IRTR.");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to find.");
  cmdp.setOption("which",&which,"Targetted eigenvalues (SM,LM,SR,or LR).");
  cmdp.setOption("tol",&tol,"Solver convergence tolerance.");
  cmdp.setOption("blocksize",&blockSize,"Block size to use in solver.");
  cmdp.setOption("numblocks",&numblocks,"Number of blocks to allocate.");
  cmdp.setOption("extrablocks",&extrablocks,
                 "Number of extra NEV blocks to allocate in BKS.");
  cmdp.setOption("maxrestarts",&maxrestarts,
                 "Maximum number of restarts in BKS or BD.");
  cmdp.setOption("maxiterations",&maxiterations,
                 "Maximum number of iterations in LOBPCG.");
  cmdp.setOption("lock","no-lock",&lock,
                 "Use Locking parameter (deflate for converged eigenvalues)");
  cmdp.setOption("initvec", &initvec,
                 "Initial vectors (random, unit, zero, random2)");
  cmdp.setOption("ortho", &ortho,
                 "Orthogonalization method (DGKS, SVQB, TSQR).");
  cmdp.setOption("relative-convergence-tol","no-relative-convergence-tol",
                 &relconvtol,
                 "Use Relative convergence tolerance "
                 "(normalized by eigenvalue)");
  cmdp.setOption("relative-lock-tol","no-relative-lock-tol",&rellocktol,
                 "Use Relative locking tolerance (normalized by eigenvalue)");
  if (cmdp.parse(narg,arg)!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    FINALIZE;
    return -1;
  }

  // Print the most essential options (not in the MyPL parameters later)
  verbose = (verbosity>0);
  if (verbose && MyPID==0){
    cout << "verbosity = " << verbosity << endl;
    cout << "method = " << method << endl;
    cout << "initvec = " << initvec << endl;
    cout << "nev = " << nev << endl;
  }

  // We need blockSize to be set so we can allocate memory with it.
  // If it wasn't set on the command line, set it to the Anasazi defaults here.
  // Defaults are those given in the documentation.
  if (blockSize < 0) {
    if (method == "BKS")
      blockSize = 1;
    else // other methods: LOBPCG, BD, IRTR
      blockSize = nev;
  }

  // Make sure Epetra was built with 64-bit global indices enabled.
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  if (testEpetra64)
    testEpetra64 = false;
#endif

  Epetra_CrsMatrix *K = NULL;

  // Read matrix from file or generate a matrix
  if ((gensize > 0 && testEpetra64)) {
    // Generate the matrix using long long for global indices
    build_simple_matrix<long long>(Comm, K, (long long)gensize, true, verbose);
  }
  else if (gensize) {
    // Generate the matrix using int for global indices
    build_simple_matrix<int>(Comm, K, gensize, false, verbose);
  }
  else {
    printf("YOU SHOULDN'T BE HERE \n");
    exit(-1);
  }

  if (verbose && (K->NumGlobalRows64() < TINYMATRIX)) {
    if (MyPID == 0) cout << "Input matrix:  " << endl;
    K->Print(cout);
  }
  Teuchos::RCP<Epetra_CrsMatrix> rcpK = Teuchos::rcp( K );

  // Set Anasazi verbosity level
  if (MyPID == 0) cout << "Setting up the problem..." << endl;

  int anasazi_verbosity = Anasazi::Errors + Anasazi::Warnings;
  if (verbosity >= 1)  // low
    anasazi_verbosity += Anasazi::FinalSummary + Anasazi::TimingDetails;
  if (verbosity >= 2)  // medium
    anasazi_verbosity += Anasazi::IterationDetails;
  if (verbosity >= 3)  // high
    anasazi_verbosity += Anasazi::StatusTestDetails
                       + Anasazi::OrthoDetails
                       + Anasazi::Debug;

  // Create parameter list to pass into solver
  Teuchos::ParameterList MyPL;
  MyPL.set("Verbosity", anasazi_verbosity);
  MyPL.set("Which", which);
  MyPL.set("Convergence Tolerance", tol);
  MyPL.set("Relative Convergence Tolerance", relconvtol);
  MyPL.set("Orthogonalization", ortho);

  // For the following, use Anasazi's defaults unless explicitly specified.
  if (numblocks > 0) MyPL.set( "Num Blocks", numblocks);
  if (maxrestarts > 0) MyPL.set( "Maximum Restarts", maxrestarts);
  if (maxiterations > 0) MyPL.set( "Maximum Iterations", maxiterations);
  if (blockSize > 0) MyPL.set( "Block Size", blockSize );

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Anasazi::MultiVecTraits<double, MV> MVT;
  typedef Anasazi::OperatorTraits<double, MV, OP> OPT;


  // Create the eigenproblem to be solved.

  // Dummy initial vectors - will be set later.
  Teuchos::RCP<Epetra_MultiVector> ivec =
    Teuchos::rcp(new Epetra_MultiVector(K->Map(), blockSize));

  Teuchos::RCP<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem;
  MyProblem =
    Teuchos::rcp(new Anasazi::BasicEigenproblem<double, MV, OP>(rcpK, ivec) );

  // Inform the eigenproblem whether K is Hermitian

  MyProblem->setHermitian(isHermitian);

  // Set the number of eigenvalues requested

  MyProblem->setNEV(nev);

  // Loop to solve the same eigenproblem numtrial times (different initial vectors)
  int numfailed = 0;
  int iter = 0;
  double solvetime = 0;

  // Set random seed to have consistent initial vectors between
  // experiments.  Different seed in each loop iteration.
  ivec->SetSeed(2*(MyPID) +1); // Odd seed

  // Set up initial vectors
  // Using random values as the initial guess.
  if (initvec == "random"){
    MVT::MvRandom(*ivec);
  }
  else if (initvec == "zero"){
    // All zero initial vector should be essentially the same,
    // but appears slightly worse in practice.
    ivec->PutScalar(0.);
  }
  else if (initvec == "unit"){
    // Orthogonal unit initial vectors.
    ivec->PutScalar(0.);
    for (int i = 0; i < blockSize; i++)
      ivec->ReplaceGlobalValue(i,i,1.);
  }
  else if (initvec == "random2"){
    // Partially random but orthogonal (0,1) initial vectors.
    // ivec(i,*) is zero in all but one column (for each i)
    // Inefficient implementation but this is only done once...
    double rowmax;
    int col;
    ivec->Random();
    for (int i = 0; i < ivec->MyLength(); i++){
      rowmax = -1;
      col = -1;
      for (int j = 0; j < blockSize; j++){
        // Make ivec(i,j) = 1 for largest random value in row i
        if ((*ivec)[j][i] > rowmax){
          rowmax = (*ivec)[j][i];
          col = j;
        }
        ivec->ReplaceMyValue(i,j,0.);
      }
      ivec->ReplaceMyValue(i,col,1.);
    }
  }
  else
    cout << "ERROR: Unknown value for initial vectors." << endl;

  if (verbose && (ivec->GlobalLength64() < TINYMATRIX))
    ivec->Print(std::cout);

  // Inform the eigenproblem that you are finished passing it information

  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error."
           << endl;
    }
    FINALIZE;
    return -1;
  }

  Teuchos::RCP<Anasazi::SolverManager<double, MV, OP> > MySolverMgr;

  if (method == "BKS") {
    // Initialize the Block Arnoldi solver
    MyPL.set("Extra NEV Blocks", extrablocks);
    MySolverMgr = Teuchos::rcp( new Anasazi::BlockKrylovSchurSolMgr<double, MV, OP>(MyProblem,MyPL) );
  }
  else if (method == "BD") {
    // Initialize the Block Davidson solver
    MyPL.set("Use Locking", lock);
    MyPL.set("Relative Locking Tolerance", rellocktol);
    MySolverMgr = Teuchos::rcp( new Anasazi::BlockDavidsonSolMgr<double, MV, OP>(MyProblem, MyPL) );
  }
  else if (method == "LOBPCG") {
    // Initialize the LOBPCG solver
    MyPL.set("Use Locking", lock);
    MyPL.set("Relative Locking Tolerance", rellocktol);
    MySolverMgr = Teuchos::rcp( new Anasazi::LOBPCGSolMgr<double, MV, OP>(MyProblem, MyPL) );
  }
  else if (method == "IRTR") {
    // Initialize the IRTR solver
    MySolverMgr = Teuchos::rcp( new Anasazi::RTRSolMgr<double, MV, OP>(MyProblem, MyPL) );
  }
  else
    cout << "Unknown solver method!" << endl;

  if (verbose && MyPID==0) MyPL.print(cout);

  // Solve the problem to the specified tolerances or length
  if (MyPID == 0) cout << "Beginning the " << method << " solve..." << endl;

  Anasazi::ReturnType returnCode = MySolverMgr->solve();
  if (returnCode != Anasazi::Converged && MyPID==0) {
    ++numfailed;
    cout << "Anasazi::SolverManager::solve() returned unconverged." << endl;
  }
  iter = MySolverMgr->getNumIters();
  solvetime = (MySolverMgr->getTimers()[0])->totalElapsedTime();

  if (MyPID == 0) {
    cout << "Iterations in this solve: " << iter << endl;
    cout << "Solve complete; beginning post-processing..."<< endl;
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem

  Anasazi::Eigensolution<double,MV> sol = MyProblem->getSolution();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  Teuchos::RCP<MV> evecs = sol.Evecs;
  std::vector<int> index = sol.index;
  int numev = sol.numVecs;

  // Compute residuals.

  if (numev > 0) {
    Teuchos::LAPACK<int,double> lapack;
    std::vector<double> normR(numev);

    if (MyProblem->isHermitian()) {
      // Get storage
      Epetra_MultiVector Kevecs(K->Map(),numev);
      Teuchos::RCP<Epetra_MultiVector> Mevecs;
      Teuchos::SerialDenseMatrix<int,double> B(numev,numev);
      B.putScalar(0.0);
      for (int i=0; i<numev; i++) {B(i,i) = evals[i].realpart;}

      // Compute A*evecs
      OPT::Apply( *rcpK, *evecs, Kevecs );
      Mevecs = evecs;

      // Compute A*evecs - lambda*evecs and its norm
      MVT::MvTimesMatAddMv( -1.0, *Mevecs, B, 1.0, Kevecs );
      MVT::MvNorm( Kevecs, normR );

      // Scale the norms by the eigenvalue if relative convergence tol was used
      if (relconvtol) {
        for (int i=0; i<numev; i++)
          normR[i] /= Teuchos::ScalarTraits<double>::magnitude(evals[i].realpart);
      }

    } else {
      printf("The problem isn't non-Hermitian; sorry.\n");
      exit(-1);
    }


    if (verbose && MyPID==0) {
      cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      cout<<endl<< "Actual Results"<<endl;
      if (MyProblem->isHermitian()) {
        cout<< std::setw(16) << "Eigenvalue "
            << std::setw(20) << "Direct Residual"
            << (relconvtol?" (normalized by eigenvalue)":" (no normalization)")
            << endl;
        cout<<"--------------------------------------------------------"<<endl;
        for (int i=0; i<numev; i++) {
          cout<< "EV" << i << std::setw(16) << evals[i].realpart
              << std::setw(20) << normR[i] << endl;
        }
        cout<<"--------------------------------------------------------"<<endl;
      }
      else {
        cout<< std::setw(16) << "Real Part"
            << std::setw(16) << "Imag Part"
            << std::setw(20) << "Direct Residual"<< endl;
        cout<<"--------------------------------------------------------"<<endl;
        for (int i=0; i<numev; i++) {
          cout<< std::setw(16) << evals[i].realpart
              << std::setw(16) << evals[i].imagpart
              << std::setw(20) << normR[i] << endl;
        }
        cout<<"--------------------------------------------------------"<<endl;
      }
    }
  }

  // Summarize iteration counts and solve time
  if (MyPID == 0) {
    cout << endl;
    cout << "DRIVER SUMMARY" << endl;
    cout << "Failed to converge: " << numfailed << endl;
    cout << "Solve time:           " << solvetime << endl;
  }

  FINALIZE;

  if (numfailed) {
    if (MyPID == 0) {
      cout << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }
  //
  // Default return value
  //
  if (MyPID == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;
}
