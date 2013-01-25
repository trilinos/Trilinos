//  This program computes the linear solution using Belos and a laplacian matrix.
//  NOTE:  The laplacian matrix generated has row and column ids larger than INT_MAX 
//         to test Epetra64 functionality with Belos solvers.
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

#include "BelosBlockCGSolMgr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosEpetraAdapter.hpp"

#include "Epetra_Time.h"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_as.hpp"

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

  // Linear solver properties

  std::string method = "BlockCG";
  std::string ortho = "DGKS";
  int blockSize = 1;
  int numrhs = 1;
  int numblocks = -1;
  int maxrestarts = 0;
  int maxiterations = -1;
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
                 "Solver method to use:  BlockCG.");
  cmdp.setOption("tol",&tol,"Solver convergence tolerance.");
  cmdp.setOption("blocksize",&blockSize,"Block size to use in solver.");
  cmdp.setOption("numrhs",&numrhs,"Number of right-hand sides.");
  cmdp.setOption("numblocks",&numblocks,"Number of blocks to allocate.");
  cmdp.setOption("maxrestarts",&maxrestarts,
                 "Maximum number of restarts.");
  cmdp.setOption("maxiterations",&maxiterations,
                 "Maximum number of iterations.");
  cmdp.setOption("ortho", &ortho,
                 "Orthogonalization method (DGKS, ICGS, IMGS).");
  if (cmdp.parse(narg,arg)!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    FINALIZE;
    return -1;
  }

  // Print the most essential options (not in the MyPL parameters later)
  verbose = (verbosity>0);
  if (verbose && MyPID==0){
    cout << "verbosity = " << verbosity << endl;
    cout << "method = " << method << endl;
    cout << "numrhs = " << numrhs << endl;
  }

  // Make sure Epetra was built with 64-bit global indices enabled.
#ifdef EPETRA_NO_64BIT_GLOBAL_INDICES
  if (testEpetra64)
    testEpetra64 = false;
#endif

  Epetra_CrsMatrix *A = NULL;
  
  // Read matrix from file or generate a matrix
  if ((gensize > 0 && testEpetra64)) {
    // Generate the matrix using long long for global indices
    build_simple_matrix<long long>(Comm, A, (long long)gensize, true, verbose);
  }
  else if (gensize) {
    // Generate the matrix using int for global indices
    build_simple_matrix<int>(Comm, A, gensize, false, verbose);
  }
  else {
    printf("YOU SHOULDN'T BE HERE \n");
    exit(-1);
  }

  if (verbose && (A->NumGlobalRows64() < TINYMATRIX)) {
    if (MyPID == 0) cout << "Input matrix:  " << endl;
    A->Print(cout);
  }

  // Set Belos verbosity level
  if (MyPID == 0) cout << "Setting up the problem..." << endl;

  int belos_verbosity = Belos::Errors + Belos::Warnings;
  if (verbosity >= 1)  // low
    belos_verbosity += Belos::FinalSummary + Belos::TimingDetails;
  if (verbosity >= 2)  // medium
    belos_verbosity += Belos::IterationDetails;
  if (verbosity >= 3)  // high
    belos_verbosity += Belos::StatusTestDetails
                       + Belos::OrthoDetails
                       + Belos::Debug;
  
  // Create parameter list to pass into solver

  Teuchos::ParameterList MyPL;
  MyPL.set("Verbosity", belos_verbosity);
  MyPL.set("Convergence Tolerance", tol);
  MyPL.set("Orthogonalization", ortho);

  // For the following, use Belos's defaults unless explicitly specified.
  if (numblocks > 0) MyPL.set( "Num Blocks", numblocks);
  if (maxrestarts > 0) MyPL.set( "Maximum Restarts", maxrestarts);
  if (maxiterations <= 0) 
    maxiterations = Teuchos::asSafe<int>( gensize );
  MyPL.set( "Maximum Iterations", maxiterations);
  if (blockSize > 0) MyPL.set( "Block Size", blockSize );

  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Belos::MultiVecTraits<double, MV> MVT;
  typedef Belos::OperatorTraits<double, MV, OP> OPT;
    

  // Create the linear problem to be solved.
    
  Teuchos::RCP<Epetra_MultiVector> X = 
    Teuchos::rcp(new Epetra_MultiVector(A->Map(), numrhs));
  Teuchos::RCP<Epetra_MultiVector> B = 
    Teuchos::rcp(new Epetra_MultiVector(A->Map(), numrhs));

  Teuchos::RCP<Belos::LinearProblem<double, MV, OP> > MyProblem;
  MyProblem = 
    Teuchos::rcp(new Belos::LinearProblem<double, MV, OP>(Teuchos::rcp(A,false), X, B) );

  // Inform the linear problem whether A is Hermitian

  MyProblem->setHermitian( isHermitian );

  int iter = 0;

  // Set random seed to have consistent initial vectors between experiments.  
  X->SetSeed(2*(MyPID) +1); // Odd seed
  MVT::MvRandom( *X );
  OPT::Apply( *A, *X, *B );
  MVT::MvInit( *X, 0. );

  // Inform the linear problem that you are finished passing it information
  bool boolret = MyProblem->setProblem();
  if (boolret != true) {
    if (verbose && MyPID == 0) {
      cout << "Belos::LinearProblem::setProblem() returned with error." 
           << endl;
    }
    FINALIZE;
    return -1;
  }
 
  Teuchos::RCP<Belos::SolverManager<double, MV, OP> > MySolverMgr;
 
  if (method == "BlockCG") {
    // Initialize the Block CG solver
    MySolverMgr = Teuchos::rcp( new Belos::BlockCGSolMgr<double, MV, OP>(MyProblem,Teuchos::rcp(&MyPL, false)) );
  }
  else
    cout << "Unknown solver method!" << endl;

  if (verbose && MyPID==0) MyPL.print(cout);
      
  // Solve the problem to the specified tolerances or length
  if (MyPID == 0) cout << "Beginning the " << method << " solve..." << endl;

  int numfailed = 0; 
  Belos::ReturnType returnCode = MySolverMgr->solve();
  if (returnCode != Belos::Converged && MyPID==0) {
    ++numfailed;
    cout << "Belos::SolverManager::solve() returned unconverged." << endl;
  }
  iter = MySolverMgr->getNumIters();
  
  if (MyPID == 0) {
    cout << "Iterations in this solve: " << iter << endl; 
    cout << "Solve complete; beginning post-processing..."<< endl;
  }
  
  // Compute residuals.
  bool badRes = false;
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Epetra_MultiVector resid(A->Map(), numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  if (MyPID==0) 
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  for ( int i=0; i<numrhs; i++) {
    double actRes = actual_resids[i]/rhs_norm[i];
    if (MyPID==0)
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
    if (actRes > tol) badRes = true;
  }
 
  FINALIZE;

  if (badRes || numfailed) {
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
