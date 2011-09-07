//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// ************************************************************************
//@HEADER
//
// This driver reads a problem from a file, which can be in Harwell-Boeing (*.hb),
// Matrix Market (*.mtx), or triplet format (*.triU, *.triS).  The right-hand side
// from the problem, if it exists, will be used instead of multiple
// random right-hand sides.  The initial guesses are all set to zero. 
//
// NOTE: No preconditioner is used in this example. 
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosLSQRSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"
#include "Epetra_Map.h"
#ifdef EPETRA_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

int main(int argc, char *argv[]) {
  //
  int MyPID = 0;
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  MyPID = Comm.MyPID();
#else
  Epetra_SerialComm Comm;
#endif
  //
  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool verbose = false, debug = false, proc_verbose = false;
  int frequency = -1;        // frequency of status test output.
  int blocksize = 1;         // blocksize
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  std::string filename("orsirr1_scaled.hb");
  MT relResTol = 2.0e-4;     // relative residual tolerance
  // Like CG, LSQR is a short recurrence method that 
  // does not have the "n" step convergence property in finite precision arithmetic.
  MT resGrowthFactor = 2.0;   // In this example, warn if |resid| > resGrowthFactor * relResTol
  // With no preconditioner, this is only the difference between the "implicit" and the "explict
  // residual.  

  MT relMatTol = 1.4e-7;     // relative Matrix error, default value sqrt(eps)
  MT maxCond  = 1.e+16;      // maximum condition number default value 1/eps
  MT damp = 0.;              // regularization (or damping) parameter 

  Teuchos::CommandLineProcessor cmdp(false,true); // e.g. ./a.out --tol=.1 --filename=foo.hb

  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nondebug",&debug,"Print debugging information from solver.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filename,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("lambda",&damp,"Regularization parameter");
  cmdp.setOption("tol",&relResTol,"Relative residual tolerance");
  cmdp.setOption("matrixTol",&relMatTol,"Relative error in Matrix");
  cmdp.setOption("max-cond",&maxCond,"Maximum condition number");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blocksize,"Block size used by LSQR."); // must be one at this point
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");



  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (!verbose)
    frequency = -1;  // reset frequency if test is not verbose
  //
  // Get the problem
  //
  RCP<Epetra_Map> Map;
  RCP<Epetra_CrsMatrix> A;
  RCP<Epetra_MultiVector> B, X;
  RCP<Epetra_Vector> vecB, vecX;
  // Trilinos_Util_ReadMatrixMarket2Epetra vecX := 0,  vecB = A*randVec
  // It is odd that the exact solution is discarded.
  EpetraExt::readEpetraLinearSystem(filename, Comm, &A, &Map, &vecX, &vecB);
  A->OptimizeStorage();
  proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

  // Check to see if the number of right-hand sides is the same as requested.
  if (numrhs>1) {
    X = rcp( new Epetra_MultiVector( *Map, numrhs ) );
    B = rcp( new Epetra_MultiVector( *Map, numrhs ) );
    X->Seed();
    X->Random();
    OPT::Apply( *A, *X, *B ); // B := AX
    X->PutScalar( 0.0 );   // annihilate X
  }
  else {
    int locNumCol = Map->MaxLID() + 1;
    int globNumCol = Map->MaxAllGID() + 1;
    for( int li = 0; li < locNumCol; li++){   // assume consecutive lid
      int gid = Map->GID(li);
      double value = (double) ( globNumCol -1 - gid ); 
      int numEntries = 1;
      vecX->ReplaceGlobalValues( numEntries, &value, &gid );
    }
    bool Trans = false;
    A->Multiply( Trans, *vecX, *vecB );
    bool goodInitGuess = false;
    if( goodInitGuess )
      {
        double value = 1.e-2;  // solution norm ~ n^(3/2) ~ 10^6
        int numEntries = 1;   
        int index = 0;
        vecX->SumIntoMyValues(  numEntries, &value, &index); 
      }
    // no initial guess
    vecX->PutScalar( 0.0 ); // code breaks with nonzero initial guess.  ok .. that's my job.
    // Is the relative residual threshold reset when there is an initial guess? 
    // If so, then where?
    X = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecX);
    B = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecB);
  }
  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->GlobalLength();
  if (maxiters == -1)
    maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
  ParameterList belosList; // mechanism for configuring specific linear solver
  belosList.set( "Block Size", blocksize );       // LSQR blocksize, must be one
  belosList.set( "Lambda", damp );                // Regularization parameter
  belosList.set( "Rel RHS Err", relResTol );      // Relative convergence tolerance requested
  belosList.set( "Rel Mat Err", relMatTol );      // Maximum number of restarts allowed
  belosList.set( "Condition limit", maxCond);     // upper bound for cond(A)
  belosList.set( "Maximum Iterations", maxiters );// Maximum number of iterations allowed
  int verbosity = Belos::Errors + Belos::Warnings;
  if (verbose) {
    verbosity += Belos::TimingDetails + Belos::StatusTestDetails;
    if (frequency > 0)
      belosList.set( "Output Frequency", frequency );
  }
  if (debug) {
    verbosity += Belos::Debug;
  }
  belosList.set( "Verbosity", verbosity );
  //
  // Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<double,MV,OP> problem( A, X, B );
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
  // *******************************************************************
  // ******************* Apply Single Vector LSQR **********************
  // *******************************************************************
  // Create an iterative solver manager.
  RCP< Belos::LSQRSolMgr<double,MV,OP> > newSolver
    = rcp( new Belos::LSQRSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)));

  if (proc_verbose) { // ******** Print a problem description *********
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Max number of iterations for a linear system: " << maxiters << std::endl; 
    std::cout << "Relative residual tolerance: " << relResTol << std::endl;
    std::cout << std::endl;
    std::cout << "Solver's Description: " << std::endl;
    std::cout << newSolver->description() << std::endl; // visually verify the parameter list
  }
  Belos::ReturnType ret = newSolver->solve(); // Perform solve
  std::vector<double> solNorm( numrhs );      // get solution norm
  MVT::MvNorm( *X, solNorm );
  int numIters = newSolver->getNumIters();    // get number of solver iterations
  MT condNum = newSolver->getMatCondNum(); 
  MT matrixNorm= newSolver->getMatNorm();
  MT resNorm = newSolver->getResNorm(); 
  MT lsResNorm = newSolver->getMatResNorm();

  if (proc_verbose)
    std::cout << "Number of iterations performed for this solve: " << numIters << std::endl
     << "matrix condition number: " << condNum << std::endl
     << "matrix norm: " << matrixNorm << std::endl
     << "residual norm: " << resNorm << std::endl
     << "solution norm: " << solNorm[0] << std::endl
     << "least squares residual Norm: " << lsResNorm << std::endl;
  bool badRes = false;                     // Compute the actual residuals.
  std::vector<double> actual_resids( numrhs );
  std::vector<double> rhs_norm( numrhs );
  Epetra_MultiVector resid(*Map, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -1.0, resid, 1.0, *B, resid ); 
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  if (proc_verbose) {
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    for ( int i=0; i<numrhs; i++) {
      double actRes = actual_resids[i]/rhs_norm[i];
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      if (actRes > relResTol * resGrowthFactor)
        { 
          badRes = true;
          if (verbose) std::cout << "residual norm > " << relResTol * resGrowthFactor <<  std::endl;	
        }
    }
  }
  if (ret!=Belos::Converged || badRes) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;	
    return -1;
  }
  if (proc_verbose)
    std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
#ifdef EPETRA_MPI
  MPI_Finalize();
#endif
  return 0; // Default return value
} 
