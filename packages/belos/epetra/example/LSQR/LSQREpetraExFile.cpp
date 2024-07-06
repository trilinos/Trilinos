// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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
#include "EpetraExt_MultiVectorIn.h"

#include "Epetra_Map.h"
#ifdef EPETRA_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

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

bool verbose = false;
bool success = true;
try {
bool proc_verbose = false;
  bool debug = false;
  int frequency = -1;        // frequency of status test output.
  int blocksize = 1;         // blocksize
  int numrhs = 1;            // number of right-hand sides to solve for
  int maxiters = -1;         // maximum number of iterations allowed per linear system
  std::string filenameMatrix("orsirr1_scaled.hb");
  std::string filenameRHS;   // blank mean unset
  MT relResTol = 3.0e-4;     // relative residual tolerance
  // Like CG, LSQR is a short recurrence method that
  // does not have the "n" step convergence property in finite precision arithmetic.
  MT resGrowthFactor = 4.0;   // In this example, warn if |resid| > resGrowthFactor * relResTol
  // With no preconditioner, this is only the difference between the "implicit" and the "explict
  // residual.

  MT relMatTol = 1.e-4;     // relative Matrix error, default value sqrt(eps)
  MT maxCond  = 1.e+8;      // maximum condition number default value 1/eps
  MT damp = 0.;             // regularization (or damping) parameter

  Teuchos::CommandLineProcessor cmdp(false,true); // e.g. ./a.out --tol=.1 --filename=foo.hb

  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nondebug",&debug,"Print debugging information from solver.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("filename",&filenameMatrix,"Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
  cmdp.setOption("rhsFilename",&filenameRHS,"Filename for right-hand side.  Acceptable file extension: *.mtx");
  cmdp.setOption("lambda",&damp,"Regularization parameter");
  cmdp.setOption("tol",&relResTol,"Relative residual tolerance");
  cmdp.setOption("matrixTol",&relMatTol,"Relative error in Matrix");
  cmdp.setOption("max-cond",&maxCond,"Maximum condition number");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  cmdp.setOption("block-size",&blocksize,"Block size used by LSQR."); // must be one at this point
  cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");

  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef EPETRA_MPI
    MPI_Finalize();
#endif
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
  EpetraExt::readEpetraLinearSystem(filenameMatrix, Comm, &A, &Map, &vecX, &vecB);
  // Rectangular matrices are embedded in square matrices.  vecX := 0,  vecB = A*randVec
  A->OptimizeStorage();
  proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */

  bool isRHS = false;
  if ( filenameRHS != "" )
    {
      isRHS = true;
    }

  // Check to see if the number of right-hand sides is the same as requested.
  if (numrhs>1) {
    isRHS = false; // numrhs > 1 not yet supported
    X = rcp( new Epetra_MultiVector( *Map, numrhs ) );
    B = rcp( new Epetra_MultiVector( *Map, numrhs ) );
    X->Random();
    OPT::Apply( *A, *X, *B ); // B := AX
    X->PutScalar( 0.0 );   // annihilate X
  }
  else {
    if ( isRHS )
      {
        Epetra_MultiVector * BmustDelete;
        int mmRHSioflag = 0;
        const char * charPtrRHSfn = filenameRHS.c_str();
        mmRHSioflag = EpetraExt::MatrixMarketFileToMultiVector(charPtrRHSfn, *Map, BmustDelete);
        //std::cout << "rhs from input file " << std::endl;
        //BmustDelete->Print(std::cout);

        if( mmRHSioflag )
          {
            if (proc_verbose)
              std::cout << "Error " <<  mmRHSioflag << " occured while attempting to read file " << filenameRHS << std::endl;
#ifdef EPETRA_MPI
            MPI_Finalize();
#endif
            return -1;
          }
        X = rcp( new Epetra_MultiVector( *Map, numrhs ) );
        X->Scale( 0.0 );
        B = rcp( new MV(*BmustDelete));
        delete BmustDelete;
      }
    else
      {
        int locNumCol = Map->MaxLID() + 1; // Create a known solution
        int globNumCol = Map->MaxAllGID() + 1;
        for( int li = 0; li < locNumCol; li++){   // assume consecutive lid
          int gid = Map->GID(li);
          double value = (double) ( globNumCol -1 - gid );
          int numEntries = 1;
          vecX->ReplaceGlobalValues( numEntries, &value, &gid );
        }
        bool Trans = false;
        A->Multiply( Trans, *vecX, *vecB ); // Create a consistent linear system
        // At this point, the initial guess is exact.
        bool goodInitGuess = true; // perturb initial guess
        bool zeroInitGuess = false; // annihilate initial guess
        if( goodInitGuess )
          {
            double value = 1.e-2; // "Rel RHS Err" and "Rel Mat Err" apply to the residual equation,
            int numEntries = 1;   // norm( b - A x_k ) ?<? relResTol norm( b- Axo).
            int index = 0;        // norm(b) is inaccessible to LSQR.
            vecX->SumIntoMyValues(  numEntries, &value, &index);
          }
        if( zeroInitGuess )
          {
            vecX->PutScalar( 0.0 ); //
          }
        X = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecX);
        B = Teuchos::rcp_implicit_cast<Epetra_MultiVector>(vecB);
      }
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
  belosList.set( "Condition Limit", maxCond);     // upper bound for cond(A)
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
      {
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      }
#ifdef EPETRA_MPI
      MPI_Finalize();
#endif
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
  success = false;
  if (proc_verbose)
    std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
} else {
  success = true;
  if (proc_verbose)
    std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
}
}
TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef EPETRA_MPI
MPI_Finalize();
#endif

return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
