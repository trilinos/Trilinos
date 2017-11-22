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
// This driver reads a problem from a Harwell-Boeing (HB) file.
// The right-hand-side from the HB file is used instead of random vectors.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosRCGSolMgr.hpp"
#include "BelosEpetraUtils.h"
#include "Trilinos_Util.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Ifpack.h"

int main(int argc, char *argv[]) {
  //
#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
#endif
  //
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool success = false;
  bool verbose = false;
  try {
    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool leftprec = true;                // left preconditioning or right.
    int frequency = -1;                  // frequency of status test output.
    std::string filename("gr_30_30.hb"); // default input filename
    double tol = 1.0e-10;                // relative residual tolerance
    int numBlocks = 30;                  // maximum number of blocks the solver can use for the Krylov subspace
    int recycleBlocks = 3;               // maximum number of blocks the solver can use for the recycle space
    int numrhs = 1;                      // number of right-hand sides to solve for
    int maxiters = -1;                   // maximum number of iterations allowed per linear system

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("left-prec","right-prec",&leftprec,"Left preconditioning or right.");
    cmdp.setOption("filename",&filename,"Filename for test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by the RCG solver.");
    cmdp.setOption("max-subspace",&numBlocks,"Maximum number of vectors in search space (not including recycle space).");
    cmdp.setOption("recycle",&recycleBlocks,"Number of vectors in recycle space.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose
    //
    // Get the problem
    //
    int MyPID;
    RCP<Epetra_CrsMatrix> A;
    RCP<Epetra_MultiVector> B, X;
    int return_val =Belos::Util::createEpetraProblem(filename,NULL,&A,&B,&X,&MyPID);
    if(return_val != 0) return return_val;
    proc_verbose = ( verbose && (MyPID==0) );
    //
    // Solve using Belos
    //
    typedef double                          ST;
    typedef Epetra_Operator                 OP;
    typedef Epetra_MultiVector              MV;
    typedef Belos::OperatorTraits<ST,MV,OP> OPT;
    typedef Belos::MultiVecTraits<ST,MV>    MVT;
    //
    // *****Construct initial guess and right-hand sides *****
    //
    if (numrhs != 1) {
      X = rcp( new Epetra_MultiVector( A->Map(), numrhs ) );
      X->Random();
      B = rcp( new Epetra_MultiVector( A->Map(), numrhs ) );
      OPT::Apply( *A, *X, *B );
      MVT::MvInit( *X, 0.0 );
    }
    else { // initialize exact solution to be vector of ones
      MVT::MvInit( *X, 1.0 );
      OPT::Apply( *A, *X, *B );
      MVT::MvInit( *X, 0.0 );
    }
    //
    // ************Construct preconditioner*************
    //
    ParameterList ifpackList;

    // allocates an IFPACK factory. No data is associated
    // to this object (only method Create()).
    Ifpack Factory;

    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    std::string PrecType = "ICT"; // incomplete Cholesky
    int OverlapLevel = 0; // must be >= 0. If Comm.NumProc() == 1,
    // it is ignored.

    RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create(PrecType, &*A, OverlapLevel) );
    assert(Prec != Teuchos::null);

    // specify parameters for ICT
    ifpackList.set("fact: drop tolerance", 1e-9);
    ifpackList.set("fact: ict level-of-fill", 1.0);
    // the combine mode is on the following:
    // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
    // Their meaning is as defined in file Epetra_CombineMode.h
    ifpackList.set("schwarz: combine mode", "Add");
    // sets the parameters
    IFPACK_CHK_ERR(Prec->SetParameters(ifpackList));

    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    IFPACK_CHK_ERR(Prec->Initialize());

    // Builds the preconditioners, by looking for the values of
    // the matrix.
    IFPACK_CHK_ERR(Prec->Compute());

    // Create the Belos preconditioned operator from the Ifpack preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    RCP<Belos::EpetraPrecOp> belosPrec = rcp( new Belos::EpetraPrecOp( Prec ) );

    //
    // ************Other information used by solver**************
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->GlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Num Blocks", numBlocks);               // Maximum number of blocks in Krylov space
    belosList.set( "Num Recycled Blocks", recycleBlocks ); // Number of vectors in recycle space
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    //
    // Construct a preconditioned linear problem instance.
    //
    Belos::LinearProblem<double,MV,OP> problem( A, X, B );
    if (leftprec) {
      problem.setLeftPrec( belosPrec );
    }
    else {
      problem.setRightPrec( belosPrec );
    }

    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // Create an iterative solver manager.
    //
    RCP< Belos::SolverManager<double,MV,OP> > newSolver
      = rcp( new Belos::RCGSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)) );
    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of RCG iterations: " << maxiters << std::endl;
      std::cout << "Max number of vectors in Krylov space: " << numBlocks << std::endl;
      std::cout << "Number of vectors in recycle space: " << recycleBlocks << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = newSolver->solve();
    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<double> actual_resids( numrhs );
    std::vector<double> rhs_norm( numrhs );
    Epetra_MultiVector resid(A->Map(), numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) badRes = true;
      }
    }

    success = ret==Belos::Converged && !badRes;

    if (success) {
      if (proc_verbose)
        std::cout << std::endl << "End Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << std::endl << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_rcg_hb.cpp
