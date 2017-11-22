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
// The right-hand-side from the problem is being used instead of multiple
// random right-hand-sides.  The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosEpetraUtils.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Epetra_Map.h"

int main(int argc, char *argv[]) {
  //
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
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
    bool pseudo = false;   // use pseudo block GMRES to solve this linear system.
    int frequency = -1;
    int blocksize = 1;
    int numrhs = 1;
    int maxrestarts = 15; // number of restarts allowed
    int maxiters = -1;    // maximum number of iterations allowed per linear system
    std::string filename("orsirr1.hb");
    MT tol = 1.0e-5;  // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("pseudo","regular",&pseudo,"Use pseudo-block GMRES to solve the linear systems.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("max-restarts",&maxrestarts,"Maximum number of restarts allowed for GMRES solver.");
    cmdp.setOption("blocksize",&blocksize,"Block size used by GMRES.");
    cmdp.setOption("maxiters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
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
    const Epetra_Map &Map = A->RowMap();
    proc_verbose = verbose && (MyPID==0);  /* Only print on the zero processor */
    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->GlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements/blocksize - 1; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Num Blocks", maxiters );               // Maximum number of blocks in Krylov factorization
    belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Maximum Restarts", maxrestarts );      // Maximum number of restarts allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    //
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
    //
    // *******************************************************************
    // *************Start the block Gmres iteration*************************
    // *******************************************************************
    //
    Teuchos::RCP< Belos::SolverManager<double,MV,OP> > solver;
    if (pseudo)
      solver = Teuchos::rcp( new Belos::PseudoBlockGmresSolMgr<double,MV,OP>( rcp(&problem,false), rcp(&belosList,false) ) );
    else
      solver = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,MV,OP>( rcp(&problem,false), rcp(&belosList,false) ) );
    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Block size used by solver: " << blocksize << std::endl;
      std::cout << "Number of restarts allowed: " << maxrestarts << std::endl;
      std::cout << "Max number of Gmres iterations per restart cycle: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = solver->solve();

    //
    // Compute actual residuals.
    //
    std::vector<double> actual_resids( numrhs );
    std::vector<double> rhs_norm( numrhs );
    Epetra_MultiVector resid(Map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        std::cout<<"Problem "<<i<<" : \t"<< actual_resids[i]/rhs_norm[i] <<std::endl;
      }
    }

    success = ret==Belos::Converged;

    if (success) {
      if (proc_verbose)
        std::cout << "End Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_bl_gmres_hb.cpp
