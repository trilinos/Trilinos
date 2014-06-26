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
// Multiple right-hand-sides are created randomly.
// The initial guesses are all set to zero.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosTFQMRSolMgr.hpp"
#include "createEpetraProblem.hpp"
#include "Ifpack_IlukGraph.h"
#include "Ifpack_CrsRiluk.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

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
    bool leftprec = true; // use left preconditioning to solve these linear systems
    int frequency = -1;  // how often residuals are printed by solver
    int numrhs = 1;
    int maxiters = -1;    // maximum iterations allowed
    std::string filename("orsirr1.hb");
    MT tol = 1.0e-5;  // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("left-prec","right-prec",&leftprec,"Left preconditioning or right.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
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
    int return_val =Belos::createEpetraProblem(filename,NULL,&A,NULL,NULL,&MyPID);
    const Epetra_Map &Map = A->RowMap();
    if(return_val != 0) return return_val;
    proc_verbose = verbose && (MyPID==0); /* Only print on zero processor */
    //
    // *****Construct the Preconditioner*****
    //
    if (proc_verbose) std::cout << std::endl << std::endl;
    if (proc_verbose) std::cout << "Constructing ILU preconditioner" << std::endl;
    int Lfill = 2;
    // if (argc > 2) Lfill = atoi(argv[2]);
    if (proc_verbose) std::cout << "Using Lfill = " << Lfill << std::endl;
    int Overlap = 2;
    // if (argc > 3) Overlap = atoi(argv[3]);
    if (proc_verbose) std::cout << "Using Level Overlap = " << Overlap << std::endl;
    double Athresh = 0.0;
    // if (argc > 4) Athresh = atof(argv[4]);
    if (proc_verbose) std::cout << "Using Absolute Threshold Value of " << Athresh << std::endl;
    double Rthresh = 1.0;
    // if (argc >5) Rthresh = atof(argv[5]);
    if (proc_verbose) std::cout << "Using Relative Threshold Value of " << Rthresh << std::endl;
    //
    Teuchos::RCP<Ifpack_IlukGraph> ilukGraph;
    Teuchos::RCP<Ifpack_CrsRiluk> ilukFactors;
    //
    if (Lfill > -1) {
      ilukGraph = Teuchos::rcp(new Ifpack_IlukGraph(A->Graph(), Lfill, Overlap));
      int info = ilukGraph->ConstructFilledGraph();
      assert( info == 0 );
      ilukFactors = Teuchos::rcp(new Ifpack_CrsRiluk(*ilukGraph));
      int initerr = ilukFactors->InitValues(*A);
      if (initerr != 0) std::cout << "InitValues error = " << initerr;
      info = ilukFactors->Factor();
      assert( info == 0 );
    }
    //
    bool transA = false;
    double Cond_Est;
    ilukFactors->Condest(transA, Cond_Est);
    if (proc_verbose) {
      std::cout << "Condition number estimate for this preconditoner = " << Cond_Est << std::endl;
      std::cout << std::endl;
    }

    //
    // Create the Belos preconditioned operator from the Ifpack preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    RCP<Belos::EpetraPrecOp> Prec = rcp( new Belos::EpetraPrecOp( ilukFactors ) );

    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = Map.NumGlobalElements();
    if (maxiters == -1)
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    if (leftprec)
      belosList.set( "Explicit Residual Test", true );     // Need to check for the explicit residual before returning
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    //
    // *****Construct solution std::vector and random right-hand-sides *****
    //
    RCP<Epetra_MultiVector> X = rcp( new Epetra_MultiVector(Map, numrhs) );
    X->PutScalar( 0.0 );
    RCP<Epetra_MultiVector> B = rcp( new Epetra_MultiVector(Map, numrhs) );
    B->Random();
    Belos::LinearProblem<double,MV,OP> problem( A, X, B );
    if (leftprec)
      problem.setLeftPrec( Prec );
    else
      problem.setRightPrec( Prec );

    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // *******************************************************************
    // *****************Start the TFQMR iteration*************************
    // *******************************************************************
    //
    Teuchos::RCP< Belos::TFQMRSolMgr<double,MV,OP> > solver =
      Teuchos::rcp( new Belos::TFQMRSolMgr<double,MV,OP>( rcp(&problem,false), rcp(&belosList,false) ) );
    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of TFQMR iterations: " << maxiters << std::endl;
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
    bool badRes = false;
    std::vector<double> actual_resids( numrhs );
    std::vector<double> rhs_norm( numrhs );
    Epetra_MultiVector R(Map, numrhs);
    OPT::Apply( *A, *X, R );
    MVT::MvAddMv( -1.0, R, 1.0, *B, R );
    MVT::MvNorm( R, actual_resids );
    MVT::MvNorm( *B, rhs_norm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol ) badRes = true;
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

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
