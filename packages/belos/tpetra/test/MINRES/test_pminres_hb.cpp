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

/* Originally convert test here: belos/epetra/test/MINRES/test_pminres_hb.cpp */

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosTpetraOperator.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MatrixIO.hpp>

// Teuchos
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

template<class ScalarType>
int run(int argc, char *argv[]){
  using ST = ScalarType;
  using SCT = typename Teuchos::ScalarTraits<ST>;
  using MT = typename SCT::magnitudeType;
  using LO = typename Tpetra::Vector<>::local_ordinal_type;
  using GO = typename Tpetra::Vector<>::global_ordinal_type;
  using NT = typename Tpetra::Vector<>::node_type;
  
  using OP = typename Tpetra::Operator<ST,LO,GO,NT>;
  using MV = typename Tpetra::MultiVector<ST,LO,GO,NT>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;

  using tcrsmatrix_t = typename Tpetra::CrsMatrix<ST,LO,GO,NT>;

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();
    
  bool success = false;
  bool verbose = false;

  try {
    int MyPID = rank(*comm);
    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool leftprec = true;      // left preconditioning or right.
    int frequency = -1;        // how often residuals are printed by solver
    int numrhs = 5;            // total number of right-hand sides to solve for
    int maxiters = -1;         // maximum number of iterations for the solver to use
    std::string filename("bcsstk14.hb");
    ST tol = sqrt(std::numeric_limits<ST>::epsilon()); // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("left-prec","right-prec",&leftprec,"Left preconditioning or right.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by Minres solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");

    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // Reset frequency if verbosity is off

    proc_verbose = ( verbose && (MyPID==0) );
    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }

    //
    // Get the problem
    //
    RCP<tcrsmatrix_t> A;
    Tpetra::Utils::readHBMatrix(filename,comm,A);
    RCP<const Tpetra::Map<> > map = A->getDomainMap();

    // Create initial vectors
    RCP<MV> X, B;
    X = rcp(new MV(map, numrhs));
    MVT::MvRandom( *X );
    B = rcp(new MV(map, numrhs));
    OPT::Apply(*A, *X, *B);
    MVT::MvInit(*X, 0.0);
    
    //
    // *****Create parameter list for the Minres solver manager*****
    //
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1)
      maxiters = NumGlobalElements - 1; // maximum number of iterations to run
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    belosList.set( "Flexible Gmres", true );               // use FGMRES
                                                           
    int verbLevel = Belos::Errors + Belos::Warnings;
    if (verbose) {
      verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
    }
    belosList.set( "Verbosity", verbLevel );
    if (verbose) {
      if (frequency > 0) {
        belosList.set( "Output Frequency", frequency );
      }
    }

    // Set parameters for the inner GMRES (preconditioning) iteration.
    ParameterList innerBelosList;
    innerBelosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    innerBelosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
    innerBelosList.set( "Verbosity", Belos::Errors + Belos::Warnings );

    // *****Construct linear problem for the inner iteration using A *****
    Belos::LinearProblem<ST,MV,OP> innerProblem;
    innerProblem.setOperator( A );
    innerProblem.setLabel( "Belos Preconditioner Solve" );

    //  
    // *****Create the inner Minres iteration********
    //  
    RCP<Belos::TpetraOperator<ST>> innerSolver;
    innerSolver = rcp( new Belos::TpetraOperator<ST>( rcpFromRef(innerProblem), rcpFromRef(innerBelosList), "Minres", true ) );

    //
    // Construct a linear problem instance with GMRES as preconditoner.
    //
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    problem.setInitResVec(B);
    problem.setRightPrec(innerSolver);
    problem.setLabel( "Belos Flexible Gmres Solve" );
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // *******************************************************************
    // *************Start the block Gmres iteration***********************
    // *******************************************************************
    //
    Belos::BlockGmresSolMgr<ST,MV,OP> solver( rcpFromRef(problem), rcpFromRef(belosList) );
    
    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of Minres iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = solver.solve();
    //
    // Compute actual residuals.
    //
    bool badRes = false;
    std::vector<MT> actualResids( numrhs );
    std::vector<MT> rhsNorm( numrhs );
    MV resid(map, numrhs);
    OPT::Apply( *A, *X, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *B, resid );
    MVT::MvNorm( resid, actualResids );
    MVT::MvNorm( *B, rhsNorm );
    if (proc_verbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        double actRes = actualResids[i]/rhsNorm[i];
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

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
  // run with different scalar types
  return run<double>(argc, argv);
  // return run<float>(argc, argv);
}
