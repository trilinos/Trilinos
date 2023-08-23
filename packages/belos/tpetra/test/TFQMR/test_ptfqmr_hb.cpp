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

// Tpetra
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>

// Teuchos
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_Assert.hpp>

// Ifpack2
#include <Ifpack2_IlukGraph.hpp>
//#include <Ifpack_CrsRiluk.h> // AM: todo, finds an equivalent

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosPseudoBlockTFQMRSolMgr.hpp"
#include "BelosTpetraTestFramework.hpp"

template<typename ScalarType>
int run(int argc, char *argv[]) {
  // Teuchos
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Tpetra
  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::MultiVector<>::local_ordinal_type;
  using GO = typename Tpetra::MultiVector<>::global_ordinal_type;
  using NT = typename Tpetra::MultiVector<>::node_type;

  using OP = Tpetra::Operator<ST,LO,GO,NT>;
  using MV = Tpetra::MultiVector<ST,LO,GO,NT>;

  using tcrsmatrix_t = Tpetra::CrsMatrix<ST,LO,GO,NT>;
  using tmap_t = Tpetra::Map<LO,GO,NT>;

  // Belos
  using OPT = typename Belos::OperatorTraits<ST,MV,OP>;
  using MVT = typename Belos::MultiVecTraits<ST,MV>;

  // MPI session
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);
  
  bool success = false;
  bool verbose = false;
  
  try {
    bool procVerbose = false;
    bool leftprec = true; // use left preconditioning to solve these linear systems
    bool explicitTest = true;
    bool pseudo = false;
    int frequency = -1;  // how often residuals are printed by solver
    int numrhs = 1;
    int maxiters = -1;    // maximum iterations allowed
    std::string filename("orsirr1.hb");
    ST tol = sqrt(std::numeric_limits<ST>::epsilon()); // relative residual tolerance

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("left-prec","right-prec",&leftprec,"Left preconditioning or right.");
    cmdp.setOption("explicit","implicit-only",&explicitTest,"Compute explicit residuals.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("pseudo","not-pseudo",&pseudo,"Use pseudo-block TFQMR solver.");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by GMRES solver.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("maxiters",&maxiters,"Maximum number of iterations per linear system (-1 = adapted to problem size).");
    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose
    
    // Get the problem
    Belos::Tpetra::HarwellBoeingReader<tcrsmatrix_t> reader( comm );
    RCP<tcrsmatrix_t> A = reader.readFromFile( filename );
    RCP<const tmap_t> map = A->getDomainMap();
    procVerbose = verbose && (MyPID==0); /* Only print on zero processor */
    
    // *****Construct the Preconditioner*****
    if (procVerbose) std::cout << std::endl << std::endl;
    if (procVerbose) std::cout << "Constructing ILU preconditioner" << std::endl;
    int Lfill = 2;
    // if (argc > 2) Lfill = atoi(argv[2]);
    if (procVerbose) std::cout << "Using Lfill = " << Lfill << std::endl;
    int Overlap = 2;
    // if (argc > 3) Overlap = atoi(argv[3]);
    if (procVerbose) std::cout << "Using Level Overlap = " << Overlap << std::endl;
    double Athresh = 0.0;
    // if (argc > 4) Athresh = atof(argv[4]);
    if (procVerbose) std::cout << "Using Absolute Threshold Value of " << Athresh << std::endl;
    double Rthresh = 1.0;
    // if (argc >5) Rthresh = atof(argv[5]);
    if (procVerbose) std::cout << "Using Relative Threshold Value of " << Rthresh << std::endl;
    
    // Init Ifpack2 classes
    RCP<Ifpack2::IlukGraph> ilukGraph;
    // AM: TODO, IFPACK ISSUE
    //RCP<Ifpack::CrsRiluk> ilukFactors;
    
    // AM: TODO, IFPACK ISSUE 
    if (Lfill > -1) {
      ilukGraph = rcp(new Ifpack2::IlukGraph(A->getGraph(), Lfill, Overlap));
      /*int info = ilukGraph->ConstructFilledGraph();
      TEUCHOS_ASSERT( info == 0 );
      ilukFactors = rcp(new Ifpack_CrsRiluk(*ilukGraph));
      int initerr = ilukFactors->InitValues(*A);
      if (initerr != 0) std::cout << "InitValues error = " << initerr;
      info = ilukFactors->Factor();
      TEUCHOS_ASSERT( info == 0 );*/
    }
    
    // AM: FROM HERE, CONVERSION TO DO BELOW 

    //
    bool transA = false;
    double Cond_Est;
    ilukFactors->Condest(transA, Cond_Est);
    if (procVerbose) {
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
    belosList.set( "Explicit Residual Test", explicit_test );       // Need to check for the explicit residual before returning
    if (verbose) {
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
          Belos::TimingDetails + Belos::StatusTestDetails );
      if (frequency > 0)
        belosList.set( "Output Frequency", frequency );
    }
    else
      belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
    
    // *****Construct solution std::vector and random right-hand-sides *****
    RCP<MV> X = rcp( new MV(Map, numrhs) );
    RCP<MV> B = rcp( new MV(Map, numrhs) );
    MVT::MvRandom( *X );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );

    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    if (leftprec)
      problem.setLeftPrec( Prec );
    else
      problem.setRightPrec( Prec );

    bool set = problem.setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    
    // Start the TFQMR or Pseudo Block TFQMR iteration
    RCP< Belos::SolverManager<ST,MV,OP> > solver;
    if (pseudo)
      solver = rcp( new Belos::PseudoBlockTFQMRSolMgr<ST,MV,OP>(rcp(&problem,false), rcp(&belosList,false)) );
    else
      solver = rcp( new Belos::TFQMRSolMgr<ST,MV,OP>(rcp(&problem, false), rcp(&belosList, false)));
    
    // Print out information about problem
    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of Pseudo Block TFQMR iterations: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << tol << std::endl;
      std::cout << std::endl;
    }
    
    // Perform solve
    Belos::ReturnType ret = solver->solve();

    // Compute actual residuals.
    bool badRes = false;
    std::vector<ST> actualResids( numrhs );
    std::vector<ST> rhsNorm( numrhs );
    MV R(Map, numrhs);
    OPT::Apply( *A, *X, R );
    MVT::MvAddMv( -1.0, R, 1.0, *B, R );
    MVT::MvNorm( R, actualResids );
    MVT::MvNorm( *B, rhsNorm );
    if (procVerbose) {
      std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        ST actRes = actualResids[i]/rhsNorm[i];
        std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol ) badRes = true;
      }
    }
    
    success = ret==Belos::Converged && !badRes;

    if (success) {
      if (procVerbose)
        std::cout << "End Result: TEST PASSED" << std::endl;
    } else {
      if (procVerbose)
        std::cout << "End Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

int main(int argc, char *argv[]) {
  return run<double>(argc,argv);

  // wrapped with a check: CMake option Trilinos_ENABLE_FLOAT=ON
  // run<float>(argc,argv);
}
