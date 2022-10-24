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
// The right-hand-side corresponds to a randomly generated solution.
// The initial guesses are all set to zero.
//
// NOTE: No preconditioner is used in this case.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosTpetraTestFramework.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

using namespace Teuchos;
using Tpetra::Operator;
using Tpetra::CrsMatrix;
using Tpetra::MultiVector;
using std::endl;
using std::cout;
using std::vector;
using Teuchos::tuple;

int main(int argc, char *argv[]) {
  typedef Tpetra::MultiVector<>::scalar_type BSC;
  typedef std::complex<BSC>                ST;
  typedef ScalarTraits<ST>                SCT;
  typedef SCT::magnitudeType               MT;
  typedef Tpetra::Operator<ST>             OP;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;

  GlobalMPISession mpisess(&argc,&argv,&cout);

  const ST one  = SCT::one();

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = rank(*comm);

  //
  // Get test parameters from command-line processor
  //
  bool verbose = false, proc_verbose = false, debug = false;
  int frequency = -1;  // how often residuals are printed by solver
  int numrhs = 1;      // total number of right-hand sides to solve for
  std::string filename("mhd1280b.cua");
  MT tol = 1.0e-5;     // relative residual tolerance

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
  cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
  cmdp.setOption("tol",&tol,"Relative residual tolerance used by GCRODR solver.");
  cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
  cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if (debug) {
    verbose = true;
  }
  if (!verbose) {
    frequency = -1;  // reset frequency if test is not verbose
  }

  proc_verbose = ( verbose && (MyPID==0) );

  if (proc_verbose) {
    std::cout << Belos::Belos_Version() << std::endl << std::endl;
  }

  Belos::Tpetra::HarwellBoeingReader<Tpetra::CrsMatrix<ST> > reader( comm );
  RCP<Tpetra::CrsMatrix<ST> > A = reader.readFromFile( filename );
  RCP<const Tpetra::Map<> > map = A->getDomainMap();

  // Create initial vectors
  RCP<MV> B, X;
  X = rcp( new MV(map,numrhs) );
  MVT::MvRandom( *X );
  B = rcp( new MV(map,numrhs) );
  OPT::Apply( *A, *X, *B );
  MVT::MvInit( *X, 0.0 );

  //
  // ********Other information used by block solver***********
  // *****************(can be user specified)******************
  //
  const int NumGlobalElements = B->getGlobalLength();
  int numIters1, numIters2, numIters3;
  int maxits = NumGlobalElements; // maximum number of iterations to run
  int numBlocks = 100;
  int numRecycledBlocks = 20;
  ParameterList belosList;
  belosList.set( "Maximum Iterations", maxits );         // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
  belosList.set( "Num Blocks", numBlocks );
  belosList.set( "Num Recycled Blocks", numRecycledBlocks );

  int verbLevel = Belos::Errors + Belos::Warnings;
  if (debug) {
    verbLevel += Belos::Debug;
  }
  if (verbose) {
    verbLevel += Belos::TimingDetails + Belos::FinalSummary + Belos::StatusTestDetails;
  }
  belosList.set( "Verbosity", verbLevel );
  if (verbose) {
    if (frequency > 0) {
      belosList.set( "Output Frequency", frequency );
    }
  }
  //
  // Construct an unpreconditioned linear problem instance.
  //
  Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
  bool set = problem.setProblem();
  if (set == false) {
    if (proc_verbose)
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
  }
  //
  // *******************************************************************
  // ******************Start the GCRODR iteration***********************
  // *******************************************************************
  //
  Belos::GCRODRSolMgr<ST,MV,OP> solver( rcpFromRef(problem), rcpFromRef(belosList) );

  //
  // **********Print out information about problem*******************
  //
  if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Max number of GCRODR iterations: " << maxits << std::endl;
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
  }
  //
  // Perform solve
  //
  Belos::ReturnType ret = solver.solve();
  numIters1=solver.getNumIters();
  //
  // Compute actual residuals.
  //
  bool badRes = false;
  std::vector<MT> actual_resids( numrhs );
  std::vector<MT> rhs_norm( numrhs );
  MV resid(map, numrhs);
  OPT::Apply( *A, *X, resid );
  MVT::MvAddMv( -one, resid, one, *B, resid );
  MVT::MvNorm( resid, actual_resids );
  MVT::MvNorm( *B, rhs_norm );
  if (proc_verbose) {
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
  }
  for ( int i=0; i<numrhs; i++) {
    MT actRes = actual_resids[i]/rhs_norm[i];
    if (proc_verbose) {
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
    }
    if (actRes > tol) badRes = true;
  }

  if (proc_verbose) { std::cout << "First solve took " << numIters1 << " iterations." << std::endl; }

  // Resolve linear system with same rhs and recycled space
  MVT::MvInit( *X, SCT::zero() );
  solver.reset(Belos::Problem);
  ret = solver.solve();
  numIters2=solver.getNumIters();

  if (proc_verbose) { std::cout << "Second solve took " << numIters2 << " iterations." << std::endl; }

  // Resolve linear system (again) with same rhs and recycled space
  MVT::MvInit( *X, SCT::zero() );
  solver.reset(Belos::Problem);
  ret = solver.solve();
  numIters3=solver.getNumIters();

  if (proc_verbose) { std::cout << "Third solve took " << numIters3 << " iterations." << std::endl; }

  if ( ret!=Belos::Converged || badRes || numIters1 < numIters2 || numIters2 < numIters3 ) {
    if (proc_verbose) {
      std::cout << "\nEnd Result: TEST FAILED" << std::endl;
    }
    return -1;
  }
  //
  // Default return value
  //
  if (proc_verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
  //
} // end test_gcrodr_complex_hb.cpp
