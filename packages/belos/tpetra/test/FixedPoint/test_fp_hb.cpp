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
#include "BelosFixedPointSolMgr.hpp"

// I/O for Harwell-Boeing files
#define HIDE_TPETRA_INOUT_IMPLEMENTATIONS
#include <Tpetra_MatrixIO.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
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

  typedef double                           ST;
  typedef ScalarTraits<ST>                SCT;
  typedef SCT::magnitudeType               MT;
  typedef Tpetra::Operator<ST>             OP;
  typedef Tpetra::MultiVector<ST>          MV;
  typedef Tpetra::Vector<ST>               VV;
  typedef Belos::OperatorTraits<ST,MV,OP> OPT;
  typedef Belos::MultiVecTraits<ST,MV>    MVT;
  using GO = MV::global_ordinal_type;

  GlobalMPISession mpisess(&argc,&argv,&cout);

  bool success = false;
  bool verbose = false;
  try {
    const ST one  = SCT::one();

    int MyPID = 0;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    //
    // Get test parameters from command-line processor
    //
    bool proc_verbose = false;
    bool debug = false;
    int frequency = -1;  // how often residuals are printed by solver
    int numrhs = 1;      // total number of right-hand sides to solve for
    int maxiters = -1;   // maximum number of iterations for solver to use
    std::string filename("bcsstk14.hb");
    MT tol = 1.0e-5;     // relative residual tolerance
    bool precond = false; // use diagonal preconditioner
    bool leftPrecond = false; // if preconditioner is used, left or right?

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");
    cmdp.setOption("frequency",&frequency,"Solvers frequency for printing residuals (#iters).");
    cmdp.setOption("tol",&tol,"Relative residual tolerance used by fixed point solver.");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("num-rhs",&numrhs,"Number of right-hand sides to be solved for.");
    cmdp.setOption("max-iters",&maxiters,"Maximum number of iterations per linear system (-1 := adapted to problem/block size).");
    cmdp.setOption("use-precond","no-precond",&precond,"Use a diagonal preconditioner.");
    cmdp.setOption("left","right",&leftPrecond,"Use a left/right preconditioner.");
    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (debug) {
      verbose = true;
    }
    if (!verbose) {
      frequency = -1;  // reset frequency if test is not verbose
    }

    MyPID = rank(*comm);
    proc_verbose = ( verbose && (MyPID==0) );

    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }

    //
    // Get the data from the HB file and build the Map,Matrix
    //
    RCP<CrsMatrix<ST> > A;
    Tpetra::Utils::readHBMatrix(filename,comm,A);
    RCP<const Tpetra::Map<> > map = A->getDomainMap();

    // Create initial vectors
    RCP<MV> B, X;
    X = rcp( new MV(map,numrhs) );
    MVT::MvRandom( *X );
    B = rcp( new MV(map,numrhs) );
    OPT::Apply( *A, *X, *B );
    MVT::MvInit( *X, 0.0 );


    // Scale down the rhs
    std::vector<ST> norm( MVT::GetNumberVecs(*B));
    MVT::MvNorm(*B,norm);
    for(int i=0; i< MVT::GetNumberVecs(*B); i++)
      norm[i] = 1.0 / norm[i];
    MVT::MvScale(*B,norm);

    // Drastically bump the diagonal and then rescale so FP can converge
    {
      VV diag(A->getRowMap());
      A->getLocalDiagCopy(diag);
      {
        Teuchos::ArrayRCP<ST> dd=diag.getDataNonConst();

        auto GlobalElements = A->getRowMap()->getNodeElementList();
        A->resumeFill();
        for(int i=0; i<(int)dd.size(); i++) {
          dd[i]=dd[i]*1e4;
          A->replaceGlobalValues(GlobalElements[i],
                                 Teuchos::tuple<GO>(GlobalElements[i]),
                                 Teuchos::tuple<ST>(dd[i]) );
        }
        A->fillComplete();

        for(int i=0; i<(int)dd.size(); i++)
          dd[i] = 1.0/sqrt(dd[i]);
      }
      A->leftScale(diag);
      A->rightScale(diag);
    }

    //
    // ********Other information used by block solver***********
    // *****************(can be user specified)******************
    //
    const int NumGlobalElements = B->getGlobalLength();
    if (maxiters == -1) {
      maxiters = NumGlobalElements- 1; // maximum number of iterations to run
    }
    //
    ParameterList belosList;
    belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
    belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
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
    // Construct an linear problem instance.
    //
    Belos::LinearProblem<ST,MV,OP> problem( A, X, B );
    // diagonal preconditioner
    if (precond) {
      VV diagonal(A->getRowMap());
      A->getLocalDiagCopy(diagonal);

      int NumMyElements    = diagonal.getMap()->getNodeNumElements();
      auto MyGlobalElements = diagonal.getMap()->getNodeElementList();
      Teuchos::ArrayRCP<ST> dd=diagonal.getDataNonConst();
      RCP<CrsMatrix<ST> > invDiagMatrix = Teuchos::rcp(new CrsMatrix<ST>(A->getRowMap(), 1, Tpetra::StaticProfile));

      for (Teuchos_Ordinal i=0; i<NumMyElements; ++i) {
        invDiagMatrix->insertGlobalValues(MyGlobalElements[i],
                                          Teuchos::tuple<GO>(MyGlobalElements[i]),
                                          Teuchos::tuple<ST>(SCT::one() / dd[i]) );
      }
      invDiagMatrix->fillComplete();

      if (leftPrecond)
        problem.setLeftPrec(invDiagMatrix);
      else
        problem.setRightPrec(invDiagMatrix);
    }
    bool set = problem.setProblem();
    if (set == false) {
      if (proc_verbose)
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }
    //
    // *******************************************************************
    // *************Start the fixed point iteration***********************
    // *******************************************************************
    //
    Belos::FixedPointSolMgr<ST,MV,OP> solver( rcpFromRef(problem), rcpFromRef(belosList) );

    //
    // **********Print out information about problem*******************
    //
    if (proc_verbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numrhs << std::endl;
      std::cout << "Max number of fixed point iterations: " << maxiters << std::endl;
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

    success = (ret==Belos::Converged && !badRes);

    if (success) {
      if (proc_verbose)
        std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    } else {
      if (proc_verbose)
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // end test_fp_hb.cpp
