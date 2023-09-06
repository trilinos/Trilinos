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
// from the problem, if it exists, will be used instead of multiple random
// right-hand-sides.  The initial guesses are all set to zero.  An ILU preconditioner
// is constructed using the Ifpack factory.
//

// Ifpack2
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Parameters.hpp>
#include <Ifpack2_Preconditioner.hpp>

// Teuchos
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Tpetra
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLSQRSolMgr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosTpetraTestFramework.hpp"

template <typename ScalarType>
int run(int argc, char *argv[]) {
  using Teuchos::CommandLineProcessor;
  using Teuchos::GlobalMPISession;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::tuple;

  using ST = typename Tpetra::MultiVector<ScalarType>::scalar_type;
  using LO = typename Tpetra::Vector<>::local_ordinal_type;
  using GO = typename Tpetra::Vector<>::global_ordinal_type;
  using NT = typename Tpetra::Vector<>::node_type;

  using V = typename Tpetra::Vector<ST, LO, GO, NT>;
  using MV = typename Tpetra::MultiVector<ST, LO, GO, NT>;
  using OP = typename Tpetra::Operator<ST, LO, GO, NT>;
  using MAP = typename Tpetra::Map<LO, GO, NT>;
  using MAT = typename Tpetra::CrsMatrix<ST, LO, GO, NT>;

  using MVT = typename Belos::MultiVecTraits<ST, MV>;
  using OPT = typename Belos::OperatorTraits<ST, MV, OP>;

  using MT = typename Teuchos::ScalarTraits<ST>::magnitudeType;

  using LinearProblem = typename Belos::LinearProblem<ST, MV, OP>;
  using Preconditioner = typename Ifpack2::Preconditioner<ST, LO, GO, NT>;
  using Solver = ::Belos::LSQRSolMgr<ST, MV, OP>;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

  bool verbose = false;
  bool success = true;

  try {
    bool procVerbose = false;
    bool leftprec = true;  // left preconditioning or right.
    // LSQR applies the operator and the transposed operator.
    // A preconditioner must support transpose multiply.
    int frequency = -1;  // frequency of status test output.
    int blocksize = 1;   // blocksize
    // LSQR as currently implemented is a single vector algorithm.
    // However some of the parameters that would be used by a block version
    // have not been removed from this file.
    int numRHS = 1;     // number of right-hand sides to solve for
    int maxiters = -1;  // maximum number of iterations allowed per linear system
    std::string filename("orsirr1_scaled.hb");
    MT relResTol = 3.0e-4;    // relative residual tolerance for the preconditioned linear system
    MT resGrowthFactor = 4.;  // In this example, warn if |resid| > resGrowthFactor * relResTol

    MT relMatTol = 1.e-4;  // relative Matrix error, default value sqrt(eps)
    MT maxCond = 1.e+8;    // maximum condition number default value 1/eps
    MT damp = 0.;          // regularization (or damping) parameter

    Teuchos::CommandLineProcessor cmdp(false, true);
    cmdp.setOption("verbose", "quiet", &verbose, "Print messages and results.");
    cmdp.setOption("left-prec", "right-prec", &leftprec, "Left preconditioning or right.");
    cmdp.setOption("frequency", &frequency, "Solvers frequency for printing residuals (#iters).");
    cmdp.setOption(
        "filename", &filename,
        "Filename for test matrix.  Acceptable file extensions: *.hb,*.mtx,*.triU,*.triS");
    cmdp.setOption("lambda", &damp, "Regularization parameter");
    cmdp.setOption("tol", &relResTol, "Relative residual tolerance");
    cmdp.setOption("matrixTol", &relMatTol, "Relative error in Matrix");
    cmdp.setOption("max-cond", &maxCond, "Maximum condition number");
    cmdp.setOption("num-rhs", &numRHS, "Number of right-hand sides to be solved for.");
    cmdp.setOption("block-size", &blocksize, "Block size used by LSQR.");
    cmdp.setOption(
        "max-iters", &maxiters,
        "Maximum number of iterations per linear system (-1 = adapted to problem/block size).");
    if (cmdp.parse(argc, argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }
    if (!verbose)
      frequency = -1;  // reset frequency if test is not verbose

    // Get the problem
    Belos::Tpetra::HarwellBoeingReader<MAT> reader(comm);
    RCP<MAT> A = reader.readFromFile(filename);
    RCP<const MAP> map = A->getDomainMap();

    // Initialize vectors
    RCP<MV> vecB = rcp(new MV(map, numRHS));
    RCP<MV> vecX = rcp(new MV(map, numRHS));
    RCP<MV> B, X;

    procVerbose = verbose && (comm->getRank() == 0);  // Only print on the zero processor

    // Check to see if the number of right-hand sides is the same as requested.
    if (numRHS > 1) {
      // numRHS > 1 not yet supported
      X = rcp(new MV(map, numRHS));
      B = rcp(new MV(map, numRHS));
      X->randomize();
      OPT::Apply(*A, *X, *B);  // B := AX
      X->putScalar(0.0);       // annihilate X
    } else {
      LO locNumCol = map->getMaxLocalIndex() + 1;  // Create a known solution
      GO globNumCol = map->getMaxGlobalIndex() + 1;
      for (LO li = 0; li < locNumCol; li++) {  // assume consecutive lid
        const auto gid = map->getGlobalElement(li);
        ST value = (ST)(globNumCol - 1 - gid);
        int numEntries = 1;
        vecX->replaceGlobalValue(numEntries, 0, value);
      }
      A->apply(*vecX, *vecB);  // Create a consistent linear system

      // At this point, the initial guess is exact.
      bool goodInitGuess = true;   // initial guess near solution
      bool zeroInitGuess = false;  // annihilate initial guess
      if (goodInitGuess) {
        ST value = 1.e-2;  // "Rel RHS Err" and "Rel Mat Err" apply to the residual equation,
        // LO numEntries = 1;   // norm( b - A x_k ) ?<? relResTol norm( b- Axo).
        LO index = 0;  // norm(b) is inaccessible to LSQR.
        vecX->sumIntoLocalValue(index, 0, value);
      }

      if (zeroInitGuess) {
        vecX->putScalar(0.0);
      }

      X = vecX;
      B = vecB;
    }

    // Construct preconditioner
    ParameterList precParams;

    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    std::string precType = "ILUT";  // incomplete LU
    int overlapLevel = 1;           // nonnegative

    RCP<Preconditioner> prec = Ifpack2::Factory::create<MAT>(precType, A, overlapLevel);
    assert(prec != Teuchos::null);

    // specify parameters for ILU
    precParams.set("fact: level-of-fill", 1);
    // the combine mode is on the following:
    // "ADD", "INSERT", "REPLACE", "ABSMAX", and "ZERO"
    // Their meaning is as defined in Tpetra::CombineMode
    precParams.set("schwarz: combine mode", "ADD");
    // sets the parameters
    prec->setParameters(precParams);

    // initialize the preconditioner. At this point the matrix must have been FillComplete()'d, but
    // actual values are ignored.
    prec->initialize();

    // Builds the preconditioners, by looking for the values of the matrix.
    prec->compute();

    // Create parameter list for the LSQR solver manager
    const int numGlobalElements = B->getGlobalLength();
    if (maxiters == -1)
      maxiters = numGlobalElements / blocksize - 1;
    RCP<ParameterList> belosList = rcp(new ParameterList());
    belosList->set("Block Size", blocksize);         // Blocksize to be used by iterative solver
    belosList->set("Lambda", damp);                  // Regularization parameter
    belosList->set("Rel RHS Err", relResTol);        // Relative convergence tolerance requested
    belosList->set("Rel Mat Err", relMatTol);        // Maximum number of restarts allowed
    belosList->set("Condition Limit", maxCond);      // upper bound for cond(A)
    belosList->set("Maximum Iterations", maxiters);  // Maximum number of iterations allowed
    if (numRHS > 1) {
      belosList->set("Show Maximum Residual Norm Only",
                     true);  // Show only the maximum residual norm
    }
    if (verbose) {
      belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails +
                                      Belos::StatusTestDetails);
      if (frequency > 0)
        belosList->set("Output Frequency", frequency);
    } else
      belosList->set("Verbosity", Belos::Errors + Belos::Warnings);

    // Construct a preconditioned linear problem
    RCP<LinearProblem> problem = rcp(new LinearProblem(A, X, B));
    if (leftprec) {
      problem->setLeftPrec(prec);
    } else {
      problem->setRightPrec(prec);
    }
    bool set = problem->setProblem();
    if (set == false) {
      if (procVerbose)
        std::cout << std::endl
                  << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return -1;
    }

    // Create an iterative solver manager.
    RCP<Solver> solver = rcp(new Solver(problem, belosList));

    // Start the LSQR iteration
    if (procVerbose) {
      std::cout << std::endl << std::endl;
      std::cout << "Dimension of matrix: " << numGlobalElements << std::endl;
      std::cout << "Number of right-hand sides: " << numRHS << std::endl;
      std::cout << "Block size used by solver: " << blocksize << std::endl;
      std::cout << "Max number of Gmres iterations per restart cycle: " << maxiters << std::endl;
      std::cout << "Relative residual tolerance: " << relResTol << std::endl;
      std::cout << std::endl;
      std::cout << "Solver's Description: " << std::endl;
      std::cout << solver->description() << std::endl;  // visually verify the parameter list
    }

    // Perform solve
    Belos::ReturnType ret = solver->solve();

    std::vector<ST> solNorm(numRHS);  // get solution norm
    MVT::MvNorm(*X, solNorm);
    int numIters = solver->getNumIters();  // get number of solver iterations
    MT condNum = solver->getMatCondNum();
    MT matrixNorm = solver->getMatNorm();
    MT resNorm = solver->getResNorm();
    MT lsResNorm = solver->getMatResNorm();

    if (procVerbose)
      std::cout << "Number of iterations performed for this solve: " << numIters << std::endl
                << "matrix condition number: " << condNum << std::endl
                << "matrix norm: " << matrixNorm << std::endl
                << "residual norm: " << resNorm << std::endl
                << "solution norm: " << solNorm[0] << std::endl
                << "least squares residual Norm: " << lsResNorm << std::endl;

    // Compute actual residuals.
    bool badRes = false;
    std::vector<ST> actualResids(numRHS);
    std::vector<ST> rhsNorm(numRHS);
    MV resid(map, numRHS);
    OPT::Apply(*A, *X, resid);
    MVT::MvAddMv(-1.0, resid, 1.0, *B, resid);
    MVT::MvNorm(resid, actualResids);
    MVT::MvNorm(*B, rhsNorm);
    if (procVerbose) {
      std::cout << "---------- Actual Residuals (normalized) ----------" << std::endl << std::endl;
      for (int i = 0; i < numRHS; i++) {
        ST actRes = actualResids[i] / rhsNorm[i];
        std::cout << "Problem " << i << " : \t" << actRes << std::endl;
        if (actRes > relResTol * resGrowthFactor) {
          badRes = true;
          if (verbose)
            std::cout << "residual norm > " << relResTol * resGrowthFactor << std::endl;
        }
      }
    }

    if (ret != Belos::Converged || badRes) {
      success = false;
      if (procVerbose)
        std::cout << std::endl << "ERROR:  Belos did not converge!" << std::endl;
    } else {
      success = true;
      if (procVerbose)
        std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(int argc, char *argv[]) {
  return run<double>(argc, argv);
  // return run<float>(argc, argv);
}  // end PrecLSQRTpetraExFile.cpp
