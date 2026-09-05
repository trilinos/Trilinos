// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver tests Belos::SolutionProjection using a complex
// Harwell-Boeing matrix.
//
// The test does the following:
//
//   1. Generate two historical exact solutions u1Exact and u2Exact.
//   2. Form b1 = A u1Exact and b2 = A u2Exact.
//   3. Solve these two systems with one persistent GCRODR solver, so that
//      GCRODR accumulates a recycle space.
//   4. Add the computed historical solutions to Belos::SolutionProjection.
//   5. Construct a target exact solution mostly in the span of the computed
//      historical solutions.
//   6. Compare:
//        a. trained GCRODR using its internal recycle space,
//        b. fresh GCRODR using a SolutionProjection initial guess.
//
// The projection utility implements Fischer's Method 1:
//
//      A U = C, C^H C = I,
//      x <- x + U C^H r,
//      r <- r - C C^H r.
//
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosGCRODRSolMgr.hpp"
#include "BelosSolutionProjection.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// I/O for Harwell-Boeing files
#include "Tpetra_Util_iohb.h"

#include "MyMultiVec.hpp"
#include "MyBetterOperator.hpp"
#include "MyOperator.hpp"

#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>

#include <stdexcept>
#include <vector>

using namespace Teuchos;

int main(int argc, char *argv[]) {
  typedef std::complex<double>              ST;
  typedef ScalarTraits<ST>                  SCT;
  typedef SCT::magnitudeType                MT;
  typedef Belos::MultiVec<ST>               MV;
  typedef Belos::Operator<ST>               OP;
  typedef Belos::MultiVecTraits<ST,MV>      MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>   OPT;

  const ST one  = SCT::one();
  const ST zero = SCT::zero();

  bool success = false;
  bool verbose = false;
  bool debug = false;

  int info = 0;

  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  const int MyPID = session.getRank();

  using Teuchos::RCP;
  using Teuchos::rcp;

  try {
    bool proc_verbose = false;
    bool requireProjectionBeatsRecycled = false;

    std::string filename("mhd1280b.cua");

    MT solverTol = 1.0e-6;
    MT projectionReductionTol = 5.0e-2;
    MT epsilon = 1.0e-2;

    int maxBasisSize = 4;
    int numBlocks = 100;
    int numRecycledBlocks = 20;
    int maxIters = 1000;
    int minIterationImprovement = 1;

    CommandLineProcessor cmdp(false,true);
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("debug","no-debug",&debug,"Print debug messages.");
    cmdp.setOption("filename",&filename,"Filename for Harwell-Boeing test matrix.");
    cmdp.setOption("solver-tol",&solverTol,
                   "Relative residual tolerance used by GCRODR solver.");
    cmdp.setOption("projection-reduction-tol",&projectionReductionTol,
                   "Required residual reduction from SolutionProjection.");
    cmdp.setOption("epsilon",&epsilon,
                   "Size of component outside the stored solution projection space.");
    cmdp.setOption("max-basis-size",&maxBasisSize,
                   "Maximum number of vectors stored by SolutionProjection.");
    cmdp.setOption("num-blocks",&numBlocks,"Number of GCRODR blocks.");
    cmdp.setOption("num-recycled-blocks",&numRecycledBlocks,
                   "Number of GCRODR recycled blocks.");
    cmdp.setOption("max-iters",&maxIters,"Maximum GCRODR iterations.");
    cmdp.setOption("require-projection-beats-recycled",
                   "no-require-projection-beats-recycled",
                   &requireProjectionBeatsRecycled,
                   "Require fresh GCRODR with SolutionProjection initial guess "
                   "to use fewer iterations than trained recycled GCRODR.");
    cmdp.setOption("min-iteration-improvement",&minIterationImprovement,
                   "Minimum required iteration reduction when "
                   "--require-projection-beats-recycled is enabled.");

    if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
      return EXIT_FAILURE;
    }

    proc_verbose = verbose && (MyPID == 0);

    if (proc_verbose) {
      std::cout << Belos::Belos_Version() << std::endl << std::endl;
    }

    // ------------------------------------------------------------------
    // Read HB matrix.
    // ------------------------------------------------------------------

    int dim = 0;
    int dim2 = 0;
    int nnz = -1;

    MT *dvals = NULL;
    int *colptr = NULL;
    int *rowind = NULL;
    ST *cvals = NULL;

    info = Tpetra::HB::readHB_newmat_double(filename.c_str(), &dim, &dim2, &nnz,
                                            &colptr, &rowind, &dvals);

    if (info == 0 || nnz < 0) {
      if (MyPID == 0) {
        std::cout << "Error reading '" << filename << "'" << std::endl;
        std::cout << "End Result: TEST FAILED" << std::endl;
      }
      return EXIT_FAILURE;
    }

    // Convert interleaved doubles to std::complex values.
    cvals = new ST[nnz];
    for (int ii = 0; ii < nnz; ++ii) {
      cvals[ii] = ST(dvals[ii*2], dvals[ii*2+1]);
    }

    RCP<MyBetterOperator<ST> > A =
      rcp(new MyBetterOperator<ST>(dim, colptr, nnz, rowind, cvals));

    // ------------------------------------------------------------------
    // Solver parameter list.
    // ------------------------------------------------------------------

    ParameterList belosList;
    belosList.set("Maximum Iterations", maxIters);
    belosList.set("Convergence Tolerance", solverTol);
    belosList.set("Num Blocks", numBlocks);
    belosList.set("Num Recycled Blocks", numRecycledBlocks);

    // Use RHS scaling because SolutionProjection changes the initial guess
    // before LinearProblem::setProblem().
    belosList.set("Implicit Residual Scaling", "Norm of RHS");
    belosList.set("Explicit Residual Scaling", "Norm of RHS");

    if (verbose) {
      int verbosity = Belos::Errors + Belos::Warnings +
        Belos::TimingDetails + Belos::StatusTestDetails;
      if (debug) {
        verbosity += Belos::OrthoDetails + Belos::Debug;
      }
      belosList.set("Verbosity", verbosity);
    }
    else {
      belosList.set("Verbosity", Belos::Errors + Belos::Warnings);
    }

    if (proc_verbose) {
      std::cout << "Dimension of matrix: " << dim << std::endl;
      std::cout << "GCRODR Num Blocks: " << numBlocks << std::endl;
      std::cout << "GCRODR Num Recycled Blocks: " << numRecycledBlocks << std::endl;
      std::cout << "GCRODR tolerance: " << solverTol << std::endl;
      std::cout << "Projection epsilon: " << epsilon << std::endl;
      std::cout << std::endl;
    }

    // ------------------------------------------------------------------
    // Helper for forming b = A x.
    // ------------------------------------------------------------------

    auto makeRHS =
      [&](const RCP<MyMultiVec<ST> >& xExact) -> RCP<MyMultiVec<ST> >
      {
        RCP<MyMultiVec<ST> > b = rcp(new MyMultiVec<ST>(dim, 1));
        OPT::Apply(*A, *xExact, *b);
        return b;
      };

    // ------------------------------------------------------------------
    // Helper for computing ||b - A x|| / ||b||.
    // ------------------------------------------------------------------

    auto computeRelativeResidual =
      [&](const RCP<MyMultiVec<ST> >& x,
          const RCP<MyMultiVec<ST> >& b) -> MT
      {
        RCP<MyMultiVec<ST> > residual = rcp(new MyMultiVec<ST>(dim, 1));

        OPT::Apply(*A, *x, *residual);
        MVT::MvAddMv(one, *b, -one, *residual, *residual);

        std::vector<MT> residualNorm(1);
        std::vector<MT> rhsNorm(1);

        MVT::MvNorm(*residual, residualNorm);
        MVT::MvNorm(*b, rhsNorm);

        return residualNorm[0] / rhsNorm[0];
      };

    // ------------------------------------------------------------------
    // Helper result type.
    // ------------------------------------------------------------------

    struct SolveResult {
      int iters;
      MT finalRelResidual;
      Belos::ReturnType ret;
      RCP<MyMultiVec<ST> > solution;
    };

    // ------------------------------------------------------------------
    // Helper for solving with an existing GCRODR solver and LinearProblem.
    //
    // This is used to preserve GCRODR's recycle space across historical
    // solves and the target solve.
    // ------------------------------------------------------------------

    auto solveWithExistingSolver =
      [&](Belos::GCRODRSolMgr<ST,MV,OP>& solver,
          const RCP<Belos::LinearProblem<ST,MV,OP> >& problem,
          const RCP<MyMultiVec<ST> >& x,
          const RCP<MyMultiVec<ST> >& b,
          const std::string& label) -> SolveResult
      {
        const bool set = problem->setProblem(x, b);

        TEUCHOS_TEST_FOR_EXCEPTION(
          !set,
          std::runtime_error,
          "Belos::LinearProblem failed to set up correctly.");

        Belos::ReturnType ret = solver.solve();

        SolveResult result;
        result.iters = solver.getNumIters();
        result.finalRelResidual = computeRelativeResidual(x, b);
        result.ret = ret;
        result.solution = x;

        if (proc_verbose) {
          std::cout << label << " iterations: "
                    << result.iters << std::endl;
          std::cout << label << " final relative true residual: "
                    << result.finalRelResidual << std::endl;
        }

        return result;
      };

    // ------------------------------------------------------------------
    // Helper for solving with a fresh GCRODR solver.
    // ------------------------------------------------------------------

    auto solveWithFreshSolver =
      [&](const RCP<MyMultiVec<ST> >& x,
          const RCP<MyMultiVec<ST> >& b,
          const std::string& label) -> SolveResult
      {
        RCP<Belos::LinearProblem<ST,MV,OP> > problem =
          rcp(new Belos::LinearProblem<ST,MV,OP>(A, x, b));

        const bool set = problem->setProblem();

        TEUCHOS_TEST_FOR_EXCEPTION(
          !set,
          std::runtime_error,
          "Belos::LinearProblem failed to set up correctly.");

        Belos::GCRODRSolMgr<ST,MV,OP> solver(problem, rcp(&belosList, false));

        Belos::ReturnType ret = solver.solve();

        SolveResult result;
        result.iters = solver.getNumIters();
        result.finalRelResidual = computeRelativeResidual(x, b);
        result.ret = ret;
        result.solution = x;

        if (proc_verbose) {
          std::cout << label << " iterations: "
                    << result.iters << std::endl;
          std::cout << label << " final relative true residual: "
                    << result.finalRelResidual << std::endl;
        }

        return result;
      };

    // ------------------------------------------------------------------
    // Build historical exact solutions and RHS vectors.
    // ------------------------------------------------------------------

    RCP<MyMultiVec<ST> > u1Exact = rcp(new MyMultiVec<ST>(dim, 1));
    RCP<MyMultiVec<ST> > u2Exact = rcp(new MyMultiVec<ST>(dim, 1));
    RCP<MyMultiVec<ST> > noise = rcp(new MyMultiVec<ST>(dim, 1));

    MVT::MvRandom(*u1Exact);
    MVT::MvRandom(*u2Exact);
    MVT::MvRandom(*noise);

    RCP<MyMultiVec<ST> > b1 = makeRHS(u1Exact);
    RCP<MyMultiVec<ST> > b2 = makeRHS(u2Exact);

    // ------------------------------------------------------------------
    // Persistent GCRODR solver for historical solves.
    // ------------------------------------------------------------------

    RCP<MyMultiVec<ST> > dummyX = rcp(new MyMultiVec<ST>(dim, 1));
    RCP<MyMultiVec<ST> > dummyB = rcp(new MyMultiVec<ST>(dim, 1));
    MVT::MvInit(*dummyX, zero);
    MVT::MvInit(*dummyB, zero);

    RCP<Belos::LinearProblem<ST,MV,OP> > recycledProblem =
      rcp(new Belos::LinearProblem<ST,MV,OP>(A, dummyX, dummyB));

    recycledProblem->setProblem();

    Belos::GCRODRSolMgr<ST,MV,OP> recycledSolver(
      recycledProblem, rcp(&belosList, false));

    // ------------------------------------------------------------------
    // SolutionProjection database.
    // ------------------------------------------------------------------

    Belos::SolutionProjection<ST,MV,OP> projection(A, maxBasisSize);

    // ------------------------------------------------------------------
    // Historical solve 1.
    // ------------------------------------------------------------------

    RCP<MyMultiVec<ST> > xHist1 = rcp(new MyMultiVec<ST>(dim, 1));
    MVT::MvInit(*xHist1, zero);

    SolveResult hist1 =
      solveWithExistingSolver(recycledSolver, recycledProblem,
                              xHist1, b1, "History solve 1");

    bool historyFailure = false;
    if (hist1.ret != Belos::Converged ||
        hist1.finalRelResidual > solverTol) {
      historyFailure = true;
      if (proc_verbose) {
        std::cout << "ERROR: History solve 1 failed." << std::endl;
      }
    }

    const int accepted1 = projection.addSolution(*hist1.solution);

    // ------------------------------------------------------------------
    // Historical solve 2.
    // ------------------------------------------------------------------

    RCP<MyMultiVec<ST> > xHist2 = rcp(new MyMultiVec<ST>(dim, 1));
    MVT::MvInit(*xHist2, zero);

    SolveResult hist2 =
      solveWithExistingSolver(recycledSolver, recycledProblem,
                              xHist2, b2, "History solve 2");

    if (hist2.ret != Belos::Converged ||
        hist2.finalRelResidual > solverTol) {
      historyFailure = true;
      if (proc_verbose) {
        std::cout << "ERROR: History solve 2 failed." << std::endl;
      }
    }

    const int accepted2 = projection.addSolution(*hist2.solution);

    if (proc_verbose) {
      std::cout << "SolutionProjection accepted historical vectors: "
                << accepted1 << " + " << accepted2 << std::endl;
      std::cout << "Projection basis size after history: "
                << projection.getBasisSize() << std::endl;
    }

    bool projectionBasisFailure = false;
    if (projection.getBasisSize() < 2) {
      projectionBasisFailure = true;
      if (proc_verbose) {
        std::cout << "ERROR: SolutionProjection failed to build a basis of size 2."
                  << std::endl;
      }
    }

    // ------------------------------------------------------------------
    // Construct target exact solution from computed historical solutions.
    //
    //   xTarget = 0.7 xHist1 - 0.25 xHist2 + epsilon * noise.
    //
    // Then targetB = A xTarget.
    // ------------------------------------------------------------------

    RCP<MyMultiVec<ST> > xTargetExact = rcp(new MyMultiVec<ST>(dim, 1));
    RCP<MyMultiVec<ST> > tmp = rcp(new MyMultiVec<ST>(dim, 1));

    MVT::MvAddMv(ST(0.7), *hist1.solution,
                 ST(-0.25), *hist2.solution,
                 *xTargetExact);

    MVT::MvAddMv(one, *xTargetExact,
                 ST(epsilon), *noise,
                 *tmp);

    MVT::Assign(*tmp, *xTargetExact);

    RCP<MyMultiVec<ST> > targetB = makeRHS(xTargetExact);

    // ------------------------------------------------------------------
    // Apply SolutionProjection to zero target initial guess.
    // ------------------------------------------------------------------

    RCP<MyMultiVec<ST> > targetXProjected =
      rcp(new MyMultiVec<ST>(dim, 1));
    RCP<MyMultiVec<ST> > targetXRecycled =
      rcp(new MyMultiVec<ST>(dim, 1));

    MVT::MvInit(*targetXProjected, zero);
    MVT::MvInit(*targetXRecycled, zero);

    RCP<MyMultiVec<ST> > projectedResidual =
      rcp(new MyMultiVec<ST>(dim, 1));

    projection.project(*targetXProjected, *targetB, projectedResidual);

    std::vector<MT> projectedResidualNorm(1);
    std::vector<MT> targetRhsNorm(1);

    MVT::MvNorm(*projectedResidual, projectedResidualNorm);
    MVT::MvNorm(*targetB, targetRhsNorm);

    const MT projectionReduction =
      projectedResidualNorm[0] / targetRhsNorm[0];

    if (proc_verbose) {
      std::cout << "Target projection residual ratio: "
                << projectionReduction << std::endl;
    }

    bool projectionFailure = false;
    if (projectionReduction > projectionReductionTol) {
      projectionFailure = true;
      if (proc_verbose) {
        std::cout << "ERROR: SolutionProjection did not sufficiently reduce "
                  << "the target residual." << std::endl;
      }
    }

    // ------------------------------------------------------------------
    // Target solve A: trained GCRODR with its internal recycle space.
    // ------------------------------------------------------------------

    SolveResult recycledTarget =
      solveWithExistingSolver(recycledSolver, recycledProblem,
                              targetXRecycled, targetB,
                              "Target trained GCRODR");

    bool recycledTargetFailure = false;
    if (recycledTarget.ret != Belos::Converged ||
        recycledTarget.finalRelResidual > solverTol) {
      recycledTargetFailure = true;
      if (proc_verbose) {
        std::cout << "ERROR: Trained recycled GCRODR target solve failed."
                  << std::endl;
      }
    }

    // ------------------------------------------------------------------
    // Target solve B: fresh GCRODR with SolutionProjection initial guess.
    // ------------------------------------------------------------------

    SolveResult projectedTarget =
      solveWithFreshSolver(targetXProjected, targetB,
                           "Target projected fresh GCRODR");

    bool projectedTargetFailure = false;
    if (projectedTarget.ret != Belos::Converged ||
        projectedTarget.finalRelResidual > solverTol) {
      projectedTargetFailure = true;
      if (proc_verbose) {
        std::cout << "ERROR: Projected fresh GCRODR target solve failed."
                  << std::endl;
      }
    }

    // ------------------------------------------------------------------
    // Optional strict comparison.
    // ------------------------------------------------------------------

    bool comparisonFailure = false;

    const int observedImprovement =
      recycledTarget.iters - projectedTarget.iters;

    if (proc_verbose) {
      std::cout << "Target trained GCRODR iterations: "
                << recycledTarget.iters << std::endl;
      std::cout << "Target projected fresh GCRODR iterations: "
                << projectedTarget.iters << std::endl;
      std::cout << "Projected-vs-recycled iteration improvement: "
                << observedImprovement << std::endl;
    }

    if (requireProjectionBeatsRecycled &&
        observedImprovement < minIterationImprovement) {
      comparisonFailure = true;
      if (proc_verbose) {
        std::cout << "ERROR: Projection initial guess did not beat trained "
                  << "GCRODR recycle space by the required amount.  Required "
                  << "improvement: " << minIterationImprovement
                  << ", observed improvement: " << observedImprovement
                  << std::endl;
      }
    }

    // ------------------------------------------------------------------
    // Add the target projected solution to the projection basis.
    // This exercises post-solve basis update.
    // ------------------------------------------------------------------

    const int acceptedTarget =
      projection.addSolution(*projectedTarget.solution);

    if (proc_verbose) {
      std::cout << "Accepted target solution vector(s): "
                << acceptedTarget << std::endl;
      std::cout << "Final projection basis size: "
                << projection.getBasisSize() << std::endl;
    }

    // Clean up.
    delete [] dvals;
    delete [] colptr;
    delete [] rowind;
    delete [] cvals;

    success =
      !historyFailure &&
      !projectionBasisFailure &&
      !projectionFailure &&
      !recycledTargetFailure &&
      !projectedTargetFailure &&
      !comparisonFailure;

    if (proc_verbose) {
      if (success) {
        std::cout << "End Result: TEST PASSED" << std::endl;
      }
      else {
        std::cout << "End Result: TEST FAILED" << std::endl;
      }
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
