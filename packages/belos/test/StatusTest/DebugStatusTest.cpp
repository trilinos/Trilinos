// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Verify that a *debug status test* installed via
// SolverManager::setDebugStatusTest() which trips mid-solve cleanly stops the
// Krylov iteration and yields Belos::Unconverged WITHOUT throwing -- for every
// canonical solver manager supported by this factory, and that this still holds
// when the same manager instance is reused for a second solve.
//
// The debug test used here (TripAtIter) is fully deterministic: it trips
// exactly when the solver reaches a fixed iteration count, independent of
// hardware/wall-clock.  This mirrors how a real wall-clock time limit would
// stop a solve, but without any timing nondeterminism.
//

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "BelosConfigDefs.hpp"
#include "BelosIteration.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosSolverManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosTypes.hpp"
#include "Belos_Details_EBelosSolverType.hpp"

#include "MyBetterOperator.hpp"
#include "MyMultiVec.hpp"

#include <string>
#include <vector>

namespace {

//
// A deterministic debug status test that "trips" (returns Belos::Passed) as
// soon as the underlying iteration has performed at least tripIter iterations.
// Before that it returns Belos::Failed (i.e. "keep going").  This stands in for
// a wall-clock time limit, but is reproducible on any hardware.
//
template <class ScalarType, class MV, class OP>
class TripAtIter : public Belos::StatusTest<ScalarType, MV, OP> {
public:
  explicit TripAtIter(const int tripIter)
      : tripIter_(tripIter), status_(Belos::Undefined) {}

  ~TripAtIter() override = default;

  Belos::StatusType
  checkStatus(Belos::Iteration<ScalarType, MV, OP> *iSolver) override {
    const int iters = iSolver->getNumIters();
    status_ = (iters >= tripIter_) ? Belos::Passed : Belos::Failed;
    return status_;
  }

  Belos::StatusType getStatus() const override { return status_; }

  void reset() override { status_ = Belos::Undefined; }

  void print(std::ostream &os, int indent = 0) const override {
    for (int j = 0; j < indent; ++j) {
      os << ' ';
    }
    this->printStatus(os, status_);
    os << "TripAtIter( trip at iteration >= " << tripIter_ << " )" << std::endl;
  }

  int tripIter() const { return tripIter_; }

private:
  int tripIter_;
  Belos::StatusType status_;
};

} // namespace

// The heart of the test: for a given solver name, create the manager via the
// factory, install a debug test that trips at 'tripIter', solve, and verify the
// solve stopped cleanly with Belos::Unconverged.  Then reuse the SAME manager
// for a second solve with a fresh debug test and verify the same behavior.
TEUCHOS_UNIT_TEST(StatusTest, DebugStatusTestStopsSolve) {
  using std::endl;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef double ST;
  typedef Belos::MultiVec<ST> MV;
  typedef Belos::Operator<ST> OP;
  typedef Belos::LinearProblem<ST, MV, OP> problem_type;
  typedef Belos::SolverManager<ST, MV, OP> solver_type;
  typedef Belos::SolverFactory<ST, MV, OP> factory_type;
  typedef TripAtIter<ST, MV, OP> trip_test_type;

  Teuchos::OSTab tab0(out);
  out << "Verify setDebugStatusTest() stop => Belos::Unconverged (no throw), "
         "for every supported canonical solver, including manager reuse"
      << endl;
  Teuchos::OSTab tab1(out);

  // -----------------------------------------------------------------------
  // Build a small, real, symmetric-positive-definite tridiagonal system
  // (the classic 1D Laplacian: diag = 2, off-diag = -1).  It is SPD (so the
  // symmetric solvers such as CG / MINRES are well posed) and its condition
  // number grows with the dimension, so with a tight tolerance every iterative
  // solver needs many more than 'tripIter' iterations to converge.  That
  // guarantees the debug test -- not natural convergence -- is what stops the
  // solve.  MyBetterOperator stores the matrix in 1-based CSC format; since the
  // matrix is symmetric, CSC and CSR coincide.
  // -----------------------------------------------------------------------
  const int dim = 100;

  std::vector<int> colptr(dim + 1);
  std::vector<int> rowind;
  std::vector<ST> vals;
  rowind.reserve(3 * dim);
  vals.reserve(3 * dim);

  int count = 0;
  colptr[0] = 1; // 1-based
  for (int j = 0; j < dim; ++j) {
    if (j - 1 >= 0) {
      rowind.push_back((j - 1) + 1);
      vals.push_back(-1.0);
      ++count;
    }
    /* diagonal */ rowind.push_back(j + 1);
    vals.push_back(2.0);
    ++count;
    if (j + 1 < dim) {
      rowind.push_back((j + 1) + 1);
      vals.push_back(-1.0);
      ++count;
    }
    colptr[j + 1] = count + 1; // 1-based pointer to start of next column
  }
  const int nnz = count;

  RCP<MyBetterOperator<ST>> A = rcp(new MyBetterOperator<ST>(
      dim, colptr.data(), nnz, rowind.data(), vals.data()));

  const int numrhs = 1;

  // Build a consistent right-hand side B = A * xTrue, then use a zero initial
  // guess for the left-hand side X.
  RCP<MyMultiVec<ST>> xTrue = rcp(new MyMultiVec<ST>(dim, numrhs));
  RCP<MyMultiVec<ST>> B = rcp(new MyMultiVec<ST>(dim, numrhs));
  RCP<MyMultiVec<ST>> X = rcp(new MyMultiVec<ST>(dim, numrhs));
  xTrue->MvRandom();
  A->Apply(*xTrue, *B);
  X->MvInit(0.0);

  RCP<problem_type> problem = rcp(new problem_type(A, X, B));
  TEST_ASSERT(problem->setProblem());

  // 'tripIter' is chosen well below "Maximum Iterations" so the debug test is
  // always what stops the solve, never maxiter or natural convergence.
  const int tripIter = 3;
  const int maxIters = 200;
  const double tol = 1.0e-12;

  // NOTE: We intentionally do NOT set "Verbosity" or "Output Style" here: some
  // managers (e.g. LSQR) build a residual-norm output test that throws when
  // those are set (see Belos test/Factory/Factory.cpp, Bug 6383 comments).

  const std::vector<std::string> names = Belos::Details::canonicalSolverNames();
  // NOTE: BlockGCRODRSolMgr is intentionally absent -- it is not registered in
  // canonicalSolverNames()/the SolverFactory, so it is not exercised here.

  for (size_t k = 0; k < names.size(); ++k) {
    const std::string &name = names[k];
    RCP<factory_type> factory = rcp(new factory_type());
    if (!factory->isSupported(name)) {
      out << "=== Solver: \"" << name << "\" (not registered; skipped) ==="
          << endl;
      continue;
    }

    out << "=== Solver: \"" << name << "\" ===" << endl;
    Teuchos::OSTab tabN(out);

    // Report failures per-solver but keep going, so one bad manager does not
    // hide the results for the others.
    try {
      RCP<ParameterList> params = parameterList("Belos");
      params->set("Maximum Iterations", maxIters);
      params->set("Convergence Tolerance", tol);
      if (name == "HYBRID BLOCK GMRES") {
        // Exercise GmresPolySolMgr's outer-solver path, where the underlying
        // Krylov solver can honor the debug status test.
        params->set("Maximum Degree", 0);
        params->set("Outer Solver", "PSEUDOBLOCK GMRES");
        params->sublist("Outer Solver Params").set("Maximum Iterations", maxIters);
        params->sublist("Outer Solver Params").set("Convergence Tolerance", tol);
      }

      RCP<solver_type> solver;
      TEST_NOTHROW(solver = factory->create(name, params));
      TEST_ASSERT(!solver.is_null());
      if (solver.is_null()) {
        success = false;
        continue;
      }

      // -------- First solve (fresh manager) --------
      X->MvInit(0.0);
      TEST_ASSERT(problem->setProblem());
      solver->setProblem(problem);

      RCP<trip_test_type> trip1 = rcp(new trip_test_type(tripIter));
      solver->setDebugStatusTest(trip1);

      Belos::ReturnType ret1 = Belos::Converged;
      TEST_NOTHROW(ret1 = solver->solve());

      const int iters1 = solver->getNumIters();
      out << "first solve: ret="
          << (ret1 == Belos::Unconverged ? "Unconverged" : "Converged")
          << ", numIters=" << iters1 << endl;

      TEST_EQUALITY(ret1, Belos::Unconverged);
      TEST_EQUALITY(trip1->getStatus(), Belos::Passed);
      // The debug test -- not maxiter -- stopped the solve.
      TEST_ASSERT(iters1 >= tripIter);
      TEST_ASSERT(iters1 < maxIters);

      // -------- Second solve (REUSE the same manager instance) --------
      // A fresh debug test is installed and setProblem() is called again on the
      // SAME manager, verifying the manager rebuilds/re-wires its status-test
      // tree correctly on reuse.
      X->MvInit(0.0);
      TEST_ASSERT(problem->setProblem());
      solver->setProblem(problem);

      RCP<trip_test_type> trip2 = rcp(new trip_test_type(tripIter));
      solver->setDebugStatusTest(trip2);

      Belos::ReturnType ret2 = Belos::Converged;
      TEST_NOTHROW(ret2 = solver->solve());

      const int iters2 = solver->getNumIters();
      out << "reused solve: ret="
          << (ret2 == Belos::Unconverged ? "Unconverged" : "Converged")
          << ", numIters=" << iters2 << endl;

      TEST_EQUALITY(ret2, Belos::Unconverged);
      TEST_EQUALITY(trip2->getStatus(), Belos::Passed);
      TEST_ASSERT(iters2 >= tripIter);
      TEST_ASSERT(iters2 < maxIters);
    } catch (std::exception &e) {
      out << "*** Solver \"" << name << "\" threw an exception: " << e.what()
          << endl;
      success = false;
    } catch (...) {
      out << "*** Solver \"" << name << "\" threw a non-std::exception" << endl;
      success = false;
    }
  }
}
