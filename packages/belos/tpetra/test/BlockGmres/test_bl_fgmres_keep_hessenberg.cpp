// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// Test for "Keep Hessenberg" feature in BlockFGmresIter / BlockGmresSolMgr.
//
// When "Flexible Gmres" = true and "Keep Hessenberg" = true, getState().H
// must hold the raw (pre-QR) upper Hessenberg matrix and getState().R must
// hold the QR-rotated upper triangular factor.  They must be distinct objects
// with distinct values: H has nonzero subdiagonal entries; R does not.
//
// When "Keep Hessenberg" = false (the default), H and R must agree — there
// is only one matrix and it contains the QR-factored form.
//
// The test uses a simple n x n tridiagonal matrix so no external data files
// are needed.
//

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBlockFGmresIter.hpp"
#include "BelosGmresIteration.hpp"
#include "BelosStatusTest.hpp"
#include "BelosTypes.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

// -----------------------------------------------------------------------
// StatusTest that snapshots H and R from BlockFGmresIter after each
// iteration.  Returns Undefined so it never affects convergence decisions.
// -----------------------------------------------------------------------
template <class SC, class MV, class OP>
class HessenbergCapture : public Belos::StatusTest<SC, MV, OP> {
 public:
  using State = Belos::GmresIterationState<SC, MV>;

  int                                     curDim = 0;
  bool                                    sawFGmresIter = false;
  RCP<Teuchos::SerialDenseMatrix<int,SC>> H;   // deep copy of raw Hessenberg
  RCP<Teuchos::SerialDenseMatrix<int,SC>> R;   // deep copy of QR-rotated R

  Belos::StatusType checkStatus(Belos::Iteration<SC, MV, OP>* it) override {
    auto* fg = dynamic_cast<Belos::BlockFGmresIter<SC, MV, OP>*>(it);
    if (fg) {
      sawFGmresIter = true;
      State s = fg->getState();
      if (s.curDim > 0) {
        curDim = s.curDim;
        if (s.H)
          H = rcp(new Teuchos::SerialDenseMatrix<int,SC>(*s.H));
        if (s.R)
          R = rcp(new Teuchos::SerialDenseMatrix<int,SC>(*s.R));
      }
    }
    return Belos::Undefined;
  }
  Belos::StatusType getStatus() const override { return Belos::Undefined; }
  void reset() override {
    curDim = 0; sawFGmresIter = false; H = Teuchos::null; R = Teuchos::null;
  }
  void print(std::ostream& os, int indent = 0) const override {
    os << std::string(indent, ' ') << "HessenbergCapture\n";
  }
};

// -----------------------------------------------------------------------
// Build an n x n tridiagonal CrsMatrix: diag = d, off-diag = o
// -----------------------------------------------------------------------
template <class SC, class LO, class GO, class NT>
RCP<Tpetra::CrsMatrix<SC,LO,GO,NT>>
buildTridiagonal(RCP<const Tpetra::Map<LO,GO,NT>> map, SC d, SC o)
{
  auto A = rcp(new Tpetra::CrsMatrix<SC,LO,GO,NT>(map, 3));
  const GO N = static_cast<GO>(map->getGlobalNumElements());
  for (LO i = 0; i < static_cast<LO>(map->getLocalNumElements()); ++i) {
    GO g = map->getGlobalElement(i);
    if (g > 0)
      A->insertGlobalValues(g, Teuchos::tuple(g-1), Teuchos::tuple(o));
    A->insertGlobalValues(g, Teuchos::tuple(g), Teuchos::tuple(d));
    if (g < N - 1)
      A->insertGlobalValues(g, Teuchos::tuple(g+1), Teuchos::tuple(o));
  }
  A->fillComplete();
  return A;
}

// -----------------------------------------------------------------------
// Run the test for one combination of keepHessenberg flag and restart size.
// numBlocks controls the Krylov subspace size before restart; setting it
// smaller than the iterations needed forces at least one restart.
// Returns true if all assertions pass.
// -----------------------------------------------------------------------
template <class SC>
bool runCase(bool keepHessenberg, int numBlocks, bool verbose)
{
  using LO  = typename Tpetra::MultiVector<SC>::local_ordinal_type;
  using GO  = typename Tpetra::MultiVector<SC>::global_ordinal_type;
  using NT  = typename Tpetra::MultiVector<SC>::node_type;
  using MV  = Tpetra::MultiVector<SC,LO,GO,NT>;
  using OP  = Tpetra::Operator<SC,LO,GO,NT>;
  using SDM = Teuchos::SerialDenseMatrix<int,SC>;
  using STS = Teuchos::ScalarTraits<SC>;
  using MT  = typename STS::magnitudeType;
  using STM = Teuchos::ScalarTraits<MT>;

  auto comm = Tpetra::getDefaultComm();
  const int me = comm->getRank();

  // Use a larger system when forcing restarts so convergence needs more
  // iterations than numBlocks.
  const GO n = (numBlocks < 20) ? 200 : 20;
  auto map = rcp(new Tpetra::Map<LO,GO,NT>(n, 0, comm));
  auto A   = buildTridiagonal<SC,LO,GO,NT>(map, SC(4), SC(-1));

  auto b = rcp(new MV(map, 1));
  auto x = rcp(new MV(map, 1));
  b->putScalar(STS::one());
  x->putScalar(STS::zero());

  // Provide a trivial (identity) right preconditioner so BlockGmresSolMgr
  // keeps isFlexible_=true and instantiates BlockFGmresIter rather than
  // falling back to BlockGmresIter when no right prec is present.
  auto I = rcp(new Tpetra::CrsMatrix<SC,LO,GO,NT>(map, 1));
  for (LO i = 0; i < static_cast<LO>(map->getLocalNumElements()); ++i) {
    GO g = map->getGlobalElement(i);
    I->insertGlobalValues(g, Teuchos::tuple(g), Teuchos::tuple(STS::one()));
  }
  I->fillComplete();

  auto problem = rcp(new Belos::LinearProblem<SC,MV,OP>(A, x, b));
  problem->setRightPrec(I);
  problem->setProblem();

  auto params = rcp(new Teuchos::ParameterList);
  params->set("Num Blocks",            numBlocks);
  params->set("Maximum Restarts",      50);
  params->set("Maximum Iterations",   5000);
  params->set("Convergence Tolerance", MT(1e-10));
  params->set("Flexible Gmres",        true);
  params->set("Keep Hessenberg",       keepHessenberg);
  params->set("Verbosity",             Belos::Errors);

  auto capture = rcp(new HessenbergCapture<SC,MV,OP>());
  Belos::BlockGmresSolMgr<SC,MV,OP> solver(problem, params);
  solver.setDebugStatusTest(capture);

  Belos::ReturnType ret = solver.solve();

  if (ret != Belos::Converged) {
    if (me == 0)
      std::cerr << "  FAIL: solver did not converge.\n";
    return false;
  }

  const int numIters = solver.getNumIters();

  // If a small numBlocks was requested, verify that the solver actually
  // restarted (total iterations must exceed the restart cycle length).
  if (numBlocks < 20) {
    if (numIters <= numBlocks) {
      if (me == 0)
        std::cerr << "  FAIL (restart test): solver converged in " << numIters
                  << " iterations without restarting (numBlocks = " << numBlocks << ")."
                     " Increase problem size or tighten tolerance.\n";
      return false;
    }
    if (verbose && me == 0)
      std::cout << "  Confirmed restart: converged in " << numIters
                << " iterations with numBlocks = " << numBlocks << ".\n";
  }

  if (!capture->sawFGmresIter) {
    if (me == 0)
      std::cerr << "  FAIL: solver did not use BlockFGmresIter (dynamic_cast never succeeded)."
                   " Check that isFlexible_=true and a right preconditioner is set.\n";
    return false;
  }

  const int dim = capture->curDim;
  if (dim <= 0) {
    if (me == 0)
      std::cerr << "  FAIL: no iteration state captured (curDim = " << dim << ").\n";
    return false;
  }
  if (!capture->H || !capture->R) {
    if (me == 0)
      std::cerr << "  FAIL: H or R is null after solve.\n";
    return false;
  }

  const SDM& H = *capture->H;
  const SDM& R = *capture->R;

  if (keepHessenberg) {
    // H must have at least one nonzero subdiagonal entry (raw Hessenberg).
    // R must have all subdiagonal entries equal to zero (QR-rotated).
    bool H_has_subdiag         = false;
    bool R_has_nonzero_subdiag = false;
    for (int j = 0; j < dim; ++j) {
      if (STS::magnitude(H(j+1, j)) > STM::zero())  H_has_subdiag = true;
      if (STS::magnitude(R(j+1, j)) > MT(1e-14))    R_has_nonzero_subdiag = true;
    }
    if (!H_has_subdiag) {
      if (me == 0)
        std::cerr << "  FAIL (Keep Hessenberg=true): H has no nonzero subdiagonal entries.\n";
      return false;
    }
    if (R_has_nonzero_subdiag) {
      if (me == 0)
        std::cerr << "  FAIL (Keep Hessenberg=true): R has nonzero subdiagonal entries.\n";
      return false;
    }
    if (verbose && me == 0)
      std::cout << "  PASS (Keep Hessenberg=true): H has subdiagonals; R is upper triangular."
                   " (" << numIters << " iters)\n";
  } else {
    // H and R must agree pointwise (same object, both hold the QR-rotated form).
    MT maxdiff = STM::zero();
    for (int j = 0; j < dim; ++j)
      for (int i = 0; i < H.numRows() && i <= j+1; ++i)
        maxdiff = std::max(maxdiff, STS::magnitude(H(i,j) - R(i,j)));

    if (maxdiff > MT(1e-14)) {
      if (me == 0)
        std::cerr << "  FAIL (Keep Hessenberg=false): H and R differ (max |H-R| = "
                  << maxdiff << ").\n";
      return false;
    }
    // R must also be upper triangular.
    bool R_has_nonzero_subdiag = false;
    for (int j = 0; j < dim; ++j)
      if (STS::magnitude(R(j+1, j)) > MT(1e-14))
        R_has_nonzero_subdiag = true;
    if (R_has_nonzero_subdiag) {
      if (me == 0)
        std::cerr << "  FAIL (Keep Hessenberg=false): R has nonzero subdiagonal entries.\n";
      return false;
    }
    if (verbose && me == 0)
      std::cout << "  PASS (Keep Hessenberg=false): H == R, both upper triangular.\n";
  }

  return true;
}

int main(int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);

  bool verbose = false;
  bool success = false;

  try {
    auto comm = Tpetra::getDefaultComm();
    const int me = comm->getRank();

    for (int i = 1; i < argc; ++i)
      if (std::string(argv[i]) == "--verbose")
        verbose = true;

    bool ok = true;

    if (verbose && me == 0)
      std::cout << "\nCase 1: Keep Hessenberg = true (no restart)\n";
    ok &= runCase<double>(true,  30, verbose);

    if (verbose && me == 0)
      std::cout << "\nCase 2: Keep Hessenberg = false (no restart)\n";
    ok &= runCase<double>(false, 30, verbose);

    if (verbose && me == 0)
      std::cout << "\nCase 3: Keep Hessenberg = true (forced restart, numBlocks=10)\n";
    ok &= runCase<double>(true,  10, verbose);

    success = ok;

    if (me == 0) {
      if (success)
        std::cout << "\nEnd Result: TEST PASSED\n";
      else
        std::cout << "\nEnd Result: TEST FAILED\n";
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
