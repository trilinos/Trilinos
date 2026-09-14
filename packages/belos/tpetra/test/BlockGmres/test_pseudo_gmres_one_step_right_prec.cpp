// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// Test the optional one-step right-preconditioned update path in
// PseudoBlockGmresSolMgr.  With one GMRES iteration, the preconditioned
// Arnoldi basis vector can be reused to form the solution update, avoiding a
// second right-preconditioner application.
//

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "BelosStatusTest.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosTypes.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include <cstdlib>
#include <iostream>
#include <string>

using Teuchos::RCP;
using Teuchos::rcp;

namespace {

template <class SC, class LO, class GO, class NT>
class CountingDiagonalOperator : public Tpetra::Operator<SC, LO, GO, NT> {
public:
  using map_type = Tpetra::Map<LO, GO, NT>;
  using mv_type = Tpetra::MultiVector<SC, LO, GO, NT>;

  CountingDiagonalOperator(const RCP<const map_type> &map, const SC evenDiag,
                           const SC oddDiag)
      : map_(map), evenDiag_(evenDiag), oddDiag_(oddDiag), numApply_(0) {}

  void apply(const mv_type &X, mv_type &Y,
             Teuchos::ETransp /* mode */ = Teuchos::NO_TRANS,
             SC alpha = Teuchos::ScalarTraits<SC>::one(),
             SC beta = Teuchos::ScalarTraits<SC>::zero()) const override {
    auto xView = X.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto yView = Y.getLocalViewHost(Tpetra::Access::ReadWrite);

    for (size_t j = 0; j < X.getNumVectors(); ++j) {
      for (LO lclRow = 0; lclRow < static_cast<LO>(X.getLocalLength());
           ++lclRow) {
        const GO gblRow = map_->getGlobalElement(lclRow);
        const SC diag = (gblRow % 2 == 0) ? evenDiag_ : oddDiag_;
        yView(lclRow, j) =
            beta * yView(lclRow, j) + alpha * diag * xView(lclRow, j);
      }
    }
    ++numApply_;
  }

  bool hasTransposeApply() const override { return false; }

  RCP<const map_type> getDomainMap() const override { return map_; }

  RCP<const map_type> getRangeMap() const override { return map_; }

  int getNumApply() const { return numApply_; }

private:
  RCP<const map_type> map_;
  SC evenDiag_;
  SC oddDiag_;
  mutable int numApply_;
};

template <class SC, class MV, class OP, class DM>
class PassAfterOneIterStatusTest : public Belos::StatusTest<SC, MV, OP, DM> {
public:
  Belos::StatusType
  checkStatus(Belos::Iteration<SC, MV, OP, DM> *iter) override {
    status_ = iter->getNumIters() >= 1 ? Belos::Passed : Belos::Failed;
    return status_;
  }

  Belos::StatusType getStatus() const override { return status_; }

  void reset() override { status_ = Belos::Undefined; }

  void print(std::ostream &os, int indent = 0) const override {
    os << std::string(indent, ' ') << "PassAfterOneIterStatusTest\n";
  }

private:
  Belos::StatusType status_ = Belos::Undefined;
};

template <class SC>
bool runCase(const bool useFlexibleOneIterUpdate, const bool useUserStatusTest,
             const bool verbose) {
  using LO = typename Tpetra::MultiVector<SC>::local_ordinal_type;
  using GO = typename Tpetra::MultiVector<SC>::global_ordinal_type;
  using NT = typename Tpetra::MultiVector<SC>::node_type;
  using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP = Tpetra::Operator<SC, LO, GO, NT>;
  using SDM = Teuchos::SerialDenseMatrix<int, SC>;
  using STS = Teuchos::ScalarTraits<SC>;
  using MT = typename STS::magnitudeType;

  auto comm = Tpetra::getDefaultComm();
  const int me = comm->getRank();

  const GO numGlobalRows = 10;
  auto map = rcp(new Tpetra::Map<LO, GO, NT>(numGlobalRows, 0, comm));

  auto A = rcp(new CountingDiagonalOperator<SC, LO, GO, NT>(map, STS::one(),
                                                            STS::one()));
  auto rightPrec =
      rcp(new CountingDiagonalOperator<SC, LO, GO, NT>(map, SC(2), SC(3)));

  auto b = rcp(new MV(map, 1));
  auto x = rcp(new MV(map, 1));
  b->putScalar(STS::one());
  x->putScalar(STS::zero());

  auto problem = rcp(new Belos::LinearProblem<SC, MV, OP, SDM>(A, x, b));
  problem->setRightPrec(rightPrec);
  problem->setProblem();

  auto params = rcp(new Teuchos::ParameterList);
  params->set("Num Blocks", 2);
  params->set("Maximum Iterations", 2);
  params->set("Maximum Restarts", 0);
  params->set("Convergence Tolerance", MT(0.25));
  params->set("Use Flexible Gmres Update for One Iteration",
              useFlexibleOneIterUpdate);
  params->set("Verbosity", Belos::Errors);

  Belos::PseudoBlockGmresSolMgr<SC, MV, OP, SDM> solver(problem, params);
  if (useUserStatusTest) {
    auto userTest = rcp(new PassAfterOneIterStatusTest<SC, MV, OP, SDM>());
    solver.setUserConvStatusTest(userTest);
  }
  Belos::ReturnType ret = solver.solve();

  bool ok = true;
  if (ret != Belos::Converged) {
    if (me == 0)
      std::cerr << "  FAIL: solver did not converge.\n";
    ok = false;
  }

  if (solver.getNumIters() != 1) {
    if (me == 0)
      std::cerr << "  FAIL: expected 1 GMRES iteration, got "
                << solver.getNumIters() << ".\n";
    ok = false;
  }

  const int expectedPrecApplies =
      useFlexibleOneIterUpdate ? 1 : (useUserStatusTest ? 3 : 2);
  if (rightPrec->getNumApply() != expectedPrecApplies) {
    if (me == 0)
      std::cerr << "  FAIL: expected " << expectedPrecApplies
                << " right preconditioner application(s), got "
                << rightPrec->getNumApply() << ".\n";
    ok = false;
  }

  // For A = I, b = [1, ...], x0 = 0, and right preconditioner
  // M^{-1} = diag(2,3,2,3,...), one GMRES step gives
  // x = (mean(diag(M^{-1})) / mean(diag(M^{-1})^2)) * diag(M^{-1}).
  auto expected = rcp(new MV(map, 1));
  auto expectedView = expected->getLocalViewHost(Tpetra::Access::ReadWrite);
  const SC coefficient = SC(2.5 / 6.5);
  for (LO lclRow = 0; lclRow < static_cast<LO>(expected->getLocalLength());
       ++lclRow) {
    const GO gblRow = map->getGlobalElement(lclRow);
    const SC diag = (gblRow % 2 == 0) ? SC(2) : SC(3);
    expectedView(lclRow, 0) = coefficient * diag;
  }

  auto error = rcp(new MV(*x, Teuchos::Copy));
  error->update(-STS::one(), *expected, STS::one());
  Teuchos::Array<MT> norms(1);
  error->norm2(norms);
  if (norms[0] > MT(1e-12)) {
    if (me == 0)
      std::cerr << "  FAIL: solution error norm is " << norms[0] << ".\n";
    ok = false;
  }

  if (verbose && me == 0) {
    std::cout << "  " << (useFlexibleOneIterUpdate ? "optimized" : "baseline")
              << (useUserStatusTest ? " user-status" : "") << " case used "
              << rightPrec->getNumApply()
              << " right preconditioner application(s).\n";
  }

  return ok;
}

} // namespace

int main(int argc, char *argv[]) {
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
      std::cout << "\nCase 1: one-step update option disabled\n";
    ok &= runCase<double>(false, false, verbose);

    if (verbose && me == 0)
      std::cout << "\nCase 2: one-step update option enabled\n";
    ok &= runCase<double>(true, false, verbose);

    if (verbose && me == 0)
      std::cout
          << "\nCase 3: one-step update option enabled with user status test\n";
    ok &= runCase<double>(true, true, verbose);

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
