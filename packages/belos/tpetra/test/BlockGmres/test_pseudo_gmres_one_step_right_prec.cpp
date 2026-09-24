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
// PseudoBlockGmresIter.  The test runs exactly one Arnoldi step directly in the
// iteration object and compares the optimized solution-space update with the
// standard update after explicitly applying the right preconditioner, all from
// the same iterator state.
//

#include "BelosConfigDefs.hpp"
#include "BelosCurrentSolutionProvider.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosPseudoBlockGmresIter.hpp"
#include "BelosStatusTestMaxIters.hpp"
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

using Teuchos::RCP;
using Teuchos::rcp;

namespace {

template <class SC, class LO, class GO, class NT>
class CountingDiagonalOperator : public Tpetra::Operator<SC, LO, GO, NT> {
public:
  using map_type = Tpetra::Map<LO, GO, NT>;
  using mv_type = Tpetra::MultiVector<SC, LO, GO, NT>;

  CountingDiagonalOperator(const RCP<const map_type>& map,
                           const SC evenDiag,
                           const SC oddDiag)
    : map_(map), evenDiag_(evenDiag), oddDiag_(oddDiag), numApply_(0)
  {}

  void apply(const mv_type& X,
             mv_type& Y,
             Teuchos::ETransp /* mode */ = Teuchos::NO_TRANS,
             SC alpha = Teuchos::ScalarTraits<SC>::one(),
             SC beta = Teuchos::ScalarTraits<SC>::zero()) const override
  {
    auto xView = X.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto yView = Y.getLocalViewHost(Tpetra::Access::ReadWrite);

    for (size_t j = 0; j < X.getNumVectors(); ++j) {
      for (LO lclRow = 0; lclRow < static_cast<LO>(X.getLocalLength()); ++lclRow) {
        const GO gblRow = map_->getGlobalElement(lclRow);
        const SC diag = (gblRow % 2 == 0) ? evenDiag_ : oddDiag_;
        yView(lclRow, j) = beta * yView(lclRow, j) + alpha * diag * xView(lclRow, j);
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

template <class SC, class MV>
bool compareUpdates(const RCP<const MV>& expected,
                    const RCP<const MV>& actual,
                    const char label[],
                    const int me)
{
  using MT = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  using STS = Teuchos::ScalarTraits<SC>;

  auto diff = rcp(new MV(*actual, Teuchos::Copy));
  diff->update(-STS::one(), *expected, STS::one());
  Teuchos::Array<MT> norms(1);
  diff->norm2(norms);
  if (norms[0] > MT(1e-12)) {
    if (me == 0)
      std::cerr << "  FAIL: " << label << " differs from standard solution-space update; norm = "
                << norms[0] << ".\n";
    return false;
  }
  return true;
}

template <class SC>
bool runCase(const bool useFlexibleOneIterUpdate, const bool verbose)
{
  using LO  = typename Tpetra::MultiVector<SC>::local_ordinal_type;
  using GO  = typename Tpetra::MultiVector<SC>::global_ordinal_type;
  using NT  = typename Tpetra::MultiVector<SC>::node_type;
  using MV  = Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP  = Tpetra::Operator<SC, LO, GO, NT>;
  using SDM = Teuchos::SerialDenseMatrix<int, SC>;
  using STS = Teuchos::ScalarTraits<SC>;
  using MVT = Belos::MultiVecTraits<SC, MV, SDM>;

  auto comm = Tpetra::getDefaultComm();
  const int me = comm->getRank();

  const GO numGlobalRows = 10;
  auto map = rcp(new Tpetra::Map<LO, GO, NT>(numGlobalRows, 0, comm));

  // A non-scalar diagonal operator avoids a first-step lucky breakdown; the
  // right preconditioner is scalar and counted so the optimized path can be
  // compared directly with the standard update conversion.
  auto A = rcp(new CountingDiagonalOperator<SC, LO, GO, NT>(map, SC(1), SC(2)));
  auto rightPrec = rcp(new CountingDiagonalOperator<SC, LO, GO, NT>(map, SC(2), SC(2)));

  auto b = rcp(new MV(map, 1));
  auto x = rcp(new MV(map, 1));
  b->putScalar(STS::one());
  x->putScalar(STS::zero());

  auto problem = rcp(new Belos::LinearProblem<SC, MV, OP, SDM>(A, x, b));
  problem->setRightPrec(rightPrec);
  problem->setProblem();
  std::vector<int> currIdx(1, 0);
  problem->setLSIndex(currIdx);

  auto printer = rcp(new Belos::OutputManager<SC>(Belos::Errors));
  auto statusTest = rcp(new Belos::StatusTestMaxIters<SC, MV, OP, SDM>(1));
  auto ortho = rcp(new Belos::ICGSOrthoManager<SC, MV, OP, SDM>());

  Teuchos::ParameterList params;
  params.set("Num Blocks", 2);
  params.set("Use Flexible Gmres Update for One Iteration", useFlexibleOneIterUpdate);

  Belos::PseudoBlockGmresIter<SC, MV, OP, SDM> iter(problem, printer, statusTest, ortho, params);

  Belos::PseudoBlockGmresIterState<SC, MV, SDM> state;
  Teuchos::RCP<MV> R0 = MVT::CloneCopy(*problem->getInitPrecResVec(), currIdx);
  Teuchos::RCP<SDM> z0 = Belos::DenseMatTraits<SC, SDM>::Create(1, 1);
  const int rank = ortho->normalize(*R0, z0);
  if (rank != 1) {
    if (me == 0)
      std::cerr << "  FAIL: initial basis normalization failed.\n";
    return false;
  }
  state.V.resize(1);
  state.Z.resize(1);
  state.V[0] = R0;
  state.Z[0] = z0;
  state.curDim = 0;
  iter.initialize(state);

  iter.iterate();

  bool ok = true;
  if (iter.getNumIters() != 1 || iter.getCurSubspaceDim() != 1) {
    if (me == 0)
      std::cerr << "  FAIL: expected one Arnoldi iteration, got "
                << iter.getNumIters() << " iterations and subspace dimension "
                << iter.getCurSubspaceDim() << ".\n";
    ok = false;
  }

  if (rightPrec->getNumApply() != 1) {
    if (me == 0)
      std::cerr << "  FAIL: expected one preconditioner application during iteration, got "
                << rightPrec->getNumApply() << ".\n";
    ok = false;
  }

  RCP<const MV> optimizedUpdate;
  if (useFlexibleOneIterUpdate) {
    auto* provider = dynamic_cast<Belos::CurrentSolutionProvider<SC, MV, OP, SDM>*>(&iter);
    if (provider == nullptr || !provider->hasCurrentSolution()) {
      if (me == 0)
        std::cerr << "  FAIL: optimized iteration did not provide current solution update.\n";
      ok = false;
    }
    else {
      optimizedUpdate = provider->getCurrentSolutionUpdate();
      if (rightPrec->getNumApply() != 1) {
        if (me == 0)
          std::cerr << "  FAIL: optimized provider applied the right preconditioner again.\n";
        ok = false;
      }
    }
  }
  else {
    auto* provider = dynamic_cast<Belos::CurrentSolutionProvider<SC, MV, OP, SDM>*>(&iter);
    if (provider != nullptr && provider->hasCurrentSolution()) {
      if (me == 0)
        std::cerr << "  FAIL: disabled optimization unexpectedly provided a solution update.\n";
      ok = false;
    }
  }

  auto standardUpdate = iter.getCurrentUpdate();
  auto standardSolutionUpdate = MVT::Clone(*standardUpdate, MVT::GetNumberVecs(*standardUpdate));
  problem->applyRightPrec(*standardUpdate, *standardSolutionUpdate);
  if (rightPrec->getNumApply() != 2) {
    if (me == 0)
      std::cerr << "  FAIL: expected two total preconditioner applications after standard conversion, got "
                << rightPrec->getNumApply() << ".\n";
    ok = false;
  }

  if (useFlexibleOneIterUpdate && optimizedUpdate != Teuchos::null) {
    ok = compareUpdates<SC, MV>(standardSolutionUpdate, optimizedUpdate,
                                "optimized solution-space update", me) && ok;
  }

  if (verbose && me == 0) {
    std::cout << "  "
              << (useFlexibleOneIterUpdate ? "optimized" : "baseline")
              << " case used " << rightPrec->getNumApply()
              << " total right preconditioner application(s).\n";
  }

  return ok;
}

} // namespace

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
      std::cout << "\nCase 1: one-step update option disabled\n";
    ok &= runCase<double>(false, verbose);

    if (verbose && me == 0)
      std::cout << "\nCase 2: one-step update option enabled\n";
    ok &= runCase<double>(true, verbose);

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
