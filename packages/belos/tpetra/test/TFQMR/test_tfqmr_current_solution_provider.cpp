// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// Test CurrentSolutionProvider with a non-GMRES iteration.
//

#include "BelosConfigDefs.hpp"
#include "BelosCurrentSolutionProvider.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosOutputManager.hpp"
#include "BelosStatusTest.hpp"
#include "BelosTFQMRIter.hpp"
#include "BelosTpetraAdapter.hpp"

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
class DiagonalOperator : public Tpetra::Operator<SC, LO, GO, NT> {
public:
  using map_type = Tpetra::Map<LO, GO, NT>;
  using mv_type = Tpetra::MultiVector<SC, LO, GO, NT>;

  DiagonalOperator(const RCP<const map_type>& map, const SC evenDiag, const SC oddDiag)
    : map_(map), evenDiag_(evenDiag), oddDiag_(oddDiag)
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
  }

  bool hasTransposeApply() const override { return false; }

  RCP<const map_type> getDomainMap() const override { return map_; }

  RCP<const map_type> getRangeMap() const override { return map_; }

private:
  RCP<const map_type> map_;
  SC evenDiag_;
  SC oddDiag_;
};

template <class SC, class MV, class OP, class DM>
class ProviderCheckStatusTest : public Belos::StatusTest<SC, MV, OP, DM> {
public:
  using MVT = Belos::MultiVecTraits<SC, MV, DM>;
  using MT = typename Teuchos::ScalarTraits<SC>::magnitudeType;

  ProviderCheckStatusTest(const MT tolerance)
    : tolerance_(tolerance), status_(Belos::Undefined), sawProvider_(false), solutionChanged_(false)
  {}

  Belos::StatusType checkStatus(Belos::Iteration<SC, MV, OP, DM>* iter) override
  {
    if (iter->getNumIters() == 0) {
      status_ = Belos::Failed;
      return status_;
    }

    auto* provider = dynamic_cast<Belos::CurrentSolutionProvider<SC, MV, OP, DM>*>(iter);
    if (provider == nullptr || !provider->hasCurrentSolution()) {
      status_ = Belos::Failed;
      return status_;
    }

    sawProvider_ = true;
    const auto& problem = iter->getProblem();
    auto before = problem.updateSolution(Teuchos::null);
    auto beforeCopy = MVT::CloneCopy(*before);
    auto currentSolution = provider->getCurrentSolution();
    auto after = problem.updateSolution(Teuchos::null);

    auto mutationCheck = MVT::CloneCopy(*after);
    MVT::MvAddMv(SC(-1), *beforeCopy, SC(1), *mutationCheck, *mutationCheck);
    std::vector<MT> mutationNorms(MVT::GetNumberVecs(*mutationCheck));
    MVT::MvNorm(*mutationCheck, mutationNorms);
    solutionChanged_ = mutationNorms[0] > tolerance_;

    auto residual = MVT::Clone(*currentSolution, MVT::GetNumberVecs(*currentSolution));
    problem.computeCurrResVec(&*residual, &*currentSolution);
    std::vector<MT> residualNorms(MVT::GetNumberVecs(*residual));
    MVT::MvNorm(*residual, residualNorms);
    status_ = residualNorms[0] <= tolerance_ ? Belos::Passed : Belos::Failed;
    return status_;
  }

  Belos::StatusType getStatus() const override { return status_; }

  void reset() override
  {
    status_ = Belos::Undefined;
    sawProvider_ = false;
    solutionChanged_ = false;
  }

  void print(std::ostream& os, int indent = 0) const override
  {
    os << std::string(indent, ' ') << "ProviderCheckStatusTest\n";
  }

  bool sawProvider() const { return sawProvider_; }

  bool solutionChanged() const { return solutionChanged_; }

private:
  MT tolerance_;
  Belos::StatusType status_;
  bool sawProvider_;
  bool solutionChanged_;
};

template <class SC>
bool runCase(const bool verbose)
{
  using LO  = typename Tpetra::MultiVector<SC>::local_ordinal_type;
  using GO  = typename Tpetra::MultiVector<SC>::global_ordinal_type;
  using NT  = typename Tpetra::MultiVector<SC>::node_type;
  using MV  = Tpetra::MultiVector<SC, LO, GO, NT>;
  using OP  = Tpetra::Operator<SC, LO, GO, NT>;
  using SDM = Teuchos::SerialDenseMatrix<int, SC>;
  using STS = Teuchos::ScalarTraits<SC>;
  using MT  = typename STS::magnitudeType;

  auto comm = Tpetra::getDefaultComm();
  const int me = comm->getRank();

  const GO numGlobalRows = 10;
  auto map = rcp(new Tpetra::Map<LO, GO, NT>(numGlobalRows, 0, comm));
  auto A = rcp(new DiagonalOperator<SC, LO, GO, NT>(map, STS::one(), STS::one()));

  auto b = rcp(new MV(map, 1));
  auto x = rcp(new MV(map, 1));
  b->putScalar(STS::one());
  x->putScalar(STS::zero());

  auto problem = rcp(new Belos::LinearProblem<SC, MV, OP, SDM>(A, x, b));
  problem->setProblem();
  std::vector<int> currIdx(1, 0);
  problem->setLSIndex(currIdx);

  auto params = rcp(new Teuchos::ParameterList);
  auto printer = rcp(new Belos::OutputManager<SC>(Belos::Errors));
  auto statusTest = rcp(new ProviderCheckStatusTest<SC, MV, OP, SDM>(MT(1e-12)));
  Belos::TFQMRIter<SC, MV, OP, SDM> iter(problem, printer, statusTest, *params);

  Belos::TFQMRIterState<SC, MV> state;
  state.R = problem->getInitResVec();
  iter.initializeTFQMR(state);
  iter.iterate();

  bool ok = true;
  if (iter.getNumIters() != 1) {
    if (me == 0)
      std::cerr << "  FAIL: expected TFQMR to stop after 1 iteration, got "
                << iter.getNumIters() << ".\n";
    ok = false;
  }
  if (!statusTest->sawProvider()) {
    if (me == 0)
      std::cerr << "  FAIL: status test did not see CurrentSolutionProvider.\n";
    ok = false;
  }
  if (statusTest->solutionChanged()) {
    if (me == 0)
      std::cerr << "  FAIL: getCurrentSolution() mutated LinearProblem solution.\n";
    ok = false;
  }

  if (verbose && me == 0 && ok)
    std::cout << "  TFQMR CurrentSolutionProvider case passed.\n";

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

    success = runCase<double>(verbose);

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
