// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <ios>
#include <utility>
#include "Teko_AdaptivePreconditionerFactory.hpp"
#include "Teko_ImplicitLinearOp.hpp"
#include "Teko_PreconditionerState.hpp"
#include "Teko_Utilities.hpp"
#include "Teuchos_ENull.hpp"

namespace Teko {

using SizeOrdinalType = AdaptivePreconditionerFactory::SizeOrdinalType;

namespace {

LinearOp buildInverse(const InverseFactory &invFact, const Teuchos::RCP<InverseFactory> &precFact,
                      const LinearOp &matrix, PreconditionerState &state,
                      const std::string &opPrefix, size_t i) {
  std::stringstream ss;
  ss << opPrefix << "_" << i;

  ModifiableLinearOp &invOp  = state.getModifiableOp(ss.str());
  ModifiableLinearOp &precOp = state.getModifiableOp("prec_" + ss.str());

  if (precFact != Teuchos::null) {
    if (precOp == Teuchos::null) {
      precOp = precFact->buildInverse(matrix);
      state.addModifiableOp("prec_" + ss.str(), precOp);
    } else {
      rebuildInverse(*precFact, matrix, precOp);
    }
  }

  if (invOp == Teuchos::null) {
    if (precOp.is_null()) {
      invOp = buildInverse(invFact, matrix);
    } else {
      invOp = buildInverse(invFact, matrix, precOp);
    }
  } else {
    if (precOp.is_null()) {
      rebuildInverse(invFact, matrix, invOp);
    } else {
      rebuildInverse(invFact, matrix, precOp, invOp);
    }
  }

  return invOp;
}

class AdaptiveLinearOp : public ImplicitLinearOp {
 public:
  AdaptiveLinearOp(LinearOp A_, const std::vector<Teuchos::RCP<InverseFactory>> &inverses_,
                   const std::vector<Teuchos::RCP<InverseFactory>> &preconditioners_,
                   const std::vector<SizeOrdinalType> &maximumSizes_, PreconditionerState &state_,
                   AdaptiveSolverSettings settings_)
      : A(std::move(A_)),
        inverses(inverses_),
        preconditioners(preconditioners_),
        maximumSizes(maximumSizes_),
        state(state_),
        settings(settings_) {
    subblockSize = A->range()->dim();
  }

  ~AdaptiveLinearOp() override = default;

  VectorSpace range() const override { return A->range(); }
  VectorSpace domain() const override { return A->domain(); }

  void implicitApply(const MultiVector &x, MultiVector &y, const double alpha = 1.0,
                     const double beta = 0.0) const override {
    bool converged                    = true;
    bool additionalIterationsRequired = false;
    double residualReduction          = 0.0;
    bool try_next_solver              = true;

    do {
      applyOp(A_inv, x, y, alpha, beta);

      residualReduction = computeMaxRelativeNorm(x, y);

      converged = residualReduction <= settings.targetResidualReduction;
      if (converged) continue;

      if (!has_another_solver()) break;

      try_next_solver = setup_next_solver();

      additionalIterationsRequired = true;

    } while (!converged && try_next_solver);

    successfulApplications++;
    if (additionalIterationsRequired || !converged) successfulApplications = 0;

    // revert back to initial solver in chain
    if (successfulApplications >= settings.numAppliesCycle && index != 0) {
      index                  = 0;
      successfulApplications = 0;
      setup_solver();
    }
  }

  void describe(Teuchos::FancyOStream &out_arg,
                const Teuchos::EVerbosityLevel verbLevel) const override {
    A->describe(out_arg, verbLevel);
  }

  void initialize_step() const { setup_solver(); }

 private:
  bool has_another_solver() const { return (index + 1) < inverses.size(); }

  bool setup_solver() const {
    if (subblockSize > maximumSizes.at(index)) {
      return false;
    }

    try {
      A_inv = buildInverse(*inverses.at(index), preconditioners.at(index), A, state, prefix, index);
      return true;
    } catch (std::exception &exc) {
      return false;
    }
  }

  bool setup_next_solver() const {
    if (!has_another_solver()) return false;

    const auto currentSolverIndex = index;
    while (has_another_solver()) {
      index++;
      const auto success = setup_solver();
      if (success) return true;
    }

    index = currentSolverIndex;
    return false;
  }

  double computeMaxRelativeNorm(const MultiVector &x, MultiVector &y) const {
    const double alpha = 1.0;
    std::vector<double> norms(y->domain()->dim());
    std::vector<double> rhsNorms(x->domain()->dim());

    auto residual = deepcopy(x);
    applyOp(A, y, residual, -1.0, alpha);

    Thyra::norms_2<double>(*residual, Teuchos::arrayViewFromVector(norms));
    Thyra::norms_2<double>(*x, Teuchos::arrayViewFromVector(rhsNorms));

    double maxRelRes = 0.0;
    for (auto i = 0U; i < norms.size(); ++i) {
      const auto relRes = rhsNorms[i] == 0 ? norms[i] : norms[i] / rhsNorms[i];
      maxRelRes         = std::max(relRes, maxRelRes);
    }

    return maxRelRes;
  }

  LinearOp A;
  std::vector<Teuchos::RCP<InverseFactory>> inverses;
  std::vector<Teuchos::RCP<InverseFactory>> preconditioners;
  std::vector<SizeOrdinalType> maximumSizes;
  PreconditionerState &state;
  AdaptiveSolverSettings settings;
  SizeOrdinalType subblockSize;
  mutable size_t index               = 0;
  mutable int successfulApplications = 0;
  mutable LinearOp A_inv;

  const std::string prefix = "adaptive";
};

LinearOp create_adaptive_linear_operator(
    const LinearOp &A, const std::vector<Teuchos::RCP<InverseFactory>> &inverses,
    const std::vector<Teuchos::RCP<InverseFactory>> &preconditioners,
    const std::vector<SizeOrdinalType> &maximumSizes, PreconditionerState &state,
    const AdaptiveSolverSettings &settings) {
  ModifiableLinearOp &adaptiveOp = state.getModifiableOp("adaptive_linear_op");
  if (adaptiveOp == Teuchos::null) {
    adaptiveOp = Teuchos::rcp(
        new AdaptiveLinearOp(A, inverses, preconditioners, maximumSizes, state, settings));
  }

  {
    auto adaptiveSolver = Teuchos::rcp_dynamic_cast<AdaptiveLinearOp>(adaptiveOp);
    adaptiveSolver->initialize_step();
  }

  return adaptiveOp;
}

}  // namespace

LinearOp AdaptivePreconditionerFactory::buildPreconditionerOperator(
    LinearOp &lo, PreconditionerState &state) const {
  return create_adaptive_linear_operator(lo, inverses, preconditioners, maximumSizes, state,
                                         settings);
}

void AdaptivePreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList &pl) {
  settings.numAppliesCycle         = 100;
  settings.targetResidualReduction = 0.1;

  if (pl.isParameter("Target Residual Reduction")) {
    settings.targetResidualReduction = pl.get<double>("Target Residual Reduction");
  }

  if (pl.isParameter("Number of Successful Applications Before Resetting")) {
    settings.numAppliesCycle = pl.get<int>("Number of Successful Applications Before Resetting");
  }

  const std::string inverse_type          = "Inverse Type";
  const std::string preconditioner_type   = "Preconditioner Type";
  const std::string maximum_size_subblock = "Maximum Size";
  auto invLib                             = getInverseLibrary();

  std::set<int> positions;
  for (const auto &entry : pl) {
    auto fieldName = entry.first;

    // figure out what the integer is
    bool isInverse = fieldName.find(inverse_type) != std::string::npos;
    if (!isInverse) continue;

    int position = -1;
    std::string inverse, type;

    // figure out position
    std::stringstream ss(fieldName);
    ss >> inverse >> type >> position;

    if (position <= 0) {
      Teko_DEBUG_MSG("Adaptive \"Inverse Type\" must be a (strictly) positive integer", 1);
    }

    positions.insert(position);
  }

  inverses.resize(positions.size());
  preconditioners.resize(positions.size());
  maximumSizes.resize(positions.size());
  std::fill(maximumSizes.begin(), maximumSizes.end(), std::numeric_limits<SizeOrdinalType>::max());

  // check for individual solvers
  for (const auto &position : positions) {
    auto inverseName        = inverse_type + std::string(" ") + std::to_string(position);
    auto preconditionerName = preconditioner_type + std::string(" ") + std::to_string(position);
    auto maximumSizeName    = maximum_size_subblock + std::string(" ") + std::to_string(position);

    // inverse or preconditioner
    const auto &invStr     = pl.get<std::string>(inverseName);
    inverses[position - 1] = invLib->getInverseFactory(invStr);

    if (pl.isParameter(preconditionerName)) {
      const auto &precStr           = pl.get<std::string>(preconditionerName);
      preconditioners[position - 1] = invLib->getInverseFactory(precStr);
    }

    if (pl.isParameter(maximumSizeName)) {
      const auto maxSizeSubblock = pl.get<SizeOrdinalType>(maximumSizeName);
      maximumSizes[position - 1] = maxSizeSubblock;
    }
  }
}

}  // end namespace Teko