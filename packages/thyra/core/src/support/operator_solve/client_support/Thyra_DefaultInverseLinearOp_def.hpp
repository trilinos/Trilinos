// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_INVERSE_LINEAR_OP_DEF_HPP
#define THYRA_DEFAULT_INVERSE_LINEAR_OP_DEF_HPP

#include "Thyra_DefaultInverseLinearOp_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_TypeNameTraits.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template <class Scalar>
DefaultInverseLinearOp<Scalar>::DefaultInverseLinearOp() {}

template <class Scalar>
DefaultInverseLinearOp<Scalar>::DefaultInverseLinearOp(
    const Teuchos::RCP<LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria,
    const EThrowOnSolveFailure throwOnFwdSolveFailure,
    const SolveCriteria<Scalar> *adjSolveCriteria,
    const EThrowOnSolveFailure throwOnAdjSolveFailure) {
  initializeImpl(
      lows, fwdSolveCriteria, throwOnFwdSolveFailure, adjSolveCriteria, throwOnAdjSolveFailure);
}

template <class Scalar>
DefaultInverseLinearOp<Scalar>::DefaultInverseLinearOp(
    const Teuchos::RCP<const LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria,
    const EThrowOnSolveFailure throwOnFwdSolveFailure,
    const SolveCriteria<Scalar> *adjSolveCriteria,
    const EThrowOnSolveFailure throwOnAdjSolveFailure) {
  initializeImpl(
      lows, fwdSolveCriteria, throwOnFwdSolveFailure, adjSolveCriteria, throwOnAdjSolveFailure);
}

template <class Scalar>
void DefaultInverseLinearOp<Scalar>::initialize(
    const Teuchos::RCP<LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria,
    const EThrowOnSolveFailure throwOnFwdSolveFailure,
    const SolveCriteria<Scalar> *adjSolveCriteria,
    const EThrowOnSolveFailure throwOnAdjSolveFailure) {
  initializeImpl(
      lows, fwdSolveCriteria, throwOnFwdSolveFailure, adjSolveCriteria, throwOnAdjSolveFailure);
}

template <class Scalar>
void DefaultInverseLinearOp<Scalar>::initialize(
    const Teuchos::RCP<const LinearOpWithSolveBase<Scalar> > &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria,
    const EThrowOnSolveFailure throwOnFwdSolveFailure,
    const SolveCriteria<Scalar> *adjSolveCriteria,
    const EThrowOnSolveFailure throwOnAdjSolveFailure) {
  initializeImpl(
      lows, fwdSolveCriteria, throwOnFwdSolveFailure, adjSolveCriteria, throwOnAdjSolveFailure);
}

template <class Scalar>
void DefaultInverseLinearOp<Scalar>::uninitialize() {
  lows_.uninitialize();
  fwdSolveCriteria_ = Teuchos::null;
  adjSolveCriteria_ = Teuchos::null;
}

// Overridden form InverseLinearOpBase

template <class Scalar>
bool DefaultInverseLinearOp<Scalar>::isLowsConst() const {
  return lows_.isConst();
}

template <class Scalar>
Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
DefaultInverseLinearOp<Scalar>::getNonconstLows() {
  return lows_.getNonconstObj();
}

template <class Scalar>
Teuchos::RCP<const LinearOpWithSolveBase<Scalar> >
DefaultInverseLinearOp<Scalar>::getLows() const {
  return lows_.getConstObj();
}

// Overridden from LinearOpBase

template <class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DefaultInverseLinearOp<Scalar>::range() const {
  assertInitialized();
  return lows_.getConstObj()->domain();
}

template <class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
DefaultInverseLinearOp<Scalar>::domain() const {
  assertInitialized();
  return lows_.getConstObj()->range();
}

template <class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultInverseLinearOp<Scalar>::clone() const {
  return Teuchos::null;  // Not supported yet but could be!
}

// Overridden from Teuchos::Describable

template <class Scalar>
std::string DefaultInverseLinearOp<Scalar>::description() const {
  assertInitialized();
  std::ostringstream oss;
  oss
      << Teuchos::Describable::description() << "{"
      << "lows=" << lows_.getConstObj()->description()
      << ",fwdSolveCriteria=" << (fwdSolveCriteria_.get() ? "..." : "DEFAULT")
      << ",adjSolveCriteria=" << (adjSolveCriteria_.get() ? "..." : "DEFAULT")
      << "}";
  return oss.str();
}

template <class Scalar>
void DefaultInverseLinearOp<Scalar>::describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel) const {
  using Teuchos::OSTab;
  using Teuchos::RCP;
  assertInitialized();
  OSTab tab(out);
  switch (verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME: {
      out
          << Teuchos::Describable::description() << "{"
          << "rangeDim=" << this->range()->dim()
          << ",domainDim=" << this->domain()->dim() << "}:\n";
      OSTab tab2(out);
      out << "lows = ";
      if (!lows_.getConstObj().get()) {
        out << " NULL\n";
      } else {
        out << Teuchos::describe(*lows_.getConstObj(), verbLevel);
      }
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);  // Should never be called!
  }
}

// protected

// Overridden from LinearOpBase

template <class Scalar>
bool DefaultInverseLinearOp<Scalar>::opSupportedImpl(EOpTransp M_trans) const {
  if (nonnull(lows_)) {
    return solveSupports(*lows_.getConstObj(), M_trans);
  }
  return false;
}

template <class Scalar>
void DefaultInverseLinearOp<Scalar>::applyImpl(
    const EOpTransp M_trans,
    const MultiVectorBase<Scalar> &X,
    const Ptr<MultiVectorBase<Scalar> > &Y,
    const Scalar alpha,
    const Scalar beta) const {
  typedef Teuchos::ScalarTraits<Scalar> ST;
  assertInitialized();
  // ToDo: Put in hooks for propogating verbosity level
  //
  // Y = alpha*op(M)*X + beta*Y
  //
  //   =>
  //
  // Y = beta*Y
  // Y += alpha*inv(op(lows))*X
  //
  Teuchos::RCP<MultiVectorBase<Scalar> > T;
  if (beta == ST::zero()) {
    T = Teuchos::rcpFromPtr(Y);
  } else {
    T = createMembers(Y->range(), Y->domain()->dim());
    scale(beta, Y);
  }
  //
  const Ptr<const SolveCriteria<Scalar> > solveCriteria =
      (real_trans(M_trans) == NOTRANS
           ? fwdSolveCriteria_.ptr()
           : adjSolveCriteria_.ptr());
  assign(T.ptr(), ST::zero());  // Have to initialize before solve!
  SolveStatus<Scalar> solveStatus =
      Thyra::solve<Scalar>(*lows_.getConstObj(), M_trans, X, T.ptr(), solveCriteria);

  TEUCHOS_TEST_FOR_EXCEPTION(
      nonnull(solveCriteria) && solveStatus.solveStatus != SOLVE_STATUS_CONVERGED && (real_trans(M_trans) == NOTRANS ? throwOnFwdSolveFailure_ == THROW_ON_SOLVE_FAILURE : throwOnAdjSolveFailure_ == THROW_ON_SOLVE_FAILURE), CatastrophicSolveFailure, "Error, the LOWS object " << lows_.getConstObj()->description() << " returned an unconverged"
                                                                                                                                                                                                                                                                                                                            "status of "
                                                                                                                                                                                                                                                                                   << toString(solveStatus.solveStatus) << " with the message " << solveStatus.message << ".");
  //
  if (beta == ST::zero()) {
    scale(alpha, Y);
  } else {
    update(alpha, *T, Y);
  }
}

// private

template <class Scalar>
template <class LOWS>
void DefaultInverseLinearOp<Scalar>::initializeImpl(
    const Teuchos::RCP<LOWS> &lows,
    const SolveCriteria<Scalar> *fwdSolveCriteria,
    const EThrowOnSolveFailure throwOnFwdSolveFailure,
    const SolveCriteria<Scalar> *adjSolveCriteria,
    const EThrowOnSolveFailure throwOnAdjSolveFailure) {
  lows_.initialize(lows);
  if (fwdSolveCriteria)
    fwdSolveCriteria_ = Teuchos::rcp(new SolveCriteria<Scalar>(*fwdSolveCriteria));
  else
    fwdSolveCriteria_ = Teuchos::null;
  if (adjSolveCriteria)
    adjSolveCriteria_ = Teuchos::rcp(new SolveCriteria<Scalar>(*adjSolveCriteria));
  else
    adjSolveCriteria_ = Teuchos::null;
  throwOnFwdSolveFailure_     = throwOnFwdSolveFailure;
  throwOnAdjSolveFailure_     = throwOnAdjSolveFailure;
  const std::string lowsLabel = lows_.getConstObj()->getObjectLabel();
  if (lowsLabel.length())
    this->setObjectLabel("inv(" + lowsLabel + ")");
}

}  // end namespace Thyra

// Related non-member functions

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstInverse(
    const RCP<LinearOpWithSolveBase<Scalar> > &A,
    const Ptr<const SolveCriteria<Scalar> > &fwdSolveCriteria,
    const EThrowOnSolveFailure throwOnFwdSolveFailure,
    const Ptr<const SolveCriteria<Scalar> > &adjSolveCriteria,
    const EThrowOnSolveFailure throwOnAdjSolveFailure) {
  return Teuchos::rcp(
      new DefaultInverseLinearOp<Scalar>(
          A, fwdSolveCriteria.get(), throwOnFwdSolveFailure,
          adjSolveCriteria.get(), throwOnAdjSolveFailure));
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::inverse(
    const RCP<const LinearOpWithSolveBase<Scalar> > &A,
    const Ptr<const SolveCriteria<Scalar> > &fwdSolveCriteria,
    const EThrowOnSolveFailure throwOnFwdSolveFailure,
    const Ptr<const SolveCriteria<Scalar> > &adjSolveCriteria,
    const EThrowOnSolveFailure throwOnAdjSolveFailure) {
  return Teuchos::rcp(
      new DefaultInverseLinearOp<Scalar>(
          A, fwdSolveCriteria.get(), throwOnFwdSolveFailure,
          adjSolveCriteria.get(), throwOnAdjSolveFailure));
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Thyra namespace!
//

#define THYRA_DEFAULT_INVERSE_LINEAR_OP_INSTANT(SCALAR)          \
                                                                 \
  template class DefaultInverseLinearOp<SCALAR>;                 \
                                                                 \
  template RCP<LinearOpBase<SCALAR> >                            \
  nonconstInverse(                                               \
      const RCP<LinearOpWithSolveBase<SCALAR> > &A,              \
      const Ptr<const SolveCriteria<SCALAR> > &fwdSolveCriteria, \
      const EThrowOnSolveFailure throwOnFwdSolveFailure,         \
      const Ptr<const SolveCriteria<SCALAR> > &adjSolveCriteria, \
      const EThrowOnSolveFailure throwOnAdjSolveFailure);        \
                                                                 \
  template RCP<LinearOpBase<SCALAR> >                            \
  inverse(                                                       \
      const RCP<const LinearOpWithSolveBase<SCALAR> > &A,        \
      const Ptr<const SolveCriteria<SCALAR> > &fwdSolveCriteria, \
      const EThrowOnSolveFailure throwOnFwdSolveFailure,         \
      const Ptr<const SolveCriteria<SCALAR> > &adjSolveCriteria, \
      const EThrowOnSolveFailure throwOnAdjSolveFailure);

#endif  // THYRA_DEFAULT_INVERSE_LINEAR_OP_DEF_HPP
