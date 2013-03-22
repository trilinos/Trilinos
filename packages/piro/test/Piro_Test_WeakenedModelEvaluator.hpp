/*
// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER
*/

#ifndef PIRO_TEST_WEAKENEDMODELEVALUATOR_HPP
#define PIRO_TEST_WEAKENEDMODELEVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

namespace Piro {

namespace Test {

/** \brief Simple ModelEvaluator wrapper with MultiVector-based DgDx derivative disabled */
class WeakenedModelEvaluator_NoDgDxMv : public Thyra::ModelEvaluatorDelegatorBase<double> {
public:
  explicit WeakenedModelEvaluator_NoDgDxMv(const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model)
  {}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const {
    const Thyra::ModelEvaluatorBase::DerivativeSupport expected_support =
      Thyra::ModelEvaluatorBase::DerivativeSupport(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
    TEUCHOS_ASSERT(expected_support.isSameSupport(outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0)));
    ModelEvaluatorBase::OutArgs<double> forwardedOutArgs = getUnderlyingModel()->createOutArgs();
    forwardedOutArgs.setArgs(outArgs);
    getUnderlyingModel()->evalModel(inArgs, forwardedOutArgs);
  }
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  virtual ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const {
    ModelEvaluatorBase::OutArgsSetup<double> outArgs = getUnderlyingModel()->createOutArgs();
    outArgs.setModelEvalDescription(this->description());
    const Thyra::ModelEvaluatorBase::DerivativeSupport newSupport(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0, newSupport);
    return outArgs;
  }
  //@}
};

/** \brief Simple ModelEvaluator wrapper with DgDp derivative disabled */
class WeakenedModelEvaluator_NoDgDp : public Thyra::ModelEvaluatorDelegatorBase<double> {
public:
  explicit WeakenedModelEvaluator_NoDgDp(const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model),
    j_(0),
    l_(0)
  {}

  WeakenedModelEvaluator_NoDgDp(
      const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model,
      int j,
      int l) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model),
    j_(j),
    l_(l)
  {}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const {
    const Thyra::ModelEvaluatorBase::DerivativeSupport expected_support;
    TEUCHOS_ASSERT(expected_support.isSameSupport(outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j_, l_)));
    ModelEvaluatorBase::OutArgs<double> forwardedOutArgs = getUnderlyingModel()->createOutArgs();
    forwardedOutArgs.setArgs(outArgs);
    getUnderlyingModel()->evalModel(inArgs, forwardedOutArgs);
  }
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  virtual ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const {
    ModelEvaluatorBase::OutArgsSetup<double> outArgs = getUnderlyingModel()->createOutArgs();
    outArgs.setModelEvalDescription(this->description());
    const Thyra::ModelEvaluatorBase::DerivativeSupport noSupport;
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j_, l_, noSupport);
    return outArgs;
  }
  //@}

  int j_, l_;
};

/** \brief Simple ModelEvaluator wrapper with MultiVector-based DgDp derivative disabled */
class WeakenedModelEvaluator_NoDgDpMv : public Thyra::ModelEvaluatorDelegatorBase<double> {
public:
  explicit WeakenedModelEvaluator_NoDgDpMv(const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model),
    j_(0),
    l_(0)
  {}

  WeakenedModelEvaluator_NoDgDpMv(
      const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model,
      int j,
      int l) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model),
    j_(j),
    l_(l)
  {}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const {
    const Thyra::ModelEvaluatorBase::DerivativeSupport expected_support =
      Thyra::ModelEvaluatorBase::DerivativeSupport(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
    TEUCHOS_ASSERT(expected_support.isSameSupport(outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j_, l_)));
    ModelEvaluatorBase::OutArgs<double> forwardedOutArgs = getUnderlyingModel()->createOutArgs();
    forwardedOutArgs.setArgs(outArgs);
    getUnderlyingModel()->evalModel(inArgs, forwardedOutArgs);
  }
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  virtual ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const {
    ModelEvaluatorBase::OutArgsSetup<double> outArgs = getUnderlyingModel()->createOutArgs();
    outArgs.setModelEvalDescription(this->description());
    const Thyra::ModelEvaluatorBase::DerivativeSupport newSupport(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j_, l_, newSupport);
    return outArgs;
  }
  //@}

  int j_, l_;
};

/** \brief Simple ModelEvaluator wrapper with Jacobian-form MultiVector-based DgDp derivative disabled */
class WeakenedModelEvaluator_NoDgDpMvJac : public Thyra::ModelEvaluatorDelegatorBase<double> {
public:
  explicit WeakenedModelEvaluator_NoDgDpMvJac(const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model)
  {}

  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  virtual ModelEvaluatorBase::OutArgs<double> createOutArgs() const {
    return WeakenedModelEvaluator_NoDgDpMvJac::createOutArgsImpl();
  }
  //@}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const {
    const Thyra::ModelEvaluatorBase::DerivativeSupport expected_support =
      Thyra::ModelEvaluatorBase::DerivativeSupport(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
    TEUCHOS_ASSERT(expected_support.isSameSupport(outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0)));
    ModelEvaluatorBase::OutArgs<double> forwardedOutArgs = getUnderlyingModel()->createOutArgs();
    forwardedOutArgs.setArgs(outArgs);
    getUnderlyingModel()->evalModel(inArgs, forwardedOutArgs);
  }
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  virtual ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const {
    ModelEvaluatorBase::OutArgsSetup<double> outArgs = getUnderlyingModel()->createOutArgs();
    outArgs.setModelEvalDescription(this->description());
    const Thyra::ModelEvaluatorBase::DerivativeSupport newSupport(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0, newSupport);
    return outArgs;
  }
  //@}
};

/** \brief Simple ModelEvaluator wrapper with MultiVector-based DfDp derivative disabled */
class WeakenedModelEvaluator_NoDfDpMv : public Thyra::ModelEvaluatorDelegatorBase<double> {
public:
  explicit WeakenedModelEvaluator_NoDfDpMv(const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model)
  {}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const {
    const Thyra::ModelEvaluatorBase::DerivativeSupport expected_support =
      Thyra::ModelEvaluatorBase::DerivativeSupport(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
    TEUCHOS_ASSERT(expected_support.isSameSupport(outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 0)));
    ModelEvaluatorBase::OutArgs<double> forwardedOutArgs = getUnderlyingModel()->createOutArgs();
    forwardedOutArgs.setArgs(outArgs);
    getUnderlyingModel()->evalModel(inArgs, forwardedOutArgs);
  }
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  virtual ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const {
    ModelEvaluatorBase::OutArgsSetup<double> outArgs = getUnderlyingModel()->createOutArgs();
    outArgs.setModelEvalDescription(this->description());
    const Thyra::ModelEvaluatorBase::DerivativeSupport newSupport(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 0, newSupport);
    return outArgs;
  }
  //@}
};

/** \brief Simple ModelEvaluator wrapper with Jacobian operator disabled */
class WeakenedModelEvaluator_NoW : public Thyra::ModelEvaluatorDelegatorBase<double> {
public:
  explicit WeakenedModelEvaluator_NoW(const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model)
  {}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const {
    TEUCHOS_ASSERT(!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W));
    TEUCHOS_ASSERT(!outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op));

    ModelEvaluatorBase::OutArgs<double> forwardedOutArgs = getUnderlyingModel()->createOutArgs();
    forwardedOutArgs.setArgs(outArgs);
    getUnderlyingModel()->evalModel(inArgs, forwardedOutArgs);
  }
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  virtual ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const {
    ModelEvaluatorBase::OutArgsSetup<double> outArgs = getUnderlyingModel()->createOutArgs();
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W, false);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op, false);
    return outArgs;
  }
  //@}
};

class NoAdjointLows : public Thyra::LinearOpWithSolveBase<double> {
public:
  explicit NoAdjointLows(const Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > &underlyingLows) :
    underlyingLows_(underlyingLows)
  {}

  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > getUnderlyingLows() {
    return underlyingLows_;
  }

  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domain() const { return underlyingLows_->domain(); }
  virtual Teuchos::RCP<const Thyra::VectorSpaceBase<double> > range() const { return underlyingLows_->range(); }

protected:
  virtual bool opSupportedImpl(Thyra::EOpTransp M_trans) const {
    return underlyingLows_->opSupported(M_trans);
  }

  virtual void applyImpl(
      const Thyra::EOpTransp M_trans,
      const Thyra::MultiVectorBase<double> &X,
      const Teuchos::Ptr<Thyra::MultiVectorBase<double> > &Y,
      const double alpha,
      const double beta) const {
    underlyingLows_->apply(M_trans, X, Y, alpha, beta);
  }

  virtual bool solveSupports(const Thyra::EOpTransp transp) const {
    return transp == Thyra::NOTRANS;
  }

  virtual Thyra::SolveStatus<double> solveImpl(
      const Thyra::EOpTransp transp,
      const Thyra::MultiVectorBase<double> &B,
      const Teuchos::Ptr<Thyra::MultiVectorBase<double> > &X,
      const Teuchos::Ptr<const Thyra::SolveCriteria<double> > solveCriteria) const {
    TEUCHOS_ASSERT(this->solveSupports(transp));
    return underlyingLows_->solve(transp, B, X, solveCriteria);
  }

private:
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > underlyingLows_;
};

/** \brief Simple ModelEvaluator wrapper with adjoint Jacobian solver disabled */
class WeakenedModelEvaluator_NoAdjointW : public Thyra::ModelEvaluatorDelegatorBase<double> {
public:
  explicit WeakenedModelEvaluator_NoAdjointW(const Teuchos::RCP<Thyra::ModelEvaluator<double> > &model) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model)
  {}

  virtual Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > create_W() const {
    const Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > W_underlying = getUnderlyingModel()->create_W();
    return Teuchos::rcp(new NoAdjointLows(W_underlying));
  }

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  virtual void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const {
    ModelEvaluatorBase::OutArgs<double> forwardedOutArgs = getUnderlyingModel()->createOutArgs();
    forwardedOutArgs.setArgs(outArgs);
    if (Teuchos::nonnull(outArgs.get_W())) {
      const Teuchos::RCP<NoAdjointLows> W_downcasted = Teuchos::rcp_dynamic_cast<NoAdjointLows>(outArgs.get_W());
      W_downcasted.assert_not_null();
      forwardedOutArgs.set_W(W_downcasted->getUnderlyingLows());
    }
    getUnderlyingModel()->evalModel(inArgs, forwardedOutArgs);
  }
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorDelegatorBase . */
  //@{
  virtual ModelEvaluatorBase::OutArgs<double> createOutArgsImpl() const {
    ModelEvaluatorBase::OutArgsSetup<double> outArgs = getUnderlyingModel()->createOutArgs();
    outArgs.setModelEvalDescription(this->description());
    {
      Thyra::ModelEvaluatorBase::DerivativeProperties W_properties = outArgs.get_W_properties();
      W_properties.supportsAdjoint = false;
      outArgs.set_W_properties(W_properties);
    }
    return outArgs;
  }
  //@}
};

} // namespace Test

} // namespace Piro

#endif /* PIRO_TEST_WEAKENEDMODELEVALUATOR_HPP */
