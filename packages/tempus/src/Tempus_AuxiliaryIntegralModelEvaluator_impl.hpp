//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_AuxiliaryIntegralModelEvaluator_impl_hpp
#define Tempus_AuxiliaryIntegralModelEvaluator_impl_hpp

#include "Thyra_ScaledIdentityLinearOpWithSolve.hpp"
#include "Thyra_ScaledIdentityLinearOpWithSolveFactory.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Tempus {

template <typename Scalar>
AuxiliaryIntegralModelEvaluator<Scalar>::AuxiliaryIntegralModelEvaluator(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
    const int g_index)
  : model_(model),
    g_index_(g_index),
    t_interp_(Teuchos::ScalarTraits<Scalar>::rmax())
{
  typedef Thyra::ModelEvaluatorBase MEB;

  space_ = model_->get_g_space(g_index);

  MEB::InArgs<Scalar> me_inArgs = model_->createInArgs();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  inArgs.setSupports(MEB::IN_ARG_t);
  if (me_inArgs.supports(MEB::IN_ARG_x_dot))
    inArgs.setSupports(MEB::IN_ARG_x_dot);
  inArgs.setSupports(MEB::IN_ARG_alpha);
  inArgs.setSupports(MEB::IN_ARG_beta);
  inArgs.set_Np(me_inArgs.Np());
  prototypeInArgs_ = inArgs;

  MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(me_inArgs.Np(), 0);
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  prototypeOutArgs_ = outArgs;
}

template <typename Scalar>
void AuxiliaryIntegralModelEvaluator<Scalar>::setForwardSolutionHistory(
    const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > &sh)
{
  sh_            = sh;
  t_interp_      = Teuchos::ScalarTraits<Scalar>::rmax();
  forward_state_ = Teuchos::null;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
AuxiliaryIntegralModelEvaluator<Scalar>::get_p_space(int p) const
{
  TEUCHOS_ASSERT(p < model_->Np());
  return model_->get_p_space(p);
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
AuxiliaryIntegralModelEvaluator<Scalar>::get_p_names(int p) const
{
  TEUCHOS_ASSERT(p < model_->Np());
  return model_->get_p_names(p);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
AuxiliaryIntegralModelEvaluator<Scalar>::get_x_space() const
{
  return space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
AuxiliaryIntegralModelEvaluator<Scalar>::get_f_space() const
{
  return space_;
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
AuxiliaryIntegralModelEvaluator<Scalar>::create_W_op() const
{
  return Thyra::scaledIdentity(space_, Scalar(1.0));
}

template <typename Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
AuxiliaryIntegralModelEvaluator<Scalar>::get_W_factory() const
{
  return Thyra::scaledIdentitySolveFactory(space_, Scalar(1.0));
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
AuxiliaryIntegralModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
AuxiliaryIntegralModelEvaluator<Scalar>::getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  MEB::InArgs<Scalar> me_nominal = model_->getNominalValues();
  MEB::InArgs<Scalar> nominal    = this->createInArgs();

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  // Set initial x, x_dot
  RCP<Thyra::VectorBase<Scalar> > x = Thyra::createMember(*space_);
  Thyra::assign(x.ptr(), zero);
  nominal.set_x(x);

  if (me_nominal.supports(MEB::IN_ARG_x_dot)) {
    RCP<Thyra::VectorBase<Scalar> > x_dot = Thyra::createMember(*space_);
    Thyra::assign(x_dot.ptr(), zero);
    nominal.set_x_dot(x_dot);
  }

  const int np = model_->Np();
  for (int i = 0; i < np; ++i) nominal.set_p(i, me_nominal.get_p(i));

  return nominal;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
AuxiliaryIntegralModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}

template <typename Scalar>
void AuxiliaryIntegralModelEvaluator<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // Interpolate forward solution at supplied time, reusing previous
  // interpolation if possible
  TEUCHOS_ASSERT(sh_ != Teuchos::null);
  const Scalar t = inArgs.get_t();
  if (forward_state_ == Teuchos::null || t_interp_ != t) {
    if (forward_state_ == Teuchos::null)
      forward_state_ = sh_->interpolateState(t);
    else
      sh_->interpolateState(t, forward_state_.get());
    t_interp_ = t;
  }

  // Residual
  RCP<Thyra::VectorBase<Scalar> > f = outArgs.get_f();
  if (f != Teuchos::null) {
    MEB::InArgs<Scalar> me_inArgs = model_->getNominalValues();
    me_inArgs.set_x(forward_state_->getX());
    if (me_inArgs.supports(MEB::IN_ARG_x_dot))
      me_inArgs.set_x_dot(forward_state_->getXDot());
    if (me_inArgs.supports(MEB::IN_ARG_t)) me_inArgs.set_t(t);
    const int np = me_inArgs.Np();
    for (int i = 0; i < np; ++i) me_inArgs.set_p(i, inArgs.get_p(i));

    MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
    me_outArgs.set_g(g_index_, f);

    model_->evalModel(me_inArgs, me_outArgs);

    // For explicit form, f = g and we are done
    // For implicit form, f = dz/dt - g
    if (inArgs.supports(MEB::IN_ARG_x_dot)) {
      RCP<const Thyra::VectorBase<Scalar> > x_dot = inArgs.get_x_dot();
      if (x_dot != Teuchos::null) Thyra::V_VmV(f.ptr(), *x_dot, *f);
    }
  }

  // W = alpha*df/z_dot + beta*df/dz = alpha * I
  RCP<Thyra::LinearOpBase<Scalar> > op = outArgs.get_W_op();
  if (op != Teuchos::null) {
    const Scalar alpha = inArgs.get_alpha();
    RCP<Thyra::ScaledIdentityLinearOpWithSolve<Scalar> > si_op =
        rcp_dynamic_cast<Thyra::ScaledIdentityLinearOpWithSolve<Scalar> >(op);
    si_op->setScale(alpha);
  }
}

}  // namespace Tempus

#endif
