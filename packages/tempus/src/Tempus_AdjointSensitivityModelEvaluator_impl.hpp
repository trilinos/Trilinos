//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_AdjointSensitivityModelEvaluator_impl_hpp
#define Tempus_AdjointSensitivityModelEvaluator_impl_hpp

#include "Teuchos_TimeMonitor.hpp"
#include "Tempus_config.hpp"

#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_MultiVectorLinearOp.hpp"
#include "Thyra_MultiVectorLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAdjointLinearOpWithSolve.hpp"
#include "Thyra_AdjointLinearOpWithSolveFactory.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace Tempus {

template <typename Scalar>
AdjointSensitivityModelEvaluator<Scalar>::AdjointSensitivityModelEvaluator(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
        adjoint_residual_model,
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
        adjoint_solve_model,
    const Scalar& t_init, const Scalar& t_final, const bool is_pseudotransient,
    const Teuchos::RCP<const Teuchos::ParameterList>& pList)
  : model_(model),
    adjoint_residual_model_(adjoint_residual_model),
    adjoint_solve_model_(adjoint_solve_model),
    t_init_(t_init),
    t_final_(t_final),
    is_pseudotransient_(is_pseudotransient),
    mass_matrix_is_computed_(false),
    jacobian_matrix_is_computed_(false),
    response_gradient_is_computed_(false),
    t_interp_(Teuchos::ScalarTraits<Scalar>::rmax())
{
  typedef Thyra::ModelEvaluatorBase MEB;

  // Set parameters
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::rcp(new Teuchos::ParameterList);
  if (pList != Teuchos::null) *pl = *pList;
  pl->validateParametersAndSetDefaults(*this->getValidParameters());
  mass_matrix_is_constant_ = pl->get<bool>("Mass Matrix Is Constant");
  mass_matrix_is_identity_ = pl->get<bool>("Mass Matrix Is Identity");
  p_index_                 = pl->get<int>("Sensitivity Parameter Index", 0);
  g_index_                 = pl->get<int>("Response Function Index", 0);
  num_adjoint_             = model_->get_g_space(g_index_)->dim();

  // We currently do not support a non-constant mass matrix
  TEUCHOS_TEST_FOR_EXCEPTION(
      mass_matrix_is_constant_ == false, std::logic_error,
      "AdjointSensitivityModelEvaluator currently does not support "
          << "non-constant mass matrix df/dx_dot!");

  adjoint_space_ =
      Thyra::multiVectorProductVectorSpace(model_->get_f_space(), num_adjoint_);
  residual_space_ =
      Thyra::multiVectorProductVectorSpace(model_->get_x_space(), num_adjoint_);
  response_space_ = Thyra::multiVectorProductVectorSpace(
      model_->get_p_space(p_index_), num_adjoint_);

  // forward and adjoint models must support same InArgs
  MEB::InArgs<Scalar> me_inArgs = model_->createInArgs();
  me_inArgs.assertSameSupport(adjoint_residual_model_->createInArgs());
  me_inArgs.assertSameSupport(adjoint_solve_model_->createInArgs());

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  inArgs.setSupports(MEB::IN_ARG_t);
  if (me_inArgs.supports(MEB::IN_ARG_x_dot))
    inArgs.setSupports(MEB::IN_ARG_x_dot);
  if (me_inArgs.supports(MEB::IN_ARG_alpha))
    inArgs.setSupports(MEB::IN_ARG_alpha);
  if (me_inArgs.supports(MEB::IN_ARG_beta))
    inArgs.setSupports(MEB::IN_ARG_beta);

  // Support additional parameters for x and xdot
  inArgs.set_Np(me_inArgs.Np());
  prototypeInArgs_ = inArgs;

  MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
  MEB::OutArgs<Scalar> adj_mer_outArgs =
      adjoint_residual_model_->createOutArgs();
  MEB::OutArgs<Scalar> adj_mes_outArgs = adjoint_solve_model_->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(me_inArgs.Np(), 2);
  outArgs.setSupports(MEB::OUT_ARG_f);
  if (adj_mes_outArgs.supports(MEB::OUT_ARG_W_op))
    outArgs.setSupports(MEB::OUT_ARG_W_op);
  prototypeOutArgs_ = outArgs;

  // Adjoint residual ME must support W_op to define adjoint ODE/DAE.
  // Must support alpha, beta if it suports x_dot
  TEUCHOS_ASSERT(me_inArgs.supports(MEB::IN_ARG_x));
  TEUCHOS_ASSERT(adj_mer_outArgs.supports(MEB::OUT_ARG_W_op));
  if (me_inArgs.supports(MEB::IN_ARG_x_dot)) {
    TEUCHOS_ASSERT(me_inArgs.supports(MEB::IN_ARG_alpha));
    TEUCHOS_ASSERT(me_inArgs.supports(MEB::IN_ARG_beta));
  }
  MEB::DerivativeSupport dgdx_support =
      me_outArgs.supports(MEB::OUT_ARG_DgDx, g_index_);
  MEB::DerivativeSupport dgdp_support =
      me_outArgs.supports(MEB::OUT_ARG_DgDp, g_index_, p_index_);
  TEUCHOS_ASSERT(dgdx_support.supports(MEB::DERIV_MV_GRADIENT_FORM));
  TEUCHOS_ASSERT(dgdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM));
}

template <typename Scalar>
void AdjointSensitivityModelEvaluator<Scalar>::setFinalTime(
    const Scalar t_final)
{
  t_final_ = t_final;
}

template <typename Scalar>
void AdjointSensitivityModelEvaluator<Scalar>::setForwardSolutionHistory(
    const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& sh)
{
  sh_ = sh;
  if (is_pseudotransient_)
    forward_state_ = sh_->getCurrentState();
  else {
    t_interp_      = Teuchos::ScalarTraits<Scalar>::rmax();
    forward_state_ = Teuchos::null;
  }

  // Reset computation flags because we have done a new forward integration
  mass_matrix_is_computed_       = false;
  jacobian_matrix_is_computed_   = false;
  response_gradient_is_computed_ = false;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
AdjointSensitivityModelEvaluator<Scalar>::get_p_space(int p) const
{
  TEUCHOS_ASSERT(p < model_->Np());
  return model_->get_p_space(p);
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
AdjointSensitivityModelEvaluator<Scalar>::get_p_names(int p) const
{
  TEUCHOS_ASSERT(p < model_->Np());
  return model_->get_p_names(p);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
AdjointSensitivityModelEvaluator<Scalar>::get_x_space() const
{
  return adjoint_space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
AdjointSensitivityModelEvaluator<Scalar>::get_f_space() const
{
  return residual_space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
AdjointSensitivityModelEvaluator<Scalar>::get_g_space(int j) const
{
  TEUCHOS_ASSERT(j == 0 || j == 1);
  if (j == 0) return response_space_;
  return model_->get_g_space(g_index_);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
AdjointSensitivityModelEvaluator<Scalar>::create_W_op() const
{
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > adjoint_op =
      adjoint_solve_model_->create_W_op();
  if (adjoint_op == Teuchos::null) return Teuchos::null;

  return Thyra::nonconstMultiVectorLinearOp(adjoint_op, num_adjoint_);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
AdjointSensitivityModelEvaluator<Scalar>::get_W_factory() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::LinearOpWithSolveFactoryBase<Scalar> LOWSFB;

  RCP<const LOWSFB> alowsfb = adjoint_solve_model_->get_W_factory();
  if (alowsfb == Teuchos::null)
    return Teuchos::null;  // adjoint_solve_model_ doesn't support W_factory

  return Thyra::multiVectorLinearOpWithSolveFactory(alowsfb, residual_space_,
                                                    adjoint_space_);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
AdjointSensitivityModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
AdjointSensitivityModelEvaluator<Scalar>::getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  MEB::InArgs<Scalar> me_nominal = model_->getNominalValues();
  MEB::InArgs<Scalar> nominal    = this->createInArgs();

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  // Set initial x, x_dot
  RCP<Thyra::VectorBase<Scalar> > x = Thyra::createMember(*adjoint_space_);
  Thyra::assign(x.ptr(), zero);
  nominal.set_x(x);

  if (me_nominal.supports(MEB::IN_ARG_x_dot)) {
    RCP<Thyra::VectorBase<Scalar> > x_dot =
        Thyra::createMember(*adjoint_space_);
    Thyra::assign(x_dot.ptr(), zero);
    nominal.set_x_dot(x_dot);
  }

  const int np = model_->Np();
  for (int i = 0; i < np; ++i) nominal.set_p(i, me_nominal.get_p(i));

  return nominal;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
AdjointSensitivityModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}

template <typename Scalar>
void AdjointSensitivityModelEvaluator<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  TEMPUS_FUNC_TIME_MONITOR_DIFF(
      "Tempus::AdjointSensitivityModelEvaluator::evalModel()", TEMPUS_EVAL);

  // Note:  adjoint models compute the transposed W (either explicitly or
  // implicitly.  Thus we need to always call their evalModel() functions
  // whenever computing the adjoint operators, and subsequent calls to apply()
  // do not transpose them.

  // Interpolate forward solution at supplied time, reusing previous
  // interpolation if possible
  TEUCHOS_ASSERT(sh_ != Teuchos::null);
  const Scalar t = inArgs.get_t();
  Scalar forward_t;
  if (is_pseudotransient_)
    forward_t = forward_state_->getTime();
  else {
    forward_t = t_final_ + t_init_ - t;
    if (forward_state_ == Teuchos::null || t_interp_ != t) {
      if (forward_state_ == Teuchos::null)
        forward_state_ = sh_->interpolateState(forward_t);
      else
        sh_->interpolateState(forward_t, forward_state_.get());
      t_interp_ = t;
    }
  }

  // setup input arguments for model
  MEB::InArgs<Scalar> me_inArgs = model_->getNominalValues();
  me_inArgs.set_x(forward_state_->getX());
  if (me_inArgs.supports(MEB::IN_ARG_x_dot)) {
    if (inArgs.get_x_dot() != Teuchos::null)
      me_inArgs.set_x_dot(forward_state_->getXDot());
    else {
      if (is_pseudotransient_) {
        // For pseudo-transient, we have to always use the same form of the
        // residual in order to reuse df/dx, df/dx_dot,..., so we force
        // the underlying ME to always compute the implicit form with x_dot == 0
        // if it wasn't provided.
        if (my_x_dot_ == Teuchos::null) {
          my_x_dot_ = Thyra::createMember(model_->get_x_space());
          Thyra::assign(my_x_dot_.ptr(), Scalar(0.0));
        }
        me_inArgs.set_x_dot(my_x_dot_);
      }
      else {
        // clear out xdot if it was set in nominalValues to get to ensure we
        // get the explicit form
        me_inArgs.set_x_dot(Teuchos::null);
      }
    }
  }
  if (me_inArgs.supports(MEB::IN_ARG_t)) me_inArgs.set_t(forward_t);
  const int np = me_inArgs.Np();
  for (int i = 0; i < np; ++i) me_inArgs.set_p(i, inArgs.get_p(i));

  // compute adjoint W == model W
  // It would be nice to not reevaluate W in the psuedo-transient case, but
  // it isn't clear how to do this in a clean way.  Probably just need to
  // control that with the nonlinear solver.
  RCP<Thyra::LinearOpBase<Scalar> > op;
  if (outArgs.supports(MEB::OUT_ARG_W_op)) op = outArgs.get_W_op();
  if (op != Teuchos::null) {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::AdjointSensitivityModelEvaluator::evalModel::W",
        TEMPUS_EVAL_W);
    if (me_inArgs.supports(MEB::IN_ARG_alpha))
      me_inArgs.set_alpha(inArgs.get_alpha());
    if (me_inArgs.supports(MEB::IN_ARG_beta)) {
      if (me_inArgs.supports(MEB::IN_ARG_x_dot))
        me_inArgs.set_beta(inArgs.get_beta());  // Implicit form (see below)
      else
        me_inArgs.set_beta(-inArgs.get_beta());  // Explicit form (see below)
    }

    RCP<Thyra::MultiVectorLinearOp<Scalar> > mv_adjoint_op =
        rcp_dynamic_cast<Thyra::MultiVectorLinearOp<Scalar> >(op, true);
    RCP<Thyra::LinearOpBase<Scalar> > adjoint_op =
        mv_adjoint_op->getNonconstLinearOp();
    MEB::OutArgs<Scalar> adj_me_outArgs = adjoint_solve_model_->createOutArgs();
    adj_me_outArgs.set_W_op(adjoint_op);
    adjoint_solve_model_->evalModel(me_inArgs, adj_me_outArgs);
  }

  RCP<Thyra::VectorBase<Scalar> > adjoint_f = outArgs.get_f();
  RCP<Thyra::VectorBase<Scalar> > adjoint_g = outArgs.get_g(0);
  RCP<Thyra::VectorBase<Scalar> > g         = outArgs.get_g(1);
  RCP<const Thyra::MultiVectorBase<Scalar> > adjoint_x_mv;
  RCP<const Thyra::VectorBase<Scalar> > adjoint_x;
  if (adjoint_f != Teuchos::null || adjoint_g != Teuchos::null) {
    adjoint_x = inArgs.get_x().assert_not_null();
    adjoint_x_mv =
        rcp_dynamic_cast<const DMVPV>(adjoint_x, true)->getMultiVector();
  }

  // Compute adjoint residual F(y):
  //   * For implicit form,  F(y) = d/dt( df/dx_dot^T*y ) + df/dx^T*y - dg/dx^T
  //   * For explict form, F(y) = -df/dx^T*y + dg/dx^T
  // For implicit form, we assume df/dx_dot is constant w.r.t. x, x_dot, and t,
  // so the residual becomes F(y) = df/dx_dot^T*y_dot + df/dx^T*y - dg/dx^T
  if (adjoint_f != Teuchos::null) {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::AdjointSensitivityModelEvaluator::evalModel::f",
        TEMPUS_EVAL_F);

    RCP<Thyra::MultiVectorBase<Scalar> > adjoint_f_mv =
        rcp_dynamic_cast<DMVPV>(adjoint_f, true)->getNonconstMultiVector();

    MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
    MEB::OutArgs<Scalar> adj_me_outArgs =
        adjoint_residual_model_->createOutArgs();

    // dg/dx^T
    // Don't re-evaluate dg/dx for pseudotransient
    {
      TEMPUS_FUNC_TIME_MONITOR_DIFF(
          "Tempus::AdjointSensitivityModelEvaluator::evalModel::dg/dx",
          TEMPUS_EVAL_DGDX);
      if (my_dgdx_mv_ == Teuchos::null)
        my_dgdx_mv_ = Thyra::createMembers(
            model_->get_x_space(), model_->get_g_space(g_index_)->dim());
      if (!response_gradient_is_computed_) {
        me_outArgs.set_DgDx(
            g_index_,
            MEB::Derivative<Scalar>(my_dgdx_mv_, MEB::DERIV_MV_GRADIENT_FORM));
        model_->evalModel(me_inArgs, me_outArgs);
        me_outArgs.set_DgDx(g_index_, MEB::Derivative<Scalar>());
        if (is_pseudotransient_) response_gradient_is_computed_ = true;
      }
      Thyra::assign(adjoint_f_mv.ptr(), *my_dgdx_mv_);
    }

    // Explicit form of the residual F(y) = -df/dx^T*y + dg/dx^T
    // Don't re-evaluate df/dx for pseudotransient
    {
      TEMPUS_FUNC_TIME_MONITOR_DIFF(
          "Tempus::AdjointSensitivityModelEvaluator::evalModel::df/dx",
          TEMPUS_EVAL_DFDX);
      if (my_dfdx_ == Teuchos::null)
        my_dfdx_ = adjoint_residual_model_->create_W_op();
      if (!jacobian_matrix_is_computed_) {
        adj_me_outArgs.set_W_op(my_dfdx_);
        if (me_inArgs.supports(MEB::IN_ARG_alpha)) me_inArgs.set_alpha(0.0);
        if (me_inArgs.supports(MEB::IN_ARG_beta)) me_inArgs.set_beta(1.0);
        adjoint_residual_model_->evalModel(me_inArgs, adj_me_outArgs);
        if (is_pseudotransient_) jacobian_matrix_is_computed_ = true;
      }
      my_dfdx_->apply(Thyra::NOTRANS, *adjoint_x_mv, adjoint_f_mv.ptr(),
                      Scalar(-1.0), Scalar(1.0));
    }

    // Implicit form residual F(y) df/dx_dot^T*y_dot + df/dx^T*y - dg/dx^T
    // using the second scalar argument to apply() to change the explicit term
    // above.
    // Don't re-evaluate df/dx_dot for pseudotransient
    if (me_inArgs.supports(MEB::IN_ARG_x_dot)) {
      TEMPUS_FUNC_TIME_MONITOR_DIFF(
          "Tempus::AdjointSensitivityModelEvaluator::evalModel::df/dx_dot",
          TEMPUS_EVAL_DFDXDOT);
      RCP<const Thyra::VectorBase<Scalar> > adjoint_x_dot = inArgs.get_x_dot();
      if (adjoint_x_dot != Teuchos::null) {
        RCP<const Thyra::MultiVectorBase<Scalar> > adjoint_x_dot_mv =
            rcp_dynamic_cast<const DMVPV>(adjoint_x_dot, true)
                ->getMultiVector();
        if (mass_matrix_is_identity_) {
          // F = -F + y_dot
          Thyra::V_StVpV(adjoint_f_mv.ptr(), Scalar(-1.0), *adjoint_f_mv,
                         *adjoint_x_dot_mv);
        }
        else {
          if (my_dfdxdot_ == Teuchos::null)
            my_dfdxdot_ = adjoint_residual_model_->create_W_op();
          if (!mass_matrix_is_computed_) {
            adj_me_outArgs.set_W_op(my_dfdxdot_);
            me_inArgs.set_alpha(1.0);
            me_inArgs.set_beta(0.0);
            adjoint_residual_model_->evalModel(me_inArgs, adj_me_outArgs);
            if (is_pseudotransient_ || mass_matrix_is_constant_)
              mass_matrix_is_computed_ = true;
          }
          my_dfdxdot_->apply(Thyra::NOTRANS, *adjoint_x_dot_mv,
                             adjoint_f_mv.ptr(), Scalar(1.0), Scalar(-1.0));
        }
      }
    }
  }

  // Compute g = dg/dp^T - df/dp^T*y for computing the model parameter term in
  // the adjoint sensitivity formula.
  // We don't add pseudotransient logic here because this part is only
  // evaluated once in that case anyway.
  if (adjoint_g != Teuchos::null) {
    TEMPUS_FUNC_TIME_MONITOR_DIFF(
        "Tempus::AdjointSensitivityModelEvaluator::evalModel::g",
        TEMPUS_EVAL_G);
    RCP<Thyra::MultiVectorBase<Scalar> > adjoint_g_mv =
        rcp_dynamic_cast<DMVPV>(adjoint_g, true)->getNonconstMultiVector();

    MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();

    // dg/dp
    MEB::DerivativeSupport dgdp_support =
        me_outArgs.supports(MEB::OUT_ARG_DgDp, g_index_, p_index_);
    if (dgdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
      me_outArgs.set_DgDp(
          g_index_, p_index_,
          MEB::Derivative<Scalar>(adjoint_g_mv, MEB::DERIV_MV_GRADIENT_FORM));
      model_->evalModel(me_inArgs, me_outArgs);
    }
    else if (dgdp_support.supports(MEB::DERIV_MV_JACOBIAN_FORM)) {
      const int num_g = model_->get_g_space(g_index_)->dim();
      const int num_p = model_->get_p_space(p_index_)->dim();
      RCP<Thyra::MultiVectorBase<Scalar> > dgdp_trans =
          createMembers(model_->get_g_space(g_index_), num_p);
      me_outArgs.set_DgDp(
          g_index_, p_index_,
          MEB::Derivative<Scalar>(dgdp_trans, MEB::DERIV_MV_JACOBIAN_FORM));
      model_->evalModel(me_inArgs, me_outArgs);
      Thyra::DetachedMultiVectorView<Scalar> dgdp_view(*adjoint_g_mv);
      Thyra::DetachedMultiVectorView<Scalar> dgdp_trans_view(*dgdp_trans);
      for (int i = 0; i < num_p; ++i)
        for (int j = 0; j < num_g; ++j) dgdp_view(i, j) = dgdp_trans_view(j, i);
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                 "Invalid dg/dp support");
    me_outArgs.set_DgDp(g_index_, p_index_, MEB::Derivative<Scalar>());

    // dg/dp - df/dp^T*y
    MEB::DerivativeSupport dfdp_support =
        me_outArgs.supports(MEB::OUT_ARG_DfDp, p_index_);
    Thyra::EOpTransp trans = Thyra::CONJTRANS;
    if (dfdp_support.supports(MEB::DERIV_LINEAR_OP)) {
      if (my_dfdp_op_ == Teuchos::null)
        my_dfdp_op_ = model_->create_DfDp_op(p_index_);
      me_outArgs.set_DfDp(p_index_, MEB::Derivative<Scalar>(my_dfdp_op_));
      trans = Thyra::CONJTRANS;
    }
    else if (dfdp_support.supports(MEB::DERIV_MV_JACOBIAN_FORM)) {
      if (my_dfdp_mv_ == Teuchos::null)
        my_dfdp_mv_ = createMembers(model_->get_f_space(),
                                    model_->get_p_space(p_index_)->dim());
      me_outArgs.set_DfDp(
          p_index_,
          MEB::Derivative<Scalar>(my_dfdp_mv_, MEB::DERIV_MV_JACOBIAN_FORM));
      my_dfdp_op_ = my_dfdp_mv_;
      trans       = Thyra::CONJTRANS;
    }
    else if (dfdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
      if (my_dfdp_mv_ == Teuchos::null)
        my_dfdp_mv_ = createMembers(model_->get_p_space(p_index_),
                                    model_->get_f_space()->dim());
      me_outArgs.set_DfDp(
          p_index_,
          MEB::Derivative<Scalar>(my_dfdp_mv_, MEB::DERIV_MV_GRADIENT_FORM));
      my_dfdp_op_ = my_dfdp_mv_;
      trans       = Thyra::CONJ;
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                 "Invalid df/dp support");
    model_->evalModel(me_inArgs, me_outArgs);
    my_dfdp_op_->apply(trans, *adjoint_x_mv, adjoint_g_mv.ptr(), Scalar(-1.0),
                       Scalar(1.0));
  }

  if (g != Teuchos::null) {
    MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
    me_outArgs.set_g(g_index_, g);
    model_->evalModel(me_inArgs, me_outArgs);
  }
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
AdjointSensitivityModelEvaluator<Scalar>::getValidParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->set<int>("Sensitivity Parameter Index", 0);
  pl->set<int>("Response Function Index", 0);
  pl->set<bool>("Mass Matrix Is Constant", true);
  pl->set<bool>("Mass Matrix Is Identity", false);
  return pl;
}

}  // namespace Tempus

#endif
