//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StaggeredForwardSensitivityModelEvaluator_impl_hpp
#define Tempus_StaggeredForwardSensitivityModelEvaluator_impl_hpp

#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_MultiVectorLinearOp.hpp"
#include "Thyra_MultiVectorLinearOpWithSolveFactory.hpp"
#include "Thyra_ReuseLinearOpWithSolveFactory.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace Tempus {

template <typename Scalar>
StaggeredForwardSensitivityModelEvaluator<Scalar>::
    StaggeredForwardSensitivityModelEvaluator(
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
            sens_residual_model,
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
            sens_solve_model,
        const bool is_pseudotransient,
        const Teuchos::RCP<const Teuchos::ParameterList>& pList,
        const Teuchos::RCP<MultiVector>& dxdp_init,
        const Teuchos::RCP<MultiVector>& dx_dotdp_init,
        const Teuchos::RCP<MultiVector>& dx_dotdotdp_init)
  : model_(model),
    sens_residual_model_(sens_residual_model),
    sens_solve_model_(sens_solve_model),
    dxdp_init_(dxdp_init),
    dx_dotdp_init_(dx_dotdp_init),
    dx_dotdotdp_init_(dx_dotdotdp_init),
    p_index_(0),
    g_index_(-1),
    x_tangent_index_(1),
    xdot_tangent_index_(2),
    xdotdot_tangent_index_(3),
    use_dfdp_as_tangent_(false),
    use_dgdp_as_tangent_(false),
    num_param_(0),
    num_response_(0),
    g_offset_(0),
    is_pseudotransient_(is_pseudotransient),
    mass_matrix_is_computed_(false),
    jacobian_matrix_is_computed_(false),
    acceleration_matrix_is_computed_(false),
    residual_sensitivity_is_computed_(false),
    t_interp_(Teuchos::ScalarTraits<Scalar>::rmax())
{
  typedef Thyra::ModelEvaluatorBase MEB;

  // Set parameters
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::rcp(new Teuchos::ParameterList);
  if (pList != Teuchos::null) *pl = *pList;
  pl->validateParametersAndSetDefaults(*this->getValidParameters());
  use_dfdp_as_tangent_   = pl->get<bool>("Use DfDp as Tangent");
  use_dgdp_as_tangent_   = pl->get<bool>("Use DgDp as Tangent");
  p_index_               = pl->get<int>("Sensitivity Parameter Index");
  g_index_               = pl->get<int>("Response Function Index", -1);
  x_tangent_index_       = pl->get<int>("Sensitivity X Tangent Index");
  xdot_tangent_index_    = pl->get<int>("Sensitivity X-Dot Tangent Index");
  xdotdot_tangent_index_ = pl->get<int>("Sensitivity X-Dot-Dot Tangent Index");

  num_param_ = model_->get_p_space(p_index_)->dim();
  if (g_index_ >= 0) {
    num_response_ = model_->get_g_space(g_index_)->dim();
    g_offset_     = 2;
  }
  dxdp_space_ =
      Thyra::multiVectorProductVectorSpace(model_->get_x_space(), num_param_);
  dfdp_space_ =
      Thyra::multiVectorProductVectorSpace(model_->get_f_space(), num_param_);
  if (g_index_ >= 0)
    dgdp_space_ = Thyra::multiVectorProductVectorSpace(
        model_->get_g_space(g_index_), num_param_);

  // forward and sensitivity models must support same InArgs
  MEB::InArgs<Scalar> me_inArgs = model_->createInArgs();
  me_inArgs.assertSameSupport(sens_residual_model_->createInArgs());
  me_inArgs.assertSameSupport(sens_solve_model_->createInArgs());

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  if (me_inArgs.supports(MEB::IN_ARG_x_dot))
    inArgs.setSupports(MEB::IN_ARG_x_dot);
  if (me_inArgs.supports(MEB::IN_ARG_t)) inArgs.setSupports(MEB::IN_ARG_t);
  if (me_inArgs.supports(MEB::IN_ARG_alpha))
    inArgs.setSupports(MEB::IN_ARG_alpha);
  if (me_inArgs.supports(MEB::IN_ARG_beta))
    inArgs.setSupports(MEB::IN_ARG_beta);
  if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
    inArgs.setSupports(MEB::IN_ARG_W_x_dot_dot_coeff);

  // Support additional parameters for x and xdot
  inArgs.set_Np(me_inArgs.Np());
  prototypeInArgs_ = inArgs;

  MEB::OutArgs<Scalar> me_outArgs       = model_->createOutArgs();
  MEB::OutArgs<Scalar> sens_mer_outArgs = sens_residual_model_->createOutArgs();
  MEB::OutArgs<Scalar> sens_mes_outArgs = sens_solve_model_->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(me_inArgs.Np(), me_outArgs.Ng() + g_offset_);
  outArgs.setSupports(MEB::OUT_ARG_f);
  if (sens_mer_outArgs.supports(MEB::OUT_ARG_W_op) ||
      sens_mes_outArgs.supports(MEB::OUT_ARG_W_op))
    outArgs.setSupports(MEB::OUT_ARG_W_op);

  // Response 0 is the reduced dg/dp (no sensitivities supported)
  // Response 1 is the reduced g (no sensitivities supported)
  for (int j = 0; j < me_outArgs.Ng(); ++j) {
    outArgs.setSupports(MEB::OUT_ARG_DgDx_dot, j + g_offset_,
                        me_outArgs.supports(MEB::OUT_ARG_DgDx_dot, j));
    outArgs.setSupports(MEB::OUT_ARG_DgDx, j + g_offset_,
                        me_outArgs.supports(MEB::OUT_ARG_DgDx, j));
    for (int l = 0; l < me_outArgs.Np(); ++l) {
      outArgs.setSupports(MEB::OUT_ARG_DgDp, j + g_offset_, l,
                          me_outArgs.supports(MEB::OUT_ARG_DgDp, j, l));
    }
  }
  prototypeOutArgs_ = outArgs;

  // Sensitivity residual ME must support W_op to define adjoint ODE/DAE.
  // Must support alpha, beta if it suports x_dot
  TEUCHOS_ASSERT(me_inArgs.supports(MEB::IN_ARG_x));
  if (!use_dfdp_as_tangent_)
    TEUCHOS_ASSERT(sens_mer_outArgs.supports(MEB::OUT_ARG_W_op));
  if (me_inArgs.supports(MEB::IN_ARG_x_dot)) {
    TEUCHOS_ASSERT(me_inArgs.supports(MEB::IN_ARG_alpha));
    TEUCHOS_ASSERT(me_inArgs.supports(MEB::IN_ARG_beta));
  }
  TEUCHOS_ASSERT(me_outArgs.supports(MEB::OUT_ARG_DfDp, p_index_)
                     .supports(MEB::DERIV_MV_JACOBIAN_FORM));
}

template <typename Scalar>
void StaggeredForwardSensitivityModelEvaluator<Scalar>::
    setForwardSolutionHistory(
        const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& sh)
{
  sh_       = sh;
  t_interp_ = Teuchos::ScalarTraits<Scalar>::rmax();
}

template <typename Scalar>
void StaggeredForwardSensitivityModelEvaluator<Scalar>::setForwardSolutionState(
    const Teuchos::RCP<const Tempus::SolutionState<Scalar> >& s)
{
  sh_            = Teuchos::null;
  forward_state_ = s;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::get_p_space(int p) const
{
  TEUCHOS_ASSERT(p < model_->Np());
  return model_->get_p_space(p);
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::get_p_names(int p) const
{
  TEUCHOS_ASSERT(p < model_->Np());
  return model_->get_p_names(p);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::get_x_space() const
{
  return dxdp_space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::get_f_space() const
{
  return dfdp_space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::get_g_space(int j) const
{
  if (g_index_ >= 0) {
    if (j == 0)
      return dgdp_space_;
    else if (j == 1)
      return model_->get_g_space(g_index_);
  }
  return model_->get_g_space(j - g_offset_);
}

template <typename Scalar>
Teuchos::ArrayView<const std::string>
StaggeredForwardSensitivityModelEvaluator<Scalar>::get_g_names(int j) const
{
  if (g_index_ >= 0) {
    if (j == 0) {
      Teuchos::Array<std::string> names = model_->get_g_names(g_index_);
      for (int i = 0; i < names.size(); ++i)
        names[i] = names[i] + "_reduced sensitivity";
      return names();
    }
    else if (j == 1)
      return model_->get_g_names(g_index_);
  }
  return model_->get_g_names(j - g_offset_);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::create_W_op() const
{
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > op;
  if (lo_ != Teuchos::null)
    op = lo_;
  else
    op = sens_solve_model_->create_W_op();
  return Thyra::nonconstMultiVectorLinearOp(op, num_param_);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::create_DgDx_dot_op(
    int j) const
{
  return model_->create_DgDx_dot_op(j - g_offset_);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::create_DgDx_op(int j) const
{
  return model_->create_DgDx_op(j - g_offset_);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::create_DgDp_op(int j,
                                                                  int l) const
{
  return model_->create_DgDp_op(j - g_offset_, l);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::get_W_factory() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
      model_factory = sens_solve_model_->get_W_factory();
  if (model_factory == Teuchos::null)
    return Teuchos::null;  // model_ doesn't support W_factory
  if (po_ != Teuchos::null) {
    factory = Thyra::reuseLinearOpWithSolveFactory<Scalar>(model_factory, po_);
  }
  else
    factory = model_factory;
  return Thyra::multiVectorLinearOpWithSolveFactory(factory, dfdp_space_,
                                                    dxdp_space_);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StaggeredForwardSensitivityModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StaggeredForwardSensitivityModelEvaluator<Scalar>::getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  MEB::InArgs<Scalar> me_nominal = model_->getNominalValues();
  MEB::InArgs<Scalar> nominal    = this->createInArgs();

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  // Set initial x.  If dxdp_init == null, set initial dx/dp = 0
  RCP<Thyra::VectorBase<Scalar> > x = Thyra::createMember(*dxdp_space_);
  RCP<DMVPV> dxdp                   = rcp_dynamic_cast<DMVPV>(x, true);
  if (dxdp_init_ == Teuchos::null)
    Thyra::assign(dxdp->getNonconstMultiVector().ptr(), zero);
  else
    Thyra::assign(dxdp->getNonconstMultiVector().ptr(), *dxdp_init_);
  nominal.set_x(x);

  // Set initial xdot.  If dx_dotdp_init == null, set initial dxdot/dp = 0
  if (me_nominal.supports(MEB::IN_ARG_x_dot)) {
    RCP<Thyra::VectorBase<Scalar> > xdot = Thyra::createMember(*dxdp_space_);
    RCP<DMVPV> dxdotdp                   = rcp_dynamic_cast<DMVPV>(xdot, true);
    if (dx_dotdp_init_ == Teuchos::null)
      Thyra::assign(dxdotdp->getNonconstMultiVector().ptr(), zero);
    else
      Thyra::assign(dxdotdp->getNonconstMultiVector().ptr(), *dx_dotdp_init_);
    nominal.set_x_dot(xdot);
  }

  // Set initial xdotdot.  If dx_dotdotdp_init == null, set initial dxdotdot/dp
  // = 0
  if (me_nominal.supports(MEB::IN_ARG_x_dot_dot)) {
    RCP<Thyra::VectorBase<Scalar> > xdotdot = Thyra::createMember(*dxdp_space_);
    RCP<DMVPV> dxdotdotdp                   = rcp_dynamic_cast<DMVPV>(xdotdot, true);
    if (dx_dotdotdp_init_ == Teuchos::null)
      Thyra::assign(dxdotdotdp->getNonconstMultiVector().ptr(), zero);
    else
      Thyra::assign(dxdotdotdp->getNonconstMultiVector().ptr(),
                    *dx_dotdotdp_init_);
    nominal.set_x_dot_dot(xdotdot);
  }

  const int np = model_->Np();
  for (int i = 0; i < np; ++i) nominal.set_p(i, me_nominal.get_p(i));
  return nominal;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
StaggeredForwardSensitivityModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}

template <typename Scalar>
void StaggeredForwardSensitivityModelEvaluator<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  // Interpolate forward solution at supplied time, reusing previous
  // interpolation if possible
  Scalar forward_t;
  if (sh_ != Teuchos::null) {
    forward_t = inArgs.get_t();
    if (t_interp_ != forward_t) {
      if (nc_forward_state_ == Teuchos::null)
        nc_forward_state_ = sh_->interpolateState(forward_t);
      else
        sh_->interpolateState(forward_t, nc_forward_state_.get());
      forward_state_ = nc_forward_state_;
      t_interp_      = forward_t;
    }
  }
  else {
    TEUCHOS_ASSERT(forward_state_ != Teuchos::null);
    forward_t = forward_state_->getTime();
  }

  const bool use_tangent = use_dfdp_as_tangent_ || use_dgdp_as_tangent_;

  // setup input arguments for model
  RCP<const Thyra::MultiVectorBase<Scalar> > dxdp, dxdotdp, dxdotdotdp;
  MEB::InArgs<Scalar> me_inArgs = model_->getNominalValues();
  dxdp                          = rcp_dynamic_cast<const DMVPV>(inArgs.get_x(), true)->getMultiVector();
  me_inArgs.set_x(forward_state_->getX());
  if (use_dfdp_as_tangent_) me_inArgs.set_p(x_tangent_index_, inArgs.get_x());
  if (me_inArgs.supports(MEB::IN_ARG_x_dot)) {
    if (inArgs.get_x_dot() != Teuchos::null) {
      dxdotdp = rcp_dynamic_cast<const DMVPV>(inArgs.get_x_dot(), true)
                    ->getMultiVector();
      me_inArgs.set_x_dot(forward_state_->getXDot());
      if (use_dfdp_as_tangent_)
        me_inArgs.set_p(xdot_tangent_index_, inArgs.get_x_dot());
    }
    else {
      // clear out xdot if it was set in nominalValues to get to ensure we
      // get the explicit form
      me_inArgs.set_x_dot(Teuchos::null);
    }
  }
  if (me_inArgs.supports(MEB::IN_ARG_x_dot_dot)) {
    if (inArgs.get_x_dot_dot() != Teuchos::null) {
      dxdotdotdp = rcp_dynamic_cast<const DMVPV>(inArgs.get_x_dot_dot(), true)
                       ->getMultiVector();
      me_inArgs.set_x_dot_dot(forward_state_->getXDotDot());
      if (use_dfdp_as_tangent_)
        me_inArgs.set_p(xdotdot_tangent_index_, inArgs.get_x_dot_dot());
    }
    else  // clear out xdotdot if it was set in nominalValues
      me_inArgs.set_x_dot_dot(Teuchos::null);
  }
  if (me_inArgs.supports(MEB::IN_ARG_t)) me_inArgs.set_t(forward_t);
  if (me_inArgs.supports(MEB::IN_ARG_alpha))
    me_inArgs.set_alpha(inArgs.get_alpha());
  if (me_inArgs.supports(MEB::IN_ARG_beta))
    me_inArgs.set_beta(inArgs.get_beta());
  if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
    me_inArgs.set_W_x_dot_dot_coeff(inArgs.get_W_x_dot_dot_coeff());

  // Set parameters -- be careful to only set ones that were set in our
  // inArgs to not null out any specified through nominalValues or
  // dx/dp above
  const int np = me_inArgs.Np();
  for (int i = 0; i < np; ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      if (!use_tangent ||
          (use_tangent && i != x_tangent_index_ && i != xdot_tangent_index_ &&
           i != xdotdot_tangent_index_))
        me_inArgs.set_p(i, inArgs.get_p(i));
  }

  // setup output arguments for model
  RCP<Thyra::MultiVectorBase<Scalar> > dfdp;
  MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
  if (outArgs.get_f() != Teuchos::null) {
    dfdp = rcp_dynamic_cast<DMVPV>(outArgs.get_f(), true)
               ->getNonconstMultiVector();
    if (!residual_sensitivity_is_computed_) {
      if (!use_dfdp_as_tangent_ && is_pseudotransient_) {
        if (my_dfdp_ == Teuchos::null)
          my_dfdp_ = Thyra::createMembers(model_->get_f_space(),
                                          model_->get_p_space(p_index_)->dim());
        me_outArgs.set_DfDp(p_index_, my_dfdp_);
      }
      else
        me_outArgs.set_DfDp(p_index_, dfdp);
    }
  }
  if (outArgs.supports(MEB::OUT_ARG_W_op) && lo_ == Teuchos::null &&
      outArgs.get_W_op() != Teuchos::null &&
      model_.ptr() == sens_solve_model_.ptr()) {
    RCP<Thyra::LinearOpBase<Scalar> > op = outArgs.get_W_op();
    RCP<Thyra::MultiVectorLinearOp<Scalar> > mv_op =
        rcp_dynamic_cast<Thyra::MultiVectorLinearOp<Scalar> >(op, true);
    me_outArgs.set_W_op(mv_op->getNonconstLinearOp());
  }
  for (int j = g_offset_; j < outArgs.Ng(); ++j) {
    me_outArgs.set_g(j - g_offset_, outArgs.get_g(j));
    if (!me_outArgs.supports(MEB::OUT_ARG_DgDx_dot, j - g_offset_).none())
      me_outArgs.set_DgDx_dot(j - g_offset_, outArgs.get_DgDx_dot(j));
    if (!me_outArgs.supports(MEB::OUT_ARG_DgDx, j - g_offset_).none())
      me_outArgs.set_DgDx(j - g_offset_, outArgs.get_DgDx(j));
    for (int l = 0; l < outArgs.Np(); ++l)
      if (!me_outArgs.supports(MEB::OUT_ARG_DgDp, j - g_offset_, l).none())
        me_outArgs.set_DgDp(j - g_offset_, l, outArgs.get_DgDp(j, l));
  }
  if (g_index_ >= 0 && outArgs.get_g(1) != Teuchos::null)
    me_outArgs.set_g(g_index_, outArgs.get_g(1));

  // build residual and jacobian
  model_->evalModel(me_inArgs, me_outArgs);

  // Compute W_op separately if we have a separate sensitivity solve ME
  if (outArgs.supports(MEB::OUT_ARG_W_op) && lo_ == Teuchos::null &&
      outArgs.get_W_op() != Teuchos::null &&
      model_.ptr() != sens_solve_model_.ptr()) {
    MEB::OutArgs<Scalar> sens_me_outArgs = sens_solve_model_->createOutArgs();
    RCP<Thyra::LinearOpBase<Scalar> > op = outArgs.get_W_op();
    RCP<Thyra::MultiVectorLinearOp<Scalar> > mv_op =
        rcp_dynamic_cast<Thyra::MultiVectorLinearOp<Scalar> >(op, true);
    sens_me_outArgs.set_W_op(mv_op->getNonconstLinearOp());
    sens_solve_model_->evalModel(me_inArgs, sens_me_outArgs);
  }

  // Compute (df/dx) * (dx/dp) + (df/dxdot) * (dxdot/dp) + (df/dxdotdot) *
  // (dxdotdot/dp) + (df/dp) if the underlying ME doesn't already do this.  This
  // requires computing df/dx, df/dxdot, df/dxdotdot as separate operators. For
  // pseudo-transient, we would like to reuse these operators, but this is
  // complicated when steppers use both implicit and explicit forms.
  if (!use_dfdp_as_tangent_) {
    if (dfdp != Teuchos::null && is_pseudotransient_) {
      residual_sensitivity_is_computed_ = true;
      Thyra::assign(dfdp.ptr(), *my_dfdp_);
    }
    if (dxdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdx_ == Teuchos::null)
        my_dfdx_ = sens_residual_model_->create_W_op();
      if (!jacobian_matrix_is_computed_) {
        MEB::OutArgs<Scalar> meo = sens_residual_model_->createOutArgs();
        meo.set_W_op(my_dfdx_);
        if (me_inArgs.supports(MEB::IN_ARG_alpha)) me_inArgs.set_alpha(0.0);
        if (me_inArgs.supports(MEB::IN_ARG_beta)) me_inArgs.set_beta(1.0);
        if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
          me_inArgs.set_W_x_dot_dot_coeff(0.0);
        sens_residual_model_->evalModel(me_inArgs, meo);
        if (is_pseudotransient_) jacobian_matrix_is_computed_ = true;
      }
      my_dfdx_->apply(Thyra::NOTRANS, *dxdp, dfdp.ptr(), Scalar(1.0),
                      Scalar(1.0));
    }
    if (dxdotdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdxdot_ == Teuchos::null)
        my_dfdxdot_ = sens_residual_model_->create_W_op();
      if (!mass_matrix_is_computed_) {
        MEB::OutArgs<Scalar> meo = sens_residual_model_->createOutArgs();
        meo.set_W_op(my_dfdxdot_);
        if (me_inArgs.supports(MEB::IN_ARG_alpha)) me_inArgs.set_alpha(1.0);
        if (me_inArgs.supports(MEB::IN_ARG_beta)) me_inArgs.set_beta(0.0);
        if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
          me_inArgs.set_W_x_dot_dot_coeff(0.0);
        sens_residual_model_->evalModel(me_inArgs, meo);
        if (is_pseudotransient_) mass_matrix_is_computed_ = true;
      }
      my_dfdxdot_->apply(Thyra::NOTRANS, *dxdotdp, dfdp.ptr(), Scalar(1.0),
                         Scalar(1.0));
    }
    if (dxdotdotdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdxdotdot_ == Teuchos::null)
        my_dfdxdotdot_ = sens_residual_model_->create_W_op();
      if (!acceleration_matrix_is_computed_) {
        MEB::OutArgs<Scalar> meo = sens_residual_model_->createOutArgs();
        meo.set_W_op(my_dfdxdotdot_);
        if (me_inArgs.supports(MEB::IN_ARG_alpha)) me_inArgs.set_alpha(0.0);
        if (me_inArgs.supports(MEB::IN_ARG_beta)) me_inArgs.set_beta(0.0);
        if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
          me_inArgs.set_W_x_dot_dot_coeff(1.0);
        sens_residual_model_->evalModel(me_inArgs, meo);
        if (is_pseudotransient_) acceleration_matrix_is_computed_ = true;
      }
      my_dfdxdotdot_->apply(Thyra::NOTRANS, *dxdotdotdp, dfdp.ptr(),
                            Scalar(1.0), Scalar(1.0));
    }
  }

  if (g_index_ >= 0 && outArgs.get_g(0) != Teuchos::null) {
    MEB::OutArgs<Scalar> meo = model_->createOutArgs();
    RCP<Thyra::MultiVectorBase<Scalar> > dgdp =
        rcp_dynamic_cast<DMVPV>(outArgs.get_g(0), true)
            ->getNonconstMultiVector();
    RCP<Thyra::MultiVectorBase<Scalar> > dgdp_trans;
    MEB::DerivativeSupport dgdp_support =
        meo.supports(MEB::OUT_ARG_DgDp, g_index_, p_index_);
    if (dgdp_support.supports(MEB::DERIV_MV_JACOBIAN_FORM))
      meo.set_DgDp(g_index_, p_index_,
                   MEB::Derivative<Scalar>(dgdp, MEB::DERIV_MV_JACOBIAN_FORM));
    else if (dgdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
      dgdp_trans = createMembers(model_->get_p_space(p_index_), num_response_);
      meo.set_DgDp(
          g_index_, p_index_,
          MEB::Derivative<Scalar>(dgdp_trans, MEB::DERIV_MV_GRADIENT_FORM));
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Operator form of dg/dp not supported for reduced sensitivity");

    if (!use_dgdp_as_tangent_ || dxdp == Teuchos::null) {
      MEB::DerivativeSupport dgdx_support =
          meo.supports(MEB::OUT_ARG_DgDx, g_index_);
      if (dgdx_support.supports(MEB::DERIV_LINEAR_OP)) {
        if (my_dgdx_ == Teuchos::null)
          my_dgdx_ = model_->create_DgDx_op(g_index_);
        meo.set_DgDx(g_index_, my_dgdx_);
      }
      else if (dgdx_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
        if (my_dgdx_mv_ == Teuchos::null)
          my_dgdx_mv_ =
              createMembers(model_->get_g_space(g_index_), num_param_);
        meo.set_DgDx(g_index_, MEB::Derivative<Scalar>(
                                   my_dgdx_mv_, MEB::DERIV_MV_GRADIENT_FORM));
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error,
            "Jacobian form of dg/dx not supported for reduced sensitivity");
    }

    // Clear dx/dp, dxdot/dp, dxdotdot/dp from inArgs if set and we are not
    // using dg/dp as the tangent.  Note, even though dxdp, dxdotdp, and
    // dxdotdotdp are not used here, they are non-null if x/xdot/xdotdot are
    // supported and supplied.
    if (!use_dgdp_as_tangent_ && use_dfdp_as_tangent_ && dxdp != Teuchos::null)
      me_inArgs.set_p(x_tangent_index_, Teuchos::null);
    if (!use_dgdp_as_tangent_ && use_dfdp_as_tangent_ &&
        dxdotdp != Teuchos::null)
      me_inArgs.set_p(xdot_tangent_index_, Teuchos::null);
    if (!use_dgdp_as_tangent_ && use_dfdp_as_tangent_ &&
        dxdotdotdp != Teuchos::null)
      me_inArgs.set_p(xdotdot_tangent_index_, Teuchos::null);

    // Set dx/dp, dxdot/dp, dxdodot/dp if necessary
    if (use_dgdp_as_tangent_ && !use_dfdp_as_tangent_ && dxdp != Teuchos::null)
      me_inArgs.set_p(x_tangent_index_, inArgs.get_x());
    if (use_dgdp_as_tangent_ && !use_dfdp_as_tangent_ &&
        dxdotdp != Teuchos::null)
      me_inArgs.set_p(xdot_tangent_index_, inArgs.get_x_dot());
    if (use_dgdp_as_tangent_ && !use_dfdp_as_tangent_ &&
        dxdotdotdp != Teuchos::null)
      me_inArgs.set_p(xdotdot_tangent_index_, inArgs.get_x_dot_dot());

    model_->evalModel(me_inArgs, meo);

    // transpose reduced dg/dp if necessary
    if (dgdp_trans != Teuchos::null) {
      Thyra::DetachedMultiVectorView<Scalar> dgdp_view(*dgdp);
      Thyra::DetachedMultiVectorView<Scalar> dgdp_trans_view(*dgdp_trans);
      for (int i = 0; i < num_param_; ++i)
        for (int j = 0; j < num_response_; ++j)
          dgdp_view(j, i) = dgdp_trans_view(i, j);
    }

    // Compute (dg/dx) * (dx/dp) + (dg/dp) if the underlying ME doesn't already
    // do this.
    if (!use_dgdp_as_tangent_ && dxdp != Teuchos::null) {
      MEB::DerivativeSupport dgdx_support =
          me_outArgs.supports(MEB::OUT_ARG_DgDx, g_index_);
      if (dgdx_support.supports(MEB::DERIV_LINEAR_OP)) {
        my_dgdx_->apply(Thyra::NOTRANS, *dxdp, dgdp.ptr(), Scalar(1.0),
                        Scalar(1.0));
      }
      else if (dgdx_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
        my_dgdx_mv_->apply(Thyra::TRANS, *dxdp, dgdp.ptr(), Scalar(1.0),
                           Scalar(1.0));
      }
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error,
            "Jacobian form of dg/dx not supported for reduced sensitivity");
    }
  }
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StaggeredForwardSensitivityModelEvaluator<Scalar>::getValidParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->set<bool>("Use DfDp as Tangent", false);
  pl->set<bool>("Use DgDp as Tangent", false);
  pl->set<int>("Sensitivity Parameter Index", 0);
  pl->set<int>("Response Function Index", -1);
  pl->set<int>("Sensitivity X Tangent Index", 1);
  pl->set<int>("Sensitivity X-Dot Tangent Index", 2);
  pl->set<int>("Sensitivity X-Dot-Dot Tangent Index", 3);
  return pl;
}

}  // namespace Tempus

#endif
