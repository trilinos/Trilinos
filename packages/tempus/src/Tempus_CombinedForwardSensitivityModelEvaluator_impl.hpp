//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_CombinedForwardSensitivityModelEvaluator_impl_hpp
#define Tempus_CombinedForwardSensitivityModelEvaluator_impl_hpp

#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_MultiVectorLinearOp.hpp"
#include "Thyra_MultiVectorLinearOpWithSolveFactory.hpp"
#include "Thyra_ReuseLinearOpWithSolveFactory.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace Tempus {

template <typename Scalar>
CombinedForwardSensitivityModelEvaluator<Scalar>::
    CombinedForwardSensitivityModelEvaluator(
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
            sens_residual_model,
        const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >&
            sens_solve_model,
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
    g_offset_(0)
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
  g_index_               = pl->get<int>("Response Function Index");
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
  x_dxdp_space_ = Thyra::multiVectorProductVectorSpace(model_->get_x_space(),
                                                       num_param_ + 1);
  dfdp_space_ =
      Thyra::multiVectorProductVectorSpace(model_->get_f_space(), num_param_);
  f_dfdp_space_ = Thyra::multiVectorProductVectorSpace(model_->get_f_space(),
                                                       num_param_ + 1);
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
  outArgs.set_Np_Ng(me_outArgs.Np(), me_outArgs.Ng() + g_offset_);
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
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
CombinedForwardSensitivityModelEvaluator<Scalar>::get_p_space(int p) const
{
  return model_->get_p_space(p);
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
CombinedForwardSensitivityModelEvaluator<Scalar>::get_p_names(int p) const
{
  return model_->get_p_names(p);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
CombinedForwardSensitivityModelEvaluator<Scalar>::get_x_space() const
{
  return x_dxdp_space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
CombinedForwardSensitivityModelEvaluator<Scalar>::get_f_space() const
{
  return f_dfdp_space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
CombinedForwardSensitivityModelEvaluator<Scalar>::get_g_space(int j) const
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
CombinedForwardSensitivityModelEvaluator<Scalar>::get_g_names(int j) const
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
CombinedForwardSensitivityModelEvaluator<Scalar>::create_W_op() const
{
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > op =
      sens_solve_model_->create_W_op();
  return Thyra::nonconstMultiVectorLinearOp(op, num_param_ + 1);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
CombinedForwardSensitivityModelEvaluator<Scalar>::create_DgDx_dot_op(
    int j) const
{
  return model_->create_DgDx_dot_op(j - g_offset_);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
CombinedForwardSensitivityModelEvaluator<Scalar>::create_DgDx_op(int j) const
{
  return model_->create_DgDx_op(j - g_offset_);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
CombinedForwardSensitivityModelEvaluator<Scalar>::create_DgDp_op(int j,
                                                                 int l) const
{
  return model_->create_DgDp_op(j - g_offset_, l);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
CombinedForwardSensitivityModelEvaluator<Scalar>::get_W_factory() const
{
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory =
      sens_solve_model_->get_W_factory();
  if (factory == Teuchos::null)
    return Teuchos::null;  // model_ doesn't support W_factory

  return Thyra::multiVectorLinearOpWithSolveFactory(factory, f_dfdp_space_,
                                                    x_dxdp_space_);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
CombinedForwardSensitivityModelEvaluator<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
CombinedForwardSensitivityModelEvaluator<Scalar>::getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::Range1D;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  MEB::InArgs<Scalar> me_nominal = model_->getNominalValues();
  MEB::InArgs<Scalar> nominal    = this->createInArgs();

  const Teuchos::Range1D rng(1, num_param_);
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  // Set initial x.  If dxdp_init == null, set initial dx/dp = 0
  RCP<const Thyra::VectorBase<Scalar> > me_x = me_nominal.get_x();
  if (me_x != Teuchos::null) {
    RCP<Thyra::VectorBase<Scalar> > x = Thyra::createMember(*x_dxdp_space_);
    RCP<DMVPV> x_dxdp                 = rcp_dynamic_cast<DMVPV>(x, true);
    Thyra::assign(x_dxdp->getNonconstMultiVector()->col(0).ptr(), *me_x);
    if (dxdp_init_ == Teuchos::null)
      Thyra::assign(x_dxdp->getNonconstMultiVector()->subView(rng).ptr(), zero);
    else
      Thyra::assign(x_dxdp->getNonconstMultiVector()->subView(rng).ptr(),
                    *dxdp_init_);
    nominal.set_x(x);
  }

  // Set initial xdot.  If dx_dotdp_init == null, set initial dxdot/dp = 0
  RCP<const Thyra::VectorBase<Scalar> > me_xdot;
  if (me_nominal.supports(MEB::IN_ARG_x_dot)) me_xdot = me_nominal.get_x_dot();
  if (me_xdot != Teuchos::null) {
    RCP<Thyra::VectorBase<Scalar> > xdot = Thyra::createMember(*x_dxdp_space_);
    RCP<DMVPV> xdot_dxdp                 = rcp_dynamic_cast<DMVPV>(xdot, true);
    Thyra::assign(xdot_dxdp->getNonconstMultiVector()->col(0).ptr(), *me_xdot);
    if (dx_dotdp_init_ == Teuchos::null)
      Thyra::assign(xdot_dxdp->getNonconstMultiVector()->subView(rng).ptr(),
                    zero);
    else
      Thyra::assign(xdot_dxdp->getNonconstMultiVector()->subView(rng).ptr(),
                    *dx_dotdp_init_);
    nominal.set_x_dot(xdot);
  }

  // Set initial xdotdot.  If dx_dotdotdp_init == null, set initial dxdotdot/dp
  // = 0
  RCP<const Thyra::VectorBase<Scalar> > me_xdotdot;
  if (me_nominal.supports(MEB::IN_ARG_x_dot_dot))
    me_xdotdot = me_nominal.get_x_dot_dot();
  if (me_xdotdot != Teuchos::null) {
    RCP<Thyra::VectorBase<Scalar> > xdotdot =
        Thyra::createMember(*x_dxdp_space_);
    RCP<DMVPV> xdotdot_dxdp = rcp_dynamic_cast<DMVPV>(xdotdot, true);
    Thyra::assign(xdotdot_dxdp->getNonconstMultiVector()->col(0).ptr(),
                  *me_xdotdot);
    if (dx_dotdotdp_init_ == Teuchos::null)
      Thyra::assign(xdotdot_dxdp->getNonconstMultiVector()->subView(rng).ptr(),
                    zero);
    else
      Thyra::assign(xdotdot_dxdp->getNonconstMultiVector()->subView(rng).ptr(),
                    *dx_dotdotdp_init_);
    nominal.set_x_dot_dot(xdotdot);
  }

  const int np = model_->Np();
  for (int i = 0; i < np; ++i) nominal.set_p(i, me_nominal.get_p(i));

  return nominal;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
CombinedForwardSensitivityModelEvaluator<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}

template <typename Scalar>
void CombinedForwardSensitivityModelEvaluator<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::Range1D;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  const bool use_tangent = use_dfdp_as_tangent_ || use_dgdp_as_tangent_;

  // setup input arguments for model
  RCP<const Thyra::VectorBase<Scalar> > x, xdot, xdotdot;
  RCP<const Thyra::MultiVectorBase<Scalar> > dxdp, dxdotdp, dxdotdotdp;
  RCP<const Thyra::VectorBase<Scalar> > dxdp_vec, dxdotdp_vec, dxdotdotdp_vec;
  MEB::InArgs<Scalar> me_inArgs = model_->getNominalValues();
  RCP<const DMVPV> x_dxdp       = rcp_dynamic_cast<const DMVPV>(inArgs.get_x(), true);
  x                             = x_dxdp->getMultiVector()->col(0);
  dxdp                          = x_dxdp->getMultiVector()->subView(Range1D(1, num_param_));
  me_inArgs.set_x(x);
  if (use_tangent) {
    dxdp_vec = Thyra::multiVectorProductVector(dxdp_space_, dxdp);
    if (use_dfdp_as_tangent_) me_inArgs.set_p(x_tangent_index_, dxdp_vec);
  }
  if (me_inArgs.supports(MEB::IN_ARG_x_dot)) {
    if (inArgs.get_x_dot() != Teuchos::null) {
      RCP<const DMVPV> xdot_dxdotdp =
          rcp_dynamic_cast<const DMVPV>(inArgs.get_x_dot(), true);
      xdot    = xdot_dxdotdp->getMultiVector()->col(0);
      dxdotdp = xdot_dxdotdp->getMultiVector()->subView(Range1D(1, num_param_));
      me_inArgs.set_x_dot(xdot);
      if (use_tangent) {
        dxdotdp_vec = Thyra::multiVectorProductVector(dxdp_space_, dxdotdp);
        if (use_dfdp_as_tangent_)
          me_inArgs.set_p(xdot_tangent_index_, dxdotdp_vec);
      }
    }
    else  // clear out xdot if it was set in nominalValues
      me_inArgs.set_x_dot(Teuchos::null);
  }
  if (me_inArgs.supports(MEB::IN_ARG_x_dot_dot)) {
    if (inArgs.get_x_dot_dot() != Teuchos::null) {
      RCP<const DMVPV> xdotdot_dxdotdotdp =
          rcp_dynamic_cast<const DMVPV>(inArgs.get_x_dot_dot(), true);
      xdotdot = xdotdot_dxdotdotdp->getMultiVector()->col(0);
      dxdotdotdp =
          xdotdot_dxdotdotdp->getMultiVector()->subView(Range1D(1, num_param_));
      me_inArgs.set_x_dot_dot(xdotdot);
      if (use_tangent) {
        dxdotdotdp_vec =
            Thyra::multiVectorProductVector(dxdp_space_, dxdotdotdp);
        if (use_dfdp_as_tangent_)
          me_inArgs.set_p(xdotdot_tangent_index_, dxdotdotdp_vec);
      }
    }
    else  // clear out xdotdot if it was set in nominalValues
      me_inArgs.set_x_dot_dot(Teuchos::null);
  }
  if (me_inArgs.supports(MEB::IN_ARG_t)) me_inArgs.set_t(inArgs.get_t());
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
  RCP<Thyra::VectorBase<Scalar> > f;
  RCP<Thyra::MultiVectorBase<Scalar> > dfdp;
  MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
  if (outArgs.get_f() != Teuchos::null) {
    RCP<DMVPV> f_dfdp = rcp_dynamic_cast<DMVPV>(outArgs.get_f(), true);
    f                 = f_dfdp->getNonconstMultiVector()->col(0);
    dfdp              = f_dfdp->getNonconstMultiVector()->subView(Range1D(1, num_param_));
    me_outArgs.set_f(f);
    me_outArgs.set_DfDp(p_index_, dfdp);
  }
  if (outArgs.supports(MEB::OUT_ARG_W_op) &&
      outArgs.get_W_op() != Teuchos::null &&
      model_.ptr() == sens_solve_model_.ptr()) {
    RCP<Thyra::LinearOpBase<Scalar> > op = outArgs.get_W_op();
    RCP<Thyra::MultiVectorLinearOp<Scalar> > mv_op =
        rcp_dynamic_cast<Thyra::MultiVectorLinearOp<Scalar> >(op, true);
    me_outArgs.set_W_op(mv_op->getNonconstLinearOp());
  }
  for (int j = g_offset_; j < outArgs.Ng(); ++j) {
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
  if (outArgs.supports(MEB::OUT_ARG_W_op) &&
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
  // requires computing df/dx, df/dxdot, df/dxdotdot as separate operators
  if (!use_dfdp_as_tangent_) {
    if (dxdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdx_ == Teuchos::null)
        my_dfdx_ = sens_residual_model_->create_W_op();
      MEB::OutArgs<Scalar> meo = sens_residual_model_->createOutArgs();
      meo.set_W_op(my_dfdx_);
      if (me_inArgs.supports(MEB::IN_ARG_alpha)) me_inArgs.set_alpha(0.0);
      if (me_inArgs.supports(MEB::IN_ARG_beta)) me_inArgs.set_beta(1.0);
      if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
        me_inArgs.set_W_x_dot_dot_coeff(0.0);
      sens_residual_model_->evalModel(me_inArgs, meo);
      my_dfdx_->apply(Thyra::NOTRANS, *dxdp, dfdp.ptr(), Scalar(1.0),
                      Scalar(1.0));
    }
    if (dxdotdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdxdot_ == Teuchos::null)
        my_dfdxdot_ = sens_residual_model_->create_W_op();
      MEB::OutArgs<Scalar> meo = sens_residual_model_->createOutArgs();
      meo.set_W_op(my_dfdxdot_);
      if (me_inArgs.supports(MEB::IN_ARG_alpha)) me_inArgs.set_alpha(1.0);
      if (me_inArgs.supports(MEB::IN_ARG_beta)) me_inArgs.set_beta(0.0);
      if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
        me_inArgs.set_W_x_dot_dot_coeff(0.0);
      sens_residual_model_->evalModel(me_inArgs, meo);
      my_dfdxdot_->apply(Thyra::NOTRANS, *dxdotdp, dfdp.ptr(), Scalar(1.0),
                         Scalar(1.0));
    }
    if (dxdotdotdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdxdotdot_ == Teuchos::null)
        my_dfdxdotdot_ = sens_residual_model_->create_W_op();
      MEB::OutArgs<Scalar> meo = sens_residual_model_->createOutArgs();
      meo.set_W_op(my_dfdxdotdot_);
      if (me_inArgs.supports(MEB::IN_ARG_alpha)) me_inArgs.set_alpha(0.0);
      if (me_inArgs.supports(MEB::IN_ARG_beta)) me_inArgs.set_beta(0.0);
      if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
        me_inArgs.set_W_x_dot_dot_coeff(1.0);
      sens_residual_model_->evalModel(me_inArgs, meo);
      my_dfdxdotdot_->apply(Thyra::NOTRANS, *dxdotdotdp, dfdp.ptr(),
                            Scalar(1.0), Scalar(1.0));
    }
  }

  // Reduced dg/dp = dg/dx*dx/dp + dg/dp
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
    // using dg/dp as the tangent
    if (!use_dgdp_as_tangent_ && use_dfdp_as_tangent_ &&
        dxdp_vec != Teuchos::null)
      me_inArgs.set_p(x_tangent_index_, Teuchos::null);
    if (!use_dgdp_as_tangent_ && use_dfdp_as_tangent_ &&
        dxdotdp_vec != Teuchos::null)
      me_inArgs.set_p(xdot_tangent_index_, Teuchos::null);
    if (!use_dgdp_as_tangent_ && use_dfdp_as_tangent_ &&
        dxdotdotdp_vec != Teuchos::null)
      me_inArgs.set_p(xdotdot_tangent_index_, Teuchos::null);

    // Set dx/dp, dxdot/dp, dxdodot/dp if necessary
    if (use_dgdp_as_tangent_ && !use_dfdp_as_tangent_ &&
        dxdp_vec != Teuchos::null)
      me_inArgs.set_p(x_tangent_index_, dxdp_vec);
    if (use_dgdp_as_tangent_ && !use_dfdp_as_tangent_ &&
        dxdotdp_vec != Teuchos::null)
      me_inArgs.set_p(xdot_tangent_index_, dxdotdp_vec);
    if (use_dgdp_as_tangent_ && !use_dfdp_as_tangent_ &&
        dxdotdotdp_vec != Teuchos::null)
      me_inArgs.set_p(xdotdot_tangent_index_, dxdotdotdp_vec);

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
          meo.supports(MEB::OUT_ARG_DgDx, g_index_);
      if (dgdx_support.supports(MEB::DERIV_LINEAR_OP))
        my_dgdx_->apply(Thyra::NOTRANS, *dxdp, dgdp.ptr(), Scalar(1.0),
                        Scalar(1.0));
      else if (dgdx_support.supports(MEB::DERIV_MV_GRADIENT_FORM))
        my_dgdx_mv_->apply(Thyra::TRANS, *dxdp, dgdp.ptr(), Scalar(1.0),
                           Scalar(1.0));
      else
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error,
            "Jacobian form of dg/dx not supported for reduced sensitivity");
    }
  }
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
CombinedForwardSensitivityModelEvaluator<Scalar>::getValidParameters()
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
