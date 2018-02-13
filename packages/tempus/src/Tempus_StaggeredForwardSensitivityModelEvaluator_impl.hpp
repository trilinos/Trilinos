// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

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
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & model,
  const Teuchos::RCP<const Teuchos::ParameterList>& pList,
  const Teuchos::RCP<MultiVector>& dxdp_init,
  const Teuchos::RCP<MultiVector>& dx_dotdp_init,
  const Teuchos::RCP<MultiVector>& dx_dotdotdp_init) :
  model_(model),
  dxdp_init_(dxdp_init),
  dx_dotdp_init_(dx_dotdp_init),
  dx_dotdotdp_init_(dx_dotdotdp_init),
  p_index_(0),
  x_tangent_index_(1),
  xdot_tangent_index_(2),
  xdotdot_tangent_index_(3),
  use_dfdp_as_tangent_(false),
  t_interp_(Teuchos::ScalarTraits<Scalar>::rmax())
{
  typedef Thyra::ModelEvaluatorBase MEB;

  // Set parameters
  Teuchos::RCP<Teuchos::ParameterList> pl =
    Teuchos::rcp(new Teuchos::ParameterList);
  if (pList != Teuchos::null)
    *pl = *pList;
  pl->validateParametersAndSetDefaults(*this->getValidParameters());
  use_dfdp_as_tangent_ = pl->get<bool>("Use DfDp as Tangent");
  p_index_ = pl->get<int>("Sensitivity Parameter Index");
  x_tangent_index_ = pl->get<int>("Sensitivity X Tangent Index");
  xdot_tangent_index_ = pl->get<int>("Sensitivity X-Dot Tangent Index");
  xdotdot_tangent_index_ = pl->get<int>("Sensitivity X-Dot-Dot Tangent Index");

  num_param_ = model_->get_p_space(p_index_)->dim();
  dxdp_space_ =
    Thyra::multiVectorProductVectorSpace(model_->get_x_space(), num_param_);
  dfdp_space_ =
    Thyra::multiVectorProductVectorSpace(model_->get_f_space(), num_param_);

  MEB::InArgs<Scalar> me_inArgs = model_->createInArgs();
  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  if (me_inArgs.supports(MEB::IN_ARG_x_dot))
    inArgs.setSupports(MEB::IN_ARG_x_dot);
  if (me_inArgs.supports(MEB::IN_ARG_t))
    inArgs.setSupports(MEB::IN_ARG_t);
  if (me_inArgs.supports(MEB::IN_ARG_alpha))
    inArgs.setSupports(MEB::IN_ARG_alpha);
  if (me_inArgs.supports(MEB::IN_ARG_beta))
    inArgs.setSupports(MEB::IN_ARG_beta);
  if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
    inArgs.setSupports(MEB::IN_ARG_W_x_dot_dot_coeff);

  // Support additional parameters for x and xdot
  inArgs.set_Np(me_inArgs.Np());
  prototypeInArgs_ = inArgs;

  MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(me_inArgs.Np(),me_outArgs.Ng());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  for (int j=0; j<me_outArgs.Ng(); ++j) {
    outArgs.setSupports(MEB::OUT_ARG_DgDx_dot, j,
                        me_outArgs.supports(MEB::OUT_ARG_DgDx_dot, j));
    outArgs.setSupports(MEB::OUT_ARG_DgDx, j,
                        me_outArgs.supports(MEB::OUT_ARG_DgDx, j));
    for (int l=0; l<me_outArgs.Np(); ++l) {
      outArgs.setSupports(MEB::OUT_ARG_DgDp, j, l,
                          me_outArgs.supports(MEB::OUT_ARG_DgDp, j, l));
    }
  }
  prototypeOutArgs_ = outArgs;

  TEUCHOS_ASSERT(me_outArgs.supports(MEB::OUT_ARG_DfDp, p_index_).supports(MEB::DERIV_MV_JACOBIAN_FORM));
  if (!use_dfdp_as_tangent_)
    TEUCHOS_ASSERT(me_outArgs.supports(MEB::OUT_ARG_W_op));
}

template <typename Scalar>
void
StaggeredForwardSensitivityModelEvaluator<Scalar>::
setForwardSolutionHistory(
  const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& sh)
{
  sh_ = sh;
  t_interp_ = Teuchos::ScalarTraits<Scalar>::rmax();
}

template <typename Scalar>
void
StaggeredForwardSensitivityModelEvaluator<Scalar>::
setForwardSolutionState(
  const Teuchos::RCP<const Tempus::SolutionState<Scalar> >& s)
{
  sh_ = Teuchos::null;
  forward_state_ = s;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
get_p_space(int p) const
{
  TEUCHOS_ASSERT(p < model_->Np());
  return model_->get_p_space(p);
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
get_p_names(int p) const
{
  TEUCHOS_ASSERT(p < model_->Np());
  return model_->get_p_names(p);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
get_x_space() const
{
  return dxdp_space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
get_f_space() const
{
  return dfdp_space_;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
get_g_space(int j) const
{
  return model_->get_g_space(j);
}

template <typename Scalar>
Teuchos::ArrayView<const std::string>
StaggeredForwardSensitivityModelEvaluator<Scalar>::
get_g_names(int j) const
{
  return model_->get_g_names(j);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
create_W_op() const
{
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > op;
  if (lo_ != Teuchos::null)
    op = lo_;
  else
    op = model_->create_W_op();
  return Thyra::nonconstMultiVectorLinearOp(op, num_param_);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
create_DgDx_dot_op(int j) const
{
  return model_->create_DgDx_dot_op(j);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
create_DgDx_op(int j) const
{
  return model_->create_DgDx_op(j);
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
create_DgDp_op(int j, int l) const
{
  return model_->create_DgDp_op(j,l);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
StaggeredForwardSensitivityModelEvaluator<Scalar>::
get_W_factory() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > factory;
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > model_factory
    = model_->get_W_factory();
  if (model_factory == Teuchos::null)
    return Teuchos::null; // model_ doesn't support W_factory
  if (po_ != Teuchos::null) {
    factory =
      Thyra::reuseLinearOpWithSolveFactory<Scalar>(model_factory, po_);
  }
  else
    factory = model_factory;
  return Thyra::multiVectorLinearOpWithSolveFactory(factory,
                                                    dfdp_space_,
                                                    dxdp_space_);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StaggeredForwardSensitivityModelEvaluator<Scalar>::
createInArgs() const
{
  return prototypeInArgs_;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
StaggeredForwardSensitivityModelEvaluator<Scalar>::
getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  MEB::InArgs<Scalar> me_nominal = model_->getNominalValues();
  MEB::InArgs<Scalar> nominal = this->createInArgs();

  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  // Set initial x.  If dxdp_init == null, set initial dx/dp = 0
  RCP< Thyra::VectorBase<Scalar> > x = Thyra::createMember(*dxdp_space_);
  RCP<DMVPV> dxdp = rcp_dynamic_cast<DMVPV>(x,true);
  if (dxdp_init_ == Teuchos::null)
    Thyra::assign(dxdp->getNonconstMultiVector().ptr(), zero);
  else
    Thyra::assign(dxdp->getNonconstMultiVector().ptr(), *dxdp_init_);
  nominal.set_x(x);

  // Set initial xdot.  If dx_dotdp_init == null, set initial dxdot/dp = 0
  if (me_nominal.supports(MEB::IN_ARG_x_dot)) {
    RCP< Thyra::VectorBase<Scalar> > xdot = Thyra::createMember(*dxdp_space_);
    RCP<DMVPV> dxdotdp = rcp_dynamic_cast<DMVPV>(xdot,true);
    if (dx_dotdp_init_ == Teuchos::null)
      Thyra::assign(dxdotdp->getNonconstMultiVector().ptr(), zero);
    else
      Thyra::assign(dxdotdp->getNonconstMultiVector().ptr(), *dx_dotdp_init_);
    nominal.set_x_dot(xdot);
  }

  // Set initial xdotdot.  If dx_dotdotdp_init == null, set initial dxdotdot/dp = 0
  if (me_nominal.supports(MEB::IN_ARG_x_dot_dot)) {
    RCP< Thyra::VectorBase<Scalar> > xdotdot =
      Thyra::createMember(*dxdp_space_);
    RCP<DMVPV> dxdotdotdp = rcp_dynamic_cast<DMVPV>(xdotdot,true);
    if (dx_dotdotdp_init_ == Teuchos::null)
      Thyra::assign(dxdotdotdp->getNonconstMultiVector().ptr(), zero);
    else
      Thyra::assign(dxdotdotdp->getNonconstMultiVector().ptr(),
                    *dx_dotdotdp_init_);
    nominal.set_x_dot_dot(xdotdot);
  }

  const int np = model_->Np();
  for (int i=0; i<np; ++i)
    nominal.set_p(i, me_nominal.get_p(i));
  return nominal;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
StaggeredForwardSensitivityModelEvaluator<Scalar>::
createOutArgsImpl() const
{
  return prototypeOutArgs_;
}

template <typename Scalar>
void
StaggeredForwardSensitivityModelEvaluator<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
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
      t_interp_ = forward_t;
    }
  }
  else {
    TEUCHOS_ASSERT(forward_state_ != Teuchos::null);
    forward_t = forward_state_->getTime();
  }

  // setup input arguments for model
  RCP< const Thyra::MultiVectorBase<Scalar> > dxdp, dxdotdp, dxdotdotdp;
  MEB::InArgs<Scalar> me_inArgs = model_->getNominalValues();
  dxdp = rcp_dynamic_cast<const DMVPV>(inArgs.get_x(),true)->getMultiVector();
  me_inArgs.set_x(forward_state_->getX());
  if (use_dfdp_as_tangent_)
    me_inArgs.set_p(x_tangent_index_, inArgs.get_x());
  if (me_inArgs.supports(MEB::IN_ARG_x_dot)) {
    if (inArgs.get_x_dot() != Teuchos::null) {
      dxdotdp =
        rcp_dynamic_cast<const DMVPV>(inArgs.get_x_dot(),true)->getMultiVector();
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
      dxdotdotdp =
        rcp_dynamic_cast<const DMVPV>(inArgs.get_x_dot_dot(),true)->getMultiVector();
      me_inArgs.set_x_dot_dot(forward_state_->getXDotDot());
      if (use_dfdp_as_tangent_)
        me_inArgs.set_p(xdotdot_tangent_index_, inArgs.get_x_dot_dot());
    }
    else // clear out xdotdot if it was set in nominalValues
      me_inArgs.set_x_dot_dot(Teuchos::null);
  }
  if (me_inArgs.supports(MEB::IN_ARG_t))
    me_inArgs.set_t(forward_t);
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
  for (int i=0; i<np; ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      if (!use_dfdp_as_tangent_ ||
          (use_dfdp_as_tangent_ && i != x_tangent_index_ &&
                                   i != xdot_tangent_index_ &&
                                   i != xdotdot_tangent_index_ ))
        me_inArgs.set_p(i, inArgs.get_p(i));
  }

  // setup output arguments for model
  RCP< Thyra::MultiVectorBase<Scalar> > dfdp;
  MEB::OutArgs<Scalar> me_outArgs = model_->createOutArgs();
  if (outArgs.get_f() != Teuchos::null) {
    dfdp = rcp_dynamic_cast<DMVPV>(outArgs.get_f(),true)->getNonconstMultiVector();
    me_outArgs.set_DfDp(p_index_, dfdp);
  }
  if (lo_ == Teuchos::null && outArgs.get_W_op() != Teuchos::null) {
    RCP<Thyra::LinearOpBase<Scalar> > op = outArgs.get_W_op();
    RCP<Thyra::MultiVectorLinearOp<Scalar> > mv_op =
      rcp_dynamic_cast<Thyra::MultiVectorLinearOp<Scalar> >(op,true);
    me_outArgs.set_W_op(mv_op->getNonconstLinearOp());
  }
  for (int j=0; j<outArgs.Ng(); ++j) {
    if (!me_outArgs.supports(MEB::OUT_ARG_DgDx_dot,j).none())
      me_outArgs.set_DgDx_dot(j, outArgs.get_DgDx_dot(j));
    if (!me_outArgs.supports(MEB::OUT_ARG_DgDx,j).none())
      me_outArgs.set_DgDx(j, outArgs.get_DgDx(j));
    for (int l=0; l<outArgs.Np(); ++l)
      if (!me_outArgs.supports(MEB::OUT_ARG_DgDp,j,l).none())
        me_outArgs.set_DgDp(j, l, outArgs.get_DgDp(j,l));
  }

  // build residual and jacobian
  model_->evalModel(me_inArgs, me_outArgs);

  // Compute (df/dx) * (dx/dp) + (df/dxdot) * (dxdot/dp) + (df/dxdotdot) * (dxdotdot/dp) + (df/dp)
  // if the underlying ME doesn't already do this.  This requires computing
  // df/dx, df/dxdot, df/dxdotdot as separate operators.
  // For pseudo-transient, we would like to reuse these operators, but this is
  // complicated when steppers use both implicit and explicit forms.
  if (!use_dfdp_as_tangent_) {
    if (dxdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdx_ == Teuchos::null)
        my_dfdx_ = model_->create_W_op();
      MEB::OutArgs<Scalar> meo = model_->createOutArgs();
      meo.set_W_op(my_dfdx_);
      if (me_inArgs.supports(MEB::IN_ARG_alpha))
        me_inArgs.set_alpha(0.0);
      if (me_inArgs.supports(MEB::IN_ARG_beta))
        me_inArgs.set_beta(1.0);
      if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
        me_inArgs.set_W_x_dot_dot_coeff(0.0);
      model_->evalModel(me_inArgs, meo);
      my_dfdx_->apply(Thyra::NOTRANS, *dxdp, dfdp.ptr(), Scalar(1.0), Scalar(1.0));
    }
    if (dxdotdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdxdot_ == Teuchos::null)
        my_dfdxdot_ = model_->create_W_op();
      MEB::OutArgs<Scalar> meo = model_->createOutArgs();
      meo.set_W_op(my_dfdxdot_);
      if (me_inArgs.supports(MEB::IN_ARG_alpha))
        me_inArgs.set_alpha(1.0);
      if (me_inArgs.supports(MEB::IN_ARG_beta))
        me_inArgs.set_beta(0.0);
      if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
        me_inArgs.set_W_x_dot_dot_coeff(0.0);
      model_->evalModel(me_inArgs, meo);
      my_dfdxdot_->apply(Thyra::NOTRANS, *dxdotdp, dfdp.ptr(), Scalar(1.0), Scalar(1.0));
    }
    if (dxdotdotdp != Teuchos::null && dfdp != Teuchos::null) {
      if (my_dfdxdotdot_ == Teuchos::null)
        my_dfdxdotdot_ = model_->create_W_op();
      MEB::OutArgs<Scalar> meo = model_->createOutArgs();
      meo.set_W_op(my_dfdxdotdot_);
      if (me_inArgs.supports(MEB::IN_ARG_alpha))
        me_inArgs.set_alpha(0.0);
      if (me_inArgs.supports(MEB::IN_ARG_beta))
        me_inArgs.set_beta(0.0);
      if (me_inArgs.supports(MEB::IN_ARG_W_x_dot_dot_coeff))
        me_inArgs.set_W_x_dot_dot_coeff(1.0);
      model_->evalModel(me_inArgs, meo);
      my_dfdxdotdot_->apply(Thyra::NOTRANS, *dxdotdotdp, dfdp.ptr(), Scalar(1.0), Scalar(1.0));
    }
  }
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StaggeredForwardSensitivityModelEvaluator<Scalar>::
getValidParameters()
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->set<bool>("Use DfDp as Tangent", false);
  pl->set<int>("Sensitivity Parameter Index", 0);
  pl->set<int>("Sensitivity X Tangent Index", 1);
  pl->set<int>("Sensitivity X-Dot Tangent Index", 2);
  pl->set<int>("Sensitivity X-Dot-Dot Tangent Index", 3);
  return pl;
}

} // namespace Tempus

#endif
