//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_ModelEvaluatorPairPartIMEX_StaggeredFSA_impl_hpp
#define Tempus_ModelEvaluatorPairPartIMEX_StaggeredFSA_impl_hpp

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

namespace Tempus {

template <typename Scalar>
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::
    WrapperModelEvaluatorPairPartIMEX_StaggeredFSA(
        const Teuchos::RCP<const WrapperModelEvaluatorPairPartIMEX_Basic<
            Scalar> >& forwardModel,
        const bool is_pseudotransient,
        const Teuchos::RCP<const Teuchos::ParameterList>& pList)
  : forwardModel_(forwardModel),
    use_dfdp_as_tangent_(false),
    y_tangent_index_(3)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Thyra::multiVectorProductVectorSpace;

  RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  if (pList != Teuchos::null) *pl = *pList;
  pl->validateParametersAndSetDefaults(*this->getValidParameters());
  use_dfdp_as_tangent_ = pl->get<bool>("Use DfDp as Tangent");
  y_tangent_index_     = pl->get<int>("Sensitivity Y Tangent Index");
  pl->remove("Sensitivity Y Tangent Index");

  appExplicitModel_ = forwardModel_->getExplicitModel();
  appImplicitModel_ = forwardModel_->getImplicitModel();
  fsaExplicitModel_ = rcp(new FSAME(appExplicitModel_, appExplicitModel_,
                                    appExplicitModel_, is_pseudotransient, pl));
  fsaImplicitModel_ = rcp(new FSAME(appImplicitModel_, appImplicitModel_,
                                    appImplicitModel_, is_pseudotransient, pl));

  const int y_param_index    = forwardModel_->getParameterIndex();
  const int sens_param_index = pl->get<int>("Sensitivity Parameter Index");
  const int num_sens_param =
      appImplicitModel_->get_p_space(sens_param_index)->dim();
  RCP<const Thyra::VectorSpaceBase<Scalar> > explicit_y_space =
      appImplicitModel_->get_p_space(y_param_index);
  RCP<const Thyra::VectorSpaceBase<Scalar> > implicit_x_space =
      appImplicitModel_->get_x_space();
  explicit_dydp_prod_space_ =
      multiVectorProductVectorSpace(explicit_y_space, num_sens_param);
  imex_dxdp_prod_space_ =
      multiVectorProductVectorSpace(implicit_x_space, num_sens_param);

  Base::setup(fsaExplicitModel_, fsaImplicitModel_,
              forwardModel_->getNumExplicitOnlyBlocks(), y_param_index);
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::initialize()
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  this->useImplicitModel_       = true;
  this->wrapperImplicitInArgs_  = this->createInArgs();
  this->wrapperImplicitOutArgs_ = this->createOutArgs();
  this->useImplicitModel_       = false;

  RCP<const Thyra::VectorBase<Scalar> > z =
      this->explicitModel_->getNominalValues().get_x();

  // A Thyra::VectorSpace requirement
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(getIMEXVector(z)->space()->isCompatible(
          *(this->implicitModel_->get_x_space()))),
      std::logic_error,
      "Error - WrapperModelEvaluatorPairIMEX_StaggeredFSA::initialize()\n"
          << "  Explicit and Implicit vector x spaces are incompatible!\n"
          << "  Explicit vector x space = "
          << *(getIMEXVector(z)->space())
          << "\n  Implicit vector x space = "
          << *(this->implicitModel_->get_x_space()) << "\n");

  // Valid number of blocks?
  const RCP<const DMVPV> z_dmvpv = rcp_dynamic_cast<const DMVPV>(z, true);
  const RCP<const Thyra::MultiVectorBase<Scalar> > z_mv =
      z_dmvpv->getMultiVector();
  RCP<const PMVB> zPVector = rcp_dynamic_cast<const PMVB>(z_mv);
  TEUCHOS_TEST_FOR_EXCEPTION(
      zPVector == Teuchos::null, std::logic_error,
      "Error - WrapperModelEvaluatorPairPartIMEX_StaggeredFSA::initialize()\n"
      "  was given a VectorBase that could not be cast to a\n"
      "  ProductVectorBase!\n");

  int numBlocks = zPVector->productSpace()->numBlocks();

  TEUCHOS_TEST_FOR_EXCEPTION(
      !(0 <= this->numExplicitOnlyBlocks_ &&
        this->numExplicitOnlyBlocks_ < numBlocks),
      std::logic_error,
      "Error - WrapperModelEvaluatorPairPartIMEX_StaggeredFSA::initialize()\n"
          << "Invalid number of explicit-only blocks = "
          << this->numExplicitOnlyBlocks_
          << "\nShould be in the interval [0, numBlocks) = [0, "
          << numBlocks << ")\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::get_p_space(int i) const
{
  if (this->useImplicitModel_) {
    if (i == this->parameterIndex_)
      return explicit_dydp_prod_space_;
    else
      return appImplicitModel_->get_p_space(i);
  }

  return appExplicitModel_->get_p_space(i);
}

template <typename Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::getIMEXVector(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::MultiVectorBase;
  using Thyra::multiVectorProductVector;
  using Thyra::VectorBase;

  // StaggeredFSA ME stores vectors as DMVPV's.  To extract the implicit
  // part of the vector, cast it to DMVPV,  extract the multi-vector,
  // cast it to a product multi-vector, extract the IMEX block, then
  // create a DMVPV from it.

  if (full == Teuchos::null) return Teuchos::null;

  if (this->numExplicitOnlyBlocks_ == 0) return full;

  const RCP<DMVPV> full_dmvpv = rcp_dynamic_cast<DMVPV>(full, true);
  const RCP<MultiVectorBase<Scalar> > full_mv =
      full_dmvpv->getNonconstMultiVector();
  const RCP<PMVB> blk_full_mv = rcp_dynamic_cast<PMVB>(full_mv, true);

  // special case where the implicit terms are not blocked
  const int numBlocks         = blk_full_mv->productSpace()->numBlocks();
  const int numExplicitBlocks = this->numExplicitOnlyBlocks_;
  if (numBlocks == numExplicitBlocks + 1) {
    const RCP<MultiVectorBase<Scalar> > imex_mv =
        blk_full_mv->getNonconstMultiVectorBlock(numExplicitBlocks);
    return multiVectorProductVector(imex_dxdp_prod_space_, imex_mv);
  }

  // Not supposed to get here, apparently
  TEUCHOS_ASSERT(false);
  return Teuchos::null;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::getIMEXVector(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::MultiVectorBase;
  using Thyra::multiVectorProductVector;
  using Thyra::VectorBase;

  // StaggeredFSA ME stores vectors as DMVPV's.  To extract the implicit
  // part of the vector, cast it to DMVPV,  extract the multi-vector,
  // cast it to a product multi-vector, extract the IMEX block, then
  // create a DMVPV from it.

  if (full == Teuchos::null) return Teuchos::null;

  if (this->numExplicitOnlyBlocks_ == 0) return full;

  const RCP<const DMVPV> full_dmvpv = rcp_dynamic_cast<const DMVPV>(full, true);
  const RCP<const MultiVectorBase<Scalar> > full_mv =
      full_dmvpv->getMultiVector();
  const RCP<const PMVB> blk_full_mv =
      rcp_dynamic_cast<const PMVB>(full_mv, true);

  // special case where the implicit terms are not blocked
  const int numBlocks         = blk_full_mv->productSpace()->numBlocks();
  const int numExplicitBlocks = this->numExplicitOnlyBlocks_;
  if (numBlocks == numExplicitBlocks + 1) {
    const RCP<const MultiVectorBase<Scalar> > imex_mv =
        blk_full_mv->getMultiVectorBlock(numExplicitBlocks);
    return multiVectorProductVector(imex_dxdp_prod_space_, imex_mv);
  }

  // Not supposed to get here, apparently
  TEUCHOS_ASSERT(false);
  return Teuchos::null;
}

template <typename Scalar>
Teuchos::RCP<Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::getExplicitOnlyVector(
    const Teuchos::RCP<Thyra::VectorBase<Scalar> >& full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::MultiVectorBase;
  using Thyra::multiVectorProductVector;
  using Thyra::multiVectorProductVectorSpace;
  using Thyra::VectorBase;

  // StaggeredFSA ME stores vectors as DMVPV's.  To extract the explicit
  // part of the vector, cast it to DMVPV,  extract the multi-vector,
  // cast it to a product multi-vector, extract the explicit block, then
  // create a DMVPV from it.

  if (full == Teuchos::null) return Teuchos::null;

  if (this->numExplicitOnlyBlocks_ == 0) return full;

  const RCP<DMVPV> full_dmvpv = rcp_dynamic_cast<DMVPV>(full, true);
  const RCP<MultiVectorBase<Scalar> > full_mv =
      full_dmvpv->getNonconstMultiVector();
  const RCP<PMVB> blk_full_mv = rcp_dynamic_cast<PMVB>(full_mv, true);

  // special case where the explicit terms are not blocked
  const int numExplicitBlocks = this->numExplicitOnlyBlocks_;
  if (numExplicitBlocks == 1) {
    const RCP<MultiVectorBase<Scalar> > explicit_mv =
        blk_full_mv->getNonconstMultiVectorBlock(0);
    return multiVectorProductVector(explicit_dydp_prod_space_, explicit_mv);
  }

  // Not supposed to get here, apparently
  TEUCHOS_ASSERT(false);
  return Teuchos::null;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::getExplicitOnlyVector(
    const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& full) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::MultiVectorBase;
  using Thyra::multiVectorProductVector;
  using Thyra::multiVectorProductVectorSpace;
  using Thyra::VectorBase;

  // StaggeredFSA ME stores vectors as DMVPV's.  To extract the explicit
  // part of the vector, cast it to DMVPV,  extract the multi-vector,
  // cast it to a product multi-vector, extract the explicit block, then
  // create a DMVPV from it.

  if (full == Teuchos::null) return Teuchos::null;

  if (this->numExplicitOnlyBlocks_ == 0) return full;

  const RCP<const DMVPV> full_dmvpv = rcp_dynamic_cast<const DMVPV>(full, true);
  const RCP<const MultiVectorBase<Scalar> > full_mv =
      full_dmvpv->getMultiVector();
  const RCP<const PMVB> blk_full_mv =
      rcp_dynamic_cast<const PMVB>(full_mv, true);

  // special case where the explicit terms are not blocked
  const int numExplicitBlocks = this->numExplicitOnlyBlocks_;
  if (numExplicitBlocks == 1) {
    const RCP<const MultiVectorBase<Scalar> > explicit_mv =
        blk_full_mv->getMultiVectorBlock(0);
    return multiVectorProductVector(explicit_dydp_prod_space_, explicit_mv);
  }

  // Not supposed to get here, apparently
  TEUCHOS_ASSERT(false);
  return Teuchos::null;
}

template <typename Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::getForwardModel() const
{
  return forwardModel_;
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::
    setForwardSolutionHistory(
        const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& sh)
{
  sh_       = sh;
  t_interp_ = Teuchos::ScalarTraits<Scalar>::rmax();
  fsaExplicitModel_->setForwardSolutionHistory(sh);
  // We handle interpolation for fsaImplicitModel in evalModelImpl()
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::
    setForwardSolutionState(
        const Teuchos::RCP<const Tempus::SolutionState<Scalar> >& s)
{
  sh_            = Teuchos::null;
  forward_state_ = s;
  buildIMEXStates();
  fsaExplicitModel_->setForwardSolutionState(s);
  fsaImplicitModel_->setForwardSolutionState(implicit_x_state_);
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::setSolver(
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    const bool force_W_update)
{
  // Make sure we have the useImplicitModel_ flag set since
  // force_W_update = true will force a model evaluation
  bool tf = forwardModel_->getUseImplicitModel();
  Teuchos::RCP<WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >
      nc_forwardModel = Teuchos::rcp_const_cast<
          WrapperModelEvaluatorPairPartIMEX_Basic<Scalar> >(forwardModel_);
  nc_forwardModel->setUseImplicitModel(true);
  fsaImplicitModel_->setSolver(solver, force_W_update);
  nc_forwardModel->setUseImplicitModel(tf);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::createInArgs() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::createMember;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = Base::createInArgs();

  // Set p to be the correct product vector form for the explicit only vector y
  if (this->useImplicitModel_ == true) {
    if (inArgs.get_p(this->parameterIndex_) != Teuchos::null) {
      RCP<DMVPV> dydp = rcp_dynamic_cast<DMVPV>(
          createMember(*explicit_dydp_prod_space_), true);
      Thyra::assign(dydp->getNonconstMultiVector().ptr(), Scalar(0.0));
      inArgs.set_p(this->parameterIndex_, dydp);
    }
  }
  return inArgs;
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::Range1D;
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
      buildIMEXStates();
      fsaImplicitModel_->setForwardSolutionState(implicit_x_state_);
    }
  }
  else {
    TEUCHOS_ASSERT(forward_state_ != Teuchos::null);
    forward_t = forward_state_->getTime();
  }

  const int p_index = this->parameterIndex_;

  //
  // From Base::evalModelImpl()
  //
  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  RCP<Thyra::VectorBase<Scalar> > x_dot =
      Thyra::createMember(fsaImplicitModel_->get_x_space());
  this->timeDer_->compute(x, x_dot);

  MEB::InArgs<Scalar> fsaImplicitInArgs(this->wrapperImplicitInArgs_);
  MEB::OutArgs<Scalar> fsaImplicitOutArgs(this->wrapperImplicitOutArgs_);
  if (fsaImplicitInArgs.supports(MEB::IN_ARG_t))
    fsaImplicitInArgs.set_t(inArgs.get_t());
  fsaImplicitInArgs.set_x(x);
  fsaImplicitInArgs.set_x_dot(x_dot);
  for (int i = 0; i < fsaImplicitModel_->Np(); ++i) {
    // Copy over parameters except for the parameter for explicit-only vector!
    if ((inArgs.get_p(i) != Teuchos::null) && (i != p_index))
      fsaImplicitInArgs.set_p(i, inArgs.get_p(i));
  }

  // p-vector for index parameterIndex_ is part of the IMEX solution vector,
  // and therefore is an n column multi-vector where n is the number of
  // sensitivity parameters.  Pull out the sensitivity components before
  // passing along to the ME, then use them for adding in dg/dy*dy/dp term.
  TEUCHOS_ASSERT(explicit_y_state_ != Teuchos::null);
  RCP<const Thyra::VectorBase<Scalar> > y = explicit_y_state_->getX();
  RCP<const Thyra::MultiVectorBase<Scalar> > dydp;
  if (fsaImplicitInArgs.get_p(p_index) != Teuchos::null) {
    RCP<const Thyra::VectorBase<Scalar> > p = fsaImplicitInArgs.get_p(p_index);
    dydp                                    = rcp_dynamic_cast<const DMVPV>(p, true)->getMultiVector();
    fsaImplicitInArgs.set_p(p_index, y);
  }
  if (use_dfdp_as_tangent_) {
    RCP<const Thyra::VectorBase<Scalar> > dydp_vec =
        Thyra::multiVectorProductVector(explicit_dydp_prod_space_, dydp);
    fsaImplicitInArgs.set_p(y_tangent_index_, dydp_vec);
  }

  fsaImplicitOutArgs.set_f(outArgs.get_f());
  fsaImplicitOutArgs.set_W_op(outArgs.get_W_op());

  fsaImplicitModel_->evalModel(fsaImplicitInArgs, fsaImplicitOutArgs);

  // Compute derivative of implicit residual with respect to explicit only
  // vector y, which is passed as a parameter
  if (!use_dfdp_as_tangent_ && outArgs.get_f() != Teuchos::null) {
    MEB::InArgs<Scalar> appImplicitInArgs =
        appImplicitModel_->getNominalValues();
    TEUCHOS_ASSERT(forward_state_ != Teuchos::null);
    RCP<const Thyra::VectorBase<Scalar> > app_x = implicit_x_state_->getX();
    RCP<const Thyra::VectorBase<Scalar> > app_x_dot =
        implicit_x_state_->getXDot();
    appImplicitInArgs.set_x(app_x);
    appImplicitInArgs.set_x_dot(app_x_dot);
    for (int i = 0; i < appImplicitModel_->Np(); ++i) {
      if (i != p_index) appImplicitInArgs.set_p(i, inArgs.get_p(i));
    }
    appImplicitInArgs.set_p(p_index, y);
    if (appImplicitInArgs.supports(MEB::IN_ARG_t))
      appImplicitInArgs.set_t(forward_t);
    MEB::OutArgs<Scalar> appImplicitOutArgs =
        appImplicitModel_->createOutArgs();
    MEB::DerivativeSupport dfdp_support =
        appImplicitOutArgs.supports(MEB::OUT_ARG_DfDp, p_index);
    Thyra::EOpTransp trans = Thyra::NOTRANS;
    if (dfdp_support.supports(MEB::DERIV_LINEAR_OP)) {
      if (my_dfdp_op_ == Teuchos::null)
        my_dfdp_op_ = appImplicitModel_->create_DfDp_op(p_index);
      appImplicitOutArgs.set_DfDp(p_index,
                                  MEB::Derivative<Scalar>(my_dfdp_op_));
      trans = Thyra::NOTRANS;
    }
    else if (dfdp_support.supports(MEB::DERIV_MV_JACOBIAN_FORM)) {
      if (my_dfdp_mv_ == Teuchos::null)
        my_dfdp_mv_ = Thyra::createMembers(
            appImplicitModel_->get_f_space(),
            appImplicitModel_->get_p_space(p_index)->dim());
      appImplicitOutArgs.set_DfDp(
          p_index,
          MEB::Derivative<Scalar>(my_dfdp_mv_, MEB::DERIV_MV_JACOBIAN_FORM));
      my_dfdp_op_ = my_dfdp_mv_;
      trans       = Thyra::NOTRANS;
    }
    else if (dfdp_support.supports(MEB::DERIV_MV_GRADIENT_FORM)) {
      if (my_dfdp_mv_ == Teuchos::null)
        my_dfdp_mv_ =
            Thyra::createMembers(appImplicitModel_->get_p_space(p_index),
                                 appImplicitModel_->get_f_space()->dim());
      appImplicitOutArgs.set_DfDp(
          p_index,
          MEB::Derivative<Scalar>(my_dfdp_mv_, MEB::DERIV_MV_GRADIENT_FORM));
      my_dfdp_op_ = my_dfdp_mv_;
      trans       = Thyra::TRANS;
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                 "Invalid df/dp support");

    appImplicitModel_->evalModel(appImplicitInArgs, appImplicitOutArgs);

    // Add df/dy*dy/dp term to residual
    RCP<Thyra::MultiVectorBase<Scalar> > dfdp =
        rcp_dynamic_cast<DMVPV>(outArgs.get_f(), true)
            ->getNonconstMultiVector();
    my_dfdp_op_->apply(trans, *dydp, dfdp.ptr(), Scalar(1.0), Scalar(1.0));
  }
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::getValidParameters()
    const
{
  Teuchos::RCP<const Teuchos::ParameterList> fsa_pl =
      StaggeredForwardSensitivityModelEvaluator<Scalar>::getValidParameters();
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList(*fsa_pl);
  pl->set<int>("Sensitivity Y Tangent Index", 3);
  return pl;
}

template <typename Scalar>
void WrapperModelEvaluatorPairPartIMEX_StaggeredFSA<Scalar>::buildIMEXStates()
    const
{
  explicit_y_state_ = Teuchos::rcp(new SolutionState<Scalar>(
      forward_state_->getMetaData(),
      forwardModel_->getExplicitOnlyVector(forward_state_->getX()),
      forwardModel_->getExplicitOnlyVector(forward_state_->getXDot()),
      forwardModel_->getExplicitOnlyVector(forward_state_->getXDotDot()),
      forward_state_->getStepperState(), Teuchos::null));
  implicit_x_state_ = Teuchos::rcp(new SolutionState<Scalar>(
      forward_state_->getMetaData(),
      forwardModel_->getIMEXVector(forward_state_->getX()),
      forwardModel_->getIMEXVector(forward_state_->getXDot()),
      forwardModel_->getIMEXVector(forward_state_->getXDotDot()),
      forward_state_->getStepperState(), Teuchos::null));
}

}  // namespace Tempus

#endif  // Tempus_ModelEvaluatorPairPartIMEX_StaggeredFSA_impl_hpp
