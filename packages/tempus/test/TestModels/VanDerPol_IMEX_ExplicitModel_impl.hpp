//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_VANDERPOL_IMEX_EXPLICITMODEL_IMPL_HPP
#define TEMPUS_TEST_VANDERPOL_IMEX_EXPLICITMODEL_IMPL_HPP

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include <iostream>

namespace Tempus_Test {

template <class Scalar>
VanDerPol_IMEX_ExplicitModel<Scalar>::VanDerPol_IMEX_ExplicitModel(
    Teuchos::RCP<Teuchos::ParameterList> pList_, bool useProductVector)
{
  useProductVector_  = useProductVector;
  isInitialized_     = false;
  dim_               = 2;
  Np_                = 3;  // Number of parameter vectors (p, dx/dp, dx_dot/dp)
  np_                = 1;  // Number of parameters in this vector (1)
  Ng_                = 0;  // Number of observation functions (0)
  ng_                = 0;  // Number of elements in this observation function (0)
  acceptModelParams_ = false;
  haveIC_            = true;
  epsilon_           = 1.0e-06;
  x0_ic_             = 2.0;
  x1_ic_             = 0.0;
  t0_ic_             = 0.0;

  if (useProductVector_ == false) {
    x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
    f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  }
  else {
    // Create using ProductVector so we can use in partitioned IMEX-RK Stepper.
    using Teuchos::RCP;
    Teuchos::Array<RCP<const Thyra::VectorSpaceBase<Scalar> > > yxSpaceArray;
    xSpace_ = Thyra::defaultSpmdVectorSpace<Scalar>(1);
    for (int i = 0; i < dim_; ++i) yxSpaceArray.push_back(xSpace_);
    x_space_ = Thyra::productVectorSpace<Scalar>(yxSpaceArray);
    f_space_ = Thyra::productVectorSpace<Scalar>(yxSpaceArray);
  }

  // Create p_space and g_space
  p_space_    = Thyra::defaultSpmdVectorSpace<Scalar>(np_);
  g_space_    = Thyra::defaultSpmdVectorSpace<Scalar>(ng_);
  dxdp_space_ = Thyra::multiVectorProductVectorSpace(x_space_, np_);

  setParameterList(pList_);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
VanDerPol_IMEX_ExplicitModel<Scalar>::get_x_space() const
{
  return x_space_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
VanDerPol_IMEX_ExplicitModel<Scalar>::get_f_space() const
{
  return f_space_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
VanDerPol_IMEX_ExplicitModel<Scalar>::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
VanDerPol_IMEX_ExplicitModel<Scalar>::create_W_op() const
{
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > W_op;
  if (useProductVector_ == false)
    W_op = Thyra::createMembers(x_space_, dim_);
  else {
    // When using product vector formulation, we have to create a block
    // operator so that its domain space is consistent with x_space_
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > A00 =
        Thyra::createMembers(xSpace_, 1);
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > A01 =
        Thyra::createMembers(xSpace_, 1);
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > A10 =
        Thyra::createMembers(xSpace_, 1);
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > A11 =
        Thyra::createMembers(xSpace_, 1);
    W_op = Thyra::nonconstBlock2x2(A00, A01, A10, A11);
  }
  return W_op;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
VanDerPol_IMEX_ExplicitModel<Scalar>::createInArgs() const
{
  setupInOutArgs_();
  return inArgs_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
VanDerPol_IMEX_ExplicitModel<Scalar>::createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}

template <class Scalar>
void VanDerPol_IMEX_ExplicitModel<Scalar>::setupInOutArgs_() const
{
  if (isInitialized_) {
    return;
  }

  {
    // Set up prototypical InArgs
    Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta);
    if (acceptModelParams_) {
      inArgs.set_Np(Np_);
    }
    inArgs_ = inArgs;
  }

  {
    // Set up prototypical OutArgs
    Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
    if (acceptModelParams_) {
      outArgs.set_Np_Ng(Np_, Ng_);
      outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 0,
                          Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
    }
    outArgs_ = outArgs;
  }

  // Set up nominal values
  nominalValues_ = inArgs_;
  if (haveIC_) {
    using Teuchos::RCP;

    nominalValues_.set_t(t0_ic_);
    const RCP<Thyra::VectorBase<Scalar> > x_ic = createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_ic_view(*x_ic);
      x_ic_view[0] = x0_ic_;
      x_ic_view[1] = x1_ic_;
    }
    nominalValues_.set_x(x_ic);

    if (acceptModelParams_) {
      const RCP<Thyra::VectorBase<Scalar> > p_ic = createMember(p_space_);
      {
        Thyra::DetachedVectorView<Scalar> p_ic_view(*p_ic);
        p_ic_view[0] = epsilon_;
      }
      nominalValues_.set_p(0, p_ic);
    }
    const RCP<Thyra::VectorBase<Scalar> > x_dot_ic = createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_dot_ic_view(*x_dot_ic);
      x_dot_ic_view[0] = x1_ic_;
      x_dot_ic_view[1] = ((1.0 - x0_ic_ * x0_ic_) * x1_ic_ - x0_ic_) / epsilon_;
    }
    nominalValues_.set_x_dot(x_dot_ic);
  }

  isInitialized_ = true;
}

template <class Scalar>
void VanDerPol_IMEX_ExplicitModel<Scalar>::setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const &paramList)
{
  using Teuchos::get;
  using Teuchos::ParameterList;
  Teuchos::RCP<ParameterList> tmpPL =
      Teuchos::rcp(new ParameterList("VanDerPol_IMEX_ExplicitModel"));
  if (paramList != Teuchos::null) tmpPL = paramList;
  tmpPL->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(tmpPL);
  Teuchos::RCP<ParameterList> pl = this->getMyNonconstParamList();
  bool acceptModelParams         = get<bool>(*pl, "Accept model parameters");
  bool haveIC                    = get<bool>(*pl, "Provide nominal values");
  bool useDfDpAsTangent          = get<bool>(*pl, "Use DfDp as Tangent");
  if ((acceptModelParams != acceptModelParams_) || (haveIC != haveIC_)) {
    isInitialized_ = false;
  }
  acceptModelParams_ = acceptModelParams;
  haveIC_            = haveIC;
  useDfDpAsTangent_  = useDfDpAsTangent;
  epsilon_           = get<Scalar>(*pl, "Coeff epsilon");
  x0_ic_             = get<Scalar>(*pl, "IC x0");
  x1_ic_             = get<Scalar>(*pl, "IC x1");
  t0_ic_             = get<Scalar>(*pl, "IC t0");
  setupInOutArgs_();
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
VanDerPol_IMEX_ExplicitModel<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Accept model parameters", false);
    pl->set("Provide nominal values", true);
    pl->set("Use DfDp as Tangent", false);
    Teuchos::setDoubleParameter("Coeff epsilon", 1.0e-06,
                                "Coefficient a in model", &*pl);
    Teuchos::setDoubleParameter("IC x0", 2.0, "Initial Condition for x0", &*pl);
    Teuchos::setDoubleParameter("IC x1", 0.0, "Initial Condition for x1", &*pl);
    Teuchos::setDoubleParameter("IC t0", 0.0, "Initial time t0", &*pl);
    validPL = pl;
  }
  return validPL;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
VanDerPol_IMEX_ExplicitModel<Scalar>::get_p_space(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l, 0, Np_);
  if (l == 0)
    return p_space_;
  else if (l == 1 || l == 2)
    return dxdp_space_;
  return Teuchos::null;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
VanDerPol_IMEX_ExplicitModel<Scalar>::get_p_names(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l, 0, Np_);
  Teuchos::RCP<Teuchos::Array<std::string> > p_strings =
      Teuchos::rcp(new Teuchos::Array<std::string>());
  if (l == 0)
    p_strings->push_back("Model Coefficient:  epsilon");
  else if (l == 1)
    p_strings->push_back("DxDp");
  else if (l == 2)
    p_strings->push_back("Dx_dotDp");
  return p_strings;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
VanDerPol_IMEX_ExplicitModel<Scalar>::get_g_space(int j) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, Ng_);
  return g_space_;
}

template <class Scalar>
void VanDerPol_IMEX_ExplicitModel<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");

  const RCP<const Thyra::VectorBase<Scalar> > x_in =
      inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x_in_view(*x_in);

  // double t = inArgs.get_t();
  Scalar beta = inArgs.get_beta();
  Scalar eps  = epsilon_;
  RCP<const Thyra::MultiVectorBase<Scalar> > dxdp_in, dxdotdp_in;
  if (acceptModelParams_) {
    const RCP<const Thyra::VectorBase<Scalar> > p_in = inArgs.get_p(0);
    if (p_in != Teuchos::null) {
      Thyra::ConstDetachedVectorView<Scalar> p_in_view(*p_in);
      eps = p_in_view[0];
    }
    if (inArgs.get_p(1) != Teuchos::null) {
      dxdp_in = rcp_dynamic_cast<const DMVPV>(inArgs.get_p(1), true)
                    ->getMultiVector();
    }
    if (inArgs.get_p(2) != Teuchos::null) {
      dxdotdp_in = rcp_dynamic_cast<const DMVPV>(inArgs.get_p(2), true)
                       ->getMultiVector();
    }
  }

  const RCP<Thyra::VectorBase<Scalar> > f_out   = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
  RCP<Thyra::MultiVectorBase<Scalar> > DfDp_out;
  if (acceptModelParams_) {
    Thyra::ModelEvaluatorBase::Derivative<Scalar> DfDp = outArgs.get_DfDp(0);
    DfDp_out                                           = DfDp.getMultiVector();
  }

  if (inArgs.get_x_dot().is_null()) {
    // Evaluate the Explicit ODE f(x,t) [= xdot]
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      f_out_view[0] = x_in_view[1];
      f_out_view[1] = -x_in_view[0] / eps;
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDp_out_view(*DfDp_out);
      DfDp_out_view(0, 0) = 0.0;
      DfDp_out_view(1, 0) = x_in_view[0] / (eps * eps);

      // Compute df/dp + (df/dx) * (dx/dp)
      if (useDfDpAsTangent_ && !is_null(dxdp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> dxdp(*dxdp_in);
        DfDp_out_view(0, 0) += dxdp(1, 0);
        DfDp_out_view(1, 0) += -dxdp(0, 0) / eps;
      }
    }
    if (!is_null(W_out)) {
      if (useProductVector_ == false) {
        RCP<Thyra::MultiVectorBase<Scalar> > W =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,
                                                                       true);
        Thyra::DetachedMultiVectorView<Scalar> W_view(*W);
        W_view(0, 0) = 0.0;          // d(f0)/d(x0_n)
        W_view(0, 1) = beta;         // d(f0)/d(x1_n)
        W_view(1, 0) = -beta / eps;  // d(f1)/d(x0_n)
        W_view(1, 1) = 0.0;          // d(f1)/d(x1_n)
        // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
      }
      else {
        RCP<Thyra::DefaultBlockedLinearOp<Scalar> > W_block =
            Teuchos::rcp_dynamic_cast<Thyra::DefaultBlockedLinearOp<Scalar> >(
                W_out, true);
        RCP<Thyra::MultiVectorBase<Scalar> > W00 =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(
                W_block->getNonconstBlock(0, 0), true);
        RCP<Thyra::MultiVectorBase<Scalar> > W01 =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(
                W_block->getNonconstBlock(0, 1), true);
        RCP<Thyra::MultiVectorBase<Scalar> > W10 =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(
                W_block->getNonconstBlock(1, 0), true);
        RCP<Thyra::MultiVectorBase<Scalar> > W11 =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(
                W_block->getNonconstBlock(1, 1), true);
        Thyra::DetachedMultiVectorView<Scalar> W00_view(*W00);
        Thyra::DetachedMultiVectorView<Scalar> W01_view(*W01);
        Thyra::DetachedMultiVectorView<Scalar> W10_view(*W10);
        Thyra::DetachedMultiVectorView<Scalar> W11_view(*W11);
        W00_view(0, 0) = 0.0;          // d(f0)/d(x0_n)
        W01_view(0, 0) = beta;         // d(f0)/d(x1_n)
        W10_view(0, 0) = -beta / eps;  // d(f1)/d(x0_n)
        W11_view(0, 0) = 0.0;          // d(f1)/d(x1_n)
        // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
      }
    }
  }
  else {
    // Evaluate the implicit ODE f(xdot, x, t) [= 0]
    RCP<const Thyra::VectorBase<Scalar> > x_dot_in;
    x_dot_in     = inArgs.get_x_dot().assert_not_null();
    Scalar alpha = inArgs.get_alpha();

    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view(*x_dot_in);
      f_out_view[0] = x_dot_in_view[0] - x_in_view[1];
      f_out_view[1] = x_dot_in_view[1] + x_in_view[0] / eps;
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDp_out_view(*DfDp_out);
      DfDp_out_view(0, 0) = 0.0;
      DfDp_out_view(1, 0) = -x_in_view[0] / (eps * eps);

      // Compute df/dp + (df/dx_dot)*(dx_dot/dp) + (df/dx) * (dx/dp)
      if (useDfDpAsTangent_ && !is_null(dxdotdp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> dxdotdp(*dxdotdp_in);
        DfDp_out_view(0, 0) += dxdotdp(0, 0);
        DfDp_out_view(1, 0) += dxdotdp(1, 0);
      }
      if (useDfDpAsTangent_ && !is_null(dxdp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> dxdp(*dxdp_in);
        DfDp_out_view(0, 0) += -dxdp(1, 0);
        DfDp_out_view(1, 0) += dxdp(0, 0) / eps;
      }
    }
    if (!is_null(W_out)) {
      if (useProductVector_ == false) {
        RCP<Thyra::MultiVectorBase<Scalar> > W =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,
                                                                       true);
        Thyra::DetachedMultiVectorView<Scalar> W_view(*W);
        W_view(0, 0) = alpha;       // d(f0)/d(x0_n)
        W_view(0, 1) = -beta;       // d(f0)/d(x1_n)
        W_view(1, 0) = beta / eps;  // d(f1)/d(x0_n)
        W_view(1, 1) = alpha;       // d(f1)/d(x1_n)
        // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
      }
      else {
        RCP<Thyra::DefaultBlockedLinearOp<Scalar> > W_block =
            Teuchos::rcp_dynamic_cast<Thyra::DefaultBlockedLinearOp<Scalar> >(
                W_out, true);
        RCP<Thyra::MultiVectorBase<Scalar> > W00 =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(
                W_block->getNonconstBlock(0, 0), true);
        RCP<Thyra::MultiVectorBase<Scalar> > W01 =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(
                W_block->getNonconstBlock(0, 1), true);
        RCP<Thyra::MultiVectorBase<Scalar> > W10 =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(
                W_block->getNonconstBlock(1, 0), true);
        RCP<Thyra::MultiVectorBase<Scalar> > W11 =
            Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(
                W_block->getNonconstBlock(1, 1), true);
        Thyra::DetachedMultiVectorView<Scalar> W00_view(*W00);
        Thyra::DetachedMultiVectorView<Scalar> W01_view(*W01);
        Thyra::DetachedMultiVectorView<Scalar> W10_view(*W10);
        Thyra::DetachedMultiVectorView<Scalar> W11_view(*W11);
        W00_view(0, 0) = alpha;       // d(f0)/d(x0_n)
        W01_view(0, 0) = -beta;       // d(f0)/d(x1_n)
        W10_view(0, 0) = beta / eps;  // d(f1)/d(x0_n)
        W11_view(0, 0) = alpha;       // d(f1)/d(x1_n)
        // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
      }
    }
  }
}

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_VANDERPOL_IMEX_EXPLICITMODEL_IMPL_HPP
