//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_VANDERPOL_IMEXPart_IMPLICITMODEL_IMPL_HPP
#define TEMPUS_TEST_VANDERPOL_IMEXPart_IMPLICITMODEL_IMPL_HPP

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

#include <iostream>

namespace Tempus_Test {

template <class Scalar>
VanDerPol_IMEXPart_ImplicitModel<Scalar>::VanDerPol_IMEXPart_ImplicitModel(
    Teuchos::RCP<Teuchos::ParameterList> pList_)
{
  isInitialized_ = false;
  dim_           = 1;
  Np_            = 5;     // Number of parameter vectors (p, dx/dp, dx_dot/dp, dydp, y)
                          // Due to deficiencies in WrapperModelEvaluatorPairPartIMEX, y must
                          // be last
  np_               = 1;  // Number of parameters in this vector (for p and y) (1)
  Ng_               = 0;  // Number of observation functions (0)
  ng_               = 0;  // Number of elements in this observation function (1)
  useDfDpAsTangent_ = false;
  haveIC_           = true;
  epsilon_          = 1.0e-06;
  x0_ic_            = 2.0;
  x1_ic_            = 0.0;
  t0_ic_            = 0.0;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  // Create parameter spaces
  y_space_    = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  p_space_    = Thyra::defaultSpmdVectorSpace<Scalar>(np_);
  dxdp_space_ = Thyra::multiVectorProductVectorSpace(x_space_, np_);
  dydp_space_ = Thyra::multiVectorProductVectorSpace(y_space_, np_);

  // g_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(ng_);

  setParameterList(pList_);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
VanDerPol_IMEXPart_ImplicitModel<Scalar>::get_x_space() const
{
  return x_space_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
VanDerPol_IMEXPart_ImplicitModel<Scalar>::get_f_space() const
{
  return f_space_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
VanDerPol_IMEXPart_ImplicitModel<Scalar>::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
VanDerPol_IMEXPart_ImplicitModel<Scalar>::create_W() const
{
  using Teuchos::RCP;
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      this->get_W_factory();
  RCP<Thyra::LinearOpBase<Scalar> > matrix = this->create_W_op();
  {
    // 01/20/09 tscoffe:  This is a total hack to provide a full rank matrix to
    // linearOpWithSolve because it ends up factoring the matrix during
    // initialization, which it really shouldn't do, or I'm doing something
    // wrong here.   The net effect is that I get exceptions thrown in
    // optimized mode due to the matrix being rank deficient unless I do this.
    RCP<Thyra::MultiVectorBase<Scalar> > multivec =
        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(matrix,
                                                                   true);
    {
      RCP<Thyra::VectorBase<Scalar> > vec = Thyra::createMember(x_space_);
      {
        Thyra::DetachedVectorView<Scalar> vec_view(*vec);
        vec_view[0] = 1.0;
      }
      V_V(multivec->col(0).ptr(), *vec);
    }
  }
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > W =
      Thyra::linearOpWithSolve<Scalar>(*W_factory, matrix);

  return W;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
VanDerPol_IMEXPart_ImplicitModel<Scalar>::create_W_op() const
{
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Thyra::createMembers(x_space_, dim_);
  return (matrix);
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
VanDerPol_IMEXPart_ImplicitModel<Scalar>::get_W_factory() const
{
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  return W_factory;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
VanDerPol_IMEXPart_ImplicitModel<Scalar>::createInArgs() const
{
  setupInOutArgs_();
  return inArgs_;
}

// Private functions overridden from ModelEvaluatorDefaultBase

template <class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
VanDerPol_IMEXPart_ImplicitModel<Scalar>::createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}

template <class Scalar>
void VanDerPol_IMEXPart_ImplicitModel<Scalar>::evalModelImpl(
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
  Scalar x1 = x_in_view[0];

  // double t = inArgs.get_t();
  Scalar beta = inArgs.get_beta();
  Scalar eps  = epsilon_;

  const RCP<const Thyra::VectorBase<Scalar> > p_in = inArgs.get_p(0);
  if (p_in != Teuchos::null) {
    Thyra::ConstDetachedVectorView<Scalar> p_in_view(*p_in);
    eps = p_in_view[0];
  }

  RCP<const Thyra::MultiVectorBase<Scalar> > dxdp_in, dxdotdp_in, dydp_in;
  if (inArgs.get_p(1) != Teuchos::null) {
    dxdp_in =
        rcp_dynamic_cast<const DMVPV>(inArgs.get_p(1), true)->getMultiVector();
  }
  if (inArgs.get_p(2) != Teuchos::null) {
    dxdotdp_in =
        rcp_dynamic_cast<const DMVPV>(inArgs.get_p(2), true)->getMultiVector();
  }
  if (inArgs.get_p(3) != Teuchos::null) {
    dydp_in =
        rcp_dynamic_cast<const DMVPV>(inArgs.get_p(3), true)->getMultiVector();
  }

  const RCP<const Thyra::VectorBase<Scalar> > y_in =
      inArgs.get_p(4).assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> y_in_view(*y_in);
  Scalar x0 = y_in_view[0];

  const RCP<Thyra::VectorBase<Scalar> > f_out   = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
  const RCP<Thyra::MultiVectorBase<Scalar> > DfDp_out =
      outArgs.get_DfDp(0).getMultiVector();
  const RCP<Thyra::MultiVectorBase<Scalar> > DfDx0_out =
      outArgs.get_DfDp(4).getMultiVector();

  if (inArgs.get_x_dot().is_null()) {
    // Evaluate the Explicit ODE f(x,t) [= xdot]
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      f_out_view[0] = (1.0 - x0 * x0) * x1 / eps;
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDp_out_view(*DfDp_out);
      DfDp_out_view(0, 0) = -(1.0 - x0 * x0) * x1 / (eps * eps);

      // Compute df/dp + (df/dx) * (dx/dp) + (df/dy)*dy/dp
      if (useDfDpAsTangent_ && !is_null(dxdp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> dxdp(*dxdp_in);
        DfDp_out_view(0, 0) += (1.0 - x0 * x0) / eps * dxdp(0, 0);
      }
      if (useDfDpAsTangent_ && !is_null(dydp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> dydp(*dydp_in);
        DfDp_out_view(0, 0) += -2.0 * x0 * x1 / eps * dydp(0, 0);
      }
    }
    if (!is_null(DfDx0_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDx0_out_view(*DfDx0_out);
      DfDx0_out_view(0, 0) = -2.0 * x0 * x1 / eps;  // d(f0)/d(x0_n);
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > W =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,
                                                                     true);
      Thyra::DetachedMultiVectorView<Scalar> W_view(*W);
      W_view(0, 0) = beta * (1.0 - x0 * x0) / eps;  // d(f0)/d(x1_n)
      // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
    }
  }
  else {
    // Evaluate the implicit ODE f(xdot, x, t) [= 0]
    RCP<const Thyra::VectorBase<Scalar> > x_dot_in;
    x_dot_in = inArgs.get_x_dot().assert_not_null();
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view(*x_dot_in);
      f_out_view[0] = x_dot_in_view[0] - (1.0 - x0 * x0) * x1 / eps;
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDp_out_view(*DfDp_out);
      DfDp_out_view(0, 0) = (1.0 - x0 * x0) * x1 / (eps * eps);

      // Compute df/dp + (df/dx_dot)*(dx_dot/dp) + (df/dx)*(dx/dp) +
      // (df/dy)*dy/dp
      if (useDfDpAsTangent_ && !is_null(dxdotdp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> dxdotdp(*dxdotdp_in);
        DfDp_out_view(0, 0) += dxdotdp(0, 0);
      }
      if (useDfDpAsTangent_ && !is_null(dxdp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> dxdp(*dxdp_in);
        DfDp_out_view(0, 0) += -(1.0 - x0 * x0) / eps * dxdp(0, 0);
      }
      if (useDfDpAsTangent_ && !is_null(dydp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> dydp(*dydp_in);
        DfDp_out_view(0, 0) += 2.0 * x0 * x1 / eps * dydp(0, 0);
      }
    }
    if (!is_null(DfDx0_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDx0_out_view(*DfDx0_out);
      DfDx0_out_view(0, 0) = 2.0 * x0 * x1 / eps;  // d(f0)/d(x0_n);
    }
    if (!is_null(W_out)) {
      Scalar alpha = inArgs.get_alpha();
      RCP<Thyra::MultiVectorBase<Scalar> > W =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,
                                                                     true);
      Thyra::DetachedMultiVectorView<Scalar> W_view(*W);
      W_view(0, 0) = alpha - beta * (1.0 - x0 * x0) / eps;  // d(f0)/d(x1_n)
      // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
    }
  }
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
VanDerPol_IMEXPart_ImplicitModel<Scalar>::get_p_space(int l) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l, 0, Np_);
  if (l == 0)
    return p_space_;
  else if (l == 1 || l == 2)
    return dxdp_space_;
  else if (l == 3)
    return dydp_space_;
  else if (l == 4)
    return y_space_;
  return Teuchos::null;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
VanDerPol_IMEXPart_ImplicitModel<Scalar>::get_p_names(int l) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l, 0, Np_);
  Teuchos::RCP<Teuchos::Array<std::string> > p_strings =
      Teuchos::rcp(new Teuchos::Array<std::string>());
  if (l == 0)
    p_strings->push_back("Model Coefficient:  epsilon");
  else if (l == 1)
    p_strings->push_back("DxDp");
  else if (l == 2)
    p_strings->push_back("Dx_dotDp");
  else if (l == 3)
    p_strings->push_back("DyDp");
  else if (l == 4)
    p_strings->push_back("EXPLICIT_ONLY_VECTOR");
  return p_strings;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
VanDerPol_IMEXPart_ImplicitModel<Scalar>::get_g_space(int j) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, Ng_);
  return g_space_;
}

template <class Scalar>
void VanDerPol_IMEXPart_ImplicitModel<Scalar>::setupInOutArgs_() const
{
  if (isInitialized_) return;

  {
    // Set up prototypical InArgs
    Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta);
    inArgs.set_Np(Np_);
    inArgs_ = inArgs;
  }

  {
    // Set up prototypical OutArgs
    Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
    outArgs.set_Np_Ng(Np_, Ng_);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 0,
                        Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, 4,
                        Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
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
      x_ic_view[0] = x1_ic_;
    }
    nominalValues_.set_x(x_ic);

    const RCP<Thyra::VectorBase<Scalar> > p_ic = createMember(p_space_);
    const RCP<Thyra::VectorBase<Scalar> > y_ic = createMember(y_space_);
    {
      Thyra::DetachedVectorView<Scalar> p_ic_view(*p_ic);
      Thyra::DetachedVectorView<Scalar> y_ic_view(*y_ic);
      p_ic_view[0] = epsilon_;
      y_ic_view[0] = x0_ic_;
    }
    nominalValues_.set_p(0, p_ic);
    nominalValues_.set_p(4, y_ic);

    const RCP<Thyra::VectorBase<Scalar> > x_dot_ic = createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_dot_ic_view(*x_dot_ic);
      x_dot_ic_view[0] = ((1.0 - x0_ic_ * x0_ic_) * x1_ic_ - x0_ic_) / epsilon_;
    }
    nominalValues_.set_x_dot(x_dot_ic);
  }

  isInitialized_ = true;
}

template <class Scalar>
void VanDerPol_IMEXPart_ImplicitModel<Scalar>::setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const &paramList)
{
  using Teuchos::get;
  using Teuchos::ParameterList;
  Teuchos::RCP<ParameterList> tmpPL =
      Teuchos::rcp(new ParameterList("VanDerPol_IMEXPart_ImplicitModel"));
  if (paramList != Teuchos::null) tmpPL = paramList;
  tmpPL->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(tmpPL);
  Teuchos::RCP<ParameterList> pl = this->getMyNonconstParamList();
  bool haveIC                    = get<bool>(*pl, "Provide nominal values");
  bool useDfDpAsTangent          = get<bool>(*pl, "Use DfDp as Tangent");
  if (haveIC != haveIC_) isInitialized_ = false;
  haveIC_           = haveIC;
  useDfDpAsTangent_ = useDfDpAsTangent;
  epsilon_          = get<Scalar>(*pl, "Coeff epsilon");
  x0_ic_            = get<Scalar>(*pl, "IC x0");
  x1_ic_            = get<Scalar>(*pl, "IC x1");
  t0_ic_            = get<Scalar>(*pl, "IC t0");
  setupInOutArgs_();
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
VanDerPol_IMEXPart_ImplicitModel<Scalar>::getValidParameters() const
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

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_VANDERPOL_IMEXPart_IMPLICITMODEL_IMPL_HPP
