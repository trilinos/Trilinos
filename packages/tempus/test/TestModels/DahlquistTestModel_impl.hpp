//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_DAHLQUIST_TEST_MODEL_IMPL_HPP
#define TEMPUS_TEST_DAHLQUIST_TEST_MODEL_IMPL_HPP

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

//#include <iostream>

namespace Tempus_Test {

template <class Scalar>
DahlquistTestModel<Scalar>::DahlquistTestModel()
{
  constructDahlquistTestModel(-1.0, false);
}

template <class Scalar>
DahlquistTestModel<Scalar>::DahlquistTestModel(Scalar lambda, bool includeXDot)
{
  constructDahlquistTestModel(lambda, includeXDot);
}

template <class Scalar>
void DahlquistTestModel<Scalar>::constructDahlquistTestModel(Scalar lambda,
                                                             bool includeXDot)
{
  lambda_            = lambda;
  includeXDot_       = includeXDot;
  isInitialized_     = false;
  xIC_               = Scalar(1.0);
  xDotIC_            = Scalar(lambda);
  dim_               = 1;
  haveIC_            = true;
  Np_                = 1;
  np_                = 1;
  Ng_                = 1;
  ng_                = dim_;
  acceptModelParams_ = false;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);

  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;
  {
    // Set up prototypical InArgs
    MEB::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(MEB::IN_ARG_t);
    inArgs.setSupports(MEB::IN_ARG_x);
    inArgs.setSupports(MEB::IN_ARG_x_dot);

    inArgs.setSupports(MEB::IN_ARG_beta);
    inArgs.setSupports(MEB::IN_ARG_alpha);
    if (acceptModelParams_) inArgs.set_Np(Np_);

    inArgs_ = inArgs;
  }

  {
    // Set up prototypical OutArgs
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(MEB::OUT_ARG_f);

    outArgs.setSupports(MEB::OUT_ARG_W_op);
    if (acceptModelParams_) {
      outArgs.set_Np_Ng(Np_, Ng_);
      outArgs.setSupports(MEB::OUT_ARG_DfDp, 0, MEB::DERIV_MV_JACOBIAN_FORM);
      outArgs.setSupports(MEB::OUT_ARG_DgDp, 0, 0, MEB::DERIV_MV_JACOBIAN_FORM);
      outArgs.setSupports(MEB::OUT_ARG_DgDx, 0, MEB::DERIV_MV_GRADIENT_FORM);
    }
    outArgs_ = outArgs;
  }

  // Set up nominal values (initial conditions)
  nominalValues_ = inArgs_;
  {
    nominalValues_.set_t(Scalar(0.0));
    const RCP<Thyra::VectorBase<Scalar> > x_ic = createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_ic_view(*x_ic);
      x_ic_view[0] = xIC_;
    }
    nominalValues_.set_x(x_ic);
  }

  if (includeXDot_) {
    const RCP<Thyra::VectorBase<Scalar> > x_dot_ic = createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> x_dot_ic_view(*x_dot_ic);
      x_dot_ic_view[0] = xDotIC_;
    }
    nominalValues_.set_x_dot(x_dot_ic);
  }

  isInitialized_ = true;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DahlquistTestModel<Scalar>::getExactSolution(double t) const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = inArgs_;
  double exact_t                                   = t;
  inArgs.set_t(exact_t);

  // Set the exact solution, x.
  Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x = createMember(x_space_);
  {  // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_view(*exact_x);
    exact_x_view[0] = exp(lambda_ * exact_t);
  }
  inArgs.set_x(exact_x);

  // Set the exact solution time derivative, xDot.
  if (includeXDot_) {
    Teuchos::RCP<Thyra::VectorBase<Scalar> > exact_x_dot =
        createMember(x_space_);
    {  // scope to delete DetachedVectorView
      Thyra::DetachedVectorView<Scalar> exact_x_dot_view(*exact_x_dot);
      exact_x_dot_view[0] = lambda_ * exp(lambda_ * exact_t);
    }
    inArgs.set_x_dot(exact_x_dot);
  }

  return inArgs;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
DahlquistTestModel<Scalar>::get_x_space() const
{
  return x_space_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
DahlquistTestModel<Scalar>::get_f_space() const
{
  return f_space_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DahlquistTestModel<Scalar>::getNominalValues() const
{
  return nominalValues_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
DahlquistTestModel<Scalar>::createInArgs() const
{
  return inArgs_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
DahlquistTestModel<Scalar>::createOutArgsImpl() const
{
  return outArgs_;
}

template <class Scalar>
void DahlquistTestModel<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::MultiVectorBase;
  using Thyra::VectorBase;
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");

  const RCP<const VectorBase<Scalar> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x_in_view(*x_in);
  const RCP<VectorBase<Scalar> > f_out          = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();

  if (inArgs.get_x_dot().is_null()) {
    // Evaluate the Explicit ODE f(x,t) [= 0]
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      f_out_view[0] = lambda_ * x_in_view[0];
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Error -- Dahlquist Test Model requires f_out!\n");
    }

    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,
                                                                     true);
      Thyra::DetachedMultiVectorView<Scalar> matrix_view(*matrix);
      matrix_view(0, 0) = lambda_;
    }
  }
  else {
    // Evaluate the implicit ODE f(xdot, x, t) [=0]
    RCP<const VectorBase<Scalar> > x_dot_in;
    x_dot_in     = inArgs.get_x_dot().assert_not_null();
    Scalar alpha = inArgs.get_alpha();
    Scalar beta  = inArgs.get_beta();

    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view(*x_dot_in);
      f_out_view[0] = x_dot_in_view[0] - lambda_ * x_in_view[0];
    }

    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,
                                                                     true);
      Thyra::DetachedMultiVectorView<Scalar> matrix_view(*matrix);
      matrix_view(0, 0) = alpha - beta * lambda_;  // d(f0)/d(x0_n)
      // Note: alpha = d(xdot)/d(x_n) and beta = d(x)/d(x_n)
    }
  }
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
DahlquistTestModel<Scalar>::create_W() const
{
  using Teuchos::RCP;
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      this->get_W_factory();
  RCP<Thyra::LinearOpBase<Scalar> > matrix = this->create_W_op();
  {
    RCP<Thyra::MultiVectorBase<Scalar> > multivec =
        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(matrix,
                                                                   true);
    {
      RCP<Thyra::VectorBase<Scalar> > vec = Thyra::createMember(x_space_);
      {
        Thyra::DetachedVectorView<Scalar> vec_view(*vec);
        vec_view[0] = lambda_;
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
DahlquistTestModel<Scalar>::create_W_op() const
{
  // const int dim_ = 1;
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Thyra::createMembers(x_space_, dim_);
  return (matrix);
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
DahlquistTestModel<Scalar>::get_W_factory() const
{
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  return W_factory;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
DahlquistTestModel<Scalar>::get_p_space(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l, 0, Np_);
  if (l == 0)
    return p_space_;
  else if (l == 1 || l == 2)
    return DxDp_space_;
  return Teuchos::null;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
DahlquistTestModel<Scalar>::get_p_names(int l) const
{
  if (!acceptModelParams_) {
    return Teuchos::null;
  }
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l, 0, Np_);
  Teuchos::RCP<Teuchos::Array<std::string> > p_strings =
      Teuchos::rcp(new Teuchos::Array<std::string>());
  if (l == 0) {
    p_strings->push_back("Model Coefficient:  a");
    // p_strings->push_back("Model Coefficient:  f");
    // p_strings->push_back("Model Coefficient:  L");
  }
  else if (l == 1)
    p_strings->push_back("DxDp");
  else if (l == 2)
    p_strings->push_back("Dx_dotDp");
  return p_strings;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
DahlquistTestModel<Scalar>::get_g_space(int j) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, Ng_);
  return g_space_;
}

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_DAHLQUIST_TEST_MODEL_IMPL_HPP
