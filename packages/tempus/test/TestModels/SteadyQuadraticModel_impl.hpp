//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_STEADY_QUADRATIC_MODEL_IMPL_HPP
#define TEMPUS_TEST_STEADY_QUADRATIC_MODEL_IMPL_HPP

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

#include <iostream>

namespace Tempus_Test {

template <class Scalar>
SteadyQuadraticModel<Scalar>::SteadyQuadraticModel(
    Teuchos::RCP<Teuchos::ParameterList> pList_)
{
  isInitialized_    = false;
  dim_              = 1;
  Np_               = 3;  // Number of parameter vectors (p, dx/dp, dx_dot/dp)
  np_               = 1;  // Number of parameters in this vector (3)
  Ng_               = 1;  // Number of observation functions (1)
  ng_               = 1;  // Number of elements in this observation function (1)
  useDfDpAsTangent_ = false;
  b_                = 1.0;

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  // Create p_space and g_space
  p_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(np_);
  g_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(ng_);

  setParameterList(pList_);

  // Create DxDp product space
  DxDp_space_ = Thyra::multiVectorProductVectorSpace(x_space_, np_);
}

template <class Scalar>
Scalar SteadyQuadraticModel<Scalar>::getSteadyStateSolution() const
{
  if (b_ < 0.0) return b_;
  return -b_;
}

template <class Scalar>
Scalar SteadyQuadraticModel<Scalar>::getSteadyStateSolutionSensitivity() const
{
  if (b_ < 0.0) return 1.0;
  return -1.0;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
SteadyQuadraticModel<Scalar>::get_x_space() const
{
  return x_space_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
SteadyQuadraticModel<Scalar>::get_f_space() const
{
  return f_space_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
SteadyQuadraticModel<Scalar>::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
SteadyQuadraticModel<Scalar>::create_W() const
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
SteadyQuadraticModel<Scalar>::create_W_op() const
{
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Thyra::createMembers(x_space_, dim_);
  return (matrix);
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
SteadyQuadraticModel<Scalar>::get_W_factory() const
{
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  return W_factory;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
SteadyQuadraticModel<Scalar>::createInArgs() const
{
  setupInOutArgs_();
  return inArgs_;
}

// Private functions overridden from ModelEvaluatorDefaultBase

template <class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
SteadyQuadraticModel<Scalar>::createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}

template <class Scalar>
void SteadyQuadraticModel<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::MultiVectorBase;
  using Thyra::VectorBase;
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");

  const RCP<const VectorBase<Scalar> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x_in_view(*x_in);

  Scalar b                                  = b_;
  const RCP<const VectorBase<Scalar> > p_in = inArgs.get_p(0);
  if (p_in != Teuchos::null) {
    Thyra::ConstDetachedVectorView<Scalar> p_in_view(*p_in);
    b = p_in_view[0];
  }

  RCP<const MultiVectorBase<Scalar> > DxDp_in, DxdotDp_in;
  if (inArgs.get_p(1) != Teuchos::null)
    DxDp_in = rcp_dynamic_cast<const DMVPV>(inArgs.get_p(1))->getMultiVector();
  if (inArgs.get_p(2) != Teuchos::null)
    DxdotDp_in =
        rcp_dynamic_cast<const DMVPV>(inArgs.get_p(2))->getMultiVector();

  Scalar beta = inArgs.get_beta();

  const RCP<VectorBase<Scalar> > f_out          = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
  RCP<Thyra::MultiVectorBase<Scalar> > DfDp_out;
  Thyra::ModelEvaluatorBase::Derivative<Scalar> DfDp = outArgs.get_DfDp(0);
  DfDp_out                                           = DfDp.getMultiVector();

  if (inArgs.get_x_dot().is_null()) {
    // Evaluate the Explicit ODE f(x,t) [= 0]
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      f_out_view[0] = x_in_view[0] * x_in_view[0] - b * b;
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,
                                                                     true);
      Thyra::DetachedMultiVectorView<Scalar> matrix_view(*matrix);
      matrix_view(0, 0) = beta * 2.0 * x_in_view[0];
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDp_out_view(*DfDp_out);
      DfDp_out_view(0, 0) = -2.0 * b;

      // Compute df/dp + (df/dx) * (dx/dp)
      if (useDfDpAsTangent_ && !is_null(DxDp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> DxDp(*DxDp_in);
        DfDp_out_view(0, 0) += 2.0 * x_in_view[0] * DxDp(0, 0);
      }
    }
  }
  else {
    // Evaluate the implicit ODE f(xdot, x, t) [= 0]
    RCP<const VectorBase<Scalar> > x_dot_in;
    x_dot_in     = inArgs.get_x_dot().assert_not_null();
    Scalar alpha = inArgs.get_alpha();
    if (!is_null(f_out)) {
      Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
      Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view(*x_dot_in);
      f_out_view[0] = x_dot_in_view[0] - (x_in_view[0] * x_in_view[0] - b * b);
    }
    if (!is_null(W_out)) {
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out,
                                                                     true);
      Thyra::DetachedMultiVectorView<Scalar> matrix_view(*matrix);
      matrix_view(0, 0) = alpha - beta * 2.0 * x_in_view[0];
    }
    if (!is_null(DfDp_out)) {
      Thyra::DetachedMultiVectorView<Scalar> DfDp_out_view(*DfDp_out);
      DfDp_out_view(0, 0) = 2.0 * b;

      // Compute df/dp + (df/dx_dot) * (dx_dot/dp) + (df/dx) * (dx/dp)
      if (useDfDpAsTangent_ && !is_null(DxdotDp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> DxdotDp(*DxdotDp_in);
        DfDp_out_view(0, 0) += DxdotDp(0, 0);
      }
      if (useDfDpAsTangent_ && !is_null(DxDp_in)) {
        Thyra::ConstDetachedMultiVectorView<Scalar> DxDp(*DxDp_in);
        DfDp_out_view(0, 0) += -2.0 * x_in_view[0] * DxDp(0, 0);
      }
    }
  }

  // Responses:  g = x
  RCP<VectorBase<Scalar> > g_out = outArgs.get_g(0);
  if (g_out != Teuchos::null) Thyra::assign(g_out.ptr(), *x_in);

  RCP<Thyra::MultiVectorBase<Scalar> > DgDp_out =
      outArgs.get_DgDp(0, 0).getMultiVector();
  if (DgDp_out != Teuchos::null) Thyra::assign(DgDp_out.ptr(), Scalar(0.0));

  RCP<Thyra::MultiVectorBase<Scalar> > DgDx_out =
      outArgs.get_DgDx(0).getMultiVector();
  if (DgDx_out != Teuchos::null) {
    Thyra::DetachedMultiVectorView<Scalar> DgDx_out_view(*DgDx_out);
    DgDx_out_view(0, 0) = 1.0;
  }
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
SteadyQuadraticModel<Scalar>::get_p_space(int l) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l, 0, Np_);
  if (l == 0)
    return p_space_;
  else if (l == 1 || l == 2)
    return DxDp_space_;
  return Teuchos::null;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
SteadyQuadraticModel<Scalar>::get_p_names(int l) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(l, 0, Np_);
  Teuchos::RCP<Teuchos::Array<std::string> > p_strings =
      Teuchos::rcp(new Teuchos::Array<std::string>());
  if (l == 0) {
    p_strings->push_back("Model Coefficient:  b");
  }
  else if (l == 1)
    p_strings->push_back("DxDp");
  else if (l == 2)
    p_strings->push_back("Dx_dotDp");
  return p_strings;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
SteadyQuadraticModel<Scalar>::get_g_space(int j) const
{
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, Ng_);
  return g_space_;
}

// private

template <class Scalar>
void SteadyQuadraticModel<Scalar>::setupInOutArgs_() const
{
  if (isInitialized_) {
    return;
  }

  {
    // Set up prototypical InArgs
    Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot);
    inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha);
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
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, 0,
                        Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
    outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0,
                        Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM);
    outArgs_ = outArgs;
  }

  // Set up nominal values
  nominalValues_ = inArgs_;
  nominalValues_.set_t(0.0);
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > x_ic = createMember(x_space_);
  {  // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> x_ic_view(*x_ic);
    x_ic_view[0] = 0.0;
  }
  nominalValues_.set_x(x_ic);
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > p_ic = createMember(p_space_);
  {
    Thyra::DetachedVectorView<Scalar> p_ic_view(*p_ic);
    p_ic_view[0] = b_;
  }
  nominalValues_.set_p(0, p_ic);

  const Teuchos::RCP<Thyra::VectorBase<Scalar> > x_dot_ic =
      createMember(x_space_);
  {  // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> x_dot_ic_view(*x_dot_ic);
    x_dot_ic_view[0] = 0.0;
  }
  nominalValues_.set_x_dot(x_dot_ic);

  isInitialized_ = true;
}

template <class Scalar>
void SteadyQuadraticModel<Scalar>::setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const &paramList)
{
  using Teuchos::get;
  using Teuchos::ParameterList;
  Teuchos::RCP<ParameterList> tmpPL =
      Teuchos::rcp(new ParameterList("SteadyQuadraticModel"));
  if (paramList != Teuchos::null) tmpPL = paramList;
  tmpPL->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(tmpPL);
  Teuchos::RCP<ParameterList> pl = this->getMyNonconstParamList();
  useDfDpAsTangent_              = get<bool>(*pl, "Use DfDp as Tangent");
  b_                             = get<Scalar>(*pl, "Coeff b");
  isInitialized_                 = false;  // For setup of new in/out args
  setupInOutArgs_();
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
SteadyQuadraticModel<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Use DfDp as Tangent", false);
    Teuchos::setDoubleParameter("Coeff b", 1.0, "Coefficient b in model", &*pl);
    validPL = pl;
  }
  return validPL;
}

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_STEADY_QUADRATIC_MODEL_IMPL_HPP
