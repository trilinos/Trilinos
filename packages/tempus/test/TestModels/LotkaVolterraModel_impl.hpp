//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_LOTKAVOLTERRA_MODEL_IMPL_HPP
#define TEMPUS_TEST_LOTKAVOLTERRA_MODEL_IMPL_HPP

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"

#include <cmath>
#include <iostream>

namespace Tempus_Test {

// ============================================================
// Constructor
// ============================================================
template <class Scalar>
LotkaVolterraModel<Scalar>::LotkaVolterraModel(
    Teuchos::RCP<Teuchos::ParameterList> pList_)
{
  isInitialized_ = false;
  dim_           = 2;
  Np_            = 0;   // No parameter vectors (no sensitivity support)
  Ng_            = 0;   // No observation functions
  haveIC_        = true;

  // Default Lotka-Volterra parameters
  alpha_  = Scalar(1.5);
  beta_   = Scalar(1.0);
  delta_  = Scalar(1.0);
  gamma_  = Scalar(3.0);
  k_      = Scalar(0.0);   // no forcing by default
  x0_ic_  = Scalar(10.0);
  y0_ic_  = Scalar(5.0);
  t0_ic_  = Scalar(0.0);

  // Create x_space and f_space
  x_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);
  f_space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim_);

  setParameterList(pList_);
}

// ============================================================
// get_x_space / get_f_space
// ============================================================
template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
LotkaVolterraModel<Scalar>::get_x_space() const
{
  return x_space_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
LotkaVolterraModel<Scalar>::get_f_space() const
{
  return f_space_;
}

// ============================================================
// getNominalValues
// ============================================================
template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
LotkaVolterraModel<Scalar>::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
    "Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

// ============================================================
// create_W (full LOWS — needed for implicit steppers)
// ============================================================
template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
LotkaVolterraModel<Scalar>::create_W() const
{
  using Teuchos::RCP;
  RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      this->get_W_factory();
  RCP<Thyra::LinearOpBase<Scalar> > matrix = this->create_W_op();
  {
    // Initialise to a full-rank identity-like matrix so that the factory
    // can factorise it during construction without hitting a rank-deficient
    // matrix error (same pattern as SinCosModel / VanDerPolModel).
    RCP<Thyra::MultiVectorBase<Scalar> > multivec =
        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(matrix, true);
    {
      RCP<Thyra::VectorBase<Scalar> > vec = Thyra::createMember(x_space_);
      {
        Thyra::DetachedVectorView<Scalar> vv(*vec);
        vv[0] = Scalar(1.0);
        vv[1] = Scalar(0.0);
      }
      V_V(multivec->col(0).ptr(), *vec);
      {
        Thyra::DetachedVectorView<Scalar> vv(*vec);
        vv[0] = Scalar(0.0);
        vv[1] = Scalar(1.0);
      }
      V_V(multivec->col(1).ptr(), *vec);
    }
  }
  RCP<Thyra::LinearOpWithSolveBase<Scalar> > W =
      Thyra::linearOpWithSolve<Scalar>(*W_factory, matrix);
  return W;
}

// ============================================================
// create_W_op (operator only)
// ============================================================
template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
LotkaVolterraModel<Scalar>::create_W_op() const
{
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Thyra::createMembers(x_space_, dim_);
  return matrix;
}

// ============================================================
// get_W_factory
// ============================================================
template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
LotkaVolterraModel<Scalar>::get_W_factory() const
{
  return Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
}

// ============================================================
// createInArgs
// ============================================================
template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
LotkaVolterraModel<Scalar>::createInArgs() const
{
  setupInOutArgs_();
  return inArgs_;
}

// ============================================================
// createOutArgsImpl (private)
// ============================================================
template <class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
LotkaVolterraModel<Scalar>::createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}

// ============================================================
// evalModelImpl (private)
//
// Handles two modes depending on whether x_dot is provided:
//
//  Explicit (x_dot is null):
//    f[0] = alpha*x[0] - beta*x[0]*x[1] - k*sin(t)
//    f[1] = delta*x[0]*x[1] - gamma*x[1]
//
//  Implicit (x_dot provided):
//    F[0] = xdot[0] - alpha*x[0] + beta*x[0]*x[1] + k*sin(t)  = 0
//    F[1] = xdot[1] - delta*x[0]*x[1] + gamma*x[1]            = 0
//
// Because the forcing -k*sin(t) is state-independent the Jacobian
// W = alpha_c * dF/d(xdot) + beta_c * dF/dx is unchanged from
// the unforced case.
// ============================================================
template <class Scalar>
void LotkaVolterraModel<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>  &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  using Teuchos::RCP;
  using Thyra::VectorBase;

  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
    "Error, setupInOutArgs_ must be called first!\n");

  // Extract state and time
  const RCP<const VectorBase<Scalar> > x_in = inArgs.get_x().assert_not_null();
  Thyra::ConstDetachedVectorView<Scalar> x(* x_in);

  const Scalar x0 = x[0];
  const Scalar x1 = x[1];
  const Scalar t  = inArgs.get_t();

  const RCP<VectorBase<Scalar> >          f_out  = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out  = outArgs.get_W_op();

  if (inArgs.get_x_dot().is_null()) {
    // -------------------------------------------------------
    // EXPLICIT ODE:  f(x, t)  [= dx/dt]
    // -------------------------------------------------------
    if (!is_null(f_out)) {
      using std::sin;
      Thyra::DetachedVectorView<Scalar> f(*f_out);
      f[0] = alpha_ * x0 - beta_  * x0 * x1 - k_ * sin(t);  // prey
      f[1] = delta_ * x0 * x1 - gamma_ * x1;                 // predator
    }

    if (!is_null(W_out)) {
      // W = beta_c * df/dx  (forcing is state-independent, no extra terms)
      const Scalar beta_c = inArgs.get_beta();
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out, true);
      Thyra::DetachedMultiVectorView<Scalar> J(*matrix);
      J(0, 0) =  beta_c * (alpha_  - beta_  * x1);   // dF0/dx0
      J(0, 1) =  beta_c * (-beta_  * x0);             // dF0/dx1
      J(1, 0) =  beta_c * ( delta_ * x1);             // dF1/dx0
      J(1, 1) =  beta_c * ( delta_ * x0 - gamma_);   // dF1/dx1
    }
  }
  else {
    // -------------------------------------------------------
    // IMPLICIT ODE:  F(xdot, x, t) = 0
    // F[0] = xdot[0] - alpha*x0 + beta*x0*x1 + k*sin(t)
    // F[1] = xdot[1] - delta*x0*x1 + gamma*x1
    // -------------------------------------------------------
    const RCP<const VectorBase<Scalar> > x_dot_in =
        inArgs.get_x_dot().assert_not_null();
    Thyra::ConstDetachedVectorView<Scalar> xdot(*x_dot_in);

    const Scalar alpha_c = inArgs.get_alpha();
    const Scalar beta_c  = inArgs.get_beta();

    if (!is_null(f_out)) {
      using std::sin;
      Thyra::DetachedVectorView<Scalar> F(*f_out);
      F[0] = xdot[0] - alpha_ * x0 + beta_  * x0 * x1 + k_ * sin(t);
      F[1] = xdot[1] - delta_ * x0 * x1 + gamma_ * x1;
    }

    if (!is_null(W_out)) {
      // W = alpha_c * dF/d(xdot) + beta_c * dF/dx
      // Forcing is state-independent: same Jacobian as the unforced case.
      //
      // dF/d(xdot):  identity
      // dF/dx:
      //   dF0/dx0 = -alpha + beta*x1
      //   dF0/dx1 =  beta*x0
      //   dF1/dx0 = -delta*x1
      //   dF1/dx1 = -delta*x0 + gamma
      RCP<Thyra::MultiVectorBase<Scalar> > matrix =
          Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out, true);
      Thyra::DetachedMultiVectorView<Scalar> J(*matrix);
      J(0, 0) = alpha_c + beta_c * (-alpha_  + beta_  * x1);
      J(0, 1) =           beta_c * ( beta_   * x0);
      J(1, 0) =           beta_c * (-delta_  * x1);
      J(1, 1) = alpha_c + beta_c * (-delta_  * x0 + gamma_);
    }
  }
}

// ============================================================
// get_p_space / get_p_names / get_g_space
// (no parameters / observations — return null)
// ============================================================
template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
LotkaVolterraModel<Scalar>::get_p_space(int /*l*/) const
{
  return Teuchos::null;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
LotkaVolterraModel<Scalar>::get_p_names(int /*l*/) const
{
  return Teuchos::null;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
LotkaVolterraModel<Scalar>::get_g_space(int /*j*/) const
{
  return Teuchos::null;
}

// ============================================================
// setupInOutArgs_ (private, called lazily)
// ============================================================
template <class Scalar>
void LotkaVolterraModel<Scalar>::setupInOutArgs_() const
{
  if (isInitialized_) return;

  using Teuchos::RCP;
  typedef Thyra::ModelEvaluatorBase MEB;

  // InArgs
  {
    MEB::InArgsSetup<Scalar> inArgs;
    inArgs.setModelEvalDescription(this->description());
    inArgs.setSupports(MEB::IN_ARG_t);
    inArgs.setSupports(MEB::IN_ARG_x);
    inArgs.setSupports(MEB::IN_ARG_beta);
    inArgs.setSupports(MEB::IN_ARG_x_dot);
    inArgs.setSupports(MEB::IN_ARG_alpha);
    inArgs_ = inArgs;
  }

  // OutArgs
  {
    MEB::OutArgsSetup<Scalar> outArgs;
    outArgs.setModelEvalDescription(this->description());
    outArgs.setSupports(MEB::OUT_ARG_f);
    outArgs.setSupports(MEB::OUT_ARG_W_op);
    outArgs_ = outArgs;
  }

  // Nominal values / initial conditions
  nominalValues_ = inArgs_;
  if (haveIC_) {
    nominalValues_.set_t(t0_ic_);

    const RCP<Thyra::VectorBase<Scalar> > x_ic = createMember(x_space_);
    {
      Thyra::DetachedVectorView<Scalar> x_ic_view(*x_ic);
      x_ic_view[0] = x0_ic_;
      x_ic_view[1] = y0_ic_;
    }
    nominalValues_.set_x(x_ic);

    // x_dot IC: use the explicit RHS evaluated at the IC
    const RCP<Thyra::VectorBase<Scalar> > x_dot_ic = createMember(x_space_);
    {
      Thyra::DetachedVectorView<Scalar> xd(*x_dot_ic);
      xd[0] = alpha_ * x0_ic_ - beta_  * x0_ic_ * y0_ic_;
      xd[1] = delta_ * x0_ic_ * y0_ic_ - gamma_ * y0_ic_;
    }
    nominalValues_.set_x_dot(x_dot_ic);
  }

  isInitialized_ = true;
}

// ============================================================
// setParameterList
// ============================================================
template <class Scalar>
void LotkaVolterraModel<Scalar>::setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const &paramList)
{
  using Teuchos::get;
  using Teuchos::ParameterList;

  Teuchos::RCP<ParameterList> tmpPL =
      Teuchos::rcp(new ParameterList("LotkaVolterraModel"));
  if (paramList != Teuchos::null) tmpPL = paramList;
  tmpPL->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(tmpPL);
  Teuchos::RCP<ParameterList> pl = this->getMyNonconstParamList();

  bool haveIC = get<bool>(*pl, "Provide nominal values");
  if (haveIC != haveIC_) isInitialized_ = false;
  haveIC_ = haveIC;

  alpha_  = get<Scalar>(*pl, "Coeff alpha");
  beta_   = get<Scalar>(*pl, "Coeff beta");
  delta_  = get<Scalar>(*pl, "Coeff delta");
  gamma_  = get<Scalar>(*pl, "Coeff gamma");
  k_      = get<Scalar>(*pl, "Coeff k");
  x0_ic_  = get<Scalar>(*pl, "IC x0");
  y0_ic_  = get<Scalar>(*pl, "IC y0");
  t0_ic_  = get<Scalar>(*pl, "IC t0");

  setupInOutArgs_();
}

// ============================================================
// getValidParameters
// ============================================================
template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
LotkaVolterraModel<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->set("Provide nominal values", true);
    Teuchos::setDoubleParameter("Coeff alpha", 1.5, "Prey growth rate",              &*pl);
    Teuchos::setDoubleParameter("Coeff beta",  1.0, "Predation rate",                &*pl);
    Teuchos::setDoubleParameter("Coeff delta", 1.0, "Predator growth rate",          &*pl);
    Teuchos::setDoubleParameter("Coeff gamma", 3.0, "Predator death rate",           &*pl);
    Teuchos::setDoubleParameter("Coeff k",     0.0, "Prey forcing amplitude k*sin(t)",&*pl);
    Teuchos::setDoubleParameter("IC x0",      10.0, "Initial prey population",       &*pl);
    Teuchos::setDoubleParameter("IC y0",       5.0, "Initial pred population",       &*pl);
    Teuchos::setDoubleParameter("IC t0",       0.0, "Initial time t0",               &*pl);
    validPL = pl;
  }
  return validPL;
}

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_LOTKAVOLTERRA_MODEL_IMPL_HPP
