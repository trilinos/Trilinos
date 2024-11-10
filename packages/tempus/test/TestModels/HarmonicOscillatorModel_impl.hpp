//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_TEST_HARMONIC_OSCILLATOR_MODEL_IMPL_HPP
#define TEMPUS_TEST_HARMONIC_OSCILLATOR_MODEL_IMPL_HPP

#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_VectorStdOps.hpp"
#include <iostream>

namespace Tempus_Test {

template <class Scalar>
HarmonicOscillatorModel<Scalar>::HarmonicOscillatorModel(
    Teuchos::RCP<Teuchos::ParameterList> pList_, const bool use_accel_IC)
  : out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  isInitialized_ = false;
  setParameterList(pList_);
  *out_ << "\n\nDamping coeff c = " << c_ << "\n";
  *out_ << "Forcing coeff f = " << f_ << "\n";
  *out_ << "x coeff k = " << k_ << "\n";
  *out_ << "Mass coeff m = " << m_ << "\n";
  // Divide all coefficients by m_
  k_ /= m_;
  f_ /= m_;
  c_ /= m_;
  m_ = 1.0;
  // Set up space and initial guess for solution vector
  vecLength_ = 1;
  x_space_   = Thyra::defaultSpmdVectorSpace<Scalar>(vecLength_);
  x_vec_     = createMember(x_space_);
  Thyra::put_scalar(0.0, x_vec_.ptr());
  x_dot_vec_ = createMember(x_space_);
  Thyra::put_scalar(1.0, x_dot_vec_.ptr());
  x_dot_dot_vec_ = createMember(x_space_);
  // The following is the initial condition for the acceleration
  // Commenting this out to check that IC for acceleration
  // is computed correctly using displacement and velocity ICs
  // inside 2nd order steppers.
  // Thyra::put_scalar(f_-c_, x_dot_dot_vec_.ptr());
  if (use_accel_IC == true) {
    Thyra::put_scalar(-2.0, x_dot_dot_vec_.ptr());
  }
  else {
    // Instead of real IC, putting arbitrary, incorrect IC to check correctness
    // in stepper involving calculation of a IC.
    Thyra::put_scalar(7.0, x_dot_dot_vec_.ptr());
  }

  // Set up responses
  numResponses_ = 1;
  g_space_      = Thyra::defaultSpmdVectorSpace<Scalar>(numResponses_);

  setupInOutArgs_();
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
HarmonicOscillatorModel<Scalar>::getExactSolution(double t) const
{
  using Thyra::VectorBase;
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");
  Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = inArgs_;
  double exact_t                                   = t;
  inArgs.set_t(exact_t);
  Teuchos::RCP<VectorBase<Scalar> > exact_x = createMember(x_space_);
  {  // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_view(*exact_x);
    if (k_ == 0) {
      if (c_ == 0)
        exact_x_view[0] = t * (1.0 + 0.5 * f_ * t);
      else
        exact_x_view[0] =
            (c_ - f_) / (c_ * c_) * (1.0 - exp(-c_ * t)) + f_ * t / c_;
    }
    else {
      exact_x_view[0] = 1.0 / sqrt(k_) * sin(sqrt(k_) * t) +
                        f_ / k_ * (1.0 - cos(sqrt(k_) * t));
    }
  }
  inArgs.set_x(exact_x);
  Teuchos::RCP<VectorBase<Scalar> > exact_x_dot = createMember(x_space_);
  {  // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_dot_view(*exact_x_dot);
    if (k_ == 0) {
      if (c_ == 0)
        exact_x_dot_view[0] = 1.0 + f_ * t;
      else
        exact_x_dot_view[0] = (c_ - f_) / c_ * exp(-c_ * t) + f_ / c_;
    }
    else {
      exact_x_dot_view[0] =
          cos(sqrt(k_) * t) + f_ / sqrt(k_) * sin(sqrt(k_) * t);
    }
  }
  inArgs.set_x_dot(exact_x_dot);
  Teuchos::RCP<VectorBase<Scalar> > exact_x_dot_dot = createMember(x_space_);
  {  // scope to delete DetachedVectorView
    Thyra::DetachedVectorView<Scalar> exact_x_dot_dot_view(*exact_x_dot_dot);
    if (k_ == 0) {
      if (c_ == 0)
        exact_x_dot_dot_view[0] = f_;
      else
        exact_x_dot_dot_view[0] = (f_ - c_) * exp(-c_ * t);
    }
    else {
      exact_x_dot_dot_view[0] =
          f_ * cos(sqrt(k_) * t) - sqrt(k_) * sin(sqrt(k_) * t);
    }
  }
  inArgs.set_x_dot_dot(exact_x_dot_dot);
  return (inArgs);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
HarmonicOscillatorModel<Scalar>::get_x_space() const
{
  return x_space_;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
HarmonicOscillatorModel<Scalar>::get_f_space() const
{
  return x_space_;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
HarmonicOscillatorModel<Scalar>::getNominalValues() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");
  return nominalValues_;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> >
HarmonicOscillatorModel<Scalar>::create_W() const
{
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      this->get_W_factory();
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > matrix = this->create_W_op();
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix_mv =
      Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(matrix, true);
  Thyra::DetachedMultiVectorView<Scalar> matrix_view(*matrix_mv);
  // IKT: is it necessary for W to be non-singular when initialized?
  matrix_view(0, 0) = 1.0;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > W =
      Thyra::linearOpWithSolve<Scalar>(*W_factory, matrix);
  return W;
}

template <class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
HarmonicOscillatorModel<Scalar>::create_W_op() const
{
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > matrix =
      Thyra::createMembers(x_space_, vecLength_);
  return (matrix);
}

template <class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
HarmonicOscillatorModel<Scalar>::get_W_factory() const
{
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory =
      Thyra::defaultSerialDenseLinearOpWithSolveFactory<Scalar>();
  return W_factory;
}

template <class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
HarmonicOscillatorModel<Scalar>::createInArgs() const
{
  setupInOutArgs_();
  return inArgs_;
}

// Private functions overridden from ModelEvaluatorDefaultBase

template <class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
HarmonicOscillatorModel<Scalar>::createOutArgsImpl() const
{
  setupInOutArgs_();
  return outArgs_;
}

template <class Scalar>
void HarmonicOscillatorModel<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  using Teuchos::RCP;
  using Thyra::VectorBase;
  TEUCHOS_TEST_FOR_EXCEPTION(!isInitialized_, std::logic_error,
                             "Error, setupInOutArgs_ must be called first!\n");

  RCP<const VectorBase<Scalar> > x_in = inArgs.get_x();
  double beta                         = inArgs.get_beta();
  if (!x_in.get()) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "\n ERROR: HarmonicOscillatorModel requires x as InArgs.\n");
  }
  Thyra::ConstDetachedVectorView<Scalar> x_in_view(*x_in);
  // IKT, FIXME: check that subDim() is the write routine to get local length of
  // a Thyra::ConstDetachedVectorView
  auto myVecLength = x_in_view.subDim();

  RCP<const VectorBase<Scalar> > x_dot_in = inArgs.get_x_dot();
  double alpha                            = inArgs.get_alpha();

  RCP<const VectorBase<Scalar> > x_dotdot_in = inArgs.get_x_dot_dot();
  double omega                               = inArgs.get_W_x_dot_dot_coeff();

  // Parse OutArgs
  RCP<VectorBase<Scalar> > f_out                = outArgs.get_f();
  RCP<VectorBase<Scalar> > g_out                = outArgs.get_g(0);
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();

  Scalar neg_sign = 1.0;
  // Explicit ODE
  if (inArgs.get_x_dot_dot().is_null()) neg_sign = -1.0;

  // Populate residual and Jacobian
  if (f_out != Teuchos::null) {
    Thyra::DetachedVectorView<Scalar> f_out_view(*f_out);
    for (int i = 0; i < myVecLength; i++) {
      f_out_view[i] = f_;
    }
    if (x_dotdot_in != Teuchos::null) {
      Thyra::ConstDetachedVectorView<Scalar> x_dotdot_in_view(*x_dotdot_in);
      for (int i = 0; i < myVecLength; i++) {
        f_out_view[i] = x_dotdot_in_view[i] - f_out_view[i];
      }
    }
    if (x_dot_in != Teuchos::null) {
      Thyra::ConstDetachedVectorView<Scalar> x_dot_in_view(*x_dot_in);
      for (int i = 0; i < myVecLength; i++) {
        f_out_view[i] += neg_sign * c_ * x_dot_in_view[i];
      }
    }
    if (x_in != Teuchos::null) {
      for (int i = 0; i < myVecLength; i++) {
        f_out_view[i] += neg_sign * k_ * x_in_view[i];
      }
    }
  }

  // Note: W = alpha*df/dxdot + beta*df/dx + omega*df/dxdotdot
  if (W_out != Teuchos::null) {
    RCP<Thyra::MultiVectorBase<Scalar> > matrix =
        Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase<Scalar> >(W_out, true);
    Thyra::DetachedMultiVectorView<Scalar> matrix_view(*matrix);
    if (omega == 0.0) {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "\n ERROR: omega = 0 in HarmonicOscillatorModel!\n");
    }
    matrix_view(0, 0) = omega;
    if (x_dot_in != Teuchos::null) {
      matrix_view(0, 0) += neg_sign * c_ * alpha;
    }
    if (x_in != Teuchos::null) {
      matrix_view(0, 0) += neg_sign * k_ * beta;
    }
  }

  // Calculated response(s) g
  // g = mean value of x
  if (g_out != Teuchos::null) {
    Thyra::DetachedVectorView<Scalar> g_out_view(*g_out);
    g_out_view[0] = Thyra::sum(*x_in) / vecLength_;
  }
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
HarmonicOscillatorModel<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "\n Error!  HarmonicOscillatorModel::get_p_space() is not supported!\n");
  return Teuchos::null;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
HarmonicOscillatorModel<Scalar>::get_p_names(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "\n Error!  HarmonicOscillatorModel::get_p_names() is not supported!\n");
  return Teuchos::null;
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
HarmonicOscillatorModel<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j != 0, std::logic_error,
      "\n Error!  HarmonicOscillatorModel::get_g_space() only "
          << " supports 1 parameter vector.  Supplied index l = " << j << "\n");
  return g_space_;
}

// private

template <class Scalar>
void HarmonicOscillatorModel<Scalar>::setupInOutArgs_() const
{
  if (isInitialized_) return;

  // Set up InArgs
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(0);
  inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x);
  inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot);
  inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x_dot_dot);
  inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_t);
  inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_W_x_dot_dot_coeff);
  inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_alpha);
  inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_beta);
  inArgs_ = inArgs;

  // Set up OutArgs
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(0, numResponses_);

  outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
  outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W_op);
  // outArgs.setSupports(OUT_ARG_W,true);
  // IKT, what is the following supposed to do??
  // outArgs.set_W_properties( DerivativeProperties(
  //     DERIV_LINEARITY_UNKNOWN, DERIV_RANK_FULL, true));
  outArgs_ = outArgs;

  // Set up nominal values
  nominalValues_ = inArgs_;
  nominalValues_.set_t(0.0);
  nominalValues_.set_x(x_vec_);
  nominalValues_.set_x_dot(x_dot_vec_);
  nominalValues_.set_x_dot_dot(x_dot_dot_vec_);

  isInitialized_ = true;
}

template <class Scalar>
void HarmonicOscillatorModel<Scalar>::setParameterList(
    Teuchos::RCP<Teuchos::ParameterList> const &paramList)
{
  using Teuchos::get;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  RCP<ParameterList> tmpPL =
      Teuchos::rcp(new ParameterList("HarmonicOscillatorModel"));
  if (paramList != Teuchos::null) tmpPL = paramList;
  tmpPL->validateParametersAndSetDefaults(*this->getValidParameters());
  this->setMyParamList(tmpPL);
  RCP<ParameterList> pl = this->getMyNonconstParamList();
  c_                    = get<Scalar>(*pl, "Damping coeff c");
  f_                    = get<Scalar>(*pl, "Forcing coeff f");
  k_                    = get<Scalar>(*pl, "x coeff k");
  m_                    = get<Scalar>(*pl, "Mass coeff m");
  if (m_ <= 0.0) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Error: invalid value of Mass coeff m = "
                                   << m_ << "!  Mass coeff m must be > 0.\n");
  }
  if (k_ < 0.0) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "Error: invalid value of x coeff k = "
                                   << k_ << "!  x coeff k must be >= 0.\n");
  }
  if ((k_ > 0.0) && (c_ != 0.0)) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error: HarmonicOscillator model only supports x coeff k > 0 when "
        "Damping coeff c = 0.  You have "
            << "specified x coeff k = " << k_ << " and Damping coeff c = " << c_
            << ".\n");
  }
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
HarmonicOscillatorModel<Scalar>::getValidParameters() const
{
  static Teuchos::RCP<const Teuchos::ParameterList> validPL;
  if (is_null(validPL)) {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    validPL                                 = pl;
    Teuchos::setDoubleParameter("Damping coeff c", 0.0,
                                "Damping coefficient in model", &*pl);
    Teuchos::setDoubleParameter("Forcing coeff f", -1.0,
                                "Forcing coefficient in model", &*pl);
    Teuchos::setDoubleParameter("x coeff k", 0.0, "x coefficient in model",
                                &*pl);
    Teuchos::setDoubleParameter("Mass coeff m", 1.0,
                                "Mass coefficient in model", &*pl);
  }
  return validPL;
}

}  // namespace Tempus_Test
#endif  // TEMPUS_TEST_HARMONIC_OSCILLATOR_MODEL_IMPL_HPP
