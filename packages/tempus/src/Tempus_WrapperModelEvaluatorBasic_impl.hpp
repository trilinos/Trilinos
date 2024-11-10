//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_WrapperModelEvaluatorBasic_impl_hpp
#define Tempus_WrapperModelEvaluatorBasic_impl_hpp

namespace Tempus {

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorBasic<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  // MEB::InArgsSetup<Scalar> inArgs(appModel_->createInArgs());
  MEB::InArgsSetup<Scalar> inArgs(appModel_->getNominalValues());

  inArgs.set_x(x_);
  if (y_ != Teuchos::null) inArgs.set_p(index_, y_);
  if (inArgs.supports(MEB::IN_ARG_x_dot)) inArgs.set_x_dot(xDot_);
  if (inArgs.supports(MEB::IN_ARG_t)) inArgs.set_t(time_);
  if (inArgs.supports(MEB::IN_ARG_step_size))
    inArgs.set_step_size(p_->timeStepSize_);
  if (inArgs.supports(MEB::IN_ARG_alpha)) inArgs.set_alpha(p_->alpha_);
  if (inArgs.supports(MEB::IN_ARG_beta)) inArgs.set_beta(p_->beta_);
  if (inArgs.supports(MEB::IN_ARG_stage_number))
    inArgs.set_stage_number(p_->stageNumber_);

  inArgs.setModelEvalDescription(this->description());
  return std::move(inArgs);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorBasic<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs(appModel_->createOutArgs());
  outArgs.setModelEvalDescription(this->description());
  return std::move(outArgs);
}

template <typename Scalar>
void WrapperModelEvaluatorBasic<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  using Teuchos::RCP;

  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> appInArgs(inArgs);
  MEB::OutArgs<Scalar> appOutArgs(outArgs);

  // Setup input and output arguments for application ME
  switch (evaluationType_) {
    case EVALUATE_RESIDUAL: {
      // Setup input arguments
      appInArgs.set_x(inArgs.get_x());
      appInArgs.set_x_dot(inArgs.get_x_dot());

      // Setup output arguments
      appOutArgs.set_f(outArgs.get_f());

      break;
    }
    case SOLVE_FOR_X: {
      // This is the normal solution scheme where we are solving for x,
      // and xDot is dependent on x through the definition of the time
      // derivative for the Stepper.

      // Setup input arguments
      RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
      appInArgs.set_x(x);  // x is from solver.
      RCP<Thyra::VectorBase<Scalar> > xDot =
          Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
              appInArgs.get_x_dot());  // xDot is from the Stepper
      timeDer_->compute(x, xDot);
      appInArgs.set_x_dot(xDot);

      // Setup output arguments
      // Note: For the use that Tempus does of this class, these three args
      //       *should* be enough. However, keep in mind that it *may* be
      //       necessary to add more outArgs in the future. The idea would
      //       be the same: if the underlying model supports the arg, then
      //       set it in the appOutArgs.
      appOutArgs.set_f(outArgs.get_f());
      appOutArgs.set_W_op(outArgs.get_W_op());
      if (outArgs.supports(MEB::OUT_ARG_W_prec)) {
        appOutArgs.set_W_prec(outArgs.get_W_prec());
      }

      break;
    }

    case SOLVE_FOR_XDOT_CONST_X: {
      // This solution scheme is solving for xDot while keeping x constant.
      // This is similar to evaluating an explicit ODE, xDot = f(x,t).  The
      // upside is the application does not need to write f(x,t) or worry
      // about inverting a mass matrix.  This solution scheme is mostly
      // used for initial conditions to make the x and xDot consistent
      // with the governing equations.

      // Setup input arguments
      appInArgs.set_x(x_);                  // x is from the Stepper.
      appInArgs.set_x_dot(inArgs.get_x());  // xDot is from the solver.

      // Setup output arguments
      appOutArgs.set_f(outArgs.get_f());
      appOutArgs.set_W_op(outArgs.get_W_op());
      if (outArgs.supports(MEB::OUT_ARG_W_prec)) {
        appOutArgs.set_W_prec(outArgs.get_W_prec());
      }

      break;
    }

    default: {
      TEUCHOS_TEST_FOR_EXCEPT("Invalid EVALUATION_TYPE!");
    }
  }

  appModel_->evalModel(appInArgs, appOutArgs);
}

}  // namespace Tempus

#endif  // Tempus_WrapperModelEvaluatorBasic_impl_hpp
