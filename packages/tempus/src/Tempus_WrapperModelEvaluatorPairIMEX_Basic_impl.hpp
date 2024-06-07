//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp
#define Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"

namespace Tempus {

template <typename Scalar>
void WrapperModelEvaluatorPairIMEX_Basic<Scalar>::initialize()
{
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;

  p_     = Teuchos::rcp(new ImplicitODEParameters<Scalar>());
  index_ = -1;

  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgs<Scalar> inArgs = implicitModel_->getNominalValues();
  x_                         = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());

  if (inArgs.supports(MEB::IN_ARG_x_dot)) {
    xDot_ = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x_dot());
  }
  else {
    xDot_ = Teuchos::null;
  }

  // A Thyra::VectorSpace requirement
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(explicitModel_->get_x_space()->isCompatible(
          *(implicitModel_->get_x_space()))),
      std::logic_error,
      "Error - WrapperModelEvaluatorPairIMEX_Basic::initialize()\n"
          << "  Explicit and Implicit vector x spaces are incompatible!\n"
          << "  Explicit vector x space = "
          << *(explicitModel_->get_x_space())
          << "\n  Implicit vector x space = "
          << *(implicitModel_->get_x_space()) << "\n");

  // A Thyra::VectorSpace requirement
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(explicitModel_->get_f_space()->isCompatible(
          *(implicitModel_->get_f_space()))),
      std::logic_error,
      "Error - WrapperModelEvaluatorPairIMEX_Basic::initialize()\n"
          << "  Explicit and Implicit vector f spaces are incompatible!\n"
          << "  Explicit vector f space = "
          << *(explicitModel_->get_f_space())
          << "\n  Implicit vector f space = "
          << *(implicitModel_->get_f_space()) << "\n");
}

template <typename Scalar>
void WrapperModelEvaluatorPairIMEX_Basic<Scalar>::setAppModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& /* me */)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error - WrapperModelEvaluatorPairIMEX_Basic<Scalar>::setAppModel\n"
      "  should not be used.  One should instead use setExplicitModel,\n"
      "  setImplicitModel, or create a new WrapperModelEvaluatorPairIMEX.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::getAppModel() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error,
      "Error - WrapperModelEvaluatorPairIMEX_Basic<Scalar>::getAppModel\n"
      "  should not be used.  One should instead use getExplicitModel,\n"
      "  getImplicitModel, or directly use this WrapperModelEvaluatorPairIMEX\n"
      "  object.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::get_x_space() const
{
  return this->implicitModel_->get_x_space();
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::get_g_space(int i) const
{
  return this->implicitModel_->get_g_space(i);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::get_p_space(int i) const
{
  return this->implicitModel_->get_p_space(i);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::InArgsSetup<Scalar> inArgs = this->createInArgs();
  return std::move(inArgs);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  // MEB::InArgsSetup<Scalar> inArgs(implicitModel_->createInArgs());
  MEB::InArgsSetup<Scalar> inArgs(implicitModel_->getNominalValues());

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
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs(implicitModel_->createOutArgs());
  outArgs.setModelEvalDescription(this->description());
  return std::move(outArgs);
}

template <typename Scalar>
void WrapperModelEvaluatorPairIMEX_Basic<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  RCP<Thyra::VectorBase<Scalar> > x_dot   = Thyra::createMember(get_x_space());
  timeDer_->compute(x, x_dot);

  MEB::InArgs<Scalar> appImplicitInArgs(inArgs);
  MEB::OutArgs<Scalar> appImplicitOutArgs(outArgs);
  appImplicitInArgs.set_x_dot(x_dot);

  implicitModel_->evalModel(appImplicitInArgs, appImplicitOutArgs);
}

}  // end namespace Tempus

#endif  // Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp
