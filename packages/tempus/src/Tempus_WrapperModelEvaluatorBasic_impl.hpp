// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_WrapperModelEvaluatorBasic_impl_hpp
#define Tempus_WrapperModelEvaluatorBasic_impl_hpp

namespace Tempus {


template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorBasic<Scalar>::
createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(appModel_->Np());
  inArgs.setSupports(MEB::IN_ARG_x);

  return inArgs;
}


template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorBasic<Scalar>::
createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(appModel_->Np(),0);
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);

  return outArgs;
}


template <typename Scalar>
void
WrapperModelEvaluatorBasic<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  // Setup the evaluation of the ImplicitODE_DAE
  MEB::InArgs<Scalar> appInArgs = appModel_->createInArgs();
  MEB::OutArgs<Scalar> appOutArgs = appModel_->createOutArgs();
  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();

  // setup the temporal input condition for application ME
  appInArgs.set_x(x);
  if (appInArgs.supports(MEB::IN_ARG_t)) appInArgs.set_t(t_);
  for (int i=0; i<appModel_->Np(); ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      appInArgs.set_p(i, inArgs.get_p(i));
  }

  // setup output condition
  appOutArgs.set_f(outArgs.get_f());

  RCP<Thyra::VectorBase<Scalar> > x_dot = Thyra::createMember(get_x_space());
  computeXDot_(*x,*x_dot);
  appInArgs.set_x_dot(x_dot);
  appInArgs.set_alpha(alpha_);
  appInArgs.set_beta(beta_);
  appOutArgs.set_W_op(outArgs.get_W_op());
  appModel_->evalModel(appInArgs,appOutArgs);
}


} // namespace Tempus

#endif  // Tempus_WrapperModelEvaluatorBasic_impl_hpp
