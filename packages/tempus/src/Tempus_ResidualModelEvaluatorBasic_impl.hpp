// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ResidualModelEvaluatorBasic_impl_hpp
#define Tempus_ResidualModelEvaluatorBasic_impl_hpp

namespace Tempus {


template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ResidualModelEvaluatorBasic<Scalar>::
createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(transientModel_->Np());
  inArgs.setSupports(MEB::IN_ARG_x);

  return inArgs;
}


template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ResidualModelEvaluatorBasic<Scalar>::
createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(transientModel_->Np(),0);
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);

  return outArgs;
}


template <typename Scalar>
void
ResidualModelEvaluatorBasic<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  // Setup the evaluation of the ImplicitODE_DAE
  MEB::InArgs<Scalar> transientInArgs = transientModel_->createInArgs();
  MEB::OutArgs<Scalar> transientOutArgs = transientModel_->createOutArgs();
  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();

  // setup the temporal input condition for application ME
  transientInArgs.set_x(x);
  if (transientInArgs.supports(MEB::IN_ARG_t)) transientInArgs.set_t(t_);
  for (int i=0; i<transientModel_->Np(); ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      transientInArgs.set_p(i, inArgs.get_p(i));
  }

  // setup output condition
  transientOutArgs.set_f(outArgs.get_f());

  RCP<Thyra::VectorBase<Scalar> > x_dot = Thyra::createMember(get_x_space());
  computeXDot_(*x,*x_dot);
  transientInArgs.set_x_dot(x_dot);
  transientInArgs.set_alpha(alpha_);
  transientInArgs.set_beta(beta_);
  transientOutArgs.set_W_op(outArgs.get_W_op());
  transientModel_->evalModel(transientInArgs,transientOutArgs);
}


} // namespace Tempus

#endif  // Tempus_ResidualModelEvaluatorBasic_impl_hpp
