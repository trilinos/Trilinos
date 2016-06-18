#ifndef TEMPUS_MODELEVALUATOR_IMPL_HPP
#define TEMPUS_MODELEVALUATOR_IMPL_HPP

namespace Tempus {


template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ResidualModelEvaluator<Scalar>::
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
ResidualModelEvaluator<Scalar>::
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
ResidualModelEvaluator<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  using Teuchos::RCP;

  RCP<const Thyra::VectorBase<Scalar> > x     = inArgs.get_x();
  RCP<Thyra::VectorBase<Scalar> >       x_dot = Thyra::createMember(get_x_space());

  // call functor to compute x dot
  computeXDot_(*x,*x_dot);

  // setup input condition for nonlinear solve
  MEB::InArgs<Scalar> me_inArgs = transientModel_->createInArgs();
  me_inArgs.set_x(x);
  me_inArgs.set_x_dot(x_dot);
  me_inArgs.set_alpha(alpha_);
  me_inArgs.set_beta(beta_);
  for (int i=0; i<transientModel_->Np(); ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      me_inArgs.set_p(i, inArgs.get_p(i));
  }

  // setup output condition
  MEB::OutArgs<Scalar> me_outArgs = transientModel_->createOutArgs();
  me_outArgs.set_f(outArgs.get_f());
  me_outArgs.set_W_op(outArgs.get_W_op());

  // build residual and jacobian
  transientModel_->evalModel(me_inArgs,me_outArgs);
}

} // namespace Tempus

#endif
