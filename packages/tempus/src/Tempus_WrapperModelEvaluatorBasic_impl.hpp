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
  //MEB::InArgsSetup<Scalar> inArgs(appModel_->createInArgs());
  MEB::InArgsSetup<Scalar> inArgs(appModel_->getNominalValues());
  inArgs.setModelEvalDescription(this->description());
  return inArgs;
}


template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorBasic<Scalar>::
createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgsSetup<Scalar> outArgs(appModel_->createOutArgs());
  outArgs.setModelEvalDescription(this->description());
  return outArgs;
}


template <typename Scalar>
void
WrapperModelEvaluatorBasic<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar>  &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  MEB::InArgs<Scalar>  appInArgs (wrapperInArgs_);
  MEB::OutArgs<Scalar> appOutArgs(wrapperOutArgs_);

  // Setup input condition for application ME
  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  appInArgs.set_x(x);
  for (int i=0; i<appModel_->Np(); ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      appInArgs.set_p(i, inArgs.get_p(i));
  }

  RCP<Thyra::VectorBase<Scalar> > x_dot = Thyra::createMember(get_x_space());
  timeDer_->compute(x, x_dot);
  appInArgs.set_x_dot(x_dot);

  // Setup output condition
  // Note: For the use that Tempus does of this class, these three args *should*
  //       be enough. However, keep in mind that it *may* be necessary to add more
  //       out args in the future. The idea would be the same: if the underlying
  //       modele supports the arg, then set it in the appOutArgs.
  appOutArgs.set_f(outArgs.get_f());
  appOutArgs.set_W_op(outArgs.get_W_op());
  if (outArgs.supports(MEB::OUT_ARG_W_prec)) {
    appOutArgs.set_W_prec(outArgs.get_W_prec());
  }

  appModel_->evalModel(appInArgs,appOutArgs);
}


} // namespace Tempus

#endif  // Tempus_WrapperModelEvaluatorBasic_impl_hpp
