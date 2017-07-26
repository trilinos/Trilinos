// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp
#define Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp

namespace Tempus {


template <typename Scalar>
void
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
setAppModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me)
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
    "Error - WrapperModelEvaluatorPairIMEX_Basic<Scalar>::setAppModel\n"
    "  should not be used.  One should instead use setExplicitModel,\n"
    "  setImplicitModel, or create a new WrapperModelEvaluatorPairIMEX.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
getAppModel() const
{
  TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
    "Error - WrapperModelEvaluatorPairIMEX_Basic<Scalar>::getAppModel\n"
    "  should not be used.  One should instead use getExplicitModel,\n"
    "  getImplicitModel, or directly use this WrapperModelEvaluatorPairIMEX\n"
    "  object.\n");
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
get_x_space() const
{
  return this->getImplicitModel()->get_x_space();
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
get_g_space(int i) const
{
  return this->getImplicitModel()->get_g_space(i);
}

template <typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
get_p_space(int i) const
{
  return this->getImplicitModel()->get_p_space(i);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
getNominalValues() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs = this->createInArgs();

  // Set x, xdot, parameters, ... as needed.

  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
createInArgs() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs = getImplicitModel()->getNominalValues();
  inArgs.setModelEvalDescription(this->description());

  inArgs.setSupports(MEB::IN_ARG_x,true);
  inArgs.setSupports(MEB::IN_ARG_x_dot,true);
  inArgs.setSupports(MEB::IN_ARG_t,true);
  inArgs.setSupports(MEB::IN_ARG_alpha,true);
  inArgs.setSupports(MEB::IN_ARG_beta,true);
  inArgs.setSupports(MEB::IN_ARG_step_size,true);
  inArgs.setSupports(MEB::IN_ARG_stage_number,true);

  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
createOutArgsImpl() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::OutArgsSetup<Scalar> outArgs = getImplicitModel()->createOutArgs();
  outArgs.setModelEvalDescription(this->description());

  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);

  return outArgs;
}

template <typename Scalar>
void
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
evalExplicitModel(
  const Teuchos::RCP<const Thyra::VectorBase<Scalar> > & X, Scalar time,
  const Teuchos::RCP<Thyra::VectorBase<Scalar> > & F) const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgs<Scalar> explicitInArgs = getExplicitModel()->createInArgs();
  // Required inArgs
  explicitInArgs.set_x(X);
  // Optional inArgs
  if (explicitInArgs.supports(MEB::IN_ARG_t))
    explicitInArgs.set_t(time);
  if (explicitInArgs.supports(MEB::IN_ARG_step_size))
    explicitInArgs.set_step_size(stepSize_);
  if (explicitInArgs.supports(MEB::IN_ARG_stage_number))
    explicitInArgs.set_stage_number(stageNumber_);

  // For model evaluators whose state function f(x, x_dot, t) describes
  // an implicit ODE, and which accept an optional x_dot input argument,
  // make sure the latter is set to null in order to request the evaluation
  // of a state function corresponding to the explicit ODE formulation
  // x_dot = f(x, t)
  if (explicitInArgs.supports(MEB::IN_ARG_x_dot))
    explicitInArgs.set_x_dot(Teuchos::null);


  MEB::OutArgs<Scalar> explicitOutArgs = getExplicitModel()->createOutArgs();
  // Required outArgs
  explicitOutArgs.set_f(F);


  getExplicitModel()->evalModel(explicitInArgs,explicitOutArgs);
}

template <typename Scalar>
void
WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
evalImplicitModel(
            const Thyra::ModelEvaluatorBase::InArgs<Scalar>  &inArgs,
            const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  RCP<const Thyra::VectorBase<Scalar> > x = inArgs.get_x();
  RCP<Thyra::VectorBase<Scalar> >   x_dot = Thyra::createMember(get_x_space());

  // call functor to compute x dot
  computeXDot_(*x,*x_dot);

  MEB::InArgs<Scalar> implicitInArgs = getImplicitModel()->createInArgs();
  // Required inArgs
  implicitInArgs.set_x(x);
  implicitInArgs.set_x_dot(x_dot);
  implicitInArgs.set_alpha(alpha_);
  implicitInArgs.set_beta(beta_);
  // Optional inArgs
  if (implicitInArgs.supports(MEB::IN_ARG_t))
    implicitInArgs.set_t(ts_);
  if (implicitInArgs.supports(MEB::IN_ARG_step_size))
    implicitInArgs.set_step_size(stepSize_);
  if (implicitInArgs.supports(MEB::IN_ARG_stage_number))
    implicitInArgs.set_stage_number(stageNumber_);
  for (int i=0; i<getImplicitModel()->Np(); ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      implicitInArgs.set_p(i, inArgs.get_p(i));
  }

  MEB::OutArgs<Scalar> implicitOutArgs = getImplicitModel()->createOutArgs();
  // Required outArgs
  implicitOutArgs.set_f(outArgs.get_f());
  implicitOutArgs.set_W_op(outArgs.get_W_op());

  getImplicitModel()->evalModel(implicitInArgs,implicitOutArgs);
}

template <typename Scalar>
void WrapperModelEvaluatorPairIMEX_Basic<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> & inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> & outArgs) const
{
  evalImplicitModel(inArgs, outArgs);
}

} // end namespace Tempus

#endif // Tempus_ModelEvaluatorIMEXPair_Basic_impl_hpp
