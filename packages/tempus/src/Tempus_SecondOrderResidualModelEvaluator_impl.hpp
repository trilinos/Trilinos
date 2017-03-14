// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_SecondOrderResidualModelEvaluator_impl_hpp
#define Tempus_SecondOrderResidualModelEvaluator_impl_hpp

namespace Tempus {


template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
SecondOrderResidualModelEvaluator<Scalar>::
createInArgs() const
{
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(transientModel_->Np());
  inArgs.setSupports(MEB::IN_ARG_x);

  return inArgs;
}


template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
SecondOrderResidualModelEvaluator<Scalar>::
createOutArgsImpl() const
{
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
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
SecondOrderResidualModelEvaluator<Scalar>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  //Setup initial condition
  //Create and populate inArgs
  MEB::InArgs<Scalar> transientInArgs = transientModel_->createInArgs();
  //the solution variable in NOX is the acceleration, a_{n+1} 
  transientInArgs.set_x_dot_dot(inArgs.get_x()); 
  RCP<Thyra::VectorBase<Scalar> > velocity = Thyra::createMember(inArgs.get_x()->space());
  //compute the velocity, v_{n+1}(a_{n+1}) = velocity_{pred} + \gamma dt a_{n+1}
  Thyra::V_StVpStV(Teuchos::ptrFromRef(*velocity), 1.0, *v_pred_, delta_t_*gamma_, *inArgs.get_x());
  transientInArgs.set_x_dot(velocity); 
  RCP<Thyra::VectorBase<Scalar> > displacement = Thyra::createMember(inArgs.get_x()->space());
  //compute the displacement, d_{n+1}(a_{n+1}) = displacement_{pred} + \beta dt^2 a_{n+1}
  Thyra::V_StVpStV(Teuchos::ptrFromRef(*displacement), 1.0, *d_pred_, beta_*delta_t_*delta_t_, *inArgs.get_x()); 
  transientInArgs.set_x(displacement); 
  transientInArgs.set_W_x_dot_dot_coeff(1.0);                 // da/da
  transientInArgs.set_alpha(gamma_*delta_t_);                 // dv/da
  transientInArgs.set_beta(beta_*delta_t_*delta_t_);          // dd/da
  transientInArgs.set_t(t_);
  for (int i=0; i<transientModel_->Np(); ++i) {
    if (inArgs.get_p(i) != Teuchos::null)
      transientInArgs.set_p(i, inArgs.get_p(i));
  }


  //Setup output condition 
  //Create and populate outArgs 
  MEB::OutArgs<Scalar> transientOutArgs = transientModel_->createOutArgs();
  transientOutArgs.set_f(outArgs.get_f());
  transientOutArgs.set_W_op(outArgs.get_W_op());

  // build residual and jacobian
  transientModel_->evalModel(transientInArgs,transientOutArgs); 
}


} // namespace Tempus

#endif  // Tempus_SecondOrderResidualModelEvaluator_impl_hpp
