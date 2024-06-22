//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_WrapperModelEvaluatorSecondOrder_impl_hpp
#define Tempus_WrapperModelEvaluatorSecondOrder_impl_hpp

namespace Tempus {

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
WrapperModelEvaluatorSecondOrder<Scalar>::createInArgs() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(appModel_->Np());
  inArgs.setSupports(MEB::IN_ARG_x);

  return std::move(inArgs);
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
WrapperModelEvaluatorSecondOrder<Scalar>::createOutArgsImpl() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  typedef Thyra::ModelEvaluatorBase MEB;

  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(appModel_->Np(), 0);
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);

  return std::move(outArgs);
}

template <typename Scalar>
void WrapperModelEvaluatorSecondOrder<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  typedef Thyra::ModelEvaluatorBase MEB;
  using Teuchos::RCP;

  // Setup initial condition
  // Create and populate inArgs
  MEB::InArgs<Scalar> appInArgs   = appModel_->createInArgs();
  MEB::OutArgs<Scalar> appOutArgs = appModel_->createOutArgs();

  switch (schemeType_) {
    case NEWMARK_IMPLICIT_AFORM: {
      // Specific for the Newmark Implicit a-Form stepper.  May want to
      // redesign this for a generic second-order scheme to not have an
      // if statement here...
      // IKT, 3/14/17: this is effectively the same as the
      // Piro::NewmarkDecorator::evalModel function.

      // The solution variable in NOX is the acceleration, a_{n+1}
      appInArgs.set_x_dot_dot(inArgs.get_x());

      // Compute the velocity.
      // v_n = velocity_{pred} + \gamma dt a_n
      auto velocity = Thyra::createMember(inArgs.get_x()->space());
      Thyra::V_StVpStV(velocity.ptr(), Scalar(1.0), *v_pred_, delta_t_ * gamma_,
                       *inArgs.get_x());
      appInArgs.set_x_dot(velocity);

      // Compute the displacement.
      // d_n = displacement_{pred} + \beta dt^2 a_n
      auto displacement = Thyra::createMember(inArgs.get_x()->space());
      Thyra::V_StVpStV(displacement.ptr(), Scalar(1.0), *d_pred_,
                       beta_ * delta_t_ * delta_t_, *inArgs.get_x());
      appInArgs.set_x(displacement);

      appInArgs.set_W_x_dot_dot_coeff(Scalar(1.0));     // da/da
      appInArgs.set_alpha(gamma_ * delta_t_);           // dv/da
      appInArgs.set_beta(beta_ * delta_t_ * delta_t_);  // dd/da

      appInArgs.set_t(t_);
      for (int i = 0; i < appModel_->Np(); ++i) {
        if (inArgs.get_p(i) != Teuchos::null)
          appInArgs.set_p(i, inArgs.get_p(i));
      }

      // Setup output condition
      // Create and populate outArgs
      appOutArgs.set_f(outArgs.get_f());
      appOutArgs.set_W_op(outArgs.get_W_op());

      // build residual and jacobian
      appModel_->evalModel(appInArgs, appOutArgs);
      break;
    }

    case NEWMARK_IMPLICIT_DFORM: {
      // Setup initial condition
      // Populate inArgs
      RCP<Thyra::VectorBase<Scalar> const> d = inArgs.get_x();

      RCP<Thyra::VectorBase<Scalar>> v =
          Thyra::createMember(inArgs.get_x()->space());

      RCP<Thyra::VectorBase<Scalar>> a =
          Thyra::createMember(inArgs.get_x()->space());

#ifdef DEBUG_OUTPUT
      Teuchos::Range1D range;

      *out_ << "\n*** d_bef ***\n";
      RTOpPack::ConstSubVectorView<Scalar> dov;
      d->acquireDetachedView(range, &dov);
      auto doa = dov.values();
      for (auto i = 0; i < doa.size(); ++i) *out_ << doa[i] << " ";
      *out_ << "\n*** d_bef ***\n";

      *out_ << "\n*** v_bef ***\n";
      RTOpPack::ConstSubVectorView<Scalar> vov;
      v->acquireDetachedView(range, &vov);
      auto voa = vov.values();
      for (auto i = 0; i < voa.size(); ++i) *out_ << voa[i] << " ";
      *out_ << "\n*** v_bef ***\n";

      *out_ << "\n*** a_bef ***\n";
      RTOpPack::ConstSubVectorView<Scalar> aov;
      a->acquireDetachedView(range, &aov);
      auto aoa = aov.values();
      for (auto i = 0; i < aoa.size(); ++i) *out_ << aoa[i] << " ";
      *out_ << "\n*** a_bef ***\n";
#endif

      Scalar const c = 1.0 / beta_ / delta_t_ / delta_t_;

      // compute acceleration
      // a_{n+1} = (d_{n+1} - d_pred) / dt / dt / beta
      Thyra::V_StVpStV(Teuchos::ptrFromRef(*a), c, *d, -c, *d_pred_);

      // compute velocity
      // v_{n+1} = v_pred + \gamma dt a_{n+1}
      Thyra::V_StVpStV(Teuchos::ptrFromRef(*v), 1.0, *v_pred_,
                       delta_t_ * gamma_, *a);

      appInArgs.set_x(d);
      appInArgs.set_x_dot(v);
      appInArgs.set_x_dot_dot(a);

      appInArgs.set_W_x_dot_dot_coeff(c);              // da/dd
      appInArgs.set_alpha(gamma_ / delta_t_ / beta_);  // dv/dd
      appInArgs.set_beta(1.0);                         // dd/dd

      appInArgs.set_t(t_);
      for (int i = 0; i < appModel_->Np(); ++i) {
        if (inArgs.get_p(i) != Teuchos::null)
          appInArgs.set_p(i, inArgs.get_p(i));
      }

      // Setup output condition
      // Create and populate outArgs
      appOutArgs.set_f(outArgs.get_f());
      appOutArgs.set_W_op(outArgs.get_W_op());

      // build residual and jacobian
      appModel_->evalModel(appInArgs, appOutArgs);

      // compute acceleration
      // a_{n+1} = (d_{n+1} - d_pred) / dt / dt / beta
      Thyra::V_StVpStV(Teuchos::ptrFromRef(*a), c, *d, -c, *d_pred_);

      // compute velocity
      // v_{n+1} = v_pred + \gamma dt a_{n+1}
      Thyra::V_StVpStV(Teuchos::ptrFromRef(*v), 1.0, *v_pred_,
                       delta_t_ * gamma_, *a);

      appInArgs.set_x(d);
      appInArgs.set_x_dot(v);
      appInArgs.set_x_dot_dot(a);

#ifdef DEBUG_OUTPUT
      *out_ << "\n*** d_aft ***\n";
      RTOpPack::ConstSubVectorView<Scalar> dnv;
      d->acquireDetachedView(range, &dnv);
      auto dna = dnv.values();
      for (auto i = 0; i < dna.size(); ++i) *out_ << dna[i] << " ";
      *out_ << "\n*** d_aft ***\n";

      *out_ << "\n*** v_aft ***\n";
      RTOpPack::ConstSubVectorView<Scalar> vnv;
      v->acquireDetachedView(range, &vnv);
      auto vna = vnv.values();
      for (auto i = 0; i < vna.size(); ++i) *out_ << vna[i] << " ";
      *out_ << "\n*** v_aft ***\n";

      *out_ << "\n*** a_aft ***\n";
      RTOpPack::ConstSubVectorView<Scalar> anv;
      a->acquireDetachedView(range, &anv);
      auto ana = anv.values();
      for (auto i = 0; i < ana.size(); ++i) *out_ << ana[i] << " ";
      *out_ << "\n*** a_aft ***\n";
#endif
      break;
    }
  }
}

}  // namespace Tempus

#endif  // Tempus_WrapperModelEvaluatorSecondOrder_impl_hpp
