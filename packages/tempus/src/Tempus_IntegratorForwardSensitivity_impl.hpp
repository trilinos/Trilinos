//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_IntegratorForwardSensitivity_impl_hpp
#define Tempus_IntegratorForwardSensitivity_impl_hpp

#include "Tempus_SolutionHistory_decl.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "Tempus_CombinedForwardSensitivityModelEvaluator.hpp"
#include "Tempus_WrapCombinedFSAModelEvaluator.hpp"

namespace Tempus {

template <class Scalar>
IntegratorForwardSensitivity<Scalar>::IntegratorForwardSensitivity(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
    const Teuchos::RCP<IntegratorBasic<Scalar>> &integrator,
    const Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar>> &sens_model,
    const Teuchos::RCP<StepperStaggeredForwardSensitivity<Scalar>>
        &sens_stepper,
    bool use_combined_method)
  : model_(model),
    integrator_(integrator),
    sens_model_(sens_model),
    sens_stepper_(sens_stepper),
    use_combined_method_(use_combined_method)
{
  integrator_->initialize();
}

template <class Scalar>
IntegratorForwardSensitivity<Scalar>::IntegratorForwardSensitivity()
{
  integrator_ = createIntegratorBasic<Scalar>();
  integrator_->initialize();
}

template <class Scalar>
void IntegratorForwardSensitivity<Scalar>::setStepper(
    Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> model)
{
  if (use_combined_method_)
    integrator_->setModel(sens_model_);
  else
    integrator_->setStepper(sens_stepper_);
}

template <class Scalar>
void IntegratorForwardSensitivity<Scalar>::initializeSolutionHistory(
    Scalar t0, Teuchos::RCP<const Thyra::VectorBase<Scalar>> x0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdot0,
    Teuchos::RCP<const Thyra::VectorBase<Scalar>> xdotdot0,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxDp0,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxdotDp0,
    Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>> DxdotdotDp0)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::assign;
  using Thyra::createMember;
  using Thyra::VectorSpaceBase;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  //
  // Create and initialize product X, Xdot, Xdotdot

  RCP<const VectorSpaceBase<Scalar>> space;
  if (use_combined_method_)
    space = sens_model_->get_x_space();
  else
    space = sens_stepper_->get_x_space();
  RCP<DMVPV> X       = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Xdot    = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Xdotdot = rcp_dynamic_cast<DMVPV>(createMember(space));

  const int num_param = X->getNonconstMultiVector()->domain()->dim() - 1;
  const Scalar zero   = Teuchos::ScalarTraits<Scalar>::zero();
  const Teuchos::Range1D rng(1, num_param);

  // x
  assign(X->getNonconstMultiVector()->col(0).ptr(), *x0);
  if (DxDp0 == Teuchos::null)
    assign(X->getNonconstMultiVector()->subView(rng).ptr(), zero);
  else
    assign(X->getNonconstMultiVector()->subView(rng).ptr(), *DxDp0);

  // xdot
  if (xdot0 == Teuchos::null)
    assign(Xdot->getNonconstMultiVector()->col(0).ptr(), zero);
  else
    assign(Xdot->getNonconstMultiVector()->col(0).ptr(), *xdot0);
  if (DxdotDp0 == Teuchos::null)
    assign(Xdot->getNonconstMultiVector()->subView(rng).ptr(), zero);
  else
    assign(Xdot->getNonconstMultiVector()->subView(rng).ptr(), *DxdotDp0);

  // xdotdot
  if (xdotdot0 == Teuchos::null)
    assign(Xdotdot->getNonconstMultiVector()->col(0).ptr(), zero);
  else
    assign(Xdotdot->getNonconstMultiVector()->col(0).ptr(), *xdotdot0);
  if (DxdotDp0 == Teuchos::null)
    assign(Xdotdot->getNonconstMultiVector()->subView(rng).ptr(), zero);
  else
    assign(Xdotdot->getNonconstMultiVector()->subView(rng).ptr(), *DxdotdotDp0);

  integrator_->initializeSolutionHistory(t0, X, Xdot, Xdotdot);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorForwardSensitivity<Scalar>::getX() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> X = rcp_dynamic_cast<const DMVPV>(integrator_->getX());
  return X->getMultiVector()->col(0);
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorForwardSensitivity<Scalar>::getDxDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> X  = rcp_dynamic_cast<const DMVPV>(integrator_->getX());
  const int num_param = X->getMultiVector()->domain()->dim() - 1;
  const Teuchos::Range1D rng(1, num_param);
  return X->getMultiVector()->subView(rng);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorForwardSensitivity<Scalar>::getXDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdot = rcp_dynamic_cast<const DMVPV>(integrator_->getXDot());
  return Xdot->getMultiVector()->col(0);
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorForwardSensitivity<Scalar>::getDXDotDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdot = rcp_dynamic_cast<const DMVPV>(integrator_->getXDot());
  const int num_param   = Xdot->getMultiVector()->domain()->dim() - 1;
  const Teuchos::Range1D rng(1, num_param);
  return Xdot->getMultiVector()->subView(rng);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorForwardSensitivity<Scalar>::getXDotDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdotdot =
      rcp_dynamic_cast<const DMVPV>(integrator_->getXDotDot());
  return Xdotdot->getMultiVector()->col(0);
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorForwardSensitivity<Scalar>::getDXDotDotDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdotdot =
      rcp_dynamic_cast<const DMVPV>(integrator_->getXDotDot());
  const int num_param = Xdotdot->getMultiVector()->domain()->dim() - 1;
  const Teuchos::Range1D rng(1, num_param);
  return Xdotdot->getMultiVector()->subView(rng);
}

template <class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar>>
IntegratorForwardSensitivity<Scalar>::getG() const
{
  typedef Thyra::ModelEvaluatorBase MEB;

  // Compute g which is computed by response 1 of the
  // sensitivity model evaluator
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> smodel;
  if (use_combined_method_)
    smodel = sens_model_;
  else
    smodel = sens_stepper_->getModel();
  MEB::InArgs<Scalar> inargs   = smodel->getNominalValues();
  MEB::OutArgs<Scalar> outargs = smodel->createOutArgs();
  inargs.set_t(integrator_->getTime());
  inargs.set_x(integrator_->getX());
  if (inargs.supports(MEB::IN_ARG_x_dot))
    inargs.set_x_dot(integrator_->getXDot());
  if (inargs.supports(MEB::IN_ARG_x_dot_dot))
    inargs.set_x_dot_dot(integrator_->getXDotDot());

  Teuchos::RCP<Thyra::VectorBase<Scalar>> g =
      Thyra::createMember(smodel->get_g_space(1));
  outargs.set_g(1, g);

  smodel->evalModel(inargs, outargs);
  return g;
}

template <class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar>>
IntegratorForwardSensitivity<Scalar>::getDgDp() const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  // Compute final dg/dp  which is computed by response 0  of the
  // sensitivity model evaluator
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>> smodel;
  if (use_combined_method_)
    smodel = sens_model_;
  else
    smodel = sens_stepper_->getModel();
  MEB::InArgs<Scalar> inargs   = smodel->getNominalValues();
  MEB::OutArgs<Scalar> outargs = smodel->createOutArgs();
  inargs.set_t(integrator_->getTime());
  inargs.set_x(integrator_->getX());
  if (inargs.supports(MEB::IN_ARG_x_dot))
    inargs.set_x_dot(integrator_->getXDot());
  if (inargs.supports(MEB::IN_ARG_x_dot_dot))
    inargs.set_x_dot_dot(integrator_->getXDotDot());

  Teuchos::RCP<Thyra::VectorBase<Scalar>> G =
      Thyra::createMember(smodel->get_g_space(0));
  Teuchos::RCP<DMVPV> dgdp = Teuchos::rcp_dynamic_cast<DMVPV>(G);
  outargs.set_g(0, G);

  smodel->evalModel(inargs, outargs);
  return dgdp->getMultiVector();
}

template <class Scalar>
std::string IntegratorForwardSensitivity<Scalar>::description() const
{
  std::string name = "Tempus::IntegratorForwardSensitivity";
  return (name);
}

template <class Scalar>
void IntegratorForwardSensitivity<Scalar>::describe(
    Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << description() << "::describe" << std::endl;
  integrator_->describe(*l_out, verbLevel);
}

template <class Scalar>
SensitivityStepMode IntegratorForwardSensitivity<Scalar>::getStepMode() const
{
  if (use_combined_method_) return SensitivityStepMode::Combined;
  return sens_stepper_->getStepMode();  // Staggered case
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar>>
createIntegratorForwardSensitivity(
    Teuchos::RCP<Teuchos::ParameterList> pList,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &sens_residual_model,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar>> &sens_solve_model)
{
  Teuchos::RCP<SensitivityModelEvaluatorBase<Scalar>> sens_model;
  Teuchos::RCP<StepperStaggeredForwardSensitivity<Scalar>> sens_stepper;

  // 1. create integrator
  auto fwd_integrator = createIntegratorBasic<Scalar>(pList, model);

  // 2. set parameter list
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::rcp(new Teuchos::ParameterList);
  {
    Teuchos::RCP<const Teuchos::ParameterList> integrator_pl =
        fwd_integrator->getValidParameters();
    Teuchos::RCP<const Teuchos::ParameterList> sensitivity_pl =
        CombinedForwardSensitivityModelEvaluator<Scalar>::getValidParameters();
    pl->setParameters(*integrator_pl);
    Teuchos::ParameterList &spl = pl->sublist("Sensitivities");
    spl.setParameters(*sensitivity_pl);
    spl.set("Sensitivity Method", "Combined");
    spl.set("Reuse State Linear Solver", false);
  }
  pList->setParametersNotAlreadySet(*pl);

  auto sens_pl = Teuchos::sublist(pList, "Sensitivities", false);
  std::string integratorName =
      pList->get<std::string>("Integrator Name", "Default Integrator");
  std::string stepperName =
      pList->sublist(integratorName).get<std::string>("Stepper Name");
  auto stepper_pl = Teuchos::sublist(pList, stepperName, true);
  std::string sensitivity_method =
      sens_pl->get<std::string>("Sensitivity Method");
  bool use_combined_method = sensitivity_method == "Combined";

  // 3. create sensitivity model and stepper
  //  createSensitivityModelAndStepper
  {
    sens_pl->remove("Sensitivity Method");

    if (use_combined_method) {
      sens_pl->remove("Reuse State Linear Solver");
      sens_model     = wrapCombinedFSAModelEvaluator(model, sens_residual_model,
                                                     sens_solve_model, sens_pl);
      fwd_integrator = createIntegratorBasic<Scalar>(pList, sens_model);
    }
    else {
      sens_stepper =
          Teuchos::rcp(new StepperStaggeredForwardSensitivity<Scalar>(
              model, sens_residual_model, sens_solve_model, stepper_pl,
              sens_pl));
      auto fsa_staggered_me =
          Teuchos::rcp_const_cast<Thyra::ModelEvaluator<Scalar>>(
              sens_stepper->getModel());
      fwd_integrator = createIntegratorBasic<Scalar>(pList, fsa_staggered_me);
      fwd_integrator->setStepper(sens_stepper);
    }
  }

  // 4. initialize propoer integrator
  Teuchos::RCP<IntegratorForwardSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorForwardSensitivity<Scalar>(
          model, fwd_integrator, sens_model, sens_stepper,
          use_combined_method));
  return (integrator);
}

/// Nonmember constructor
template <class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar>>
createIntegratorForwardSensitivity()
{
  Teuchos::RCP<IntegratorForwardSensitivity<Scalar>> integrator =
      Teuchos::rcp(new IntegratorForwardSensitivity<Scalar>());
  return (integrator);
}

}  // namespace Tempus
#endif  // Tempus_IntegratorForwardSensitivity_impl_hpp
