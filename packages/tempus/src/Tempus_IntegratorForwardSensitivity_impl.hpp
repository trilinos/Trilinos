// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorForwardSensitivity_impl_hpp
#define Tempus_IntegratorForwardSensitivity_impl_hpp

#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

#include "Tempus_CombinedForwardSensitivityModelEvaluator.hpp"
#include "Tempus_WrapCombinedFSAModelEvaluator.hpp"


namespace Tempus {

template<class Scalar>
IntegratorForwardSensitivity<Scalar>::
IntegratorForwardSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                inputPL,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model)
{
  model_ = model;
  integrator_ = integratorBasic<Scalar>();
  this->setParameterList(inputPL);
  createSensitivityModelAndStepper(model);
  if (use_combined_method_)
    integrator_ = integratorBasic<Scalar>(tempus_pl_, sens_model_);
  else {
    integrator_ = integratorBasic<Scalar>();
    integrator_->setTempusParameterList(tempus_pl_);
    integrator_->setStepperWStepper(sens_stepper_);
    integrator_->initialize();
  }
}

template<class Scalar>
IntegratorForwardSensitivity<Scalar>::
IntegratorForwardSensitivity(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  std::string stepperType)
{
  model_ = model;
  integrator_ = integratorBasic<Scalar>();
  this->setParameterList(Teuchos::null);
  createSensitivityModelAndStepper(model);
  if (use_combined_method_)
    integrator_ = integratorBasic<Scalar>(sens_model_, stepperType);
  else {
    integrator_ = integratorBasic<Scalar>();
    integrator_->setParameterList(tempus_pl_);
    integrator_->setStepperWStepper(sens_stepper_);
    integrator_->initialize();
  }

}

template<class Scalar>
IntegratorForwardSensitivity<Scalar>::
IntegratorForwardSensitivity()
{
  integrator_ = integratorBasic<Scalar>();
  this->setParameterList(Teuchos::null);
}

template<class Scalar>
void IntegratorForwardSensitivity<Scalar>::
setStepper(
  Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > model)
{
  createSensitivityModelAndStepper(model);
  if (use_combined_method_)
    integrator_->setStepper(sens_model_);
  else
    integrator_->setStepperWStepper(sens_stepper_);
}

template<class Scalar>
void IntegratorForwardSensitivity<Scalar>::
initializeSolutionHistory(Scalar t0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxDp0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxdotDp0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > DxdotdotDp0)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorSpaceBase;
  using Thyra::assign;
  using Thyra::createMember;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  //
  // Create and initialize product X, Xdot, Xdotdot

  RCP< const VectorSpaceBase<Scalar> > space;
  if (use_combined_method_)
    space = sens_model_->get_x_space();
  else
    space = sens_stepper_->get_x_space();
  RCP<DMVPV> X       = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Xdot    = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Xdotdot = rcp_dynamic_cast<DMVPV>(createMember(space));

  const int num_param = X->getNonconstMultiVector()->domain()->dim()-1;
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
  const Teuchos::Range1D rng(1,num_param);

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

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorForwardSensitivity<Scalar>::
getX() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> X = rcp_dynamic_cast<const DMVPV>(integrator_->getX());
  return X->getMultiVector()->col(0);
}

template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
IntegratorForwardSensitivity<Scalar>::
getDxDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> X = rcp_dynamic_cast<const DMVPV>(integrator_->getX());
  const int num_param = X->getMultiVector()->domain()->dim()-1;
  const Teuchos::Range1D rng(1,num_param);
  return X->getMultiVector()->subView(rng);
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorForwardSensitivity<Scalar>::
getXDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdot = rcp_dynamic_cast<const DMVPV>(integrator_->getXDot());
  return Xdot->getMultiVector()->col(0);
}

template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
IntegratorForwardSensitivity<Scalar>::
getDXDotDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdot = rcp_dynamic_cast<const DMVPV>(integrator_->getXDot());
  const int num_param = Xdot->getMultiVector()->domain()->dim()-1;
  const Teuchos::Range1D rng(1,num_param);
  return Xdot->getMultiVector()->subView(rng);
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorForwardSensitivity<Scalar>::
getXDotDot() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdotdot =
    rcp_dynamic_cast<const DMVPV>(integrator_->getXDotDot());
  return Xdotdot->getMultiVector()->col(0);
}

template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
IntegratorForwardSensitivity<Scalar>::
getDXDotDotDp() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  RCP<const DMVPV> Xdotdot =
    rcp_dynamic_cast<const DMVPV>(integrator_->getXDotDot());
  const int num_param = Xdotdot->getMultiVector()->domain()->dim()-1;
  const Teuchos::Range1D rng(1,num_param);
  return Xdotdot->getMultiVector()->subView(rng);
}

template<class Scalar>
std::string
IntegratorForwardSensitivity<Scalar>::
description() const
{
  std::string name = "Tempus::IntegratorForwardSensitivity";
  return(name);
}

template<class Scalar>
void
IntegratorForwardSensitivity<Scalar>::
describe(
  Teuchos::FancyOStream          &in_out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  auto out = Teuchos::fancyOStream( in_out.getOStream() );
  out->setOutputToRootOnly(0);
  *out << description() << "::describe" << std::endl;
  integrator_->describe(in_out, verbLevel);
}

template<class Scalar>
void
IntegratorForwardSensitivity<Scalar>::
setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & inputPL)
{
  tempus_pl_ = Teuchos::parameterList();
  if (inputPL != Teuchos::null)
    *tempus_pl_ = *inputPL;
  tempus_pl_->setParametersNotAlreadySet(*this->getValidParameters());
  sens_pl_ = Teuchos::sublist(tempus_pl_, "Sensitivities", false);
  std::string integratorName =
    tempus_pl_->get<std::string>("Integrator Name", "Default Integrator");
  std::string stepperName =
    tempus_pl_->sublist(integratorName).get<std::string>("Stepper Name");
  stepper_pl_ = Teuchos::sublist(tempus_pl_, stepperName, true);
  use_combined_method_ =
    sens_pl_->get<std::string>("Sensitivity Method") == "Combined";
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorForwardSensitivity<Scalar>::
unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_param_list = tempus_pl_;
  tempus_pl_ = Teuchos::null;
  sens_pl_ = Teuchos::null;
  stepper_pl_ = Teuchos::null;
  integrator_->unsetParameterList();
  return temp_param_list;
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorForwardSensitivity<Scalar>::
getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::RCP<const Teuchos::ParameterList> integrator_pl =
    integrator_->getValidParameters();
  Teuchos::RCP<const Teuchos::ParameterList> sensitivity_pl =
    CombinedForwardSensitivityModelEvaluator<Scalar>::getValidParameters();
  pl->setParameters(*integrator_pl);
  Teuchos::ParameterList& spl = pl->sublist("Sensitivities");
  spl.setParameters(*sensitivity_pl);
  spl.set("Sensitivity Method", "Combined");
  spl.set("Reuse State Linear Solver", false);

  return pl;
}

template <class Scalar>
void
IntegratorForwardSensitivity<Scalar>::
createSensitivityModelAndStepper(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& /* model */)
{
  using Teuchos::rcp;

  Teuchos::RCP<Teuchos::ParameterList> spl = Teuchos::parameterList();
  *spl = *sens_pl_;
  spl->remove("Sensitivity Method");

  if (use_combined_method_) {
    spl->remove("Reuse State Linear Solver");
    sens_model_ =
      wrapCombinedFSAModelEvaluator(model_, spl);
  }
  else {
    sens_stepper_ =
      rcp(new StepperStaggeredForwardSensitivity<Scalar>(
            model_, stepper_pl_, spl));
  }
}

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar> >
integratorForwardSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                     pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model)
{
  Teuchos::RCP<IntegratorForwardSensitivity<Scalar> > integrator =
    Teuchos::rcp(new IntegratorForwardSensitivity<Scalar>(pList, model));
  return(integrator);
}

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar> >
integratorForwardSensitivity(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model,
  std::string stepperType)
{
  Teuchos::RCP<IntegratorForwardSensitivity<Scalar> > integrator =
    Teuchos::rcp(new IntegratorForwardSensitivity<Scalar>(model, stepperType));
  return(integrator);
}

/// Nonmember constructor
template<class Scalar>
Teuchos::RCP<IntegratorForwardSensitivity<Scalar> >
integratorForwardSensitivity()
{
  Teuchos::RCP<IntegratorForwardSensitivity<Scalar> > integrator =
    Teuchos::rcp(new IntegratorForwardSensitivity<Scalar>());
  return(integrator);
}

} // namespace Tempus
#endif // Tempus_IntegratorForwardSensitivity_impl_hpp
