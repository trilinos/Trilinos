// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_IntegratorPseudoTransientAdjointSensitivity_impl_hpp
#define Tempus_IntegratorPseudoTransientAdjointSensitivity_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "NOX_Thyra.H"

namespace Tempus {

template<class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
IntegratorPseudoTransientAdjointSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                inputPL,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model)
{
  model_ = model;
  state_integrator_ = integratorBasic<Scalar>(inputPL, model_);
  sens_model_ = createSensitivityModel(model_, inputPL);
  sens_integrator_ = integratorBasic<Scalar>(inputPL, sens_model_);
}

template<class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
IntegratorPseudoTransientAdjointSensitivity(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  std::string stepperType)
{
  model_ = model;
  state_integrator_ = integratorBasic<Scalar>(model_, stepperType);
  sens_model_ = createSensitivityModel(model, Teuchos::null);
  sens_integrator_ = integratorBasic<Scalar>(sens_model_, stepperType);
}

template<class Scalar>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
IntegratorPseudoTransientAdjointSensitivity()
{
  state_integrator_ = integratorBasic<Scalar>();
  sens_integrator_ = integratorBasic<Scalar>();
}

template<class Scalar>
bool
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
advanceTime()
{
  const Scalar tfinal =
    state_integrator_->getTimeStepControl()->getFinalTime();
  return advanceTime(tfinal);
}

template<class Scalar>
bool
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
advanceTime(const Scalar timeFinal)
{
  using Teuchos::RCP;
  using Thyra::VectorBase;
  typedef Thyra::ModelEvaluatorBase MEB;

  // Run state integrator and get solution
  bool state_status = state_integrator_->advanceTime(timeFinal);

  // Set solution in sensitivity ME
  sens_model_->setForwardSolutionHistory(
    state_integrator_->getSolutionHistory());

  // Run sensitivity integrator
  bool sens_status = sens_integrator_->advanceTime(timeFinal);

  // Compute final dg/dp which is computed by response 0 of the adjoint
  // model evaluator
  MEB::InArgs<Scalar> inargs = sens_model_->getNominalValues();
  MEB::OutArgs<Scalar> outargs = sens_model_->createOutArgs();
  inargs.set_t(sens_integrator_->getTime());
  inargs.set_x(sens_integrator_->getX());
  if (inargs.supports(MEB::IN_ARG_x_dot))
    inargs.set_x_dot(sens_integrator_->getXdot());
  if (inargs.supports(MEB::IN_ARG_x_dot_dot))
    inargs.set_x_dot_dot(sens_integrator_->getXdotdot());
  RCP<VectorBase<Scalar> > G = dgdp_;
  if (G == Teuchos::null) {
    G = Thyra::createMember(sens_model_->get_g_space(0));
    dgdp_ = Teuchos::rcp_dynamic_cast<DMVPV>(G);
  }
  outargs.set_g(0, G);
  sens_model_->evalModel(inargs, outargs);

  buildSolutionHistory();

  return state_status && sens_status;
}

template<class Scalar>
Scalar
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getTime() const
{
  return solutionHistory_->getCurrentTime();
}

template<class Scalar>
Scalar
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getIndex() const
{
  return solutionHistory_->getCurrentIndex();
}

template<class Scalar>
Status
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getStatus() const
{
  Status state_status = state_integrator_->getStatus();
  Status sens_status = sens_integrator_->getStatus();
  if (state_status == FAILED || sens_status == FAILED)
    return FAILED;
  if (state_status == WORKING || sens_status == WORKING)
    return WORKING;
  return PASSED;
}

template<class Scalar>
Teuchos::RCP<Stepper<Scalar> >
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getStepper() const
{
  return state_integrator_->getStepper();
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getTempusParameterList()
{
  return state_integrator_->getTempusParameterList();
}

template<class Scalar>
void
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
setTempusParameterList(Teuchos::RCP<Teuchos::ParameterList> pl)
{
  state_integrator_->setTempusParameterList(pl);
  sens_integrator_->setTempusParameterList(pl);
}

template<class Scalar>
Teuchos::RCP<const SolutionHistory<Scalar> >
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getSolutionHistory() const
{
  return solutionHistory_;
}

template<class Scalar>
Teuchos::RCP<const TimeStepControl<Scalar> >
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getTimeStepControl() const
{
  return state_integrator_->getTimeStepControl();
}

template<class Scalar>
void IntegratorPseudoTransientAdjointSensitivity<Scalar>::
setInitialState(Scalar t0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > x0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdot0,
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > xdotdot0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > y0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > ydot0,
  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > ydotdot0)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorSpaceBase;
  using Thyra::assign;
  using Thyra::createMember;
  typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;

  //
  // Create and initialize product X, Xdot, Xdotdot

  RCP< const VectorSpaceBase<Scalar> > space = sens_model_->get_x_space();
  RCP<DMVPV> Y       = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Ydot    = rcp_dynamic_cast<DMVPV>(createMember(space));
  RCP<DMVPV> Ydotdot = rcp_dynamic_cast<DMVPV>(createMember(space));
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  // y
  if (y0 == Teuchos::null)
    assign(Y->getNonconstMultiVector().ptr(), zero);
  else
    assign(Y->getNonconstMultiVector().ptr(), *y0);

  // ydot
  if (ydot0 == Teuchos::null)
    assign(Ydot->getNonconstMultiVector().ptr(), zero);
  else
    assign(Ydot->getNonconstMultiVector().ptr(), *ydot0);

  // ydotdot
  if (ydotdot0 == Teuchos::null)
    assign(Ydotdot->getNonconstMultiVector().ptr(), zero);
  else
    assign(Ydotdot->getNonconstMultiVector().ptr(), *ydotdot0);

  state_integrator_->setInitialState(t0, x0, xdot0, xdotdot0);
  sens_integrator_->setInitialState(t0, Y, Ydot, Ydotdot);
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getX() const
{
  return state_integrator_->getX();
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getXdot() const
{
  return state_integrator_->getXdot();
}

template<class Scalar>
Teuchos::RCP<const Thyra::VectorBase<Scalar> >
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getXdotdot() const
{
  return state_integrator_->getXdotdot();
}

template<class Scalar>
Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> >
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getDgDp() const
{
  return dgdp_->getMultiVector();
}

template<class Scalar>
std::string
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
description() const
{
  std::string name = "Tempus::IntegratorPseudoTransientAdjointSensitivity";
  return(name);
}

template<class Scalar>
void
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
describe(
  Teuchos::FancyOStream          &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  out << description() << "::describe" << std::endl;
  state_integrator_->describe(out, verbLevel);
  sens_integrator_->describe(out, verbLevel);
}

template<class Scalar>
void
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & inputPL)
{
  state_integrator_->setParameterList(inputPL);
  sens_integrator_->setParameterList(inputPL);
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
unsetParameterList()
{
  state_integrator_->unsetParameterList();
  return sens_integrator_->unsetParameterList();
}

template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::RCP<const Teuchos::ParameterList> integrator_pl =
    state_integrator_->getValidParameters();
  Teuchos::RCP<const Teuchos::ParameterList> sensitivity_pl =
    AdjointSensitivityModelEvaluator<Scalar>::getValidParameters();
  pl->setParameters(*integrator_pl);
  pl->sublist("Sensitivities").setParameters(*sensitivity_pl);

  return pl;
}

template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
getNonconstParameterList()
{
  return state_integrator_->getNonconstParameterList();
}

template <class Scalar>
Teuchos::RCP<Tempus::AdjointSensitivityModelEvaluator<Scalar> >
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
createSensitivityModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
  const Teuchos::RCP<Teuchos::ParameterList>& inputPL)
{
  using Teuchos::rcp;

  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  if (inputPL != Teuchos::null) {
    *pl = inputPL->sublist("Sensitivities");
  }
  const Scalar tfinal = state_integrator_->getTimeStepControl()->getFinalTime();
  return rcp(new AdjointSensitivityModelEvaluator<Scalar>(
               model, tfinal, true, pl));
}

template<class Scalar>
void
IntegratorPseudoTransientAdjointSensitivity<Scalar>::
buildSolutionHistory()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::ParameterList;
  using Thyra::VectorBase;
  using Thyra::VectorSpaceBase;
  using Thyra::assign;
  using Thyra::defaultProductVector;
  typedef Thyra::DefaultProductVector<Scalar> DPV;
  typedef Thyra::DefaultProductVectorSpace<Scalar> DPVS;

  // Create combined solution histories, first for the states with zero
  // adjoint and then for the adjoint with frozen states
  RCP<ParameterList> shPL =
    Teuchos::sublist(state_integrator_->getIntegratorParameterList(),
                     "Solution History", true);
  solutionHistory_ = rcp(new SolutionHistory<Scalar>(shPL));

  RCP<const VectorSpaceBase<Scalar> > x_space =
    model_->get_x_space();
  RCP<const VectorSpaceBase<Scalar> > adjoint_space =
    sens_model_->get_x_space();
  Teuchos::Array< RCP<const VectorSpaceBase<Scalar> > > spaces(2);
  spaces[0] = x_space;
  spaces[1] = adjoint_space;
  RCP<const DPVS > prod_space = Thyra::productVectorSpace(spaces());
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  RCP<const SolutionHistory<Scalar> > state_solution_history =
    state_integrator_->getSolutionHistory();
  int num_states = state_solution_history->getNumStates();
  for (int i=0; i<num_states; ++i) {
    RCP<const SolutionState<Scalar> > state = (*state_solution_history)[i];

    // X
    RCP<DPV> x = defaultProductVector(prod_space);
    assign(x->getNonconstVectorBlock(0).ptr(), *(state->getX()));
    assign(x->getNonconstVectorBlock(1).ptr(), zero);

    // X-Dot
    RCP<DPV> x_dot = defaultProductVector(prod_space);
    assign(x_dot->getNonconstVectorBlock(0).ptr(), *(state->getXDot()));
    assign(x_dot->getNonconstVectorBlock(1).ptr(), zero);

    // X-Dot-Dot
    RCP<DPV> x_dot_dot = defaultProductVector(prod_space);
    assign(x_dot_dot->getNonconstVectorBlock(0).ptr(), *(state->getXDotDot()));
    assign(x_dot_dot->getNonconstVectorBlock(1).ptr(), zero);

    RCP<SolutionState<Scalar> > prod_state =
      rcp(new SolutionState<Scalar>(state->getMetaData()->clone(),
                                    x, x_dot, x_dot_dot,
                                    state->getStepperState()->clone()));
    solutionHistory_->addState(prod_state);
  }

  RCP<const VectorBase<Scalar> > frozen_x =
    state_solution_history->getCurrentState()->getX();
  RCP<const VectorBase<Scalar> > frozen_x_dot =
    state_solution_history->getCurrentState()->getXDot();
  RCP<const VectorBase<Scalar> > frozen_x_dot_dot =
    state_solution_history->getCurrentState()->getXDotDot();
  RCP<const SolutionHistory<Scalar> > sens_solution_history =
    sens_integrator_->getSolutionHistory();
  num_states = sens_solution_history->getNumStates();
  for (int i=0; i<num_states; ++i) {
    RCP<const SolutionState<Scalar> > state = (*sens_solution_history)[i];

    // X
    RCP<DPV> x = defaultProductVector(prod_space);
    assign(x->getNonconstVectorBlock(0).ptr(), *frozen_x);
    assign(x->getNonconstVectorBlock(1).ptr(), *(state->getX()));

    // X-Dot
    RCP<DPV> x_dot = defaultProductVector(prod_space);
    assign(x_dot->getNonconstVectorBlock(0).ptr(), *frozen_x_dot);
    assign(x_dot->getNonconstVectorBlock(1).ptr(), *(state->getXDot()));

    // X-Dot-Dot
    RCP<DPV> x_dot_dot = defaultProductVector(prod_space);
    assign(x_dot_dot->getNonconstVectorBlock(0).ptr(), *frozen_x_dot_dot);
    assign(x_dot_dot->getNonconstVectorBlock(1).ptr(), *(state->getXDotDot()));

    RCP<SolutionState<Scalar> > prod_state =
      rcp(new SolutionState<Scalar>(state->getMetaData()->clone(),
                                    x, x_dot, x_dot_dot,
                                    state->getStepperState()->clone()));
    solutionHistory_->addState(prod_state);
  }
}

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
  Teuchos::RCP<Teuchos::ParameterList>                     pList,
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model)
{
  Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar>(pList, model));
  return(integrator);
}

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >&      model,
  std::string stepperType)
{
  Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar>(model, stepperType));
  return(integrator);
}

/// Non-member constructor
template<class Scalar>
Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> >
integratorPseudoTransientAdjointSensitivity()
{
  Teuchos::RCP<Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar> > integrator =
    Teuchos::rcp(new Tempus::IntegratorPseudoTransientAdjointSensitivity<Scalar>());
  return(integrator);
}

} // namespace Tempus
#endif // Tempus_IntegratorPseudoTransientAdjointSensitivity_impl_hpp
