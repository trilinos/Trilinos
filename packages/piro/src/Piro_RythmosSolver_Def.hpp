// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Piro_RythmosSolver.hpp"
#include "Piro_ValidPiroParameters.hpp"

#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_RampingIntegrationControlStrategy.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Rythmos_CompositeIntegrationObserver.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"

#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"

#ifdef Piro_ENABLE_NOX
#  include "Thyra_NonlinearSolver_NOX.hpp"
#endif

#include <string>
#include <stdexcept>
#include <iostream>

template <typename Scalar>
Piro::RythmosSolver<Scalar>::RythmosSolver() :
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  isInitialized(false)
{
}

template <typename Scalar>
Piro::RythmosSolver<Scalar>::RythmosSolver(
    Teuchos::RCP<Teuchos::ParameterList> appParams,
    Teuchos::RCP< Thyra::ModelEvaluator<Scalar> > in_model,
    Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer) :
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  isInitialized(false)
{
  initialize(appParams,in_model,observer);
}

template <typename Scalar>
void Piro::RythmosSolver<Scalar>::initialize(
    Teuchos::RCP<Teuchos::ParameterList> appParams,
    Teuchos::RCP< Thyra::ModelEvaluator<Scalar> > in_model,
    Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer)
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // set some internals
  model = in_model;
  num_p = in_model->Np();
  num_g = in_model->Ng();

  //
  *out << "\nA) Get the base parameter list ...\n";
  //

  RCP<Teuchos::ParameterList> rythmosPL = sublist(appParams, "Rythmos", true);
  rythmosPL->validateParameters(*getValidRythmosParameters(),0);

  {
    const std::string verbosity = rythmosPL->get("Verbosity Level", "VERB_DEFAULT");
    solnVerbLevel = Teuchos::VERB_DEFAULT;
    if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
    else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
    else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
    else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
    else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
  }

  t_initial = rythmosPL->get("Initial Time", 0.0);
  t_final = rythmosPL->get("Final Time", 0.1);

  const std::string stepperType = rythmosPL->get("Stepper Type", "Backward Euler");

  //
  *out << "\nC) Create and initalize the forward model ...\n";
  //

  *out << "\nD) Create the stepper and integrator for the forward problem ...\n";
  //

  if (rythmosPL->get<std::string>("Nonlinear Solver Type") == "Rythmos") {
    Teuchos::RCP<Rythmos::TimeStepNonlinearSolver<Scalar> > rythmosTimeStepSolver =
      Rythmos::timeStepNonlinearSolver<Scalar>();
    if (rythmosPL->getEntryPtr("NonLinear Solver")) {
      RCP<Teuchos::ParameterList> nonlinePL =
	sublist(rythmosPL, "NonLinear Solver", true);
      rythmosTimeStepSolver->setParameterList(nonlinePL);
    }
    fwdTimeStepSolver = rythmosTimeStepSolver;
  }
  else if (rythmosPL->get<std::string>("Nonlinear Solver Type") == "NOX") {
#ifdef Piro_ENABLE_NOX
    Teuchos::RCP<Thyra::NOXNonlinearSolver> nox_solver =  Teuchos::rcp(new Thyra::NOXNonlinearSolver);
    Teuchos::RCP<Teuchos::ParameterList> nox_params = Teuchos::rcp(new Teuchos::ParameterList);
    *nox_params = appParams->sublist("NOX");
    nox_solver->setParameterList(nox_params);
    fwdTimeStepSolver = nox_solver;
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Requested NOX solver for a Rythmos Transient solve, Trilinos was not built with NOX enabled.  Please rebuild Trilinos or use the native Rythmos nonlinear solver.");
#endif

  }

  if (stepperType == "Backward Euler") {
    fwdStateStepper = Rythmos::backwardEulerStepper<Scalar> (model, fwdTimeStepSolver);
    fwdStateStepper->setParameterList(sublist(rythmosPL, "Rythmos Stepper", true));
  }
  else if (stepperType == "Explicit RK") {
    fwdStateStepper = Rythmos::explicitRKStepper<Scalar>(model);
    fwdStateStepper->setParameterList(sublist(rythmosPL, "Rythmos Stepper", true));
  }
  else if (stepperType == "BDF") {
    Teuchos::RCP<Teuchos::ParameterList> BDFparams =
      Teuchos::sublist(rythmosPL, "Rythmos Stepper", true);
    Teuchos::RCP<Teuchos::ParameterList> BDFStepControlPL =
      Teuchos::sublist(BDFparams,"Step Control Settings");

    fwdStateStepper = Teuchos::rcp( new Rythmos::ImplicitBDFStepper<Scalar>(model,fwdTimeStepSolver,BDFparams) );
    fwdStateStepper->setInitialCondition(model->getNominalValues());

  }
  else {
    // first (before failing) check to see if the user has added stepper factory
    typename std::map<std::string,Teuchos::RCP<RythmosStepperFactory<Scalar> > >::const_iterator 
        stepFactItr = stepperFactories.find(stepperType);
    if(stepFactItr!=stepperFactories.end()) { 
      // the user has added it, hot dog lets build a new stepper!
      Teuchos::RCP<Teuchos::ParameterList> stepperParams = Teuchos::sublist(rythmosPL, "Rythmos Stepper", true);

      // build the stepper using the factory
      fwdStateStepper = stepFactItr->second->buildStepper(model,fwdTimeStepSolver,stepperParams);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Error! Piro::Epetra::RythmosSolver: Invalid Steper Type: "
          << stepperType << std::endl);
    }
  }

  // Step control strategy
  {
    // If the stepper can accept a step control strategy, then attempt to build one.
    RCP<Rythmos::StepControlStrategyAcceptingStepperBase<Scalar> > scsa_stepper =
      Teuchos::rcp_dynamic_cast<Rythmos::StepControlStrategyAcceptingStepperBase<Scalar> >(fwdStateStepper);

    if (Teuchos::nonnull(scsa_stepper)) {
      const std::string step_control_strategy = rythmosPL->get("Step Control Strategy Type", "None");

      if (step_control_strategy == "None") {
        // don't do anything, stepper will build default
      } else if (step_control_strategy == "ImplicitBDFRamping") {

        const RCP<Rythmos::ImplicitBDFStepperRampingStepControl<Scalar> > rscs =
          rcp(new Rythmos::ImplicitBDFStepperRampingStepControl<Scalar>);

        const RCP<ParameterList> p = parameterList(rythmosPL->sublist("Rythmos Step Control Strategy"));
        rscs->setParameterList(p);

        scsa_stepper->setStepControlStrategy(rscs);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error,
            "Error! Piro::Epetra::RythmosSolver: Invalid step control strategy type: "
            << step_control_strategy << std::endl);
      }
    }
  }

  {
    const RCP<Teuchos::ParameterList> integrationControlPL =
      Teuchos::sublist(rythmosPL, "Rythmos Integration Control", true);

    RCP<Rythmos::DefaultIntegrator<Scalar> > defaultIntegrator;
    if (rythmosPL->get("Rythmos Integration Control Strategy", "Simple") == "Simple") {
      defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(Rythmos::simpleIntegrationControlStrategy<Scalar>(integrationControlPL));
    }
    else if(rythmosPL->get<std::string>("Rythmos Integration Control Strategy") == "Ramping") {
      defaultIntegrator = Rythmos::controlledDefaultIntegrator<Scalar>(Rythmos::rampingIntegrationControlStrategy<Scalar>(integrationControlPL));
    }
    fwdStateIntegrator = defaultIntegrator;
  }

  fwdStateIntegrator->setParameterList(sublist(rythmosPL, "Rythmos Integrator", true));

  if (Teuchos::nonnull(observer)) {
    fwdStateIntegrator->setIntegrationObserver(observer);
  }

  isInitialized = true;
}


template <typename Scalar>
Piro::RythmosSolver<Scalar>::RythmosSolver(
    const Teuchos::RCP<Rythmos::DefaultIntegrator<Scalar> > &stateIntegrator,
    const Teuchos::RCP<Rythmos::StepperBase<Scalar> > &stateStepper,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &underlyingModel,
    Scalar finalTime,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &icModel,
    Teuchos::EVerbosityLevel verbosityLevel) :
  fwdStateIntegrator(stateIntegrator),
  fwdStateStepper(stateStepper),
  fwdTimeStepSolver(timeStepSolver),
  model(underlyingModel),
  initialConditionModel(icModel),
  num_p(model->Np()),
  num_g(model->Ng()),
  t_initial(0.0),
  t_final(finalTime),
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  solnVerbLevel(verbosityLevel),
  isInitialized(true)
{
  if (fwdStateStepper->acceptsModel() && fwdStateStepper->getModel() != underlyingModel) {
    fwdStateStepper->setNonconstModel(underlyingModel);
  }
}

template <typename Scalar>
Piro::RythmosSolver<Scalar>::RythmosSolver(
    const Teuchos::RCP<Rythmos::DefaultIntegrator<Scalar> > &stateIntegrator,
    const Teuchos::RCP<Rythmos::StepperBase<Scalar> > &stateStepper,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &timeStepSolver,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &underlyingModel,
    Scalar initialTime,
    Scalar finalTime,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &icModel,
    Teuchos::EVerbosityLevel verbosityLevel) :
  fwdStateIntegrator(stateIntegrator),
  fwdStateStepper(stateStepper),
  fwdTimeStepSolver(timeStepSolver),
  model(underlyingModel),
  initialConditionModel(icModel),
  num_p(model->Np()),
  num_g(model->Ng()),
  t_initial(initialTime),
  t_final(finalTime),
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  solnVerbLevel(verbosityLevel),
  isInitialized(true)
{
  if (fwdStateStepper->acceptsModel() && fwdStateStepper->getModel() != underlyingModel) {
    fwdStateStepper->setNonconstModel(underlyingModel);
  }
}


template <typename Scalar>
Teuchos::RCP<const Rythmos::IntegratorBase<Scalar> >
Piro::RythmosSolver<Scalar>::getRythmosIntegrator() const
{
  return fwdStateIntegrator;
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::RythmosSolver<Scalar>::get_p_space(int l) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      l >= num_p || l < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<
      "Error in Piro::RythmosSolver::get_p_map():  " <<
      "Invalid parameter index l = " <<
      l << std::endl);

  return model->get_p_space(l);
}

template<typename Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
Piro::RythmosSolver<Scalar>::get_g_space(int j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      j > num_g || j < 0,
      Teuchos::Exceptions::InvalidParameter,
      std::endl <<
      "Error in Piro::RythmosSolver::get_g_map():  " <<
      "Invalid response index j = " <<
      j << std::endl);

  if (j < num_g) {
    return model->get_g_space(j);
  } else {
    // j == num_g
    return model->get_x_space();
  }
}

template<typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::RythmosSolver<Scalar>::getNominalValues() const
{
  Thyra::ModelEvaluatorBase::InArgs<Scalar> result = this->createInArgs();
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> modelNominalValues = model->getNominalValues();
  for (int l = 0; l < num_p; ++l) {
    result.set_p(l, modelNominalValues.get_p(l));
  }
  return result;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar> Piro::RythmosSolver<Scalar>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar> Piro::RythmosSolver<Scalar>::createOutArgsImpl() const
{
  Thyra::ModelEvaluatorBase::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());

  // One additional response slot for the solution vector
  outArgs.set_Np_Ng(num_p, num_g + 1);

  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();

  if (num_p > 0) {
    // Only one parameter supported
    const int l = 0;

    if (Teuchos::nonnull(initialConditionModel)) {
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> initCondOutArgs =
        initialConditionModel->createOutArgs();
      const Thyra::ModelEvaluatorBase::DerivativeSupport init_dxdp_support =
        initCondOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, initCondOutArgs.Ng() - 1, l);
      if (!init_dxdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
        return outArgs;
      }
    }

    // Solution sensitivity
    outArgs.setSupports(
        Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
        num_g,
        l,
        Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);

    if (num_g > 0) {
      // Only one response supported
      const int j = 0;

      const Thyra::ModelEvaluatorBase::DerivativeSupport model_dgdx_support =
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, l);
      if (!model_dgdx_support.none()) {
        const Thyra::ModelEvaluatorBase::DerivativeSupport model_dgdp_support =
          modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
        // Response sensitivity
        Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support;
        if (model_dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
          dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        }
        if (model_dgdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP)) {
          dgdp_support.plus(Thyra::ModelEvaluatorBase::DERIV_LINEAR_OP);
        }
        outArgs.setSupports(
            Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,
            j,
            l,
            dgdp_support);
      }
    }
  }

  return outArgs;
}

template <typename Scalar>
void Piro::RythmosSolver<Scalar>::evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // TODO: Support more than 1 parameter and 1 response
  const int j = 0;
  const int l = 0;

  // Parse InArgs
  RCP<const Thyra::VectorBase<Scalar> > p_in;
  if (num_p > 0) {
    p_in = inArgs.get_p(l);
  }

  // Parse OutArgs
  RCP<Thyra::VectorBase<Scalar> > g_out;
  if (num_g > 0) {
    g_out = outArgs.get_g(j);
  }
  const RCP<Thyra::VectorBase<Scalar> > gx_out = outArgs.get_g(num_g);

  Thyra::ModelEvaluatorBase::InArgs<Scalar> state_ic = model->getNominalValues();

  // Set initial time in ME if needed

  if(t_initial > 0.0 && state_ic.supports(Thyra::ModelEvaluatorBase::IN_ARG_t))

    state_ic.set_t(t_initial);


  if (Teuchos::nonnull(initialConditionModel)) {
    // The initial condition depends on the parameter
    // It is found by querying the auxiliary model evaluator as the last response
    const RCP<Thyra::VectorBase<Scalar> > initialState =
      Thyra::createMember(model->get_x_space());

    {
      Thyra::ModelEvaluatorBase::InArgs<Scalar> initCondInArgs = initialConditionModel->createInArgs();
      if (num_p > 0) {
        initCondInArgs.set_p(l, inArgs.get_p(l));
      }

      Thyra::ModelEvaluatorBase::OutArgs<Scalar> initCondOutArgs = initialConditionModel->createOutArgs();
      initCondOutArgs.set_g(initCondOutArgs.Ng() - 1, initialState);

      initialConditionModel->evalModel(initCondInArgs, initCondOutArgs);
    }

    state_ic.set_x(initialState);
  }

  // Set paramters p_in as part of initial conditions
  if (num_p > 0) {
    state_ic.set_p(l, p_in);
  }

  *out << "\nstate_ic:\n" << Teuchos::describe(state_ic, solnVerbLevel);

  RCP<Thyra::MultiVectorBase<Scalar> > dgxdp_out;
  Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv_out;
  if (num_p > 0) {
    const Thyra::ModelEvaluatorBase::DerivativeSupport dgxdp_support =
      outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, num_g, l);
    if (dgxdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
      const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgxdp_deriv =
        outArgs.get_DgDp(num_g, l);
      dgxdp_out = dgxdp_deriv.getMultiVector();
    }

    if (num_g > 0) {
      const Thyra::ModelEvaluatorBase::DerivativeSupport dgdp_support =
        outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, j, l);
      if (dgxdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
        dgdp_deriv_out = outArgs.get_DgDp(j, l);
      }
    }
  }

  const bool requestedSensitivities =
    Teuchos::nonnull(dgxdp_out) || !dgdp_deriv_out.isEmpty();

  RCP<const Thyra::VectorBase<Scalar> > finalSolution;
  if (!requestedSensitivities) {
    //
    *out << "\nE) Solve the forward problem ...\n";
    //

    fwdStateStepper->setInitialCondition(state_ic);
    fwdStateIntegrator->setStepper(fwdStateStepper, t_final, true);

    Teuchos::Array<RCP<const Thyra::VectorBase<Scalar> > > x_final_array;
    fwdStateIntegrator->getFwdPoints(
        Teuchos::tuple<Scalar>(t_final), &x_final_array, NULL, NULL);
    finalSolution = x_final_array[0];

    if (Teuchos::VERB_MEDIUM <= solnVerbLevel) {
      std::cout << "Final Solution\n" << *finalSolution << std::endl;
    }
  } else { // Computing sensitivities
    //
    *out << "\nE) Solve the forward problem with Sensitivities...\n";
    //
    RCP<Rythmos::ForwardSensitivityStepper<Scalar> > stateAndSensStepper =
      Rythmos::forwardSensitivityStepper<Scalar>();
    stateAndSensStepper->initializeSyncedSteppers(
        model, l, model->getNominalValues(),
        fwdStateStepper, fwdTimeStepSolver);

    //
    // Set the initial condition for the state and forward sensitivities
    //

    const RCP<Thyra::VectorBase<Scalar> > s_bar_init =
      Thyra::createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());
    const RCP<Thyra::VectorBase<Scalar> > s_bar_dot_init =
      Thyra::createMember(stateAndSensStepper->getFwdSensModel()->get_x_space());

    if (Teuchos::is_null(initialConditionModel)) {
      // The initial condition is assumed to be independent from the parameters
      // Therefore, the initial condition for the sensitivity is zero
      Thyra::assign(s_bar_init.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
    } else {
      // Use initialConditionModel to compute initial condition for sensitivity
      Thyra::ModelEvaluatorBase::InArgs<Scalar> initCondInArgs = initialConditionModel->createInArgs();
      initCondInArgs.set_p(l, inArgs.get_p(l));

      Thyra::ModelEvaluatorBase::OutArgs<Scalar> initCondOutArgs = initialConditionModel->createOutArgs();
      typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
      const RCP<DMVPV> s_bar_init_downcasted = Teuchos::rcp_dynamic_cast<DMVPV>(s_bar_init);
      const Thyra::ModelEvaluatorBase::Derivative<Scalar> initCond_deriv(
          s_bar_init_downcasted->getNonconstMultiVector(),
          Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
      initCondOutArgs.set_DgDp(initCondOutArgs.Ng() - 1, l, initCond_deriv);

      initialConditionModel->evalModel(initCondInArgs, initCondOutArgs);
    }
    Thyra::assign(s_bar_dot_init.ptr(), Teuchos::ScalarTraits<Scalar>::zero());

    RCP<const Rythmos::StateAndForwardSensitivityModelEvaluator<Scalar> >
      stateAndSensModel = stateAndSensStepper->getStateAndFwdSensModel();

    Thyra::ModelEvaluatorBase::InArgs<Scalar>
      state_and_sens_ic = stateAndSensStepper->getModel()->createInArgs();

    // Copy time, parameters etc.
    state_and_sens_ic.setArgs(state_ic);
    // Set initial condition for x_bar = [ x; s_bar ]
    state_and_sens_ic.set_x(stateAndSensModel->create_x_bar_vec(state_ic.get_x(), s_bar_init));
    // Set initial condition for x_bar_dot = [ x_dot; s_bar_dot ]
    state_and_sens_ic.set_x_dot(stateAndSensModel->create_x_bar_vec(state_ic.get_x_dot(), s_bar_dot_init));

    stateAndSensStepper->setInitialCondition(state_and_sens_ic);

    //
    // Use a StepperAsModelEvaluator to integrate the state+sens
    //

    const RCP<Rythmos::StepperAsModelEvaluator<Scalar> >
      stateAndSensIntegratorAsModel = Rythmos::stepperAsModelEvaluator(
          Teuchos::rcp_implicit_cast<Rythmos::StepperBase<Scalar> >(stateAndSensStepper),
          Teuchos::rcp_implicit_cast<Rythmos::IntegratorBase<Scalar> >(fwdStateIntegrator),
          state_and_sens_ic);
    // StepperAsModelEvaluator outputs the solution as its last response
    const int stateAndSensModelStateResponseIndex = stateAndSensIntegratorAsModel->Ng() - 1;

    *out << "\nUse the StepperAsModelEvaluator to integrate state + sens x_bar(p,t_final) ... \n";
    Teuchos::OSTab tab(out);

    // Solution sensitivity in column-oriented (Jacobian) MultiVector form
    RCP<const Thyra::MultiVectorBase<Scalar> > dxdp;

    const RCP<Thyra::VectorBase<Scalar> > x_bar_final =
      Thyra::createMember(stateAndSensIntegratorAsModel->get_g_space(stateAndSensModelStateResponseIndex));
    // Extract pieces of x_bar_final to prepare output
    {
      const RCP<const Thyra::ProductVectorBase<Scalar> > x_bar_final_downcasted =
        Thyra::productVectorBase<Scalar>(x_bar_final);

      // Solution
      const int solutionBlockIndex = 0;
      finalSolution = x_bar_final_downcasted->getVectorBlock(solutionBlockIndex);

      // Sensitivity
      const int sensitivityBlockIndex = 1;
      const RCP<const Thyra::VectorBase<Scalar> > s_bar_final =
        x_bar_final_downcasted->getVectorBlock(sensitivityBlockIndex);
      {
        typedef Thyra::DefaultMultiVectorProductVector<Scalar> DMVPV;
        const RCP<const DMVPV> s_bar_final_downcasted = Teuchos::rcp_dynamic_cast<const DMVPV>(s_bar_final);

        dxdp = s_bar_final_downcasted->getMultiVector();
      }
    }

    Thyra::eval_g(
        *stateAndSensIntegratorAsModel,
        l, *state_ic.get_p(l),
        t_final,
        stateAndSensModelStateResponseIndex, x_bar_final.get()
        );

    *out
      << "\nx_bar_final = x_bar(p,t_final) evaluated using "
      << "stateAndSensIntegratorAsModel:\n"
      << Teuchos::describe(*x_bar_final,solnVerbLevel);

    if (Teuchos::nonnull(dgxdp_out)) {
      Thyra::assign(dgxdp_out.ptr(), *dxdp);
    }

    if (!dgdp_deriv_out.isEmpty()) {
      RCP<Thyra::DefaultAddedLinearOp<Scalar> > dgdp_op_out;
      {
        const RCP<Thyra::LinearOpBase<Scalar> > dgdp_op = dgdp_deriv_out.getLinearOp();
        if (Teuchos::nonnull(dgdp_op)) {
          dgdp_op_out = Teuchos::rcp_dynamic_cast<Thyra::DefaultAddedLinearOp<Scalar> >(dgdp_op);
          dgdp_op_out.assert_not_null();
        }
      }

      Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = model->createInArgs();
      {
        modelInArgs.set_x(finalSolution);
        if (num_p > 0) {
          modelInArgs.set_p(l, p_in);
        }
      }

      // require dgdx, dgdp from model
      Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();
      {
        const Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdx_deriv(model->create_DgDx_op(j));
        modelOutArgs.set_DgDx(j, dgdx_deriv);

        Thyra::ModelEvaluatorBase::Derivative<Scalar> dgdp_deriv;
        if (Teuchos::nonnull(dgdp_op_out)) {
          dgdp_deriv = model->create_DgDp_op(j, l);
        } else {
          dgdp_deriv = dgdp_deriv_out;
        }
        modelOutArgs.set_DgDp(j, l, dgdp_deriv);
      }

      model->evalModel(modelInArgs, modelOutArgs);

      const RCP<const Thyra::LinearOpBase<Scalar> > dgdx =
        modelOutArgs.get_DgDx(j).getLinearOp();

      // dgdp_out = dgdp + <dgdx, dxdp>
      if (Teuchos::nonnull(dgdp_op_out)) {
        Teuchos::Array<RCP<const Thyra::LinearOpBase<Scalar> > > op_args(2);
        {
          op_args[0] = modelOutArgs.get_DgDp(j, l).getLinearOp();
          op_args[1] = Thyra::multiply<Scalar>(dgdx, dxdp);
        }
        dgdp_op_out->initialize(op_args);
      } else {
        const RCP<Thyra::MultiVectorBase<Scalar> > dgdp_mv_out = dgdp_deriv_out.getMultiVector();
        Thyra::apply(
            *dgdx,
            Thyra::NOTRANS,
            *dxdp,
            dgdp_mv_out.ptr(),
            Teuchos::ScalarTraits<Scalar>::one(),
            Teuchos::ScalarTraits<Scalar>::one());
      }
    }
  }

  *out << "\nF) Check the solution to the forward problem ...\n";

  // As post-processing step, calculate responses at final solution
  {
    Thyra::ModelEvaluatorBase::InArgs<Scalar> modelInArgs = model->createInArgs();
    {
      modelInArgs.set_x(finalSolution);
      if (num_p > 0) {
        modelInArgs.set_p(l, p_in);
      }
    }

    Thyra::ModelEvaluatorBase::OutArgs<Scalar> modelOutArgs = model->createOutArgs();
    if (Teuchos::nonnull(g_out)) {
      Thyra::put_scalar(Teuchos::ScalarTraits<Scalar>::zero(), g_out.ptr());
      modelOutArgs.set_g(j, g_out);
    }

    model->evalModel(modelInArgs, modelOutArgs);
  }

  // Return the final solution as an additional g-vector, if requested
  if (Teuchos::nonnull(gx_out)) {
    Thyra::copy(*finalSolution, gx_out.ptr());
  }
}

template <typename Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Piro::RythmosSolver<Scalar>::create_DgDp_op_impl(int j, int l) const
{
  TEUCHOS_ASSERT(j != num_g);
  const Teuchos::Array<Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > > dummy =
    Teuchos::tuple(Thyra::zero<Scalar>(this->get_g_space(j), this->get_p_space(l)));
  return Teuchos::rcp(new Thyra::DefaultAddedLinearOp<Scalar>(dummy));
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::RythmosSolver<Scalar>::getValidRythmosParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
     Teuchos::rcp(new Teuchos::ParameterList("ValidRythmosParams"));

  validPL->set<std::string>("Nonlinear Solver Type", "");

  Teuchos::setStringToIntegralParameter<int>(
    "Nonlinear Solver Type",
    "Rythmos",
    "Determines which nonlinear solver to use.",
    Teuchos::tuple<std::string>("Rythmos","NOX"),
    validPL.get()
    );

  validPL->sublist("NonLinear Solver", false, "");
  validPL->sublist("Rythmos Builder", false, "");


  validPL->set<double>("Initial Time", 0.0, "");
  validPL->set<double>("Final Time", 1.0, "");
  validPL->sublist("Rythmos Stepper", false, "");
  validPL->sublist("Rythmos Integrator", false, "");
  validPL->set<std::string>("Rythmos Integration Control Strategy", "Simple", "");
  validPL->set<std::string>("Step Control Strategy Type", "None", "");
  validPL->sublist("Rythmos Step Control Strategy", false, "");
  validPL->sublist("Rythmos Integration Control", false, "");
  validPL->sublist("Stratimikos", false, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<std::string>("Stepper Type", "", "");

  validPL->set<double>("Alpha", 1.0, "");
  validPL->set<double>("Beta", 1.0, "");
  validPL->set<double>("Max State Error", 1.0, "");
  validPL->set<std::string>("Name", "", "");
  validPL->set<bool>("Invert Mass Matrix", false, "");

  return validPL;
}

template <typename Scalar>
void Piro::RythmosSolver<Scalar>::
addStepperFactory(const std::string & stepperName,const Teuchos::RCP<RythmosStepperFactory<Scalar> > & factory)
{
  stepperFactories[stepperName] = factory;
}
