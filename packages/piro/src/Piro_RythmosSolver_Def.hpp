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

#include "Piro_ObserverToRythmosIntegrationObserverAdapter.hpp"
#include "Piro_ValidPiroParameters.hpp"
#include "Piro_MatrixFreeDecorator.hpp"

#include "Rythmos_BackwardEulerStepper.hpp"
#include "Rythmos_ForwardEulerStepper.hpp"
#include "Rythmos_ExplicitRKStepper.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_ThetaStepper.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Rythmos_RampingIntegrationControlStrategy.hpp"
#include "Rythmos_ForwardSensitivityStepper.hpp"
#include "Rythmos_ImplicitBDFStepperRampingStepControl.hpp"
#include "Rythmos_StepperAsModelEvaluator.hpp"
#include "Rythmos_CompositeIntegrationObserver.hpp"
#include "Rythmos_IntegratorBuilder.hpp"

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Assert.hpp"

#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Piro_InvertMassMatrixDecorator.hpp"

#ifdef HAVE_PIRO_IFPACK2
#include "Thyra_Ifpack2PreconditionerFactory.hpp"
#include "Tpetra_CrsMatrix.hpp"
#endif

#ifdef HAVE_PIRO_MUELU
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Stratimikos_MueLuHelpers.hpp"
#endif

#ifdef HAVE_PIRO_NOX
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
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
    const Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > &observer) :
  out(Teuchos::VerboseObjectBase::getDefaultOStream()),
  isInitialized(false)
{
  std::string jacobianSource = appParams->get("Jacobian Operator", "Have Jacobian");
  if (jacobianSource == "Matrix-Free") {
    Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > mf_model;
    if (appParams->isParameter("Matrix-Free Perturbation")) {
      mf_model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model,
                           appParams->get<double>("Matrix-Free Perturbation")));
    }
    else mf_model = Teuchos::rcp(new Piro::MatrixFreeDecorator<Scalar>(in_model));
    initialize(appParams, mf_model, observer);
  }
  else
    initialize(appParams, in_model, observer);
}

template <typename Scalar>
void Piro::RythmosSolver<Scalar>::initialize(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP< Thyra::ModelEvaluator<Scalar> > &in_model,
    const Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > &observer)
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


  if (appParams->isSublist("Rythmos")) {
  RCP<Teuchos::ParameterList> rythmosPL = sublist(appParams, "Rythmos", true);
  rythmosPL->validateParameters(*getValidRythmosParameters(),0);

  {
    const std::string verbosity = rythmosPL->get("Verbosity Level", "VERB_DEFAULT");
    if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
    else if (verbosity == "VERB_DEFAULT") solnVerbLevel = Teuchos::VERB_DEFAULT;
    else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
    else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
    else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
    else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
    else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Unknown verbosity option specified in Piro_RythmosSolver.");
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
#ifdef HAVE_PIRO_NOX
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
  else if (stepperType == "Forward Euler") {
    fwdStateStepper = Rythmos::forwardEulerStepper<Scalar> (model);
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

  } else if (stepperType == "Theta Stepper") {

    Teuchos::RCP<Teuchos::ParameterList> CrankNicholsonPL =
      Teuchos::sublist(rythmosPL, "Rythmos Stepper", true);
    fwdStateStepper = Rythmos::thetaStepper<Scalar>(model, fwdTimeStepSolver, CrankNicholsonPL);
  }
  else {
    // first (before failing) check to see if the user has added stepper factory
    typename std::map<std::string,Teuchos::RCP<Piro::RythmosStepperFactory<Scalar> > >::const_iterator
        stepFactItr = stepperFactories.find(stepperType);
    if(stepFactItr!=stepperFactories.end()) {
      // the user has added it, hot dog lets build a new stepper!
      Teuchos::RCP<Teuchos::ParameterList> stepperParams = Teuchos::sublist(rythmosPL, "Rythmos Stepper", true);

      // build the stepper using the factory
      fwdStateStepper = stepFactItr->second->buildStepper(model,fwdTimeStepSolver,stepperParams);

      // the user decided to override the model being used (let them)
      if(fwdStateStepper->getModel()!=model && fwdStateStepper->getModel()!=Teuchos::null) {
        model = Teuchos::rcp_const_cast<Thyra::ModelEvaluator<Scalar> >(fwdStateStepper->getModel());

        num_p = in_model->Np();
        num_g = in_model->Ng();
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, Teuchos::Exceptions::InvalidParameter,
          std::endl << "Error! Piro::RythmosSolver: Invalid Steper Type: "
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
      }
    else {
       // first (before failing) check to see if the user has added step control factory
       typename std::map<std::string,Teuchos::RCP<Piro::RythmosStepControlFactory<Scalar> > >::const_iterator
           stepControlFactItr = stepControlFactories.find(step_control_strategy);
      if (stepControlFactItr != stepControlFactories.end())
         {

        const RCP<Rythmos::StepControlStrategyBase<Scalar> > rscs = stepControlFactItr->second->buildStepControl();

        const RCP<ParameterList> p = parameterList(rythmosPL -> sublist("Rythmos Step Control Strategy"));

        rscs->setParameterList(p);

        scsa_stepper->setStepControlStrategy(rscs);
        }
        else {
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error,
            "Error! Piro::RythmosSolver: Invalid step control strategy type: "
            << step_control_strategy << std::endl);
      }
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
 }

  else if (appParams->isSublist("Rythmos Solver")) {
    /** New parameter list format **/
     RCP<Teuchos::ParameterList> rythmosSolverPL = sublist(appParams, "Rythmos Solver", true);
    RCP<Teuchos::ParameterList> rythmosPL = sublist(rythmosSolverPL, "Rythmos", true);

    {
      const std::string verbosity = rythmosSolverPL->get("Verbosity Level", "VERB_DEFAULT");
      if      (verbosity == "VERB_NONE")    solnVerbLevel = Teuchos::VERB_NONE;
      else if (verbosity == "VERB_DEFAULT") solnVerbLevel = Teuchos::VERB_DEFAULT;
      else if (verbosity == "VERB_LOW")     solnVerbLevel = Teuchos::VERB_LOW;
      else if (verbosity == "VERB_MEDIUM")  solnVerbLevel = Teuchos::VERB_MEDIUM;
      else if (verbosity == "VERB_HIGH")    solnVerbLevel = Teuchos::VERB_HIGH;
      else if (verbosity == "VERB_EXTREME") solnVerbLevel = Teuchos::VERB_EXTREME;
      else TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
         "Unknown verbosity option specified in Piro_RythmosSolver.");
    }

    t_initial = rythmosPL->sublist("Integrator Settings").get("Initial Time", 0.0);
    t_final = rythmosPL->sublist("Integrator Settings").get("Final Time", 0.1);

    const std::string stepperType = rythmosPL->sublist("Stepper Settings")
      .sublist("Stepper Selection").get("Stepper Type", "Backward Euler");
     //
     //    *out << "\nB) Create the Stratimikos linear solver factory ...\n";
     //
     // This is the linear solve strategy that will be used to solve for the
     // linear system with the W.
     //
     Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

#ifdef HAVE_PIRO_IFPACK2
     typedef Thyra::PreconditionerFactoryBase<double> Base;
     typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<double> > Impl;
     linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
#ifdef HAVE_PIRO_MUELU
     Stratimikos::enableMueLu(linearSolverBuilder);
#endif

     linearSolverBuilder.setParameterList(sublist(rythmosSolverPL, "Stratimikos", true));
     rythmosSolverPL->validateParameters(*getValidRythmosSolverParameters(),0);
     RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
                                  createLinearSolveStrategy(linearSolverBuilder);
     //
     *out << "\nC) Create and initalize the forward model ...\n";
     //
     // C.1) Create the underlying Thyra::ModelEvaluator
     // already constructed as "model". Decorate if needed.
     if (
      stepperType == "Explicit RK" ||
      stepperType == "Forward Euler" ||
      stepperType == "Explicit Taylor Polynomial") {

      Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > origModel = model;
      model = Teuchos::rcp(new Piro::InvertMassMatrixDecorator<Scalar>(
              sublist(rythmosSolverPL,"Stratimikos", true), origModel,
              rythmosSolverPL->get("Constant Mass Matrix", false),rythmosSolverPL->get("Lump Mass Matrix", false),false));
     }
     // C.2) Create the Thyra-wrapped ModelEvaluator

      thyraModel = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<Scalar>(model, lowsFactory));

    const RCP<const Thyra::VectorSpaceBase<double> > x_space =
      thyraModel->get_x_space();

    //
    *out << "\nD) Create the stepper and integrator for the forward problem ...\n";
    //
    fwdTimeStepSolver = Rythmos::timeStepNonlinearSolver<double>();

    if (rythmosSolverPL->getEntryPtr("NonLinear Solver")) {
      const RCP<Teuchos::ParameterList> nonlinePL =
                                 sublist(rythmosSolverPL, "NonLinear Solver", true);
      fwdTimeStepSolver->setParameterList(nonlinePL);
    }
    // Force Default Integrator since this is needed for Observers
        rythmosPL->sublist("Integrator Settings").sublist("Integrator Selection").
      set("Integrator Type","Default Integrator");

    RCP<Rythmos::IntegratorBuilder<double> > ib = Rythmos::integratorBuilder<double>();
    ib->setParameterList(rythmosPL);
    Thyra::ModelEvaluatorBase::InArgs<double> ic = thyraModel->getNominalValues();
    RCP<Rythmos::IntegratorBase<double> > integrator = ib->create(thyraModel,ic,fwdTimeStepSolver);
    fwdStateIntegrator = Teuchos::rcp_dynamic_cast<Rythmos::DefaultIntegrator<double> >(integrator,true);

    fwdStateStepper = fwdStateIntegrator->getNonconstStepper();

    if (Teuchos::nonnull(observer))
      fwdStateIntegrator->setIntegrationObserver(observer);

  }
else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        appParams->isSublist("Rythmos") || appParams->isSublist("Rythmos Solver"),
        Teuchos::Exceptions::InvalidParameter, std::endl <<
        "Error! Piro::RythmosSolver: must have either Rythmos or Rythmos Solver sublist ");

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
  t_initial(0.0),
  t_final(finalTime),
  num_p(model->Np()),
  num_g(model->Ng()),
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
  t_initial(initialTime),
  t_final(finalTime),
  num_p(model->Np()),
  num_g(model->Ng()),
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


template <typename Scalar>
Teuchos::RCP<const Thyra::NonlinearSolverBase<Scalar> >
Piro::RythmosSolver<Scalar>::getTimeStepSolver() const
{
  return fwdTimeStepSolver;
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
Thyra::ModelEvaluatorBase::InArgs<Scalar>
Piro::RythmosSolver<Scalar>::createInArgs() const
{
  Thyra::ModelEvaluatorBase::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(num_p);
  return inArgs;
}

template <typename Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
Piro::RythmosSolver<Scalar>::createOutArgsImpl() const
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
        // Ok to return early since only one parameter supported
        return outArgs;
      }
    }

    // Computing the DxDp sensitivity for a transient problem currently requires the evaluation of
    // the mutilivector-based, Jacobian-oriented DfDp derivatives of the underlying transient model.
    const Thyra::ModelEvaluatorBase::DerivativeSupport model_dfdp_support =
      modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DfDp, l);
    if (!model_dfdp_support.supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM)) {
      // Ok to return early since only one parameter supported
      return outArgs;
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
        modelOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDx, j);
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
  RCP<const Thyra::VectorBase<Scalar> > p_in2;  //JF add for multipoint
  if (num_p > 1) {
    p_in2 = inArgs.get_p(l+1);
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
    if (Teuchos::nonnull(p_in)) {
      state_ic.set_p(l, p_in);
    }
  }
  if (num_p > 1) { //JF added for multipoint
    if (Teuchos::nonnull(p_in2)) {
      state_ic.set_p(l+1, p_in2);
    }
  }

  *out << "\nstate_ic:\n" << Teuchos::describe(state_ic, solnVerbLevel);

  //JF  may need a version of the following for multipoint, i.e. num_p>1, l+1, if we want sensitivities
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
      if (!dgdp_support.none()) {
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
    *out << "T final : " << t_final << " \n";

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
      if (num_p > 1) {  //JF added for multipoint
        modelInArgs.set_p(l+1, p_in2);
      }
      //Set time to be final time at which the solve occurs (< t_final in the case we don't make it to t_final).
      modelInArgs.set_t(fwdStateStepper->getTimeRange().lower());
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

  return validPL;
}

template <typename Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
Piro::RythmosSolver<Scalar>::getValidRythmosSolverParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> validPL =
    Teuchos::rcp(new Teuchos::ParameterList("ValidRythmosSolverParams"));;
  validPL->sublist("Rythmos", false, "");
  validPL->sublist("Stratimikos", false, "");
  validPL->sublist("NonLinear Solver", false, "");
  validPL->set<std::string>("Verbosity Level", "", "");
  validPL->set<bool>("Lump Mass Matrix", false, "Boolean to tell code to lump mass matrix");
  validPL->set<bool>("Constant Mass Matrix", false, "Boolean to tell code if mass matrix is constant in time");
  return validPL;
}

template <typename Scalar>
void Piro::RythmosSolver<Scalar>::
addStepperFactory(const std::string & stepperName,const Teuchos::RCP<Piro::RythmosStepperFactory<Scalar> > & factory)
{
  stepperFactories[stepperName] = factory;
}

template <typename Scalar>
void Piro::RythmosSolver<Scalar>::
addStepControlFactory(const std::string & stepControlName,
                      const Teuchos::RCP<Piro::RythmosStepControlFactory<Scalar>> & step_control_strategy)
{
  stepControlFactories[stepControlName] = step_control_strategy;
}

template <typename Scalar>
Teuchos::RCP<Piro::RythmosSolver<Scalar> >
Piro::rythmosSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &appParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &in_model,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &piroObserver)
{
  Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer;
  if (Teuchos::nonnull(piroObserver)) {
    observer = Teuchos::rcp(
        new ObserverToRythmosIntegrationObserverAdapter<Scalar>(piroObserver));
  }

  return Teuchos::rcp(new RythmosSolver<Scalar>(appParams, in_model, observer));

}

