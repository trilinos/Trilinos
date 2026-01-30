//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperEPI3_impl_hpp
#define Tempus_StepperEPI3_impl_hpp

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_StepperEPI3ModifierDefault.hpp"

#include "Tempus_PhiEvaluator.hpp"
#include "Tempus_PhiEvaluatorFactory.hpp"

namespace Tempus {

template <class Scalar>
StepperEPI3<Scalar>::StepperEPI3()
{
  this->setStepperName("EPI3");
  this->setStepperType("EPI3");
  this->setUseFSAL(true);
  this->setTaylorExpansionOrder(2);
  this->setICConsistency("Consistent");
  this->setICConsistencyCheck(false);
  this->setAppAction(Teuchos::null);
}

template <class Scalar>
StepperEPI3<Scalar>::StepperEPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
    const Teuchos::RCP<StepperEPI3AppAction<Scalar> >&
        stepperEPI3AppAction)
{
  this->setStepperName("EPI3");
  this->setStepperType("EPI3");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);

  this->setAppAction(stepperEPI3AppAction);
  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template <class Scalar>
void StepperEPI3<Scalar>::setAppAction(
    Teuchos::RCP<StepperEPI3AppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperEPI3AppAction_ =
        Teuchos::rcp(new StepperEPI3ModifierDefault<Scalar>());
  }
  else {
    stepperEPI3AppAction_ = appAction;
  }
}

template<class Scalar>
void StepperEPI3<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  StepperExplicit<Scalar>::setModel(appModel);

  auto phif = Teuchos::rcp(new PhiEvaluatorFactory<Scalar>());
  phiEvaluator_ = phif->createPhiEvaluator("Taylor", appModel);

  TEUCHOS_TEST_FOR_EXCEPTION(
  phiEvaluator_.is_null(), std::logic_error,
  "phiEvaluator_ is null in StepperEPI3::setModel");

  phiEvaluator_->setModel(appModel);
  phiEvaluator_->initialize();

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperEPI3<Scalar>::setInitialConditions(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();

  // Check if we need Stepper storage for xDot
  if (initialState->getXDot() == Teuchos::null)
    this->setStepperXDot(initialState->getX()->clone_v());
  else
    this->setStepperXDot(initialState->getXDot());

  StepperExplicit<Scalar>::setInitialConditions(solutionHistory);
}

template <class Scalar>
void StepperEPI3<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();
  phiEvaluator_->checkInitialized();
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperEPI3::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        solutionHistory->getNumStates() < 2, std::logic_error,
        "Error - StepperEPI3<Scalar>::takeStep(...)\n"
            << "Need at least two SolutionStates for EPI3.\n"
            << "  Number of States = " << solutionHistory->getNumStates()
            << "\n Try setting in \"Solution History\" \"Storage Type\" = "
            << "\"Undo\"\n or \"Storage Type\" = \"Static\" and \"Storage Limit\" = "
            << "\"2\"\n");

    RCP<StepperEPI3<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperEPI3AppAction_->execute(
        solutionHistory, thisStepper,
        StepperEPI3AppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > currentState =
        solutionHistory->getCurrentState();
    RCP<SolutionState<Scalar> > workingState =
        solutionHistory->getWorkingState();
    if (currentState->getXDot() != Teuchos::null)
      this->setStepperXDot(currentState->getXDot());
    RCP<Thyra::VectorBase<Scalar> > xDot = this->getStepperXDot();
    const Scalar dt                      = workingState->getTimeStep();

    // if (!(this->getUseFSAL()) || workingState->getNConsecutiveFailures() != 0) {
      // Need to compute XDotOld.
      stepperEPI3AppAction_->execute(
          solutionHistory, thisStepper,
          StepperEPI3AppAction<
              Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);

      auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

      // Evaluate xDot = f(x,t).
      this->evaluateExplicitODE(xDot, currentState->getX(),
                                currentState->getTime(), p);

      RCP<const Thyra::ModelEvaluator<Scalar>> appModel = this->getModel();
      Thyra::ModelEvaluatorBase::InArgs<Scalar> inArgs = appModel->createInArgs();

      // Model evaluator builds: alpha*u_dot + beta*F(u) = 0
      const Scalar time = workingState->getTime();
      inArgs.set_x(currentState->getX());
      inArgs.set_t(time);
      inArgs.set_x_dot(xDot);
      RCP<const Thyra::VectorBase<Scalar> > xDot_Exponential = phiEvaluator_->buildATildeMatrix(2, this->getTaylorExpansionOrder(), dt, xDot, inArgs);

      // For UseFSAL=false, x and xDot are now sync'ed or consistent
      // at the same time level for the currentState.
      currentState->setIsSynced(true);
    // }

    // EPI3 update, x^n = x^{n-1} + dt^n * xDot^{n-1}
    // Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
    //                *(currentState->getX()), dt, *(xDot));
    Thyra::V_VpStV(Teuchos::outArg(*(workingState->getX())),
                   *(currentState->getX()), dt, *(xDot_Exponential));

    if (workingState->getXDot() != Teuchos::null)
      this->setStepperXDot(workingState->getXDot());
    xDot = this->getStepperXDot();

    if (this->getUseFSAL()) {
      // Get consistent xDot^n.
      stepperEPI3AppAction_->execute(
          solutionHistory, thisStepper,
          StepperEPI3AppAction<
              Scalar>::ACTION_LOCATION::BEFORE_EXPLICIT_EVAL);

      auto p = Teuchos::rcp(new ExplicitODEParameters<Scalar>(dt));

      // Evaluate xDot = f(x,t).
      this->evaluateExplicitODE(xDot, workingState->getX(),
                                workingState->getTime(), p);

      // For UseFSAL=true, x and xDot are now sync'ed or consistent
      // for the workingState.
      workingState->setIsSynced(true);
    }
    else {
      assign(xDot.ptr(), Teuchos::ScalarTraits<Scalar>::zero());
      workingState->setIsSynced(false);
    }

    workingState->setSolutionStatus(Status::PASSED);
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);
    stepperEPI3AppAction_->execute(
        solutionHistory, thisStepper,
        StepperEPI3AppAction<Scalar>::ACTION_LOCATION::END_STEP);
  }
  return;
}

/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template <class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperEPI3<Scalar>::getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperEPI3<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << std::endl;
  Stepper<Scalar>::describe(*l_out, verbLevel);
  StepperExplicit<Scalar>::describe(*l_out, verbLevel);
  *l_out << "  stepperEPI3AppAction_ = " << stepperEPI3AppAction_ << std::endl
         << "----------------------------" << std::endl;
}

template <class Scalar>
bool StepperEPI3<Scalar>::isValidSetup(Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);

  bool isValidSetup = true;

  if (!Stepper<Scalar>::isValidSetup(out)) isValidSetup = false;
  if (!StepperExplicit<Scalar>::isValidSetup(out)) isValidSetup = false;
  if (stepperEPI3AppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The EPI3 AppAction is not set!\n";
  }
  return isValidSetup;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList> StepperEPI3<Scalar>::getValidParameters()
    const
{
  return this->getValidParametersBasic();
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList> StepperEPI3<Scalar>::getValidParametersBasic()
    const
{
  auto pl = Teuchos::parameterList(this->getStepperName());
  pl->template set<std::string>("Stepper Type", this->getStepperType());

  pl->template set<bool>(
      "Use FSAL", this->getUseFSAL(),
      "The First-Same-As-Last (FSAL) principle is the situation where the\n"
      "last function evaluation, f(x^{n-1},t^{n-1}) [a.k.a. xDot^{n-1}],\n"
      "can be used for the first function evaluation, f(x^n,t^n)\n"
      "[a.k.a. xDot^n].  For RK methods, this applies to the stages.\n"
      "\n"
      "Often the FSAL priniciple can be used to save an evaluation.\n"
      "However there are cases when it cannot be used, e.g., operator\n"
      "splitting where other steppers/operators have modified the solution,\n"
      "x^*, and thus require the function evaluation, f(x^*, t^{n-1}).\n"
      "\n"
      "It should be noted that when the FSAL priniciple can be used\n"
      "(can set useFSAL=true), setting useFSAL=false will give the\n"
      "same solution but at additional expense.  However, the reverse\n"
      "is not true.  When the FSAL priniciple can not be used\n"
      "(need to set useFSAL=false), setting useFSAL=true will produce\n"
      "incorrect solutions.\n"
      "\n"
      "Default in general for explicit and implicit steppers is false,\n"
      "but individual steppers can override this default.");

  pl->template set<std::string>(
      "Initial Condition Consistency", this->getICConsistency(),
      "This indicates which type of consistency should be applied to\n"
      "the initial conditions (ICs):\n"
      "\n"
      "  'None'   - Do nothing to the ICs provided in the SolutionHistory.\n"
      "  'Zero'   - Set the derivative of the SolutionState to zero in the\n"
      "             SolutionHistory provided, e.g., xDot^0 = 0, or \n"
      "             xDotDot^0 = 0.\n"
      "  'App'    - Use the application's ICs, e.g., getNominalValues().\n"
      "  'Consistent' - Make the initial conditions for x and xDot\n"
      "             consistent with the governing equations, e.g.,\n"
      "             xDot = f(x,t), and f(x, xDot, t) = 0.  For implicit\n"
      "             ODEs, this requires a solve of f(x, xDot, t) = 0 for\n"
      "             xDot, and another Jacobian and residual may be\n"
      "             needed, e.g., boundary conditions on xDot may need\n"
      "             to replace boundary conditions on x.\n"
      "\n"
      "In general for explicit steppers, the default is 'Consistent',\n"
      "because it is fairly cheap with just one residual evaluation.\n"
      "In general for implicit steppers, the default is 'None', because\n"
      "the application often knows its IC and can set it the initial\n"
      "SolutionState.  Also, as noted above, 'Consistent' may require\n"
      "another Jacobian from the application.  Individual steppers may\n"
      "override these defaults.");

  pl->template set<bool>(
      "Initial Condition Consistency Check", this->getICConsistencyCheck(),
      "Check if the initial condition, x and xDot, is consistent with the\n"
      "governing equations, xDot = f(x,t), or f(x, xDot, t) = 0.\n"
      "\n"
      "In general for explicit and implicit steppers, the default is true,\n"
      "because it is fairly cheap with just one residual evaluation.\n"
      "Individual steppers may override this default.");

  pl->template set<int>(
      "Taylor Expansion Order", 2,
      "The order of the Taylor expansion used in the EPI3 stepper.\n"
      "\n"
      "The default is 2.");

  return pl;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperEPI3<Scalar> > createStepperEPI3(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperEPI3<Scalar>());
  stepper->setStepperExplicitValues(pl);

  if (pl != Teuchos::null) {
    stepper->setTaylorExpansionOrder(pl->template get<int>("Taylor Expansion Order", 2));
  }

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

}  // namespace Tempus
#endif  // Tempus_StepperEPI3_impl_hpp
