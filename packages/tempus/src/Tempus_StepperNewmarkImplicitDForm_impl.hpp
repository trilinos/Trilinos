//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_StepperNewmarkImplicitDForm_impl_hpp
#define Tempus_StepperNewmarkImplicitDForm_impl_hpp

#include "Tempus_StepperNewmarkImplicitDFormModifierDefault.hpp"

//#define VERBOSE_DEBUG_OUTPUT
//#define DEBUG_OUTPUT

namespace Tempus {

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::predictVelocity(
    Thyra::VectorBase<Scalar>& vPred, const Thyra::VectorBase<Scalar>& v,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // vPred = v + dt*(1.0-gamma_)*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0, v, dt * (1.0 - gamma_), a);
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::predictDisplacement(
    Thyra::VectorBase<Scalar>& dPred, const Thyra::VectorBase<Scalar>& d,
    const Thyra::VectorBase<Scalar>& v, const Thyra::VectorBase<Scalar>& a,
    const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> tmp =
      Thyra::createMember<Scalar>(dPred.space());
  // dPred = dt*v + dt*dt/2.0*(1.0-2.0*beta_)*a
  Scalar aConst = dt * dt / 2.0 * (1.0 - 2.0 * beta_);
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), dt, v, aConst, a);
  // dPred += d;
  Thyra::Vp_V(Teuchos::ptrFromRef(dPred), d, 1.0);
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::correctVelocity(
    Thyra::VectorBase<Scalar>& v, const Thyra::VectorBase<Scalar>& vPred,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // v = vPred + dt*gamma_*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(v), 1.0, vPred, dt * gamma_, a);
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::correctDisplacement(
    Thyra::VectorBase<Scalar>& d, const Thyra::VectorBase<Scalar>& dPred,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // d = dPred + beta_*dt*dt*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(d), 1.0, dPred, beta_ * dt * dt, a);
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::correctAcceleration(
    Thyra::VectorBase<Scalar>& a, const Thyra::VectorBase<Scalar>& dPred,
    const Thyra::VectorBase<Scalar>& d, const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // a = (d - dPred) / (beta_*dt*dt)
  Scalar const c = 1.0 / beta_ / dt / dt;
  Thyra::V_StVpStV(Teuchos::ptrFromRef(a), c, d, -c, dPred);
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::setBeta(Scalar beta)
{
  if (schemeName_ != "User Defined") {
    *out_ << "\nWARNING: schemeName != 'User Defined' (=" << schemeName_
          << ").\n"
          << " Not setting beta, and leaving as beta = " << beta_ << "!\n";
    return;
  }

  beta_ = beta;

  if (beta_ == 0.0) {
    *out_ << "\nWARNING: Running (implicit implementation of) Newmark "
          << "Implicit a-Form Stepper with Beta = 0.0, which \n"
          << "specifies an explicit scheme.  Mass lumping is not possible, "
          << "so this will be slow!  To run explicit \n"
          << "implementation of Newmark Implicit a-Form Stepper, please "
          << "re-run with 'Stepper Type' = 'Newmark Explicit a-Form'.\n"
          << "This stepper allows for mass lumping when called through "
          << "Piro::TempusSolver.\n";
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      (beta_ > 1.0) || (beta_ < 0.0), std::logic_error,
      "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Beta = "
          << beta_ << ".  Please select Beta >= 0 and <= 1. \n");

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::setGamma(Scalar gamma)
{
  if (schemeName_ != "User Defined") {
    *out_ << "\nWARNING: schemeName != 'User Defined' (=" << schemeName_
          << ").\n"
          << " Not setting gamma, and leaving as gamma = " << gamma_ << "!\n";
    return;
  }

  gamma_ = gamma;

  TEUCHOS_TEST_FOR_EXCEPTION(
      (gamma_ > 1.0) || (gamma_ < 0.0), std::logic_error,
      "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Gamma ="
          << gamma_ << ".  Please select Gamma >= 0 and <= 1. \n");

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::setSchemeName(std::string schemeName)
{
  schemeName_ = schemeName;

  if (schemeName_ == "Average Acceleration") {
    beta_  = 0.25;
    gamma_ = 0.5;
  }
  else if (schemeName_ == "Linear Acceleration") {
    beta_  = 0.25;
    gamma_ = 1.0 / 6.0;
  }
  else if (schemeName_ == "Central Difference") {
    beta_  = 0.0;
    gamma_ = 0.5;
  }
  else if (schemeName_ == "User Defined") {
    beta_  = 0.25;
    gamma_ = 0.5;  // Use defaults until setBeta and setGamma calls.
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "\nError in Tempus::StepperNewmarkImplicitDForm!  "
            << "Invalid Scheme Name = " << schemeName_ << ".  \n"
            << "Valid Scheme Names are: 'Average Acceleration', "
            << "'Linear Acceleration', \n"
            << "'Central Difference' and 'User Defined'.\n");
  }

  this->isInitialized_ = false;
}

template <class Scalar>
StepperNewmarkImplicitDForm<Scalar>::StepperNewmarkImplicitDForm()
  : out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  this->setStepperName("Newmark Implicit d-Form");
  this->setStepperType("Newmark Implicit d-Form");
  this->setUseFSAL(false);
  this->setICConsistency("None");
  this->setICConsistencyCheck(false);
  this->setZeroInitialGuess(false);
  this->setSchemeName("Average Acceleration");

  this->setAppAction(Teuchos::null);
  this->setDefaultSolver();
}

template <class Scalar>
StepperNewmarkImplicitDForm<Scalar>::StepperNewmarkImplicitDForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar>>& solver,
    bool useFSAL, std::string ICConsistency, bool ICConsistencyCheck,
    bool zeroInitialGuess, std::string schemeName, Scalar beta, Scalar gamma,
    const Teuchos::RCP<StepperNewmarkImplicitDFormAppAction<Scalar>>&
        stepperAppAction)
  : out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  this->setStepperName("Newmark Implicit d-Form");
  this->setStepperType("Newmark Implicit d-Form");
  this->setUseFSAL(useFSAL);
  this->setICConsistency(ICConsistency);
  this->setICConsistencyCheck(ICConsistencyCheck);
  this->setZeroInitialGuess(zeroInitialGuess);
  this->setSchemeName(schemeName);
  this->setBeta(beta);
  this->setGamma(gamma);
  this->setAppAction(stepperAppAction);
  this->setSolver(solver);

  if (appModel != Teuchos::null) {
    this->setModel(appModel);
    this->initialize();
  }
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::setModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& appModel)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  validSecondOrderODE_DAE(appModel);
  this->wrapperModel_ =
      Teuchos::rcp(new WrapperModelEvaluatorSecondOrder<Scalar>(
          appModel, "Newmark Implicit d-Form"));

  TEUCHOS_TEST_FOR_EXCEPTION(this->getSolver() == Teuchos::null,
                             std::logic_error, "Error - Solver is not set!\n");
  this->getSolver()->setModel(this->wrapperModel_);

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::setAppAction(
    Teuchos::RCP<StepperNewmarkImplicitDFormAppAction<Scalar>> appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction
    stepperNewmarkImpAppAction_ =
        Teuchos::rcp(new StepperNewmarkImplicitDFormModifierDefault<Scalar>());
  }
  else {
    stepperNewmarkImpAppAction_ = appAction;
  }

  this->isInitialized_ = false;
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar>>& solutionHistory)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperNewmarkImplicitDForm::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        solutionHistory->getNumStates() < 2, std::logic_error,
        "Error - StepperNewmarkImplicitDForm<Scalar>::takeStep(...)\n"
            << "Need at least two SolutionStates for NewmarkImplicitDForm.\n"
            << "  Number of States = " << solutionHistory->getNumStates()
            << "\nTry setting in \"Solution History\" \"Storage Type\" = "
            << "\"Undo\"\n"
            << "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = "
            << "\"2\"\n");

    auto thisStepper = Teuchos::rcpFromRef(*this);
    stepperNewmarkImpAppAction_->execute(
        solutionHistory, thisStepper,
        StepperNewmarkImplicitDFormAppAction<
            Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar>> workingState =
        solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar>> currentState =
        solutionHistory->getCurrentState();

    Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar>> wrapperModel =
        Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar>>(
            this->wrapperModel_);

    // Get values of d, v and a from previous step
    RCP<const Thyra::VectorBase<Scalar>> d_old = currentState->getX();
    RCP<Thyra::VectorBase<Scalar>> v_old       = currentState->getXDot();
    RCP<Thyra::VectorBase<Scalar>> a_old       = currentState->getXDotDot();

    // Get new values of d, v and a from current workingState
    //(to be updated here)
    RCP<Thyra::VectorBase<Scalar>> d_new = workingState->getX();
    RCP<Thyra::VectorBase<Scalar>> v_new = workingState->getXDot();
    RCP<Thyra::VectorBase<Scalar>> a_new = workingState->getXDotDot();

    // Get time and dt
    const Scalar time = currentState->getTime();
    const Scalar dt   = workingState->getTimeStep();
    // Update time
    Scalar t = time + dt;

#ifdef DEBUG_OUTPUT
    Teuchos::Range1D range;

    *out_ << "\n*** d_old ***\n";
    RTOpPack::ConstSubVectorView<Scalar> dov;
    d_old->acquireDetachedView(range, &dov);
    auto doa = dov.values();
    for (auto i = 0; i < doa.size(); ++i) *out_ << doa[i] << " ";
    *out_ << "\n*** d_old ***\n";

    *out_ << "\n*** v_old ***\n";
    RTOpPack::ConstSubVectorView<Scalar> vov;
    v_old->acquireDetachedView(range, &vov);
    auto voa = vov.values();
    for (auto i = 0; i < voa.size(); ++i) *out_ << voa[i] << " ";
    *out_ << "\n*** v_old ***\n";

    *out_ << "\n*** a_old ***\n";
    RTOpPack::ConstSubVectorView<Scalar> aov;
    a_old->acquireDetachedView(range, &aov);
    auto aoa = aov.values();
    for (auto i = 0; i < aoa.size(); ++i) *out_ << aoa[i] << " ";
    *out_ << "\n*** a_old ***\n";
#endif

    // allocate d and v predictors
    RCP<Thyra::VectorBase<Scalar>> d_pred = Thyra::createMember(d_old->space());
    RCP<Thyra::VectorBase<Scalar>> v_pred = Thyra::createMember(v_old->space());

    // compute displacement and velocity predictors
    predictDisplacement(*d_pred, *d_old, *v_old, *a_old, dt);
    predictVelocity(*v_pred, *v_old, *a_old, dt);

#ifdef DEBUG_OUTPUT
    *out_ << "\n*** d_pred ***\n";
    RTOpPack::ConstSubVectorView<Scalar> dpv;
    d_pred->acquireDetachedView(range, &dpv);
    auto dpa = dpv.values();
    for (auto i = 0; i < dpa.size(); ++i) *out_ << dpa[i] << " ";
    *out_ << "\n*** d_pred ***\n";

    *out_ << "\n*** v_pred ***\n";
    RTOpPack::ConstSubVectorView<Scalar> vpv;
    v_pred->acquireDetachedView(range, &vpv);
    auto vpa = vpv.values();
    for (auto i = 0; i < vpa.size(); ++i) *out_ << vpa[i] << " ";
    *out_ << "\n*** v_pred ***\n";

#endif
    // inject d_pred, v_pred, a and other relevant data into wrapperModel
    wrapperModel->initializeNewmark(v_pred, d_pred, dt, t, beta_, gamma_);

    // create initial guess in NOX solver
    RCP<Thyra::VectorBase<Scalar>> initial_guess =
        Thyra::createMember(d_pred->space());
    if ((time == solutionHistory->minTime()) &&
        (this->initialGuess_ != Teuchos::null)) {
      // if first time step and initialGuess_ is provided, set initial_guess =
      // initialGuess_ Throw an exception if initial_guess is not compatible
      // with solution
      bool is_compatible =
          (initial_guess->space())->isCompatible(*this->initialGuess_->space());
      TEUCHOS_TEST_FOR_EXCEPTION(
          is_compatible != true, std::logic_error,
          "Error in Tempus::NemwarkImplicitDForm takeStep(): user-provided "
          "initial guess'!\n"
              << "for Newton is not compatible with solution vector!\n");
      Thyra::copy(*this->initialGuess_, initial_guess.ptr());
    }
    else {
      // Otherwise, set initial guess = diplacement predictor
      Thyra::copy(*d_pred, initial_guess.ptr());
    }

    stepperNewmarkImpAppAction_->execute(
        solutionHistory, thisStepper,
        StepperNewmarkImplicitDFormAppAction<
            Scalar>::ACTION_LOCATION::BEFORE_SOLVE);

    // Set d_pred as initial guess for NOX solver, and solve nonlinear system.
    const Thyra::SolveStatus<Scalar> sStatus =
        (*(this->solver_)).solve(&*initial_guess);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.

    stepperNewmarkImpAppAction_->execute(
        solutionHistory, thisStepper,
        StepperNewmarkImplicitDFormAppAction<
            Scalar>::ACTION_LOCATION::AFTER_SOLVE);

    // solve will return converged solution in initial_guess
    // vector.  Copy it here to d_new, to define the new displacement.
    Thyra::copy(*initial_guess, d_new.ptr());

    // correct acceleration, velocity
    correctAcceleration(*a_new, *d_pred, *d_new, dt);
    correctVelocity(*v_new, *v_pred, *a_new, dt);

#ifdef DEBUG_OUTPUT
    *out_ << "\n*** d_new ***\n";
    RTOpPack::ConstSubVectorView<Scalar> dnv;
    d_new->acquireDetachedView(range, &dnv);
    auto dna = dnv.values();
    for (auto i = 0; i < dna.size(); ++i) *out_ << dna[i] << " ";
    *out_ << "\n*** d_new ***\n";

    *out_ << "\n*** v_new ***\n";
    RTOpPack::ConstSubVectorView<Scalar> vnv;
    v_new->acquireDetachedView(range, &vnv);
    auto vna = vnv.values();
    for (auto i = 0; i < vna.size(); ++i) *out_ << vna[i] << " ";
    *out_ << "\n*** v_new ***\n";

    *out_ << "\n*** a_new ***\n";
    RTOpPack::ConstSubVectorView<Scalar> anv;
    a_new->acquireDetachedView(range, &anv);
    auto ana = anv.values();
    for (auto i = 0; i < ana.size(); ++i) *out_ << ana[i] << " ";
    *out_ << "\n*** a_new ***\n";
#endif

    workingState->setOrder(this->getOrder());
    workingState->computeNorms(currentState);

    stepperNewmarkImpAppAction_->execute(
        solutionHistory, thisStepper,
        StepperNewmarkImplicitDFormAppAction<
            Scalar>::ACTION_LOCATION::END_STEP);
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
Teuchos::RCP<Tempus::StepperState<Scalar>>
StepperNewmarkImplicitDForm<Scalar>::getDefaultStepperState()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Tempus::StepperState<Scalar>> stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  out.setOutputToRootOnly(0);
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);
  StepperImplicit<Scalar>::describe(out, verbLevel);

  out << "--- StepperNewmarkImplicitDForm ---\n";
  out << "  schemeName_ = " << schemeName_ << std::endl;
  out << "  beta_       = " << beta_ << std::endl;
  out << "  gamma_      = " << gamma_ << std::endl;
  out << "-----------------------------------" << std::endl;
}

template <class Scalar>
bool StepperNewmarkImplicitDForm<Scalar>::isValidSetup(
    Teuchos::FancyOStream& out) const
{
  out.setOutputToRootOnly(0);
  bool isValidSetup = true;

  if (!Stepper<Scalar>::isValidSetup(out)) isValidSetup = false;

  // if ( !StepperImplicit<Scalar>::isValidSetup(out) ) isValidSetup = false;
  if (this->wrapperModel_->getAppModel() == Teuchos::null) {
    isValidSetup = false;
    out << "The application ModelEvaluator is not set!\n";
  }

  if (this->wrapperModel_ == Teuchos::null) {
    isValidSetup = false;
    out << "The wrapper ModelEvaluator is not set!\n";
  }

  if (this->solver_ == Teuchos::null) {
    isValidSetup = false;
    out << "The solver is not set!\n";
  }

  return isValidSetup;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperNewmarkImplicitDForm<Scalar>::getValidParameters() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  auto pl = this->getValidParametersBasicImplicit();

  auto newmarkPL = Teuchos::parameterList("Newmark Parameters");
  newmarkPL->set<std::string>("Scheme Name", schemeName_);
  newmarkPL->set<double>("Beta", beta_);
  newmarkPL->set<double>("Gamma", gamma_);
  pl->set("Newmark Parameters", *newmarkPL);

  return pl;
}

// Nonmember constructor - ModelEvaluator and ParameterList
// ------------------------------------------------------------------------
template <class Scalar>
Teuchos::RCP<StepperNewmarkImplicitDForm<Scalar>>
createStepperNewmarkImplicitDForm(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& model,
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto stepper = Teuchos::rcp(new StepperNewmarkImplicitDForm<Scalar>());
  stepper->setStepperImplicitValues(pl);

  if (pl != Teuchos::null) {
    if (pl->isSublist("Newmark Parameters")) {
      auto newmarkPL = pl->sublist("Newmark Parameters", true);
      std::string schemeName =
          newmarkPL.get<std::string>("Scheme Name", "Average Acceleration");
      stepper->setSchemeName(schemeName);
      if (schemeName == "User Defined") {
        stepper->setBeta(newmarkPL.get<double>("Beta", 0.25));
        stepper->setGamma(newmarkPL.get<double>("Gamma", 0.5));
      }
    }
    else {
      stepper->setSchemeName("Average Acceleration");
    }
  }

  if (model != Teuchos::null) {
    stepper->setModel(model);
    stepper->initialize();
  }

  return stepper;
}

}  // namespace Tempus
#endif  // Tempus_StepperNewmarkImplicitDForm_impl_hpp
