// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkImplicitDForm_impl_hpp
#define Tempus_StepperNewmarkImplicitDForm_impl_hpp

#include "NOX_Thyra.H"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_config.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

//#define VERBOSE_DEBUG_OUTPUT
//#define DEBUG_OUTPUT

namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template <class Scalar>
class StepperFactory;

template <class Scalar>
void
StepperNewmarkImplicitDForm<Scalar>::predictVelocity(
    Thyra::VectorBase<Scalar>& vPred, const Thyra::VectorBase<Scalar>& v,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // vPred = v + dt*(1.0-gamma_)*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0, v, dt * (1.0 - gamma_), a);
}

template <class Scalar>
void
StepperNewmarkImplicitDForm<Scalar>::predictDisplacement(
    Thyra::VectorBase<Scalar>& dPred, const Thyra::VectorBase<Scalar>& d,
    const Thyra::VectorBase<Scalar>& v, const Thyra::VectorBase<Scalar>& a,
    const Scalar dt) const {
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
void
StepperNewmarkImplicitDForm<Scalar>::correctVelocity(
    Thyra::VectorBase<Scalar>& v, const Thyra::VectorBase<Scalar>& vPred,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // v = vPred + dt*gamma_*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(v), 1.0, vPred, dt * gamma_, a);
}

template <class Scalar>
void
StepperNewmarkImplicitDForm<Scalar>::correctDisplacement(
    Thyra::VectorBase<Scalar>& d, const Thyra::VectorBase<Scalar>& dPred,
    const Thyra::VectorBase<Scalar>& a, const Scalar dt) const {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // d = dPred + beta_*dt*dt*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(d), 1.0, dPred, beta_ * dt * dt, a);
}

template <class Scalar>
void
StepperNewmarkImplicitDForm<Scalar>::correctAcceleration(
    Thyra::VectorBase<Scalar>& a, const Thyra::VectorBase<Scalar>& dPred,
    const Thyra::VectorBase<Scalar>& d, const Scalar dt) const {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  // a = (d - dPred) / (beta_*dt*dt)
  Scalar const c = 1.0 / beta_ / dt / dt;
  Thyra::V_StVpStV(Teuchos::ptrFromRef(a), c, d, -c, dPred);
}

template<class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::setBeta(Scalar beta)
{
  if (schemeName_ != "User Defined") {
    *out_ << "\nWARNING: schemeName != 'User Defined' (=" <<schemeName_<< ").\n"
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

  TEUCHOS_TEST_FOR_EXCEPTION( (beta_ > 1.0) || (beta_ < 0.0),
    std::logic_error,
    "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Beta = "
    << beta_ << ".  Please select Beta >= 0 and <= 1. \n");
}


template<class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::setGamma(Scalar gamma)
{
  if (schemeName_ != "User Defined") {
    *out_ << "\nWARNING: schemeName != 'User Defined' (=" <<schemeName_<< ").\n"
          << " Not setting gamma, and leaving as gamma = " << gamma_ << "!\n";
    return;
  }

  gamma_ = gamma;

  TEUCHOS_TEST_FOR_EXCEPTION( (gamma_ > 1.0) || (gamma_ < 0.0),
    std::logic_error,
    "\nError in 'Newmark Implicit a-Form' stepper: invalid value of Gamma ="
    <<gamma_ << ".  Please select Gamma >= 0 and <= 1. \n");
}


template<class Scalar>
void StepperNewmarkImplicitDForm<Scalar>::setSchemeName(
  std::string schemeName)
{
  schemeName_ = schemeName;

  if (schemeName_ == "Average Acceleration") {
    beta_= 0.25; gamma_ = 0.5;
  }
  else if (schemeName_ == "Linear Acceleration") {
    beta_= 0.25; gamma_ = 1.0/6.0;
  }
  else if (schemeName_ == "Central Difference") {
    beta_=  0.0; gamma_ = 0.5;
  }
  else if (schemeName_ == "User Defined") {
    beta_= 0.25; gamma_ = 0.5; // Use defaults until setBeta and setGamma calls.
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,
       std::logic_error,
       "\nError in Tempus::StepperNewmarkImplicitDForm!  "
       <<"Invalid Scheme Name = " << schemeName_ <<".  \n"
       <<"Valid Scheme Names are: 'Average Acceleration', "
       <<"'Linear Acceleration', \n"
       <<"'Central Difference' and 'User Defined'.\n");
  }
}


template <class Scalar>
StepperNewmarkImplicitDForm<Scalar>::StepperNewmarkImplicitDForm()
    : out_(Teuchos::VerboseObjectBase::getDefaultOStream()) {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  this->setStepperType(        "Newmark Implicit d-Form");
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());
  this->setZeroInitialGuess(   false);
  this->setSchemeName(         "Average Acceleration");

  this->setObserver();
}

template <class Scalar>
StepperNewmarkImplicitDForm<Scalar>::StepperNewmarkImplicitDForm(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& appModel,
  const Teuchos::RCP<StepperObserver<Scalar> >& obs,
  const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  bool zeroInitialGuess,
  std::string schemeName,
  Scalar beta,
  Scalar gamma)
  : out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
  this->setStepperType(        "Newmark Implicit d-Form");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);
  this->setZeroInitialGuess(   zeroInitialGuess);
  this->setSchemeName(         schemeName);
  this->setBeta(               beta);
  this->setGamma(              gamma);

  this->setObserver(obs);

  if (appModel != Teuchos::null) {

    this->setModel(appModel);
    this->setSolver(solver);
    this->initialize();
  }
}


template <class Scalar>
void
StepperNewmarkImplicitDForm<Scalar>::setModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar>>& appModel) {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  validSecondOrderODE_DAE(appModel);
  Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
    Teuchos::rcp(new WrapperModelEvaluatorSecondOrder<Scalar>(
      appModel, "Newmark Implicit d-Form"));
  this->wrapperModel_ = wrapperModel;
  inArgs_ = this->wrapperModel_->getNominalValues();
  outArgs_ = this->wrapperModel_->createOutArgs();
}

template <class Scalar>
void
StepperNewmarkImplicitDForm<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperNewmarkImplicitDForm::initialize()\n");

#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
}

template <class Scalar>
void
StepperNewmarkImplicitDForm<Scalar>::takeStep(
    const Teuchos::RCP<SolutionHistory<Scalar>>& solutionHistory) {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperNewmarkImplicitDForm::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperNewmarkImplicitDForm<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for NewmarkImplicitDForm.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<SolutionState<Scalar>> workingState =solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar>> currentState =solutionHistory->getCurrentState();

    Teuchos::RCP<WrapperModelEvaluatorSecondOrder<Scalar> > wrapperModel =
      Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
        this->wrapperModel_);

    // Get values of d, v and a from previous step
    RCP<const Thyra::VectorBase<Scalar>> d_old = currentState->getX();
    RCP<Thyra::VectorBase<Scalar>> v_old = currentState->getXDot();
    RCP<Thyra::VectorBase<Scalar>> a_old = currentState->getXDotDot();

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
    RCP<Thyra::VectorBase<Scalar>> initial_guess = Thyra::createMember(d_pred->space());
    if ((time == solutionHistory->minTime()) && (this->initial_guess_ != Teuchos::null)) {
      //if first time step and initial_guess_ is provided, set initial_guess = initial_guess_
      //Throw an exception if initial_guess is not compatible with solution
      bool is_compatible = (initial_guess->space())->isCompatible(*this->initial_guess_->space());
      TEUCHOS_TEST_FOR_EXCEPTION(
          is_compatible != true, std::logic_error,
            "Error in Tempus::NemwarkImplicitDForm takeStep(): user-provided initial guess'!\n"
            << "for Newton is not compatible with solution vector!\n");
      Thyra::copy(*this->initial_guess_, initial_guess.ptr());
    }
    else {
      //Otherwise, set initial guess = diplacement predictor
      Thyra::copy(*d_pred, initial_guess.ptr());
    }

    //Set d_pred as initial guess for NOX solver, and solve nonlinear system.
    const Thyra::SolveStatus<Scalar> sStatus =
      this->solveImplicitODE(initial_guess);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.

    //solveImplicitODE will return converged solution in initial_guess
    //vector.  Copy it here to d_new, to define the new displacement.
    Thyra::copy(*initial_guess, d_new.ptr());

    //correct acceleration, velocity
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
StepperNewmarkImplicitDForm<Scalar>::getDefaultStepperState() {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Tempus::StepperState<Scalar>> stepperState =
      rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}

template <class Scalar>
void
StepperNewmarkImplicitDForm<Scalar>::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel /* verbLevel */) const {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  out << this->getStepperType() << "::describe:" << std::endl
      << "wrapperModel_ = " << this->wrapperModel_->description() << std::endl;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperNewmarkImplicitDForm<Scalar>::getValidParameters() const {
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<std::string>("Scheme Name", "Average Acceleration");
  pl->set<double>     ("Beta" , 0.25);
  pl->set<double>     ("Gamma", 0.5 );
  pl->set<std::string>("Solver Name", "Default Solver");
  pl->set<bool>       ("Zero Initial Guess", false);
  Teuchos::RCP<Teuchos::ParameterList> solverPL = defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}

}  // namespace Tempus
#endif  // Tempus_StepperNewmarkImplicitDForm_impl_hpp
