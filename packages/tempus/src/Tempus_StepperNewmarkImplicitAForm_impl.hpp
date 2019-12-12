// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperNewmarkImplicitAForm_impl_hpp
#define Tempus_StepperNewmarkImplicitAForm_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"

//#define VERBOSE_DEBUG_OUTPUT
//#define DEBUG_OUTPUT

namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::
predictVelocity(Thyra::VectorBase<Scalar>& vPred,
                const Thyra::VectorBase<Scalar>& v,
                const Thyra::VectorBase<Scalar>& a,
                const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //vPred = v + dt*(1.0-gamma_)*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(vPred), 1.0, v, dt*(1.0-gamma_), a);
}

template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::
predictDisplacement(Thyra::VectorBase<Scalar>& dPred,
                    const Thyra::VectorBase<Scalar>& d,
                    const Thyra::VectorBase<Scalar>& v,
                    const Thyra::VectorBase<Scalar>& a,
                    const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<const Thyra::VectorBase<Scalar> > tmp =
    Thyra::createMember<Scalar>(dPred.space());
  //dPred = dt*v + dt*dt/2.0*(1.0-2.0*beta_)*a
  Scalar aConst = dt*dt/2.0*(1.0-2.0*beta_);
  Thyra::V_StVpStV(Teuchos::ptrFromRef(dPred), dt, v, aConst, a);
  //dPred += d;
  Thyra::Vp_V(Teuchos::ptrFromRef(dPred), d, 1.0);
}

template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::
correctVelocity(Thyra::VectorBase<Scalar>& v,
                const Thyra::VectorBase<Scalar>& vPred,
                const Thyra::VectorBase<Scalar>& a,
                const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //v = vPred + dt*gamma_*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(v), 1.0, vPred, dt*gamma_, a);
}

template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::
correctDisplacement(Thyra::VectorBase<Scalar>& d,
                    const Thyra::VectorBase<Scalar>& dPred,
                    const Thyra::VectorBase<Scalar>& a,
                    const Scalar dt) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  //d = dPred + beta_*dt*dt*a
  Thyra::V_StVpStV(Teuchos::ptrFromRef(d), 1.0, dPred, beta_*dt*dt, a);
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setBeta(Scalar beta)
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
void StepperNewmarkImplicitAForm<Scalar>::setGamma(Scalar gamma)
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
void StepperNewmarkImplicitAForm<Scalar>::setSchemeName(
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
       "\nError in Tempus::StepperNewmarkImplicitAForm!  "
       <<"Invalid Scheme Name = " << schemeName_ <<".  \n"
       <<"Valid Scheme Names are: 'Average Acceleration', "
       <<"'Linear Acceleration', \n"
       <<"'Central Difference' and 'User Defined'.\n");
  }
}


template<class Scalar>
StepperNewmarkImplicitAForm<Scalar>::StepperNewmarkImplicitAForm() :
  out_(Teuchos::VerboseObjectBase::getDefaultOStream())
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif

  this->setStepperType(        "Newmark Implicit a-Form");
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());
  this->setZeroInitialGuess(   false);
  this->setSchemeName(         "Average Acceleration");

  this->setObserver();
}


template<class Scalar>
StepperNewmarkImplicitAForm<Scalar>::StepperNewmarkImplicitAForm(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
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
  this->setStepperType(        "Newmark Implicit a-Form");
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


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  validSecondOrderODE_DAE(appModel);
  auto wrapperModel =
    Teuchos::rcp(new WrapperModelEvaluatorSecondOrder<Scalar>(appModel,
                                              "Newmark Implicit a-Form"));
  this->wrapperModel_ = wrapperModel;
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( this->wrapperModel_ == Teuchos::null,
    std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperNewmarkImplicitAForm::initialize()\n");

#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(numStates < 1, std::logic_error,
    "Error - setInitialConditions() needs at least one SolutionState\n"
    "        to set the initial condition.  Number of States = " << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperNewmarkImplicitAForm::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"<<std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x    = initialState->getX();
  RCP<Thyra::VectorBase<Scalar> > xDot = initialState->getXDot();

  auto inArgs = this->wrapperModel_->getNominalValues();
  TEUCHOS_TEST_FOR_EXCEPTION(
    !((x != Teuchos::null && xDot != Teuchos::null) ||
      (inArgs.get_x() != Teuchos::null &&
       inArgs.get_x_dot() != Teuchos::null)), std::logic_error,
    "Error - We need to set the initial conditions for x and xDot from\n"
    "        either initialState or appModel_->getNominalValues::InArgs\n"
    "        (but not from a mixture of the two).\n");

  // Use x and xDot from inArgs as ICs, if needed.
  if ( x == Teuchos::null || xDot == Teuchos::null ) {
    using Teuchos::rcp_const_cast;
    TEUCHOS_TEST_FOR_EXCEPTION( (inArgs.get_x() == Teuchos::null) ||
      (inArgs.get_x_dot() == Teuchos::null), std::logic_error,
      "Error - setInitialConditions() needs the ICs from the initialState\n"
      "        or getNominalValues()!\n");
    x    = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
    initialState->setX(x);
    xDot = rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x_dot());
    initialState->setXDot(xDot);
  }

  // Check if we need Stepper storage for xDotDot
  if (initialState->getXDotDot() == Teuchos::null)
    initialState->setXDotDot(initialState->getX()->clone_v());

  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    if (initialState->getXDotDot() == Teuchos::null) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,
        "StepperNewmarkImplicitAForm::setInitialConditions()");
      *out << "Warning -- Requested IC consistency of 'None' but\n"
           << "           initialState does not have an xDot.\n"
           << "           Setting a 'Zero' xDot!\n" << std::endl;

      Thyra::assign(this->getStepperXDotDot(initialState).ptr(), Scalar(0.0));
    }
  }
  else if (icConsistency == "Zero")
    Thyra::assign(this->getStepperXDotDot(initialState).ptr(), Scalar(0.0));
  else if (icConsistency == "App") {
    auto xDotDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                inArgs.get_x_dot_dot());
    TEUCHOS_TEST_FOR_EXCEPTION(xDotDot == Teuchos::null, std::logic_error,
      "Error - setInitialConditions() requested 'App' for IC consistency,\n"
      "        but 'App' returned a null pointer for xDotDot!\n");
    Thyra::assign(this->getStepperXDotDot(initialState).ptr(), *xDotDot);
  }
  else if (icConsistency == "Consistent") {
    // Solve f(x, xDot, xDotDot, t) = 0.
    const Scalar time = initialState->getTime();
    auto xDotDot = this->getStepperXDotDot(initialState);

    // Compute initial acceleration using initial displacement
    // and initial velocity.
    if (this->initial_guess_ != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        !((xDotDot->space())->isCompatible(*this->initial_guess_->space())),
        std::logic_error,
        "Error - User-provided initial guess for Newton is not compatible\n"
        "        with solution vector!\n");
      Thyra::copy(*this->initial_guess_, xDotDot.ptr());
    }
    else {
      Thyra::put_scalar(0.0, xDotDot.ptr());
    }

    auto wrapperModel =
      Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
        this->wrapperModel_);

    wrapperModel->initializeNewmark(xDot, x, 0.0, time, beta_, gamma_);
    const Thyra::SolveStatus<Scalar> sStatus = this->solveImplicitODE(xDotDot);

    TEUCHOS_TEST_FOR_EXCEPTION(
      sStatus.solveStatus != Thyra::SOLVE_STATUS_CONVERGED, std::logic_error,
      "Error - Solver failed while determining the initial conditions.\n"
      "        Solver status is "<<Thyra::toString(sStatus.solveStatus)<<".\n");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "Error - setInitialConditions() invalid IC consistency, "
      << icConsistency << ".\n");
  }

  // At this point, x, xDot and xDotDot are sync'ed or consistent
  // at the same time level for the initialState.
  initialState->setIsSynced(true);

  // Test for consistency.
  if (this->getICConsistencyCheck()) {
    auto f       = initialState->getX()->clone_v();
    auto xDotDot = this->getStepperXDotDot(initialState);

    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<Scalar>  appInArgs =
      this->wrapperModel_->getAppModel()->createInArgs();
    MEB::OutArgs<Scalar> appOutArgs =
      this->wrapperModel_->getAppModel()->createOutArgs();

    appInArgs.set_x        (x      );
    appInArgs.set_x_dot    (xDot   );
    appInArgs.set_x_dot_dot(xDotDot);

    appOutArgs.set_f(appOutArgs.get_f());

    appInArgs.set_W_x_dot_dot_coeff(Scalar(0.0));     // da/da
    appInArgs.set_alpha            (Scalar(0.0));     // dv/da
    appInArgs.set_beta             (Scalar(0.0));     // dd/da

    appInArgs.set_t        (initialState->getTime()    );

    this->wrapperModel_->getAppModel()->evalModel(appInArgs, appOutArgs);

    Scalar reldiff = Thyra::norm(*f);
    Scalar normx = Thyra::norm(*x);
    Scalar eps = Scalar(100.0)*std::abs(Teuchos::ScalarTraits<Scalar>::eps());
    if (normx > eps*reldiff) reldiff /= normx;

    if (reldiff > eps) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,
        "StepperNewmarkImplicitAForm::setInitialConditions()");
      *out << "Warning -- Failed consistency check but continuing!\n"
         << "  ||f(x,xDot,xDotDot,t)||/||x|| > eps" << std::endl
         << "  ||f(x,xDot,xDotDot,t)||       = " << Thyra::norm(*f)<< std::endl
         << "  ||x||                         = " << Thyra::norm(*x)<< std::endl
         << "  ||f(x,xDot,xDotDot,t)||/||x|| = " << reldiff        << std::endl
         << "                            eps = " << eps            << std::endl;
    }
  }

  if (!(this->getUseFSAL())) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,
      "StepperNewmarkImplicitAForm::setInitialConditions()");
    *out << "\nWarning -- The First-Step-As-Last (FSAL) principle is "
         << "part of the Newmark Implicit A-Form.  The default is to "
         << "set useFSAL=true, and useFSAL=false will be ignored." << std::endl;
  }
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperNewmarkImplicitAForm::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperNewmarkImplicitAForm<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for NewmarkImplicitAForm.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    auto wrapperModel =
      Teuchos::rcp_dynamic_cast<WrapperModelEvaluatorSecondOrder<Scalar> >(
        this->wrapperModel_);

    // Get values of d, v and a from previous step
    RCP<const Thyra::VectorBase<Scalar> > d_old = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > v_old = currentState->getXDot();
    RCP<      Thyra::VectorBase<Scalar> > a_old = currentState->getXDotDot();

    // Get new values of d, v and a from workingState
    RCP<Thyra::VectorBase<Scalar> > d_new = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > v_new = workingState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > a_new = workingState->getXDotDot();

    // Get time and dt
    const Scalar time = currentState->getTime();
    const Scalar dt = workingState->getTimeStep();
    Scalar t = time+dt;

    // Compute acceleration, a_old, using displacement (d_old) and
    // velocity (v_old), if needed.
    if (!(this->getUseFSAL()) && workingState->getNConsecutiveFailures() == 0) {
      wrapperModel->initializeNewmark(v_old, d_old, dt, time,
                                      Scalar(0.0), Scalar(0.0));
      const Thyra::SolveStatus<Scalar> sStatus = this->solveImplicitODE(a_old);

      workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    }

    // Compute displacement and velocity predictors
    predictDisplacement(*d_new, *d_old, *v_old, *a_old, dt);
    predictVelocity(*v_new, *v_old, *a_old, dt);

    // Inject d_new, v_new, a and other relevant data into wrapperModel
    wrapperModel->initializeNewmark(v_new,d_new,dt,t,beta_,gamma_);

    // Solve nonlinear system with a_new as initial guess
    const Thyra::SolveStatus<Scalar> sStatus = this->solveImplicitODE(a_new);

    // Correct velocity, displacement.
    correctVelocity(*v_new, *v_new, *a_new, dt);
    correctDisplacement(*d_new, *d_new, *a_new, dt);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
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
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperNewmarkImplicitAForm<Scalar>::
getDefaultStepperState()
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperNewmarkImplicitAForm<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      /* verbLevel */) const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  out << this->getStepperType() << "::describe:" << std::endl
      << "wrapperModel = " << this->wrapperModel_->description() << std::endl;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperNewmarkImplicitAForm<Scalar>::getValidParameters() const
{
#ifdef VERBOSE_DEBUG_OUTPUT
  *out_ << "DEBUG: " << __PRETTY_FUNCTION__ << "\n";
#endif
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<std::string>("Scheme Name", "Average Acceleration");
  pl->set<double>     ("Beta" , 0.25);
  pl->set<double>     ("Gamma", 0.5 );
  pl->set<bool>       ("Use FSAL", this->getUseFSALDefault());
  pl->set<std::string>("Initial Condition Consistency",
                       this->getICConsistencyDefault());
  pl->set<std::string>("Solver Name", "Default Solver");
  pl->set<bool>       ("Zero Initial Guess", false);
  Teuchos::RCP<Teuchos::ParameterList> solverPL = defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperNewmarkImplicitAForm_impl_hpp
