// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperTrapezoidal_impl_hpp
#define Tempus_StepperTrapezoidal_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


// StepperTrapezoidal definitions:
template<class Scalar>
StepperTrapezoidal<Scalar>::StepperTrapezoidal(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
void StepperTrapezoidal<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (this->stepperObserver_ == Teuchos::null) {
      stepperTrapObserver_ =
        Teuchos::rcp(new StepperTrapezoidalObserver<Scalar>());
      this->stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >
          (stepperTrapObserver_);
     }
  } else {
    this->stepperObserver_ = obs;
    stepperTrapObserver_ =
      Teuchos::rcp_dynamic_cast<StepperTrapezoidalObserver<Scalar> >
        (this->stepperObserver_);
  }
}


template<class Scalar>
void StepperTrapezoidal<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->wrapperModel_ == Teuchos::null, std::logic_error,
    "Error - Need to set the model, setModel(), before calling "
    "StepperTrapezoidal::initialize()\n");

  this->setParameterList(this->stepperPL_);
  this->setSolver();
  this->setObserver();
}


template<class Scalar>
void StepperTrapezoidal<Scalar>::setInitialConditions (
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  int numStates = solutionHistory->getNumStates();

  TEUCHOS_TEST_FOR_EXCEPTION(numStates < 1, std::logic_error,
    "Error - setInitialConditions() needs at least one SolutionState\n"
    "        to set the initial condition.  Number of States = " << numStates);

  if (numStates > 1) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperTrapezoidal::setInitialConditions()");
    *out << "Warning -- SolutionHistory has more than one state!\n"
         << "Setting the initial conditions on the currentState.\n"<<std::endl;
  }

  RCP<SolutionState<Scalar> > initialState = solutionHistory->getCurrentState();
  RCP<Thyra::VectorBase<Scalar> > x = initialState->getX();

  // Use x from inArgs as ICs, if needed.
  auto inArgs = this->wrapperModel_->getNominalValues();
  if (x == Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION( (x == Teuchos::null) &&
      (inArgs.get_x() == Teuchos::null), std::logic_error,
      "Error - setInitialConditions() needs the ICs from the SolutionHistory\n"
      "        or getNominalValues()!\n");

    x = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(inArgs.get_x());
    initialState->setX(x);
  }

  // Check if we need to create an xDot.
  if (initialState->getXDot() == Teuchos::null) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperTrapezoidal::setInitialConditions()");
    *out << "Warning -- Trapezoidal needs xDot in the SolutionStates.\n"
         << "  Creating xDot on the initialState.\n"<<std::endl;
    this->setStepperXDot(initialState->getX()->clone_v());
  }

  // Perform IC Consistency
  std::string icConsistency = this->getICConsistency();
  if (icConsistency == "None") {
    if (initialState->getXDot() == Teuchos::null) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,
        "StepperBackwardEuler::setInitialConditions()");
      *out << "Warning -- Requested IC consistency of 'None' but\n"
           << "           initialState does not have an xDot.\n"
           << "           Setting a 'Zero' xDot!\n" << std::endl;

      Thyra::assign(this->getStepperXDot(initialState).ptr(), Scalar(0.0));
    }
  }
  else if (icConsistency == "Zero") {
    Thyra::assign(this->getStepperXDot(initialState).ptr(), Scalar(0.0));
  }
  else if (icConsistency == "App") {
    auto xDot = Teuchos::rcp_const_cast<Thyra::VectorBase<Scalar> >(
                  inArgs.get_x_dot());
    TEUCHOS_TEST_FOR_EXCEPTION(xDot == Teuchos::null, std::logic_error,
      "Error - setInitialConditions() requested 'App' for IC consistency,\n"
      "        but 'App' returned a null pointer for xDot!\n");
    Thyra::assign(this->getStepperXDot(initialState).ptr(), *xDot);
  }
  else if (icConsistency == "Consistent") {
    // Solve f(x, xDot,t) = 0.
    const Scalar time = initialState->getTime();
    const Scalar dt   = initialState->getTimeStep();
    RCP<TimeDerivative<Scalar> > timeDer = Teuchos::null;
    const Scalar alpha = 1.0;    // d(xDot)/d(xDot)
    const Scalar beta  = 0.0;    // d(x   )/d(xDot)
    RCP<ImplicitODEParameters<Scalar> > p =
      Teuchos::rcp(new ImplicitODEParameters<Scalar>(timeDer,dt,alpha,beta,
                                                     SOLVE_FOR_XDOT_CONST_X));

    auto xDot = this->getStepperXDot(initialState);
    const Thyra::SolveStatus<Scalar> sStatus =
      this->solveImplicitODE(x, xDot, time, p);

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

  if (initialState->getXDot() == Teuchos::null) {
    initialState->setXDot(this->getStepperXDot(initialState));
    this->setStepperXDot(Teuchos::null);
  }

  // At this point, x and xDot are sync'ed or consistent
  // at the same time level for the initialState.
  initialState->setIsSynced(true);

  TEUCHOS_TEST_FOR_EXCEPTION( !(this->getUseFSAL()), std::logic_error,
    "Error - The First-Step-As-Last (FSAL) principle is required\n"
    "        for the Trapezoidal Stepper (i.e., useFSAL=true)!\n");
//   There are at least two ways around this, but are not implemented.
//    - Do a solve for xDotOld, xDot_{n-1}, at each time step as for the
//      initial conditions.  This is expensive since you would be doing
//      two solves every time step.
//    - Use evaluateExplicitODE to get xDot_{n-1} if the application
//      provides it.  Explicit evaluations are cheaper but requires the
//      application to implement xDot = f(x,t).
}


template<class Scalar>
void StepperTrapezoidal<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperTrapezoidal::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperTrapezoidal<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for Trapezoidal.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");

    this->stepperObserver_->observeBeginTakeStep(solutionHistory, *this);
    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<const Thyra::VectorBase<Scalar> > xOld    = currentState->getX();
    RCP<const Thyra::VectorBase<Scalar> > xDotOld = currentState->getXDot();
    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > xDot = workingState->getXDot();

    const Scalar time  = workingState->getTime();
    const Scalar dt    = workingState->getTimeStep();
    const Scalar alpha = 2.0/dt;

    // Setup TimeDerivative
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
      Teuchos::rcp(new StepperTrapezoidalTimeDerivative<Scalar>(
        alpha, xOld, xDotOld));

    // Setup InArgs and OutArgs
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<Scalar>  inArgs  = this->wrapperModel_->getInArgs();
    MEB::OutArgs<Scalar> outArgs = this->wrapperModel_->getOutArgs();
    inArgs.set_x(x);
    if (inArgs.supports(MEB::IN_ARG_x_dot    )) inArgs.set_x_dot    (xDot);
    if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (time);
    if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(dt);
    if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (alpha);
    if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (1.0);

    this->wrapperModel_->setForSolve(timeDer, inArgs, outArgs);

    stepperTrapObserver_->observeBeforeSolve(solutionHistory, *this);

    const Thyra::SolveStatus<Scalar> sStatus = this->solveImplicitODE(x);

    stepperTrapObserver_->observeAfterSolve(solutionHistory, *this);

    if (workingState->getXDot() != Teuchos::null)
      timeDer->compute(x, xDot);

    workingState->setSolutionStatus(sStatus);  // Converged --> pass.
    workingState->setOrder(this->getOrder());
    this->stepperObserver_->observeEndTakeStep(solutionHistory, *this);
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
StepperTrapezoidal<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperTrapezoidal<Scalar>::description() const
{
  std::string name = "Trapezoidal Method";
  return(name);
}


template<class Scalar>
void StepperTrapezoidal<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "wrapperModel_ = " << this->wrapperModel_->description() << std::endl;
}


template <class Scalar>
void StepperTrapezoidal<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  Teuchos::RCP<Teuchos::ParameterList> stepperPL = this->stepperPL_;
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (stepperPL == Teuchos::null) stepperPL = this->getDefaultParameters();
  } else {
    stepperPL = pList;
  }
  if (!(stepperPL->isParameter("Solver Name"))) {
    stepperPL->set<std::string>("Solver Name", "Default Solver");
    Teuchos::RCP<Teuchos::ParameterList> solverPL =
      this->defaultSolverParameters();
    stepperPL->set("Default Solver", *solverPL);
  }
  // Can not validate because of optional Parameters (e.g., Solver Name).
  // stepperPL->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType = stepperPL->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Trapezoidal Method",
    std::logic_error,
       "Error - Stepper Type is not 'Trapezoidal Method'!\n"
    << "  Stepper Type = "<<stepperPL->get<std::string>("Stepper Type")<<"\n");

  this->stepperPL_ = stepperPL;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperTrapezoidal<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", this->description());
  this->getValidParametersBasic(pl);
  pl->set<bool>("Initial Condition Consistency Check",
                false);              // Default is false for this stepper.
  pl->set<bool>       ("Zero Initial Guess", false);
  pl->set<std::string>("Solver Name", "",
    "Name of ParameterList containing the solver specifications.");

  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperTrapezoidal<Scalar>::getDefaultParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;
  using Teuchos::rcp_const_cast;

  RCP<ParameterList> pl =
    rcp_const_cast<ParameterList>(this->getValidParameters());

  pl->set<std::string>("Solver Name", "Default Solver");
  RCP<ParameterList> solverPL = this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperTrapezoidal<Scalar>::getNonconstParameterList()
{
  return(this->stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperTrapezoidal<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = this->stepperPL_;
  this->stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperTrapezoidal_impl_hpp
