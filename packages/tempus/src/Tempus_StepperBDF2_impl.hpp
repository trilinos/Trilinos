// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperBDF2_impl_hpp
#define Tempus_StepperBDF2_impl_hpp

#include "Tempus_config.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_WrapperModelEvaluatorBasic.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "NOX_Thyra.H"


namespace Tempus {

// Forward Declaration for recursive includes (this Stepper <--> StepperFactory)
template<class Scalar> class StepperFactory;


// StepperBDF2 definitions:
template<class Scalar>
StepperBDF2<Scalar>::StepperBDF2(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Set all the input parameters and call initialize
  this->setParameterList(pList);
  this->setModel(appModel);
  this->initialize();
}


template<class Scalar>
void StepperBDF2<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->validImplicitODE_DAE(appModel);
  if (wrapperModel_ != Teuchos::null) wrapperModel_ = Teuchos::null;
  wrapperModel_ =
    Teuchos::rcp(new WrapperModelEvaluatorBasic<Scalar>(appModel));
}


template<class Scalar>
void StepperBDF2<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  this->setModel(appModel);
}


/** \brief Set the solver to a pre-defined solver in the ParameterList.
 *  The solver is set to solverName sublist in the Stepper's ParameterList.
 *  The solverName sublist should already be defined in the Stepper's
 *  ParameterList.  Otherwise it will fail.
 */
template<class Scalar>
void StepperBDF2<Scalar>::setSolver(std::string solverName)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> solverPL = Teuchos::sublist(stepperPL_, solverName, true);
  stepperPL_->set("Solver Name", solverName);
  solver_ = rcp(new Thyra::NOXNonlinearSolver());
  RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
  solver_->setParameterList(noxPL);
}


/** \brief Set the solver to the supplied Parameter sublist.
 *  This adds a new solver Parameter sublist to the Stepper's ParameterList.
 *  If the solver sublist is null, the solver is set to the solver name
 *  in the Stepper's ParameterList.
 */
template<class Scalar>
void StepperBDF2<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::string solverName = stepperPL_->get<std::string>("Solver Name");
  if (is_null(solverPL)) {
    // Create default solver, otherwise keep current solver.
    if (solver_ == Teuchos::null) {
      solverPL = Teuchos::sublist(stepperPL_, solverName, true);
      solver_ = rcp(new Thyra::NOXNonlinearSolver());
      RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
      solver_->setParameterList(noxPL);
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( solverName == solverPL->name(),
      std::logic_error,
         "Error - Trying to add a solver that is already in ParameterList!\n"
      << "  Stepper Type = "<< stepperPL_->get<std::string>("Stepper Type")
      << "\n" << "  Solver Name  = "<<solverName<<"\n");
    solverName = solverPL->name();
    stepperPL_->set("Solver Name", solverName);
    stepperPL_->set(solverName, solverPL);      // Add sublist
    solver_ = rcp(new Thyra::NOXNonlinearSolver());
    RCP<ParameterList> noxPL = Teuchos::sublist(solverPL, "NOX", true);
    solver_->setParameterList(noxPL);
  }
}


/** \brief Set the solver.
 *  This sets the solver to supplied solver and adds solver's ParameterList
 *  to the Stepper ParameterList.
 */
template<class Scalar>
void StepperBDF2<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> solverPL = solver->getNonconstParameterList();
  std::string solverName = solverPL->name();
  stepperPL_->set("Solver Name", solverName);
  stepperPL_->set(solverName, solverPL);      // Add sublist
  solver_ = solver;
}


/** \brief Set the startup stepper to a pre-defined stepper in the ParameterList.
 *  The startup stepper is set to stepperName sublist in the Stepper's
 *  ParameterList.  The stepperName sublist should already be defined
 *  in the Stepper's ParameterList.  Otherwise it will fail.
 */
template<class Scalar>
void StepperBDF2<Scalar>::setStartUpStepper(std::string startupStepperName)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> startupStepperPL = Teuchos::sublist(stepperPL_, startupStepperName, true);
  stepperPL_->set("Start Up Stepper Name", startupStepperName);
  if (startUpStepper_ != Teuchos::null) startUpStepper_ = Teuchos::null;
  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
}


/** \brief Set the start up stepper to the supplied Parameter sublist.
 *  This adds a new start up stepper Parameter sublist to the Stepper's ParameterList.
 *  If the start up stepper sublist is null, it tests if the stepper sublist is set in
 *  the Stepper's ParameterList.
 */
template<class Scalar>
void StepperBDF2<Scalar>::setStartUpStepper(
  Teuchos::RCP<Teuchos::ParameterList> startupStepperPL)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  std::string startupStepperName =
    stepperPL_->get<std::string>("Start Up Stepper Name","None");
  if (is_null(startupStepperPL)) {
    if (startupStepperName != "None") {
      RCP<ParameterList> startupStepperPL =
        Teuchos::sublist(stepperPL_, startupStepperName, true);
      RCP<StepperFactory<Scalar> > sf =
        Teuchos::rcp(new StepperFactory<Scalar>());
      startUpStepper_ =
        sf->createStepper(wrapperModel_->getAppModel(), startupStepperPL);
    }
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION( startupStepperName == startupStepperPL->name(),
      std::logic_error,
         "Error - Trying to add a startup stepper that is already in ParameterList!\n"
      << "  Stepper Type = "<< stepperPL_->get<std::string>("Stepper Type")
      << "\n" << "  Start Up Stepper Name  = "<<startupStepperName<<"\n");
    startupStepperName = startupStepperPL->name();
    stepperPL_->set("Start Up Stepper Name", startupStepperName);
    stepperPL_->set(startupStepperName, startupStepperPL);           // Add sublist
    if (startUpStepper_ != Teuchos::null) startUpStepper_ = Teuchos::null;
    RCP<StepperFactory<Scalar> > sf =
      Teuchos::rcp(new StepperFactory<Scalar>());
    startUpStepper_ =
      sf->createStepper(wrapperModel_->getAppModel(), startupStepperPL);
  }
}


template<class Scalar>
void StepperBDF2<Scalar>::setObserver(
  Teuchos::RCP<StepperBDF2Observer<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (stepperBDF2Observer_ == Teuchos::null) {
      stepperBDF2Observer_ =
        Teuchos::rcp(new StepperBDF2Observer<Scalar>());
    }
  } else {
    stepperBDF2Observer_ = obs;
  }
}


template<class Scalar>
void StepperBDF2<Scalar>::initialize()
{
  this->setSolver();
  this->setStartUpStepper();
  this->setObserver();
  order_ = 2.0;
}


template<class Scalar>
void StepperBDF2<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperBDF2::takeStep()");
  {
    int numStates = solutionHistory->getNumStates();

    RCP<Thyra::VectorBase<Scalar> > xOld;
    RCP<Thyra::VectorBase<Scalar> > xOldOld;

    //If we are on first time step, call startup stepper
    //There is test of exception in computeStartUp that start up stepper
    //did not fail, so no need to check for failure here.
    if (numStates < 3) {
      computeStartUp(solutionHistory);
      numStates = solutionHistory->getNumStates();
    }
    if (numStates < 3) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error in Tempus::StepperBDF2::takeStep(): numStates after \n"
          << "startup stepper must be at least 3, whereas numStates = "
          << numStates <<"!\n" << "If running with Storage Type = Static, "
          << "make sure Storage Limit > 2.\n");
    }

    //IKT, FIXME: add error checking regarding states being consecutive and
    //whether interpolated states are OK to use.

    stepperBDF2Observer_->observeBeginTakeStep(solutionHistory, *this);

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();
    RCP<SolutionState<Scalar> > currentState=solutionHistory->getCurrentState();

    RCP<Thyra::VectorBase<Scalar> > x    = workingState->getX();
    RCP<Thyra::VectorBase<Scalar> > xDot = workingState->getXDot();

    //get time, dt and dtOld
    const Scalar time  = workingState->getTime();
    const Scalar dt    = workingState->getTimeStep();
    const Scalar dtOld = currentState->getTimeStep();

    //IKT, FIXME: the following should be changed once CO adds
    //methods to obtain past SolutionStates, e.g., getNM1State().
    //get previous 2 states
    xOld = (*solutionHistory)[numStates-2]->getX();
    xOldOld = (*solutionHistory)[numStates-3]->getX();
    order_ = 2.0;

    const Scalar alpha = (1.0/(dt + dtOld))*(1.0/dt)*(2.0*dt + dtOld);
    const Scalar beta = 1.0;

    // Setup TimeDerivative
    Teuchos::RCP<TimeDerivative<Scalar> > timeDer =
      Teuchos::rcp(new StepperBDF2TimeDerivative<Scalar>(dt, dtOld, xOld, xOldOld));

    // Setup InArgs and OutArgs
    typedef Thyra::ModelEvaluatorBase MEB;
    MEB::InArgs<Scalar>  inArgs  = wrapperModel_->getInArgs();
    MEB::OutArgs<Scalar> outArgs = wrapperModel_->getOutArgs();
    inArgs.set_x(x);
    if (inArgs.supports(MEB::IN_ARG_x_dot    )) inArgs.set_x_dot    (xDot);
    if (inArgs.supports(MEB::IN_ARG_t        )) inArgs.set_t        (time+dt);
    if (inArgs.supports(MEB::IN_ARG_step_size)) inArgs.set_step_size(dt);
    if (inArgs.supports(MEB::IN_ARG_alpha    )) inArgs.set_alpha    (alpha);
    if (inArgs.supports(MEB::IN_ARG_beta     )) inArgs.set_beta     (beta);

    wrapperModel_->setForSolve(timeDer, inArgs, outArgs);

    stepperBDF2Observer_->observeBeforeSolve(solutionHistory, *this);

    const Thyra::SolveStatus<Scalar> sStatus =
      this->solveNonLinear(wrapperModel_, *solver_, x);

    stepperBDF2Observer_->observeAfterSolve(solutionHistory, *this);

    timeDer->compute(x, xDot);

    if (sStatus.solveStatus == Thyra::SOLVE_STATUS_CONVERGED )
      workingState->getStepperState()->stepperStatus_ = Status::PASSED;
    else
      workingState->getStepperState()->stepperStatus_ = Status::FAILED;
    workingState->setOrder(this->getOrder());
    stepperBDF2Observer_->observeEndTakeStep(solutionHistory, *this);
  }
  return;
}

template<class Scalar>
void StepperBDF2<Scalar>::computeStartUp(
      const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  if (startUpStepper_ == Teuchos::null) return;

  //Take one step using startUpStepper_
  startUpStepper_->takeStep(solutionHistory);

  order_ = startUpStepper_->getOrder();

  Status & stepperStatus =
    solutionHistory->getWorkingState()->getStepperState()->stepperStatus_;

  //If startUpStepper_ failed, abort; otherwise, augment solutionHistory with newly-
  //computed state
  if (stepperStatus == Status::FAILED) {
    TEUCHOS_TEST_FOR_EXCEPTION(true,
      std::logic_error,
       "Error: Start Up Stepper Failed'!\n");
  }
  else {
    solutionHistory->promoteWorkingState();
    solutionHistory->initWorkingState();
  }
}



/** \brief Provide a StepperState to the SolutionState.
 *  This Stepper does not have any special state data,
 *  so just provide the base class StepperState with the
 *  Stepper description.  This can be checked to ensure
 *  that the input StepperState can be used by this Stepper.
 */
template<class Scalar>
Teuchos::RCP<Tempus::StepperState<Scalar> >
StepperBDF2<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperBDF2<Scalar>::description() const
{
  std::string name = "BDF2";
  return(name);
}


template<class Scalar>
void StepperBDF2<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl
      << "wrapperModel_ = " << wrapperModel_->description() << std::endl;
}


template <class Scalar>
void StepperBDF2<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (stepperPL_ == Teuchos::null) stepperPL_ = this->getDefaultParameters();
  } else {
    stepperPL_ = pList;
  }
  // Can not validate because of optional Parameters (e.g., Solver Name).
  //stepperPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType = stepperPL_->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "BDF2",
    std::logic_error,
       "Error - Stepper Type is not 'BDF2'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperBDF2<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set("Stepper Type", this->description());
  pl->set("Solver Name", "",
    "Name of ParameterList containing the solver specifications.");

  return pl;
}


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBDF2<Scalar>::getDefaultParameters() const
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", this->description());
  pl->set<std::string>("Solver Name", "Default Solver");

  RCP<ParameterList> solverPL = this->defaultSolverParameters();
  pl->set("Default Solver", *solverPL);

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBDF2<Scalar>::getNonconstParameterList()
{
  return(stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperBDF2<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperBDF2_impl_hpp
