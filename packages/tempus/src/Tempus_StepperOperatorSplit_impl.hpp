// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperOperatorSplit_impl_hpp
#define Tempus_StepperOperatorSplit_impl_hpp

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Tempus_StepperFactory.hpp"
#include "Tempus_StepperOperatorSplitModifierDefault.hpp"

namespace Tempus {


template<class Scalar>
StepperOperatorSplit<Scalar>::StepperOperatorSplit()
{
  this->setStepperType(        "Operator Split");
  this->setUseFSAL(            false);
  this->setICConsistency(      "None");
  this->setICConsistencyCheck( false);

  this->setOrder   (1);
  this->setOrderMin(1);
  this->setOrderMax(1);
  this->setAppAction(Teuchos::null);

  OpSpSolnHistory_ = rcp(new SolutionHistory<Scalar>());
  OpSpSolnHistory_->setStorageLimit(2);
  OpSpSolnHistory_->setStorageType(Tempus::STORAGE_TYPE_STATIC);
}


template<class Scalar>
StepperOperatorSplit<Scalar>::StepperOperatorSplit(
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels,
  std::vector<Teuchos::RCP<Stepper<Scalar> > > subStepperList,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  int order,
  int orderMin,
  int orderMax,
  const Teuchos::RCP<StepperOperatorSplitAppAction<Scalar> >& stepperOSAppAction)
{
  this->setStepperType(        "Operator Split");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);

  this->setSubStepperList(subStepperList);
  this->setOrder   (order);
  this->setOrderMin(orderMin);
  this->setOrderMax(orderMax);

  this->setAppAction(stepperOSAppAction);
  OpSpSolnHistory_ = rcp(new SolutionHistory<Scalar>());
  OpSpSolnHistory_->setStorageLimit(2);
  OpSpSolnHistory_->setStorageType(Tempus::STORAGE_TYPE_STATIC);

  if ( !(appModels.empty()) ) {
    this->setModels(appModels);
    this->initialize();
  }
}


template<class Scalar>
void StepperOperatorSplit<Scalar>::setModel(
  const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel)
{
  if (appModel != Teuchos::null) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::setModel()");
    *out << "Warning -- No ModelEvaluator to set for StepperOperatorSplit, "
         << "because it is a Stepper of Steppers.\n" << std::endl;
  }
}

template<class Scalar>
Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
StepperOperatorSplit<Scalar>::getModel()
{
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model;
  typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::const_iterator
    subStepperIter = subStepperList_.begin();
  for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
    model = (*subStepperIter)->getModel();
    if (model != Teuchos::null) break;
  }
  if ( model == Teuchos::null ) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::getModel()");
    *out << "Warning -- StepperOperatorSplit::getModel() "
         << "Could not find a valid model!  Returning null!" << std::endl;
  }

  return model;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > /* solver */)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::setSolver()");
  *out << "Warning -- No solver to set for StepperOperatorSplit "
       << "because it is a Stepper of Steppers.\n" << std::endl;

  this->isInitialized_ = false;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setAppAction(
  Teuchos::RCP<StepperOperatorSplitAppAction<Scalar> > appAction)
{
  if (appAction == Teuchos::null) {
    // Create default appAction, otherwise keep current.
    if (stepperOSAppAction_ == Teuchos::null) {
      stepperOSAppAction_ =
        Teuchos::rcp(new StepperOperatorSplitModifierDefault<Scalar>());
    }
  } else {
  stepperOSAppAction_ =
    Teuchos::rcp_dynamic_cast<StepperOperatorSplitAppAction<Scalar> > (appAction, true);
  }
  this->isInitialized_ = false;
}


template<class Scalar>
void StepperOperatorSplit<Scalar>::addStepper(
  Teuchos::RCP<Stepper<Scalar> > stepper, bool useFSAL)
{
    stepper->setUseFSAL(useFSAL);
    stepper->initialize();
    subStepperList_.push_back(stepper);
}


template<class Scalar>
void StepperOperatorSplit<Scalar>::setSubStepperList(
  std::vector<Teuchos::RCP<Stepper<Scalar> > > subStepperList)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  subStepperList_ = subStepperList;

  typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::iterator
    subStepperIter = subStepperList_.begin();

  for (; subStepperIter<subStepperList_.end(); subStepperIter++) {
    auto subStepper = *(subStepperIter);
    bool useFSAL = subStepper->getUseFSAL();
    if (useFSAL) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::createSubSteppers()");
      *out << "Warning -- subStepper = '"
           << subStepper->getStepperType() << "' has \n"
           << "  subStepper->getUseFSAL() = " << useFSAL << ".\n"
           << "  subSteppers usually can not use the FSAL priniciple with\n"
           << "  operator splitting.  Proceeding with it set to true.\n"
           << std::endl;
    }
  }

  this->isInitialized_ = false;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setModels(
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  TEUCHOS_TEST_FOR_EXCEPTION(subStepperList_.size() != appModels.size(),
    std::logic_error, "Error - Number of models and Steppers do not match!\n"
    << "  There are " << appModels.size() << " models.\n"
    << "  There are " << subStepperList_.size() << " steppers.\n");

  typename std::vector<RCP<const Thyra::ModelEvaluator<Scalar> > >::iterator
    appModelIter = appModels.begin();

  typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::iterator
    subStepperIter = subStepperList_.begin();

  for (; appModelIter<appModels.end() || subStepperIter<subStepperList_.end();
       appModelIter++, subStepperIter++)
  {
    auto appModel = *(appModelIter);
    auto subStepper = *(subStepperIter);
    subStepper->setModel(appModel);
    subStepper->initialize();
  }

  this->isInitialized_ = false;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::initialize()
{
  if (tempState_ == Teuchos::null) {
    Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >model = this->getModel();
    TEUCHOS_TEST_FOR_EXCEPTION( model == Teuchos::null, std::logic_error,
      "Error - StepperOperatorSplit::initialize() Could not find "
      "a valid model!\n");
    tempState_ = createSolutionStateME(model, this->getDefaultStepperState());
  }

  if (!isOneStepMethod() ) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::initialize()");
    typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::const_iterator
      subStepperIter = subStepperList_.begin();
    for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
      *out << "SubStepper, " << (*subStepperIter)->getStepperType()
           << ", isOneStepMethod = " << (*subStepperIter)->isOneStepMethod()
           << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!isOneStepMethod(), std::logic_error,
    "Error - OperatorSplit only works for one-step methods!\n");
  }

  Stepper<Scalar>::initialize();
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::iterator
    subStepperIter = subStepperList_.begin();
  for (; subStepperIter < subStepperList_.end(); subStepperIter++)
    (*subStepperIter)->setInitialConditions(solutionHistory);

  Teuchos::RCP<SolutionState<Scalar> > initialState =
    solutionHistory->getCurrentState();

  // Check if we need Stepper storage for xDot
  this->setStepperXDot(initialState->getXDot());
  if (initialState->getXDot() == Teuchos::null)
    this->setStepperXDot(initialState->getX()->clone_v());

}

template<class Scalar>
void StepperOperatorSplit<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  this->checkInitialized();

  using Teuchos::RCP;

  TEMPUS_FUNC_TIME_MONITOR("Tempus::StepperOperatorSplit::takeStep()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(solutionHistory->getNumStates() < 2,
      std::logic_error,
      "Error - StepperOperatorSplit<Scalar>::takeStep(...)\n"
      "Need at least two SolutionStates for OperatorSplit.\n"
      "  Number of States = " << solutionHistory->getNumStates() << "\n"
      "Try setting in \"Solution History\" \"Storage Type\" = \"Undo\"\n"
      "  or \"Storage Type\" = \"Static\" and \"Storage Limit\" = \"2\"\n");
    RCP<StepperOperatorSplit<Scalar> > thisStepper = Teuchos::rcpFromRef(*this);
    stepperOSAppAction_->execute(solutionHistory, thisStepper,
      StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION::BEGIN_STEP);

    RCP<SolutionState<Scalar> > workingState=solutionHistory->getWorkingState();

    // Create OperatorSplit SolutionHistory to pass to subSteppers.
    tempState_->copy(solutionHistory->getCurrentState());
    OpSpSolnHistory_->clear();
    OpSpSolnHistory_->addState(tempState_);
    OpSpSolnHistory_->addWorkingState(workingState, false);

    RCP<SolutionState<Scalar> > currentSubState =
      OpSpSolnHistory_->getCurrentState();
    RCP<SolutionState<Scalar> > workingSubState =
      OpSpSolnHistory_->getWorkingState();

    bool pass = true;
    typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::iterator
      subStepperIter = subStepperList_.begin();
    for (; subStepperIter < subStepperList_.end() and pass; subStepperIter++) {

      stepperOSAppAction_->execute(solutionHistory, thisStepper,
        StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION::BEFORE_STEPPER);

      (*subStepperIter)->takeStep(OpSpSolnHistory_);

      stepperOSAppAction_->execute(solutionHistory, thisStepper,
        StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION::AFTER_STEPPER);

      if (workingSubState->getSolutionStatus() == Status::FAILED) {
        pass = false;
        Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::takeStep()");
        *out << "SubStepper, " << (*subStepperIter)->getStepperType()
             << ", failed!" << std::endl;
        break;
      }

      // "promote" workingSubState
      currentSubState = OpSpSolnHistory_->getCurrentState();
      currentSubState->copySolutionData(workingSubState);
    }

    if (pass == true) workingState->setSolutionStatus(Status::PASSED);
    else              workingState->setSolutionStatus(Status::FAILED);
    workingState->setOrder(this->getOrder());
    workingState->computeNorms(solutionHistory->getCurrentState());
    OpSpSolnHistory_->clear();

    stepperOSAppAction_->execute(solutionHistory, thisStepper,
      StepperOperatorSplitAppAction<Scalar>::ACTION_LOCATION::END_STEP);
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
Teuchos::RCP<Tempus::StepperState<Scalar> > StepperOperatorSplit<Scalar>::
getDefaultStepperState()
{
  Teuchos::RCP<Tempus::StepperState<Scalar> > stepperState =
    rcp(new StepperState<Scalar>(this->getStepperType()));
  return stepperState;
}


template<class Scalar>
void StepperOperatorSplit<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << std::endl;
  Stepper<Scalar>::describe(out, verbLevel);

  out << "--- StepperOperatorSplit ---\n";
  out << "  subStepperList_.size() = " << subStepperList_.size() << std::endl;
  typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::const_iterator
    subStepperIter = subStepperList_.begin();
  for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
    out << "    SubStepper       = " << (*subStepperIter)->getStepperType()<<std::endl;
    out << "                     = " << (*subStepperIter)->isInitialized() << std::endl;
    out << "                     = " << (*subStepperIter) << std::endl;
  }
  out << "  OpSpSolnHistory_    = " << OpSpSolnHistory_      << std::endl;
  out << "  tempState_          = " << tempState_      << std::endl;
  out << "  stepperOSAppAction_ = " << stepperOSAppAction_ << std::endl;
  out << "  order_              = " << order_      << std::endl;
  out << "  orderMin_           = " << orderMin_    << std::endl;
  out << "  orderMax_           = " << orderMax_    << std::endl;
  out << "----------------------------" << std::endl;
}


template<class Scalar>
bool StepperOperatorSplit<Scalar>::isValidSetup(Teuchos::FancyOStream & out) const
{
  bool isValidSetup = true;

  if ( !Stepper<Scalar>::isValidSetup(out) ) isValidSetup = false;

  if (subStepperList_.size() == 0) {
    isValidSetup = false;
    out << "The substepper list is empty!\n";
  }

  typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::const_iterator
    subStepperIter = subStepperList_.begin();

  for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
    auto subStepper = *(subStepperIter);
    if ( !subStepper->isInitialized() ) {
      isValidSetup = false;
      out << "The subStepper, " << subStepper->description()
          << ", is not initialized!\n";
    }
  }
  if (stepperOSAppAction_ == Teuchos::null) {
    isValidSetup = false;
    out << "The Operator-Split AppAction is not set!\n";
  }

  return isValidSetup;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperOperatorSplit<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  getValidParametersBasic(pl, this->getStepperType());
  pl->set<int>   ("Minimum Order", 1,
    "Minimum Operator-split order.  (default = 1)\n");
  pl->set<int>   ("Order", 1,
    "Operator-split order.  (default = 1)\n");
  pl->set<int>   ("Maximum Order", 1,
    "Maximum Operator-split order.  (default = 1)\n");

  pl->set<std::string>("Stepper List", "",
    "Comma deliminated list of single quoted Steppers, e.g., \"'Operator 1', 'Operator 2'\".");

  return pl;
}


} // namespace Tempus
#endif // Tempus_StepperOperatorSplit_impl_hpp
