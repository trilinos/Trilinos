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


namespace Tempus {


template<class Scalar>
StepperOperatorSplit<Scalar>::StepperOperatorSplit()
  : OpSpSolnHistory_(Teuchos::null)
{
  this->setStepperType(        "Operator Split");
  this->setUseFSAL(            this->getUseFSALDefault());
  this->setICConsistency(      this->getICConsistencyDefault());
  this->setICConsistencyCheck( this->getICConsistencyCheckDefault());

  this->setOrder   (1);
  this->setOrderMin(1);
  this->setOrderMax(1);

  this->setObserver();
}

template<class Scalar>
StepperOperatorSplit<Scalar>::StepperOperatorSplit(
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels,
  std::vector<Teuchos::RCP<Stepper<Scalar> > > subStepperList,
  const Teuchos::RCP<StepperObserver<Scalar> >& obs,
  bool useFSAL,
  std::string ICConsistency,
  bool ICConsistencyCheck,
  int order,
  int orderMin,
  int orderMax)
    : OpSpSolnHistory_(Teuchos::null)
{
  this->setStepperType(        "Operator Split");
  this->setUseFSAL(            useFSAL);
  this->setICConsistency(      ICConsistency);
  this->setICConsistencyCheck( ICConsistencyCheck);

  this->setSubStepperList(subStepperList);
  this->setOrder   (order);
  this->setOrderMin(orderMin);
  this->setOrderMax(orderMax);

  this->setObserver(obs);

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
  return;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setNonConstModel(
  const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& appModel)
{
  if (appModel != Teuchos::null) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::setModel()");
    *out << "Warning -- No ModelEvaluator to set for StepperOperatorSplit, "
         << "because it is a Stepper of Steppers.\n" << std::endl;
  }
  return;
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
  return;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setObserver(
  Teuchos::RCP<StepperObserver<Scalar> > obs)
{
  if (obs == Teuchos::null) {
    // Create default observer, otherwise keep current observer.
    if (stepperOSObserver_ == Teuchos::null) {
      stepperOSObserver_ =
        Teuchos::rcp(new StepperOperatorSplitObserver<Scalar>());
     }
  } else {
    stepperOSObserver_ =
      Teuchos::rcp_dynamic_cast<StepperOperatorSplitObserver<Scalar> > (obs, true);
  }
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
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::initialize()
{
  TEUCHOS_TEST_FOR_EXCEPTION( subStepperList_.size() == 0, std::logic_error,
    "Error - Need to set the subSteppers, createSubSteppers(), before calling "
    "StepperOperatorSplit::initialize()\n");

  OpSpSolnHistory_ = rcp(new SolutionHistory<Scalar>());
  OpSpSolnHistory_->setStorageLimit(2);
  OpSpSolnHistory_->setStorageType(Tempus::STORAGE_TYPE_STATIC);

  if (tempState_ == Teuchos::null) {
    Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >model = this->getModel();
    TEUCHOS_TEST_FOR_EXCEPTION( model == Teuchos::null, std::logic_error,
      "Error - StepperOperatorSplit::initialize() Could not find "
      "a valid model!\n");
    //tempState_ = rcp(new SolutionState<Scalar>()); Doesn't seem to work?!
    tempState_ = rcp(new SolutionState<Scalar>(
      model, this->getDefaultStepperState()));
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
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setInitialConditions(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
  typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::iterator
    subStepperIter = subStepperList_.begin();
  for (; subStepperIter < subStepperList_.end(); subStepperIter++)
    (*subStepperIter)->setInitialConditions(solutionHistory);
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::takeStep(
  const Teuchos::RCP<SolutionHistory<Scalar> >& solutionHistory)
{
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

    stepperOSObserver_->observeBeginTakeStep(solutionHistory, *this);

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
      int index = subStepperIter - subStepperList_.begin();

      stepperOSObserver_->observeBeforeStepper(index, solutionHistory, *this);

      (*subStepperIter)->takeStep(OpSpSolnHistory_);

      stepperOSObserver_->observeAfterStepper(index, solutionHistory, *this);

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
    OpSpSolnHistory_->clear();
    stepperOSObserver_->observeEndTakeStep(solutionHistory, *this);
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
   const Teuchos::EVerbosityLevel      /* verbLevel */) const
{
  out << this->getStepperType() << "::describe:" << std::endl;
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
