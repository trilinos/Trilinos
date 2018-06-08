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
StepperOperatorSplit<Scalar>::StepperOperatorSplit(
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels,
  Teuchos::RCP<Teuchos::ParameterList> pList)
  : stepperPL_(Teuchos::null), OpSpSolnHistory_(Teuchos::null),
    stepperOSObserver_(Teuchos::null)
{
  this->setParameterList(pList);
  this->createSubSteppers(appModels);
  this->initialize();
}

template<class Scalar>
StepperOperatorSplit<Scalar>::StepperOperatorSplit()
  : stepperPL_(Teuchos::null), OpSpSolnHistory_(Teuchos::null),
    stepperOSObserver_(Teuchos::null)
{
  this->setParameterList(Teuchos::null);
  // Still require
  //  * Setting models and steppers, i.e., addStepper()
  //  * Calling initialize()
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
  TEUCHOS_TEST_FOR_EXCEPTION( model == Teuchos::null, std::logic_error,
    "Error - StepperOperatorSplit::getModel() Could not find a valid model!\n");

  return model;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setSolver(std::string solverName)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::setSolver()");
  *out << "Warning -- No solver to set for StepperOperatorSplit, "
       << "because it is a Stepper of Steppers.\n" << std::endl;
  return;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setSolver(
  Teuchos::RCP<Teuchos::ParameterList> solverPL)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::setSolver()");
  *out << "Warning -- No solver to set for StepperOperatorSplit "
       << "because it is a Stepper of Steppers.\n" << std::endl;
  return;
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::setSolver(
  Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > solver)
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
    if (stepperObserver_ == Teuchos::null) {
      stepperOSObserver_ =
        Teuchos::rcp(new StepperOperatorSplitObserver<Scalar>());
      stepperObserver_ =
        Teuchos::rcp_dynamic_cast<StepperObserver<Scalar> >(stepperOSObserver_);
     }
  } else {
    stepperObserver_ = obs;
    stepperOSObserver_ =
      Teuchos::rcp_dynamic_cast<StepperOperatorSplitObserver<Scalar> >
        (stepperObserver_);
  }
}

template<class Scalar>
void StepperOperatorSplit<Scalar>::createSubSteppers(
  std::vector<Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > > appModels)
{
  using Teuchos::RCP;
  using Teuchos::ParameterList;

  // Parse Stepper List String
  std::vector<std::string> stepperListStr;
  stepperListStr.clear();
  std::string str = stepperPL_->get<std::string>("Stepper List");
  std::string delimiters(",");
  // Skip delimiters at the beginning
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find the first delimiter
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
    std::string token = str.substr(lastPos,pos-lastPos);
    // Strip single quotes
    std::string::size_type beg = token.find_first_of("'") + 1;
    std::string::size_type end = token.find_last_of ("'");
    stepperListStr.push_back(token.substr(beg,end-beg));

    lastPos = str.find_first_not_of(delimiters, pos); // Skip delimiters
    pos = str.find_first_of(delimiters, lastPos);     // Find next delimiter
  }

  TEUCHOS_TEST_FOR_EXCEPTION(stepperListStr.size() != appModels.size(),
    std::logic_error, "Error - Number of models and Steppers do not match!\n"
    << "  There are " << appModels.size() << " models.\n"
    << "  There are " << stepperListStr.size() << " steppers.\n"
    << "    " << str << "\n");

  RCP<StepperFactory<Scalar> > sf = Teuchos::rcp(new StepperFactory<Scalar>());
  typename
    std::vector<RCP<const Thyra::ModelEvaluator<Scalar> > >::iterator
      aMI = appModels.begin();
  typename std::vector<std::string>::iterator sLSI = stepperListStr.begin();

  for (; aMI<appModels.end() || sLSI<stepperListStr.end(); aMI++, sLSI++) {
    RCP<ParameterList> subStepperPL = Teuchos::sublist(stepperPL_,*sLSI,true);
    subStepperList_.push_back(sf->createStepper(*aMI, subStepperPL));
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

  tempState_ = rcp(new SolutionState<Scalar>(this->getModel(),
                                             this->getDefaultStepperState()));
  this->setObserver();

  if (!isOneStepMethod() ) {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::initialize()");
    typename std::vector<Teuchos::RCP<Stepper<Scalar> > >::const_iterator
      subStepperIter = subStepperList_.begin();
    for (; subStepperIter < subStepperList_.end(); subStepperIter++) {
      *out << "SubStepper, " << (*subStepperIter)->description()
           << ", isOneStepMethod = " << (*subStepperIter)->isOneStepMethod()
           << std::endl;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!isOneStepMethod(), std::logic_error,
    "Error - OperatorSplit only works for one-step methods!\n");
  }
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

      if (workingSubState->getStepperStatus() == Status::FAILED) {
        pass = false;
        Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"StepperOperatorSplit::takeStep()");
        *out << "SubStepper, " << (*subStepperIter)->description()
             << ", failed!" << std::endl;
        break;
      }

      // "promote" workingSubState
      currentSubState = OpSpSolnHistory_->getCurrentState();
      currentSubState->copySolutionStepperState(workingSubState);
    }

    if (pass == true) workingState->setStepperStatus(Status::PASSED);
    else              workingState->setStepperStatus(Status::FAILED);
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
    rcp(new StepperState<Scalar>(description()));
  return stepperState;
}


template<class Scalar>
std::string StepperOperatorSplit<Scalar>::description() const
{
  std::string name = "Operator Split";
  return(name);
}


template<class Scalar>
void StepperOperatorSplit<Scalar>::describe(
   Teuchos::FancyOStream               &out,
   const Teuchos::EVerbosityLevel      verbLevel) const
{
  out << description() << "::describe:" << std::endl;
}


template <class Scalar>
void StepperOperatorSplit<Scalar>::setParameterList(
  const Teuchos::RCP<Teuchos::ParameterList> & pList)
{
  Teuchos::RCP<Teuchos::ParameterList> stepperPL = this->stepperPL_;
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (stepperPL == Teuchos::null) stepperPL = this->getDefaultParameters();
  } else {
    stepperPL = pList;
  }
  // Can not validate because of optional Parameters, e.g. operators.
  //stepperPL->validateParametersAndSetDefaults(*this->getValidParameters());

  std::string stepperType = stepperPL->get<std::string>("Stepper Type");
  TEUCHOS_TEST_FOR_EXCEPTION( stepperType != "Operator Split", std::logic_error,
       "Error - Stepper Type is not 'Operator Split'!\n"
    << "  Stepper Type = "<< pList->get<std::string>("Stepper Type") << "\n");

  this->stepperPL_ = stepperPL;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
StepperOperatorSplit<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  pl->setName("Default Stepper - " + this->description());
  pl->set<std::string>("Stepper Type", "Operator Split",
    "'Stepper Type' must be 'Operator Split'.");
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


template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperOperatorSplit<Scalar>::getDefaultParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
  *pl = *(this->getValidParameters());
  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperOperatorSplit<Scalar>::getNonconstParameterList()
{
  return(stepperPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
StepperOperatorSplit<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = stepperPL_;
  stepperPL_ = Teuchos::null;
  return(temp_plist);
}


} // namespace Tempus
#endif // Tempus_StepperOperatorSplit_impl_hpp
