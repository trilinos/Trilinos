// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_SolutionHistory_impl_hpp
#define Tempus_SolutionHistory_impl_hpp

// Teuchos
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

// Tempus
#include "Tempus_SolutionStateMetaData.hpp"

//#include "Thyra_VectorStdOps.hpp"


namespace {

  static std::string Invalid_name    = "Invalid";
  static std::string KeepNewest_name = "Keep Newest";
  static std::string Undo_name       = "Undo";
  static std::string Static_name     = "Static";
  static std::string Unlimited_name  = "Unlimited";
  static std::string Storage_name    = "Storage Type";
  static std::string Storage_default = Unlimited_name;

  static std::string StorageLimit_name    = "Storage Limit";
  static int         StorageLimit_default = 2;

  std::vector<std::string> HistoryPolicies =
    {Invalid_name, KeepNewest_name, Undo_name, Static_name, Unlimited_name};

  const Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<Tempus::StorageType> >
    StorageTypeValidator = Teuchos::rcp(
        new Teuchos::StringToIntegralParameterEntryValidator<Tempus::StorageType>(
          HistoryPolicies,
          Teuchos::tuple<Tempus::StorageType>(
            Tempus::STORAGE_TYPE_INVALID,
            Tempus::STORAGE_TYPE_KEEP_NEWEST,
            Tempus::STORAGE_TYPE_UNDO,
            Tempus::STORAGE_TYPE_STATIC,
            Tempus::STORAGE_TYPE_UNLIMITED),
          Storage_name));

} // namespace


namespace Tempus {

template<class Scalar>
SolutionHistory<Scalar>::SolutionHistory(
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  using Teuchos::RCP;
  // Create history, an array of solution states.
  history_ = rcp(new std::vector<RCP<SolutionState<Scalar> > >);

  this->setParameterList(pList);

  if (Teuchos::as<int>(this->getVerbLevel()) >=
      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"SolutionHistory::SolutionHistory");
    *out << this->description() << std::endl;
  }
}


template<class Scalar>
void SolutionHistory<Scalar>::addState(
  const Teuchos::RCP<SolutionState<Scalar> >& state_)
{
  // Check that we're not going to exceed our storage limit:
  if (Teuchos::as<int>(history_->size()+1) > storageLimit) {
    switch (storageType) {
    case STORAGE_TYPE_INVALID: {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - Storage type is STORAGE_TYPE_INVALID.\n");
    }
    case STORAGE_TYPE_STATIC:
    case STORAGE_TYPE_KEEP_NEWEST:
    case STORAGE_TYPE_UNDO: {
      if (state_->getTime() > history_->front()->getTime()) {
        // Case:  State is older than the youngest state in history.
        // Remove state from the beginning of history, then add new state.
        history_->erase(history_->begin());
      } else {
        // Case:  State is younger than the youngest state in history.
        Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"SolutionHistory::addState");
        *out << "Warning, state is younger than youngest state in history.  "
             << "State not added!\n" << std::endl;
        return;
      }
      break;
    }
    case STORAGE_TYPE_UNLIMITED:
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - unknown storage type.\n");
    }
  }

  // Add new state in chronological order.
  if (history_->size() == 0) {
    history_->push_back(state_);
  } else {
    typename std::vector<Teuchos::RCP<SolutionState<Scalar> > >::iterator
      state_it = history_->begin();
    for (; state_it < history_->end(); state_it++) {
      if (state_->getTime() < (*state_it)->getTime()) break;
    }
    history_->insert(state_it, state_);
  }

  if      (getNumStates() > 1) currentState_ = (*history_)[getNumStates()-2];
  else if (getNumStates() == 1) currentState_ = (*history_)[0];
  TEUCHOS_TEST_FOR_EXCEPTION(getNumStates() <= 0, std::logic_error,
    "Error - SolutionHistory::addState() Invalid history size!\n");

  return;
}

template<class Scalar>
void SolutionHistory<Scalar>::removeState(
  const Teuchos::RCP<SolutionState<Scalar> >& state_)
{
  if (history_->size() != 0) {
    auto state_it = history_->rbegin();
    for ( ; state_it < history_->rend(); state_it++) {
      if (state_->getTime() == (*state_it)->getTime()) break;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(state_it == history_->rend(), std::logic_error,
      "Error - removeState() Could not remove state = "
      // << state_it->describe()
      );

    // Need to be careful when erasing a reverse iterator.
    history_->erase(std::next(state_it).base());
  }
  return;
}


template<class Scalar>
void SolutionHistory<Scalar>::removeState(const Scalar time)
{
  Teuchos::RCP<SolutionState<Scalar> > tmpState = findState(time);
  removeState(tmpState);
}


template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::findState(const Scalar time) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(minTime() <= time and time <= maxTime()), std::logic_error,
    "Error - SolutionHistory::findState() Requested time out of range!\n"
    "        [Min, Max] = [" << minTime() << ", " << maxTime() << "]\n"
    "        time = "<< time <<"\n");

  const Scalar relTol = 1.0e-14;
  auto state_it = history_->begin();
  // Linear search
  for ( ; state_it < history_->end(); ++state_it) {
    if (std::abs((*state_it)->getTime()-time)/((*state_it)->getTime()) < relTol)
      break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(state_it == history_->end(), std::logic_error,
    "Error - SolutionHistory::findState()!\n"
    "        Did not find a SolutionState with time = " <<time<< std::endl);

  return *state_it;
}


template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::interpolateState(const Scalar time) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(minTime() <= time and time <= maxTime()), std::logic_error,
    "Error - SolutionHistory::getTime() Requested time out of range!\n"
    "        [Min, Max] = [" << minTime() << ", " << maxTime() << "]\n"
    "        time = "<< time <<"\n");

  // Interpolate the state.
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
    "Error - Implement interpolation.\n");
  //SolutionState<Scalar> state_out;
  //interpolate<Scalar>(*interpolator, history_, time_vec, &state_out);
}


/** Initialize the working state */
template<class Scalar>
void SolutionHistory<Scalar>::initWorkingState()
{
  TEMPUS_FUNC_TIME_MONITOR("Tempus::SolutionHistory::initWorkingState()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(getCurrentState() == Teuchos::null,
      std::logic_error,
      "Error - SolutionHistory::initWorkingState()\n"
      "Can not initialize working state without a current state!\n");

    // If workingState_ has a valid pointer, we are still working on it,
    // i.e., step failed and trying again, so do not initialize it.
    if (getWorkingState() != Teuchos::null) return;

    Teuchos::RCP<SolutionState<Scalar> > newState;
    if (getNumStates() < storageLimit) {
      // Create newState which is duplicate of currentState_
      newState = getCurrentState()->clone();
    } else {
      // Recycle old state and duplicate currentState_
      newState = (*history_)[0];
      history_->erase(history_->begin());
      newState->copy(getCurrentState());
      // When using the Griewank algorithm, we will want to select which
      // older state to recycle.
    }

    // Add newState and sort
    addState(newState);

    // Set workingState_
    workingState_ = (*history_)[getNumStates()-1];

    getWorkingState()->getMetaData()->setSolutionStatus(Status::WORKING);
  }

  return;
}


template<class Scalar>
void SolutionHistory<Scalar>::promoteWorkingState()
{
  Teuchos::RCP<SolutionStateMetaData<Scalar> > md =
    getWorkingState()->getMetaData();
  md->setTime(md->getTime() + md->getDt());
  md->setIStep(md->getIStep()+1);
  md->setNFailures(std::max(0,md->getNFailures()-1));
  md->setNConsecutiveFailures(0);
  md->setSolutionStatus(Status::PASSED);
  //md->setIsSynced(true);
  md->setIsInterpolated(false);
  currentState_ = workingState_;
  workingState_ = Teuchos::null;
}


template<class Scalar>
void SolutionHistory<Scalar>::setStorageLimit(int storage_limit)
{
  storageLimit = std::max(1,storage_limit);

  TEUCHOS_TEST_FOR_EXCEPTION(
    (Teuchos::as<int>(history_->size()) > storageLimit), std::logic_error,
    "Error - requested storage limit = " << storageLimit
    << " is smaller than the current number of states stored = "
    << history_->size() << "!\n");
}


template<class Scalar>
std::string SolutionHistory<Scalar>::description() const
{
  std::string name = "Tempus::SolutionHistory";
  return(name);
}


template<class Scalar>
void SolutionHistory<Scalar>::describe(
  Teuchos::FancyOStream          &out,
  const Teuchos::EVerbosityLevel verbLevel) const
{
  if ((Teuchos::as<int>(verbLevel)==Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel)>=Teuchos::as<int>(Teuchos::VERB_LOW)    )  ){
    out << description() << "::describe" << std::endl;
    //out << "interpolator     = " << interpolator->description() << std::endl;
    out << "storageLimit     = " << storageLimit << std::endl;
    out << "storageType      = " << storageType << std::endl;
    out << "number of states = " << history_->size() << std::endl;
    out << "time range       = (" << history_->front()->getTime() << ", "
                                  << history_->back()->getTime() << ")"
                                  << std::endl;
  } else if (Teuchos::as<int>(verbLevel) >=
             Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "SolutionStates: " << std::endl;
    for (int i=0; i<(int)history_->size() ; ++i) {
      out << "SolutionState[" << i << "] = " << std::endl;
      (*history_)[i]->describe(out,this->getVerbLevel());
    }
  }
}


template <class Scalar>
void SolutionHistory<Scalar>::setParameterList(
  Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  if (pList == Teuchos::null) *shPL_ = *(this->getValidParameters());
  else shPL_ = pList;
  shPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  //interpolator  = Teuchos::null;
  //setInterpolator(interpolator);

  storageType = StorageTypeValidator->getIntegralValue(
    *shPL_, Storage_name, Storage_default);

  int storage_limit = shPL_->get(StorageLimit_name, StorageLimit_default);

  switch (storageType) {
  case STORAGE_TYPE_INVALID:
  case STORAGE_TYPE_KEEP_NEWEST: {
    storageType = STORAGE_TYPE_KEEP_NEWEST;
    if (storage_limit != 1) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"SolutionHistory::setParameterList");
      *out << "Warning - 'Storage Limit' for 'Keep Newest' is 1.\n"
           << "  (Storage Limit = "<<storage_limit<<").  Resetting to 1."
           << std::endl;
      storage_limit = 1;
    }
    setStorageLimit(storage_limit);
    break;
  }
  case STORAGE_TYPE_UNDO: {
    if (storage_limit != 2) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"SolutionHistory::setParameterList");
      *out << "Warning - 'Storage Limit' for 'Undo' is 2.\n"
           << "  (Storage Limit = "<<storage_limit<<").  Resetting to 2."
           << std::endl;
      storage_limit = 2;
    }
    setStorageLimit(storage_limit);
    break;
  }
  case STORAGE_TYPE_STATIC: {
    break;
  }
  case STORAGE_TYPE_UNLIMITED: {
    storage_limit = std::numeric_limits<int>::max();
    break;
  }
  }
  setStorageLimit(storage_limit);
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
SolutionHistory<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

  pl->set(Storage_name, Storage_default,
    "'Storage Type' sets the memory storage.  "
    "'Keep Newest' - will retain the single newest solution state.  "
    "'Undo' - will retain two solution states in order to do a single undo.  "
    "'Static' - will retain 'Storage Limit' number of solution states.  "
    "'Unlimited' - will not remove any solution states!",
    StorageTypeValidator);

  pl->set(StorageLimit_name, StorageLimit_default,
    "Storage limit for the solution history.");

  return pl;
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
SolutionHistory<Scalar>::getNonconstParameterList()
{
  return(shPL_);
}


template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
SolutionHistory<Scalar>::unsetParameterList()
{
  Teuchos::RCP<Teuchos::ParameterList> temp_plist = shPL_;
  shPL_ = Teuchos::null;
  return(temp_plist);
}

// Nonmember constructor.
template<class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > solutionHistory(
  Teuchos::RCP<Teuchos::ParameterList> pList)
{
  Teuchos::RCP<SolutionHistory<Scalar> > sh=rcp(new SolutionHistory<Scalar>(pList));
  return sh;
}

//template<class Scalar>
//void SolutionHistory<Scalar>::setInterpolator(
//  const Teuchos::RCP<InterpolatorBase<Scalar> >& interpolator_)
//{
//  if (interpolator_ == Teuchos::null) {
//    interpolator = linearInterpolator<Scalar>();
//  } else {
//    interpolator = interpolator_;
//  }
//  if (Teuchos::as<int>(this->getVerbLevel()) >=
//      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
//    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
//    Teuchos::OSTab ostab(out,1,"SolutionHistory::setInterpolator");
//    *out << "interpolator = " << interpolator_->description() << std::endl;
//  }
//}
//
//template<class Scalar>
//Teuchos::RCP<InterpolatorBase<Scalar> >
//  SolutionHistory<Scalar>::getNonconstInterpolator()
//{
//  return interpolator;
//}
//
//template<class Scalar>
//Teuchos::RCP<const InterpolatorBase<Scalar> >
//  SolutionHistory<Scalar>::getInterpolator() const
//{
//  return interpolator;
//}
//
//template<class Scalar>
//Teuchos::RCP<InterpolatorBase<Scalar> > SolutionHistory<Scalar>::unSetInterpolator()
//{
//  Teuchos::RCP<InterpolatorBase<Scalar> > old_interpolator = interpolator;
//  interpolator = linearInterpolator<Scalar>();
//  return old_interpolator;
//}


} // namespace Tempus
#endif // Tempus_SolutionHistory_impl_hpp
