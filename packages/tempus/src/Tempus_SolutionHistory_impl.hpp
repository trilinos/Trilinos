#ifndef TEMPUS_SOLUTIONHISTORY_IMPL_HPP
#define TEMPUS_SOLUTIONHISTORY_IMPL_HPP

// Teuchos
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
// Tempus
#include "Tempus_SolutionStateMetaData.hpp"

//#include "Thyra_VectorStdOps.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;

namespace {

  static std::string Invalid_name    = "Invalid";
  static std::string KeepNewest_name = "Keep Newest";
  static std::string Undo_name       = "Undo";
  static std::string Static_name     = "Static";
  static std::string Unlimited_name  = "Unlimited";
  static std::string Storage_name    = "Storage Type";
  static std::string Storage_default = KeepNewest_name;

  static std::string StorageLimit_name    = "Storage Limit";
  static int         StorageLimit_default = 1;

  Teuchos::Array<std::string>
    HistoryPolicies = Teuchos::tuple<std::string>(
      Invalid_name,
      KeepNewest_name,
      Undo_name,
      Static_name,
      Unlimited_name);

  const RCP<Teuchos::StringToIntegralParameterEntryValidator<Tempus::StorageType> >
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
SolutionHistory<Scalar>::SolutionHistory(RCP<ParameterList> pList_)
{
  // Create history, an array of solution states.
  history = rcp(new Teuchos::Array<RCP<SolutionState<Scalar> > >);

  if (pList_ == Teuchos::null) {
    pList     = Teuchos::null;
    //interpolator  = Teuchos::null;
    storageLimit = StorageLimit_default;
    storageType = STORAGE_TYPE_KEEP_NEWEST;
  } else {
    this->setParameterList(pList_);
  }

  if (Teuchos::as<int>(this->getVerbLevel()) >=
      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"SolutionHistory::SolutionHistory");
    *out << this->description() << std::endl;
  }
}


template<class Scalar>
void SolutionHistory<Scalar>::addState(
  const RCP<SolutionState<Scalar> >& state_)
{
  // Check that we're not going to exceed our storage limit:
  if (Teuchos::as<int>(history->size()+1) > storageLimit) {
    switch (storageType) {
    case STORAGE_TYPE_INVALID: {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - Storage type is STORAGE_TYPE_INVALID.\n");
      break;
    }
    case STORAGE_TYPE_STATIC:
    case STORAGE_TYPE_KEEP_NEWEST:
    case STORAGE_TYPE_UNDO: {
      if (state_->getTime() > history->front()->getTime()) {
        // Case:  State is older than the youngest state in history.
        // Remove state from the beginning of history, then add new state.
        history->erase(history->begin());
      } else {
        // Case:  State is younger than the youngest state in history.
        RCP<Teuchos::FancyOStream> out = this->getOStream();
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
  if (history->size() == 0) {
    history->push_back(state_);
  } else {
    typename Teuchos::Array<RCP<SolutionState<Scalar> > >::iterator
      state_it = history->begin();
    for (; state_it < history->end(); state_it++) {
      if (state_->getTime() < (*state_it)->getTime()) break;
    }
    history->insert(state_it, state_);
  }
  currentState = history->back();

  return;
}

template<class Scalar>
void SolutionHistory<Scalar>::removeState(
  const RCP<SolutionState<Scalar> >& state_)
{
  if (history->size() != 0) {
    typename Teuchos::Array<SolutionState<Scalar> >::reverse_iterator
      state_it = history->rbegin();
    for (state_it; state_it < history->rend(); state_it++) {
      if (state_->getTime() == (*state_it)->getTime()) break;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(state_it == history->rend(), std::logic_error,
      "Error - removeState() Could not remove state = "
      // << state_it->describe()
      );

    history->erase(state_it);
  }
  return;
}


template<class Scalar>
void SolutionHistory<Scalar>::removeState(const Scalar time)
{
  RCP<SolutionState<Scalar> > tmpState = findState(time);
  removeState(tmpState);
}


template<class Scalar>
RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::findState(const Scalar time) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(minTime() <= time and time <= maxTime()), std::logic_error,
    "Error - SolutionHistory::findState() Requested time out of range!\n"
    "        [Min, Max] = [" << minTime() << ", " << maxTime() << "]\n"
    "        time = "<< time <<"\n");

  const Scalar relTol = 1.0e-14;
  typename Teuchos::Array<SolutionState<Scalar> >::iterator
    state_it = history->begin();
  // Linear search
  for (state_it; state_it < history->end(); ++state_it) {
    if (abs((state_it.getTime()-time)/state_it.getTime()) < relTol) break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(state_it == history->end(), std::logic_error,
    "Error - SolutionHistory::findState()!\n"
    "        Did not find a SolutionState with time = " <<time<< std::endl);

  return history[state_it];
}


template<class Scalar>
RCP<SolutionState<Scalar> >
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
  //interpolate<Scalar>(*interpolator, history, time_vec, &state_out);
}


template<class Scalar>
RCP<SolutionState<Scalar> > SolutionHistory<Scalar>::getCurrentState() const
{
  return currentState;
}


template<class Scalar>
RCP<SolutionState<Scalar> > SolutionHistory<Scalar>::getWorkingState() const
{
  return workingState;
}

/** Initialize the working state */
template<class Scalar>
void SolutionHistory<Scalar>::initWorkingState()
{
  TEUCHOS_TEST_FOR_EXCEPTION(currentState == Teuchos::null, std::logic_error,
    "Error - SolutionHistory::initWorkingState()\n"
    "Can not initialize working state without a current state!\n");

  // If workingState has a valid pointer, we are still working on it,
  // i.e., step failed and trying again.  So do not initialize it.
  if (workingState != Teuchos::null) return;

  if (storageLimit == 1) {
    workingState = currentState;
    workingState->metaData->solutionStatus = Status::WORKING;
    currentState = Teuchos::null;
  } else {
    workingState = currentState->clone();
    workingState->metaData->solutionStatus = Status::WORKING;
    addState(workingState);
  }

  return;
}


template<class Scalar>
void SolutionHistory<Scalar>::promoteWorkingState()
{
  RCP<SolutionStateMetaData<Scalar> > md = workingState->metaData;
  md->time += md->dt;
  md->iStep++;
  md->nFailures = std::max(0,md->nFailures-1);
  md->nConsecutiveFailures = 0;
  md->solutionStatus = Status::PASSED;
  md->isRestartable = true;
  md->isInterpolated = false;
  currentState = workingState;
  workingState = Teuchos::null;
}


template<class Scalar>
void SolutionHistory<Scalar>::setStorageLimit(int storage_limit)
{
  storageLimit = std::max(1,storage_limit);

  TEUCHOS_TEST_FOR_EXCEPTION(
    (Teuchos::as<int>(history->size()) > storageLimit), std::logic_error,
    "Error - requested storage limit = " << storageLimit
    << " is smaller than the current number of states stored = "
    << history->size() << "!\n");
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
    out << "number of states = " << history->size() << std::endl;
    out << "time range       = (" << history->front()->getTime() << ", "
                                  << history->back()->getTime() << ")"
                                  << std::endl;
  } else if (Teuchos::as<int>(verbLevel) >=
             Teuchos::as<int>(Teuchos::VERB_HIGH)) {
    out << "SolutionStates: " << std::endl;
    for (Teuchos::Ordinal i=0; i<history->size() ; ++i) {
      out << "SolutionState[" << i << "] = " << std::endl;
      (*history)[i]->describe(out,this->getVerbLevel());
    }
  }
}


template <class Scalar>
void SolutionHistory<Scalar>::setParameterList(
  RCP<ParameterList> const& pList_)
{
  TEUCHOS_TEST_FOR_EXCEPT(is_null(pList_));
  pList_->validateParameters(*this->getValidParameters());
  pList = pList_;

  Teuchos::readVerboseObjectSublist(&*pList,this);

  //setInterpolator(interpolator);

  storageType = StorageTypeValidator->getIntegralValue(
    *pList, Storage_name, Storage_default);

  int storage_limit = pList->get(StorageLimit_name, StorageLimit_default);

  switch (storageType) {
  case STORAGE_TYPE_INVALID:
  case STORAGE_TYPE_KEEP_NEWEST: {
    storageType = STORAGE_TYPE_KEEP_NEWEST;
    if (storage_limit != 1) {
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"SoltuionHistory::setParameterList");
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
      RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out,1,"SoltuionHistory::setParameterList");
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
RCP<const ParameterList> SolutionHistory<Scalar>::getValidParameters() const
{
  static RCP<ParameterList> validPL;

  if (is_null(validPL)) {

    RCP<ParameterList> pl = Teuchos::parameterList();
    Teuchos::setupVerboseObjectSublist(&*pl);

    pl->set(Storage_name, Storage_default,
      "'Storage Type' sets the memory storage.  "
      "'Keep Newest' - will retain the single newest solution state.  "
      "'Undo' - will retain two solution states in order to do a single undo.  "
      "'Static' - will retain 'Storage Limit' number of solution states.  "
      "'Unlimited' - will not remove any solution states!",
      StorageTypeValidator);

    pl->set(StorageLimit_name, StorageLimit_default,
      "Storage limit for the solution history.");

    validPL = pl;

  }
  return validPL;
}


template <class Scalar>
RCP<ParameterList>
SolutionHistory<Scalar>::getNonconstParameterList()
{
  return(pList);
}


template <class Scalar>
RCP<ParameterList> SolutionHistory<Scalar>::unsetParameterList()
{
  RCP<ParameterList> temp_plist = pList;
  pList = Teuchos::null;
  return(temp_plist);
}


//template<class Scalar>
//void SolutionHistory<Scalar>::setInterpolator(
//  const RCP<InterpolatorBase<Scalar> >& interpolator_)
//{
//  if (interpolator_ == Teuchos::null) {
//    interpolator = linearInterpolator<Scalar>();
//  } else {
//    interpolator = interpolator_;
//  }
//  if (Teuchos::as<int>(this->getVerbLevel()) >=
//      Teuchos::as<int>(Teuchos::VERB_HIGH)) {
//    RCP<Teuchos::FancyOStream> out = this->getOStream();
//    Teuchos::OSTab ostab(out,1,"SolutionHistory::setInterpolator");
//    *out << "interpolator = " << interpolator_->description() << std::endl;
//  }
//}
//
//template<class Scalar>
//RCP<InterpolatorBase<Scalar> >
//  SolutionHistory<Scalar>::getNonconstInterpolator()
//{
//  return interpolator;
//}
//
//template<class Scalar>
//RCP<const InterpolatorBase<Scalar> >
//  SolutionHistory<Scalar>::getInterpolator() const
//{
//  return interpolator;
//}
//
//template<class Scalar>
//RCP<InterpolatorBase<Scalar> > SolutionHistory<Scalar>::unSetInterpolator()
//{
//  RCP<InterpolatorBase<Scalar> > old_interpolator = interpolator;
//  interpolator = linearInterpolator<Scalar>();
//  return old_interpolator;
//}


} // namespace Tempus
#endif // TEMPUS_SOLUTIONHISTORY_IMPL_HPP
