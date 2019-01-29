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
#include "Tempus_InterpolatorFactory.hpp"

//#include "Thyra_VectorStdOps.hpp"


namespace {

  static std::string Invalid_name    = "Invalid";
  static std::string KeepNewest_name = "Keep Newest";
  static std::string Undo_name       = "Undo";
  static std::string Static_name     = "Static";
  static std::string Unlimited_name  = "Unlimited";
  static std::string Storage_name    = "Storage Type";
  static std::string Storage_default = Undo_name;

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
  const Teuchos::RCP<SolutionState<Scalar> >& state)
{
  // Check that we're not going to exceed our storage limit:
  if (Teuchos::as<int>(history_->size()+1) > storageLimit_) {
    switch (storageType_) {
    case STORAGE_TYPE_INVALID: {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Error - Storage type is STORAGE_TYPE_INVALID.\n");
      break;
    }
    case STORAGE_TYPE_STATIC:
    case STORAGE_TYPE_KEEP_NEWEST:
    case STORAGE_TYPE_UNDO: {
      if (state->getTime() >= history_->front()->getTime()) {
        // Case:  State is older than the youngest state in history.
        // Remove state from the beginning of history, then add new state.
        history_->erase(history_->begin());
      } else {
        // Case:  State is younger than the youngest state in history.
        Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out,1,"SolutionHistory::addState");
        *out << "Warning, state is younger than youngest state in history.  "
             << "State not added!" << std::endl;
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
    history_->push_back(state);
  } else {
    typename std::vector<Teuchos::RCP<SolutionState<Scalar> > >::iterator
      state_it = history_->begin();
    for (; state_it < history_->end(); state_it++) {
      if (state->getTime() < (*state_it)->getTime()) break;
    }
    history_->insert(state_it, state);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(getNumStates() <= 0, std::logic_error,
    "Error - SolutionHistory::addState() Invalid history size!\n");

  return;
}

template<class Scalar>
void SolutionHistory<Scalar>::addWorkingState(
  const Teuchos::RCP<SolutionState<Scalar> >& state, const bool updateTime)
{
  using Teuchos::RCP;

  addState(state);
  workingState_ = (*history_)[getNumStates()-1];
  RCP<SolutionStateMetaData<Scalar> > csmd = getCurrentState()->getMetaData();
  RCP<SolutionStateMetaData<Scalar> > wsmd = workingState_    ->getMetaData();
  wsmd->setSolutionStatus(Status::WORKING);
  wsmd->setIStep(csmd->getIStep()+1);
  if (updateTime) {
    wsmd->setTime(csmd->getTime() + csmd->getDt());
    wsmd->setDt(csmd->getDt());
  }
}

template<class Scalar>
void SolutionHistory<Scalar>::removeState(
  const Teuchos::RCP<SolutionState<Scalar> >& state)
{
  if (history_->size() != 0) {
    auto state_it = history_->rbegin();
    for ( ; state_it < history_->rend(); state_it++) {
      if (state->getTime() == (*state_it)->getTime()) break;
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

  // Use last step in solution history as the scale for comparing times
  const Scalar scale =
    history_->size() > 0 ? (*history_)[history_->size()-1]->getTime() : Scalar(1.0);
  // Linear search
  auto state_it = history_->begin();
  for ( ; state_it < history_->end(); ++state_it) {
    if (floating_compare_equals((*state_it)->getTime(),time,scale))
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
  Teuchos::RCP<SolutionState<Scalar> > state_out = getCurrentState()->clone();
  interpolate<Scalar>(*interpolator_, history_, time, state_out.get());
  return state_out;
}


template<class Scalar>
void
SolutionHistory<Scalar>::interpolateState(
  const Scalar time, SolutionState<Scalar>* state_out) const
{
  interpolate<Scalar>(*interpolator_, history_, time, state_out);
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
    if (getWorkingState(false) != Teuchos::null) return;

    Teuchos::RCP<SolutionState<Scalar> > newState;
    if (getNumStates() < storageLimit_) {
      // Create newState which is duplicate of currentState
      newState = getCurrentState()->clone();
    } else {
      // Recycle old state and copy currentState
      newState = (*history_)[0];
      history_->erase(history_->begin());
      if (getNumStates() > 0) newState->copy(getCurrentState());
      // When using the Griewank algorithm, we will want to select which
      // older state to recycle.
    }

    addWorkingState(newState);

  }
  return;
}


template<class Scalar>
void SolutionHistory<Scalar>::promoteWorkingState()
{
  Teuchos::RCP<SolutionStateMetaData<Scalar> > md =
    getWorkingState()->getMetaData();

  if ( md->getSolutionStatus() == Status::PASSED ) {
    md->setNFailures(std::max(0,md->getNFailures()-1));
    md->setNConsecutiveFailures(0);
    md->setSolutionStatus(Status::PASSED);
    //md->setIsSynced(true);
    md->setIsInterpolated(false);
    workingState_ = Teuchos::null;
  } else {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out,1,"SolutionHistory::promoteWorkingState()");
    *out << "Warning - WorkingState is not passing, so not promoted!\n"
         << std::endl;
  }
}


template<class Scalar>
void SolutionHistory<Scalar>::setStorageLimit(int storage_limit)
{
  storageLimit_ = std::max(1,storage_limit);

  TEUCHOS_TEST_FOR_EXCEPTION(
    (Teuchos::as<int>(history_->size()) > storageLimit_), std::logic_error,
    "Error - requested storage limit = " << storageLimit_
    << " is smaller than the current number of states stored = "
    << history_->size() << "!\n");
}


template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::getStateTimeIndexN() const
{
  const int m = history_->size();
  TEUCHOS_TEST_FOR_EXCEPTION( (m < 1), std::out_of_range,
    "Error - getStateTimeIndexN() No states in SolutionHistory!\n");
  return (*history_)[m-1];
}


template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::getStateTimeIndexNM1() const
{
  const int m   = history_->size();
  TEUCHOS_TEST_FOR_EXCEPTION( (m < 2), std::out_of_range,
    "Error - getStateTimeIndexNM1() Not enough states in "
    << "SolutionHistory!\n");
  const int n   = (*history_)[m-1]->getIndex();
  const int nm1 = (*history_)[m-2]->getIndex();

  // No need to search SolutionHistory as states n and nm1 should be
  // next to each other.
  TEUCHOS_TEST_FOR_EXCEPTION( (nm1 != n-1), std::out_of_range,
    "Error - getStateTimeIndexNM1() Timestep index n-1 is not in "
    << "SolutionHistory!\n"
    << "    (n)th   index = " << n << "\n"
    << "    (n-1)th index = " << nm1 << "\n");

  return (*history_)[m-2];
}


template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::getStateTimeIndexNM2() const
{
  const int m   = history_->size();
  TEUCHOS_TEST_FOR_EXCEPTION( (m < 3), std::out_of_range,
    "Error - getStateTimeIndexNM1() Not enough states in "
    << "SolutionHistory!\n");
  const int n   = (*history_)[m-1]->getIndex();
  const int nm2 = (*history_)[m-3]->getIndex();

  // Assume states n and nm2 are one away from each other.
  // May need to do a search otherwise.
  TEUCHOS_TEST_FOR_EXCEPTION( (nm2 != n-2), std::out_of_range,
    "Error - getStateTimeIndexNM1() Timestep index n-2 is not in "
    << "SolutionHistory!\n"
    << "    (n)th   index = " << n << "\n"
    << "    (n-2)th index = " << nm2 << "\n");

  return (*history_)[m-3];
}


template<class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::getStateTimeIndex(int index) const
{
  typename std::vector<Teuchos::RCP<SolutionState<Scalar> > >::iterator
    state_it = history_->begin();
  for (; state_it < history_->end(); state_it++) {
    if ((*state_it)->getIndex() == index) break;
  }
  TEUCHOS_TEST_FOR_EXCEPTION( state_it==history_->end(), std::out_of_range,
    "Error - getStateTimeIndex() Timestep index is not in "
    << "SolutionHistory!\n"
    << "    index = " << index << "\n");
  return (*state_it);
}


template<class Scalar>
std::string SolutionHistory<Scalar>::description() const
{
  return ("Tempus::SolutionHistory - name = '" + name_ + "'");
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
    out << "storageLimit     = " << storageLimit_ << std::endl;
    out << "storageType      = " << storageType_ << std::endl;
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
  if (pList == Teuchos::null) {
    // Create default parameters if null, otherwise keep current parameters.
    if (shPL_ == Teuchos::null) {
      shPL_ = Teuchos::parameterList("Solution History");
      *shPL_ = *(this->getValidParameters());
    }
  } else {
    shPL_ = pList;
  }
  shPL_->validateParametersAndSetDefaults(*this->getValidParameters());

  name_ = shPL_->name();

  storageType_ = StorageTypeValidator->getIntegralValue(
    *shPL_, Storage_name, Storage_default);

  int storage_limit = shPL_->get(StorageLimit_name, StorageLimit_default);

  switch (storageType_) {
  case STORAGE_TYPE_INVALID:
  case STORAGE_TYPE_KEEP_NEWEST: {
    storageType_ = STORAGE_TYPE_KEEP_NEWEST;
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
    storage_limit = 1000000000;
    break;
  }
  }
  setStorageLimit(storage_limit);

  interpolator_ = InterpolatorFactory<Scalar>::createInterpolator(
    Teuchos::sublist(shPL_, "Interpolator"));
}


template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
SolutionHistory<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();

  pl->setName("Valid ParameterList");

  pl->set(Storage_name, Storage_default,
    "'Storage Type' sets the memory storage.  "
    "'Keep Newest' - will retain the single newest solution state.  "
    "'Undo' - will retain two solution states in order to do a single undo.  "
    "'Static' - will retain 'Storage Limit' number of solution states.  "
    "'Unlimited' - will not remove any solution states!",
    StorageTypeValidator);

  pl->set(StorageLimit_name, StorageLimit_default,
    "Storage limit for the solution history.");

  // Interpolator
  pl->sublist("Interpolator",false,"").disableRecursiveValidation();

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

template<class Scalar>
void SolutionHistory<Scalar>::setInterpolator(
 const Teuchos::RCP<Interpolator<Scalar> >& interpolator)
{
 if (interpolator == Teuchos::null) {
   interpolator_ = InterpolatorFactory<Scalar>::createInterpolator();
 } else {
   interpolator_ = interpolator;
 }
 if (Teuchos::as<int>(this->getVerbLevel()) >=
     Teuchos::as<int>(Teuchos::VERB_HIGH)) {
   Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
   Teuchos::OSTab ostab(out,1,"SolutionHistory::setInterpolator");
   *out << "interpolator = " << interpolator_->description() << std::endl;
 }
}

template<class Scalar>
Teuchos::RCP<Interpolator<Scalar> >
SolutionHistory<Scalar>::getNonconstInterpolator()
{
 return interpolator_;
}

template<class Scalar>
Teuchos::RCP<const Interpolator<Scalar> >
SolutionHistory<Scalar>::getInterpolator() const
{
 return interpolator_;
}

template<class Scalar>
Teuchos::RCP<Interpolator<Scalar> >
SolutionHistory<Scalar>::unSetInterpolator()
{
 Teuchos::RCP<Interpolator<Scalar> > old_interpolator = interpolator_;
 interpolator_ = lagrangeInterpolator<Scalar>();
 return old_interpolator;
}


} // namespace Tempus
#endif // Tempus_SolutionHistory_impl_hpp
