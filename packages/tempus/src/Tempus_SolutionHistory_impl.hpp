//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_SolutionHistory_impl_hpp
#define Tempus_SolutionHistory_impl_hpp

#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_InterpolatorFactory.hpp"
#include "Tempus_NumericalUtils.hpp"

namespace Tempus {

template <class Scalar>
SolutionHistory<Scalar>::SolutionHistory()
  : name_("Solution History"), storageType_(STORAGE_TYPE_UNDO), storageLimit_(2)
{
  using Teuchos::RCP;
  // Create history, a vector of solution states.
  history_       = rcp(new std::vector<RCP<SolutionState<Scalar> > >);
  interpolator_  = InterpolatorFactory<Scalar>::createInterpolator();
  isInitialized_ = false;
}

template <class Scalar>
SolutionHistory<Scalar>::SolutionHistory(
    std::string name,
    Teuchos::RCP<std::vector<Teuchos::RCP<SolutionState<Scalar> > > > history,
    Teuchos::RCP<Interpolator<Scalar> > interpolator, StorageType storageType,
    int storageLimit)
{
  setName(name);
  setHistory(history);
  setInterpolator(interpolator);
  setStorageType(storageType);
  setStorageLimit(storageLimit);

  isInitialized_ = false;
  if (getNumStates() > 0) isInitialized_ = true;
}

template <class Scalar>
void SolutionHistory<Scalar>::addState(
    const Teuchos::RCP<SolutionState<Scalar> >& state, bool doChecks)
{
  // Check that we're not going to exceed our storage limit:
  if (Teuchos::as<int>(history_->size() + 1) > storageLimit_) {
    switch (storageType_) {
      case STORAGE_TYPE_INVALID: {
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error,
            "Error - Storage type is STORAGE_TYPE_INVALID.\n");
        break;
      }
      case STORAGE_TYPE_STATIC:
      case STORAGE_TYPE_KEEP_NEWEST:
      case STORAGE_TYPE_UNDO: {
        if (state->getTime() >= history_->front()->getTime()) {
          // Case:  State is younger than the oldest state in history.
          // Remove state from the beginning of history, then add new state.
          history_->erase(history_->begin());
        }
        else {
          // Case:  State is older than the oldest state in history.
          Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out, 1, "SolutionHistory::addState");
          *out << "Warning, state is older than oldest state in history.  "
               << "State not added!" << std::endl;
          return;
        }
        break;
      }
      case STORAGE_TYPE_UNLIMITED: break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "Error - unknown storage type.\n");
    }
  }

  // Add new state in chronological order.
  if (history_->size() == 0) {
    history_->push_back(state);
  }
  else {
    typename std::vector<Teuchos::RCP<SolutionState<Scalar> > >::iterator
        state_it         = history_->begin();
    bool equal           = false;
    const Scalar newTime = state->getTime();
    for (; state_it < history_->end(); state_it++) {
      const Scalar shTime = (*state_it)->getTime();
      if (doChecks) {
        const Scalar denom = std::max(std::fabs(shTime), std::fabs(newTime));
        if (std::fabs(newTime - shTime) < 1.0e-14 * denom) {
          equal                                   = true;
          Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
          Teuchos::OSTab ostab(out, 1, "SolutionHistory::addState");
          *out << "Warning, state to be added matches state in history.  "
               << "State not added!" << std::endl;

          *out << "===============" << std::endl;
          *out << "Added SolutionState -- ";
          (*state_it)->describe(*out, Teuchos::VERB_MEDIUM);
          *out << "===============" << std::endl;
          this->describe(*out, Teuchos::VERB_MEDIUM);

          std::exit(-1);
          break;
        }
      }
      if (newTime < shTime) break;
    }
    if (!equal) history_->insert(state_it, state);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      getNumStates() <= 0, std::logic_error,
      "Error - SolutionHistory::addState() Invalid history size!\n");

  return;
}

template <class Scalar>
void SolutionHistory<Scalar>::addWorkingState(
    const Teuchos::RCP<SolutionState<Scalar> >& state, const bool updateTime)
{
  using Teuchos::RCP;

  auto cs = getCurrentState();
  state->setSolutionStatus(Status::WORKING);
  state->setIndex(cs->getIndex() + 1);
  if (updateTime) {
    state->setTime(cs->getTime() + cs->getTimeStep());
    state->setTimeStep(cs->getTimeStep());
  }

  addState(state);
  workingState_ = (*history_)[getNumStates() - 1];
}

template <class Scalar>
void SolutionHistory<Scalar>::removeState(
    const Teuchos::RCP<SolutionState<Scalar> >& state)
{
  if (history_->size() != 0) {
    auto state_it = history_->rbegin();
    for (; state_it < history_->rend(); state_it++) {
      if (state->getTime() == (*state_it)->getTime()) break;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(state_it == history_->rend(), std::logic_error,
                               "Error - removeState() Could not remove state = "
                                   << (*state_it)->description());

    // Need to be careful when erasing a reverse iterator.
    history_->erase(std::next(state_it).base());
  }
  return;
}

template <class Scalar>
void SolutionHistory<Scalar>::removeState(const Scalar time)
{
  Teuchos::RCP<SolutionState<Scalar> > tmpState = findState(time);
  removeState(tmpState);
}

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> > SolutionHistory<Scalar>::findState(
    const Scalar time) const
{
  // Use last step in solution history as the scale for comparing times
  Scalar scale = 1.0;
  if (history_->size() > 0)
    scale = (*history_)[history_->size() - 1]->getTime();
  if (approxZero(scale)) scale = Scalar(1.0);

  const Scalar absTol = scale * numericalTol<Scalar>();
  TEUCHOS_TEST_FOR_EXCEPTION(
      !(minTime() - absTol <= time && time <= maxTime() + absTol),
      std::logic_error,
      "Error - SolutionHistory::findState() Requested time out of range!\n"
          << "        [Min, Max] = ["
          << minTime() << ", " << maxTime()
          << "]\n        time = "
          << time << "\n");

  // Linear search
  auto state_it = history_->begin();
  for (; state_it < history_->end(); ++state_it) {
    if (approxEqualAbsTol(time, (*state_it)->getTime(), absTol)) break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(state_it == history_->end(), std::logic_error,
                             "Error - SolutionHistory::findState()!\n"
                             "        Did not find a SolutionState with time = "
                                 << time << std::endl);

  return *state_it;
}

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> > SolutionHistory<Scalar>::interpolateState(
    const Scalar time) const
{
  Teuchos::RCP<SolutionState<Scalar> > state_out = getCurrentState()->clone();
  interpolate<Scalar>(*interpolator_, history_, time, state_out.get());
  return state_out;
}

template <class Scalar>
void SolutionHistory<Scalar>::interpolateState(
    const Scalar time, SolutionState<Scalar>* state_out) const
{
  interpolate<Scalar>(*interpolator_, history_, time, state_out);
}

/** Initialize the working state */
template <class Scalar>
void SolutionHistory<Scalar>::initWorkingState()
{
  TEMPUS_FUNC_TIME_MONITOR("Tempus::SolutionHistory::initWorkingState()");
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        getCurrentState() == Teuchos::null, std::logic_error,
        "Error - SolutionHistory::initWorkingState()\n"
        "Can not initialize working state without a current state!\n");

    // If workingState_ has a valid pointer, we are still working on it,
    // i.e., step failed and trying again.  There are a couple options:
    //   1. Reuse the workingState as it might be a good guess.  This
    //      could help with performance.  This was the initial implementation.
    //   2. Reset the workingState to the last time step.  This could
    //      be more robust in the case of the workingState failing with nans.
    //      This is the current implementation.
    if (getWorkingState(false) != Teuchos::null) {
      Thyra::V_V(getWorkingState()->getX().ptr(), *(getCurrentState()->getX()));
      return;
    }

    Teuchos::RCP<SolutionState<Scalar> > newState;
    if (getNumStates() < storageLimit_) {
      // Create newState which is duplicate of currentState
      newState = getCurrentState()->clone();
    }
    else {
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

template <class Scalar>
void SolutionHistory<Scalar>::promoteWorkingState()
{
  auto ws = getWorkingState();

  if (ws->getSolutionStatus() == Status::PASSED) {
    ws->setNFailures(std::max(0, ws->getNFailures() - 1));
    ws->setNConsecutiveFailures(0);
    ws->setSolutionStatus(Status::PASSED);
    // ws->setIsSynced(true);
    ws->setIsInterpolated(false);
    workingState_ = Teuchos::null;
  }
  else {
    Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
    Teuchos::OSTab ostab(out, 1, "SolutionHistory::promoteWorkingState()");
    *out << "Warning - WorkingState is not passing, so not promoted!\n"
         << std::endl;
  }
}

template <class Scalar>
void SolutionHistory<Scalar>::copy(
    Teuchos::RCP<const SolutionHistory<Scalar> > sh)
{
  this->setName(sh->getName());

  this->clear();
  auto sh_history = sh->getHistory();
  typename std::vector<Teuchos::RCP<SolutionState<Scalar> > >::iterator
      state_it = sh_history->begin();
  for (; state_it < sh_history->end(); state_it++) this->addState(*state_it);

  auto interpolator =
      Teuchos::rcp_const_cast<Interpolator<Scalar> >(sh->getInterpolator());
  this->setInterpolator(interpolator);

  this->setStorageType(sh->getStorageType());
  this->setStorageLimit(sh->getStorageLimit());
}

template <class Scalar>
void SolutionHistory<Scalar>::setStorageLimit(int storage_limit)
{
  storageLimit_ = std::max(1, storage_limit);

  if (storageType_ == STORAGE_TYPE_INVALID ||
      storageType_ == STORAGE_TYPE_KEEP_NEWEST) {
    if (storage_limit != 1) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out, 1, "SolutionHistory::setStorageLimit");
      *out << "Warning - 'Storage Limit' for 'Keep Newest' is 1.\n"
           << "  (Storage Limit = " << storage_limit << ").  Resetting to 1."
           << std::endl;
      storageLimit_ = 1;
    }
  }
  else if (storageType_ == STORAGE_TYPE_UNDO) {
    if (storage_limit != 2) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out, 1, "SolutionHistory::setStorageLimit");
      *out << "Warning - 'Storage Limit' for 'Undo' is 2.\n"
           << "  (Storage Limit = " << storage_limit << ").  Resetting to 2."
           << std::endl;
      storageLimit_ = 2;
    }
  }
  else if (storageType_ == STORAGE_TYPE_STATIC) {
    storageLimit_ = storage_limit;
  }
  else if (storageType_ == STORAGE_TYPE_UNLIMITED) {
    storageLimit_ = std::numeric_limits<int>::max();
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      (Teuchos::as<int>(history_->size()) > storageLimit_), std::logic_error,
      "Error - requested storage limit = "
          << storageLimit_
          << " is smaller than the current number of states stored = "
          << history_->size() << "!\n");

  isInitialized_ = false;
}

template <class Scalar>
void SolutionHistory<Scalar>::setStorageType(StorageType st)
{
  storageType_ = st;
  if (storageType_ == STORAGE_TYPE_KEEP_NEWEST)
    setStorageLimit(1);
  else if (storageType_ == STORAGE_TYPE_UNDO)
    setStorageLimit(2);
  else if (storageType_ == STORAGE_TYPE_UNLIMITED)
    setStorageLimit(std::numeric_limits<int>::max());
  isInitialized_ = false;
}

template <class Scalar>
void SolutionHistory<Scalar>::setStorageTypeString(std::string s)
{
  if (s == "Keep Newest") {  // Keep the single newest state
    storageType_  = STORAGE_TYPE_KEEP_NEWEST;
    storageLimit_ = 1;
  }
  else if (s == "Undo") {  // Keep the 2 newest states for undo
    storageType_  = STORAGE_TYPE_UNDO;
    storageLimit_ = 2;
  }
  else if (s == "Static") {  // Keep a fix number of states
    storageType_ = STORAGE_TYPE_STATIC;
  }
  else if (s == "Unlimited") {  // Grow the history as needed
    storageType_ = STORAGE_TYPE_UNLIMITED;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "Error - Unknown 'Storage Type' = '" << s << "'\n");
  }
  isInitialized_ = false;
}

template <class Scalar>
std::string SolutionHistory<Scalar>::getStorageTypeString() const
{
  std::string s = "Invalid";
  if (storageType_ == STORAGE_TYPE_KEEP_NEWEST)
    s = "Keep Newest";
  else if (storageType_ == STORAGE_TYPE_UNDO)
    s = "Undo";
  else if (storageType_ == STORAGE_TYPE_STATIC)
    s = "Static";
  else if (storageType_ == STORAGE_TYPE_UNLIMITED)
    s = "Unlimited";
  return s;
}

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::getStateTimeIndexN(bool warn) const
{
  Teuchos::RCP<SolutionState<Scalar> > state = Teuchos::null;
  const int m                                = history_->size();
  if (m < 1) {
    if (warn) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out, 1, "SolutionHistory::getStateTimeIndexN");
      *out << "Warning - getStateTimeIndexN() No states in SolutionHistory!"
           << std::endl;
    }
  }
  else {
    state = (*history_)[m - 1];
  }
  return state;
}

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::getStateTimeIndexNM1(bool warn) const
{
  Teuchos::RCP<SolutionState<Scalar> > state = Teuchos::null;
  const int m                                = history_->size();
  if (m < 2) {
    if (warn) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out, 1, "SolutionHistory::getStateTimeIndexNM1");
      *out << "Warning - getStateTimeIndexNM1() Not enough states in "
           << "SolutionHistory!  Size of history = " << getNumStates()
           << std::endl;
    }
  }
  else {
    const int n   = (*history_)[m - 1]->getIndex();
    const int nm1 = (*history_)[m - 2]->getIndex();

    // No need to search SolutionHistory as states n and nm1 should be
    // next to each other.
    if (nm1 != n - 1) {
      if (warn) {
        Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out, 1, "SolutionHistory::getStateTimeIndexNM1");
        *out << "Warning - getStateTimeIndexNM1() Timestep index n-1 is not in "
             << "SolutionHistory!\n"
             << "    (n)th   index = " << n << "\n"
             << "    (n-1)th index = " << nm1 << std::endl;
      }
    }
    else {
      state = (*history_)[m - 2];
    }
  }

  return state;
}

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> >
SolutionHistory<Scalar>::getStateTimeIndexNM2(bool warn) const
{
  Teuchos::RCP<SolutionState<Scalar> > state = Teuchos::null;
  const int m                                = history_->size();
  if (m < 3) {
    if (warn) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out, 1, "SolutionHistory::getStateTimeIndexNM2");
      *out << "Warning - getStateTimeIndexNM2() Not enough states in "
           << "SolutionHistory!  Size of history = " << getNumStates()
           << std::endl;
    }
  }
  else {
    const int n   = (*history_)[m - 1]->getIndex();
    const int nm2 = (*history_)[m - 3]->getIndex();

    // Assume states n and nm2 are one away from each other.
    if (nm2 != n - 2) {
      // Check if it is at nm1
      const int nm1 = (*history_)[m - 2]->getIndex();
      if (nm1 == n - 2) {
        state = (*history_)[m - 2];
      }
      else if (warn) {
        Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
        Teuchos::OSTab ostab(out, 1, "SolutionHistory::getStateTimeIndexNM2");
        *out << "Warning - getStateTimeIndexNM2() Timestep index n-2 is not in "
             << "SolutionHistory!\n"
             << "    (n)th   index = " << n << "\n"
             << "    (n-2)th index = " << nm2 << std::endl;
      }
    }
    else {
      state = (*history_)[m - 3];
    }
  }

  return state;
}

template <class Scalar>
Teuchos::RCP<SolutionState<Scalar> > SolutionHistory<Scalar>::getStateTimeIndex(
    int index, bool warn) const
{
  typename std::vector<Teuchos::RCP<SolutionState<Scalar> > >::iterator
      state_it = history_->begin();
  for (; state_it < history_->end(); state_it++) {
    if ((*state_it)->getIndex() == index) break;
  }

  Teuchos::RCP<SolutionState<Scalar> > state = Teuchos::null;
  if (state_it == history_->end()) {
    if (warn) {
      Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
      Teuchos::OSTab ostab(out, 1, "SolutionHistory::getStateTimeIndex");
      *out << "Warning - getStateTimeIndex() Timestep index is not in "
           << "SolutionHistory!\n"
           << "    index = " << index << std::endl;
    }
  }
  else {
    state = *state_it;
  }
  return state;
}

template <class Scalar>
std::string SolutionHistory<Scalar>::description() const
{
  return ("Tempus::SolutionHistory - '" + name_ + "'");
}

template <class Scalar>
void SolutionHistory<Scalar>::describe(
    Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "\n--- " << this->description() << " ---" << std::endl;

  if ((Teuchos::as<int>(verbLevel) ==
       Teuchos::as<int>(Teuchos::VERB_DEFAULT)) ||
      (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_LOW))) {
    //*l_out << "  interpolator     = " << interpolator->description() <<
    // std::endl;
    *l_out << "  storageLimit     = " << storageLimit_ << std::endl;
    *l_out << "  storageType      = " << getStorageTypeString() << std::endl;
    *l_out << "  number of states = " << history_->size() << std::endl;
    if (history_->size() > 0) {
      *l_out << "  time range       = (" << history_->front()->getTime() << ", "
             << history_->back()->getTime() << ")" << std::endl;
    }
  }

  if (Teuchos::as<int>(verbLevel) >= Teuchos::as<int>(Teuchos::VERB_MEDIUM)) {
    for (int i = 0; i < (int)history_->size(); ++i)
      (*history_)[i]->describe(*l_out, verbLevel);
  }
  *l_out << std::string(this->description().length() + 8, '-') << std::endl;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
SolutionHistory<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Solution History");

  pl->setName(getName());

  pl->set(
      "Storage Type", getStorageTypeString(),
      "'Storage Type' sets the memory storage.  "
      "'Keep Newest' - will retain the single newest solution state.  "
      "'Undo' - will retain two solution states in order to do a single undo.  "
      "'Static' - will retain 'Storage Limit' number of solution states.  "
      "'Unlimited' - will not remove any solution states!");

  pl->set(
      "Storage Limit", getStorageLimit(),
      "Limit on the number of SolutionStates that SolutionHistory can have.");

  pl->set("Interpolator", *interpolator_->getNonconstParameterList());

  return pl;
}

template <class Scalar>
Teuchos::RCP<Teuchos::ParameterList>
SolutionHistory<Scalar>::getNonconstParameterList()
{
  return Teuchos::rcp_const_cast<Teuchos::ParameterList>(getValidParameters());
}

template <class Scalar>
void SolutionHistory<Scalar>::setInterpolator(
    const Teuchos::RCP<Interpolator<Scalar> >& interpolator)
{
  if (interpolator == Teuchos::null) {
    interpolator_ = InterpolatorFactory<Scalar>::createInterpolator();
  }
  else {
    interpolator_ = interpolator;
  }
  isInitialized_ = false;
}

template <class Scalar>
Teuchos::RCP<Interpolator<Scalar> >
SolutionHistory<Scalar>::getNonconstInterpolator()
{
  return interpolator_;
}

template <class Scalar>
Teuchos::RCP<const Interpolator<Scalar> >
SolutionHistory<Scalar>::getInterpolator() const
{
  return interpolator_;
}

template <class Scalar>
Teuchos::RCP<Interpolator<Scalar> > SolutionHistory<Scalar>::unSetInterpolator()
{
  Teuchos::RCP<Interpolator<Scalar> > old_interpolator = interpolator_;
  interpolator_                                        = lagrangeInterpolator<Scalar>();
  return old_interpolator;
}

template <class Scalar>
void SolutionHistory<Scalar>::printHistory(std::string verb) const
{
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  Teuchos::OSTab ostab(out, 1, "SolutionHistory::printHistory");
  *out << name_ << "  (size=" << history_->size() << ")"
       << "  (w - working; c - current; i - interpolated)" << std::endl;
  for (int i = 0; i < (int)history_->size(); ++i) {
    auto state = (*history_)[i];
    *out << "  ";
    if (state == getWorkingState())
      *out << "w - ";
    else if (state == getCurrentState())
      *out << "c - ";
    else if (state->getIsInterpolated() == true)
      *out << "i - ";
    else
      *out << "    ";
    *out << "[" << i << "] = " << state << std::endl;
    if (verb == "medium" || verb == "high") {
      if (state != Teuchos::null) {
        auto x = state->getX();
        *out << "      x       = " << x << std::endl
             << "      norm(x) = " << Thyra::norm(*x) << std::endl;
      }
    }
    if (verb == "high") {
      (*history_)[i]->describe(*out, this->getVerbLevel());
    }
  }
}

template <class Scalar>
void SolutionHistory<Scalar>::initialize() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
      getNumStates() <= 0, std::logic_error,
      "Error - SolutionHistory::initialize() Invalid history size!\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
      interpolator_ == Teuchos::null, std::logic_error,
      "Error - SolutionHistory::initialize() Interpolator is not set!\n");

  TEUCHOS_TEST_FOR_EXCEPTION(storageLimit_ < 1, std::logic_error,
                             "Error - SolutionHistory::initialize() Storage "
                             "Limit needs to a positive integer!\n"
                                 << "  Storage Limit = " << storageLimit_
                                 << "\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
      (storageType_ == STORAGE_TYPE_INVALID), std::logic_error,
      "Error - SolutionHistory::initialize() Storage Type is invalid!\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
      (storageType_ == STORAGE_TYPE_KEEP_NEWEST && storageLimit_ != 1),
      std::logic_error,
      "Error - SolutionHistory::initialize() \n"
          << "  For Storage Type = '" << getStorageTypeString()
          << "', Storage Limit needs to be one.\n"
          << "  Storage Limit = " << storageLimit_ << "\n");

  TEUCHOS_TEST_FOR_EXCEPTION(
      (storageType_ == STORAGE_TYPE_UNDO && storageLimit_ != 2),
      std::logic_error,
      "Error - SolutionHistory::initialize() \n"
          << "  For Storage Type = '" << getStorageTypeString()
          << "', Storage Limit needs to be two.\n"
          << "  Storage Limit = " << storageLimit_ << "\n");

  isInitialized_ = true;  // Only place where this is set to true!
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > createSolutionHistory()
{
  auto sh = rcp(new SolutionHistory<Scalar>());
  sh->setName("From createSolutionHistory");

  return sh;
}

template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > createSolutionHistoryPL(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto sh = rcp(new SolutionHistory<Scalar>());
  sh->setName("From createSolutionHistoryPL");

  if (pl == Teuchos::null || pl->numParams() == 0) return sh;

  pl->validateParametersAndSetDefaults(*sh->getValidParameters());

  sh->setName(pl->name());
  sh->setStorageTypeString(pl->get("Storage Type", "Undo"));
  sh->setStorageLimit(pl->get("Storage Limit", 2));

  sh->setInterpolator(InterpolatorFactory<Scalar>::createInterpolator(
      Teuchos::sublist(pl, "Interpolator")));

  return sh;
}

template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > createSolutionHistoryState(
    const Teuchos::RCP<SolutionState<Scalar> >& state)
{
  auto sh = rcp(new SolutionHistory<Scalar>());
  sh->setName("From createSolutionHistoryState");
  sh->addState(state);
  return sh;
}

template <class Scalar>
Teuchos::RCP<SolutionHistory<Scalar> > createSolutionHistoryME(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model)
{
  // Setup initial condition SolutionState --------------------
  auto state = createSolutionStateME(model);
  state->setTime(0.0);
  state->setIndex(0);
  state->setTimeStep(0.0);
  state->setOrder(1);

  // Setup SolutionHistory ------------------------------------
  auto sh = rcp(new SolutionHistory<Scalar>());
  sh->setName("From createSolutionHistoryME");
  sh->addState(state);

  return sh;
}

}  // namespace Tempus
#endif  // Tempus_SolutionHistory_impl_hpp
