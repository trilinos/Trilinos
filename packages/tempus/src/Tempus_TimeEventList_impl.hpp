//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventList_impl_hpp
#define Tempus_TimeEventList_impl_hpp

#include "Tempus_NumericalUtils.hpp"

namespace Tempus {

template <class Scalar>
TimeEventList<Scalar>::TimeEventList()
  : timeScale_(1.0),
    relTol_(this->getDefaultTol()),
    absTol_(this->getDefaultTol()),
    landOnExactly_(true)
{
  this->setType("List");
  this->setName("TimeEventList");
  setRelTol(this->getDefaultTol());
  setTimeScale();
  setLandOnExactly(true);
}

template <class Scalar>
TimeEventList<Scalar>::TimeEventList(std::vector<Scalar> timeList,
                                     std::string name, bool landOnExactly,
                                     Scalar relTol)
  : timeScale_(1.0),
    relTol_(this->getDefaultTol()),
    absTol_(this->getDefaultTol()),
    landOnExactly_(true)
{
  this->setType("List");
  this->setName(name);
  setRelTol(relTol);
  setTimeScale();
  setLandOnExactly(landOnExactly);
  setTimeList(timeList);
}

template <class Scalar>
void TimeEventList<Scalar>::setTimeScale()
{
  if (timeList_.size() == 0) {
    timeScale_ = 1.0;
    absTol_    = relTol_ * timeScale_;
    return;
  }

  timeScale_ =
      std::max(std::abs(timeList_.front()), std::abs(timeList_.back()));
  absTol_ = relTol_ * timeScale_;

  // Check if timeScale is near zero.
  if (approxZero(timeScale_, absTol_)) {
    timeScale_ = 1.0;
    absTol_    = relTol_ * timeScale_;
  }
}

template <class Scalar>
void TimeEventList<Scalar>::setTimeList(std::vector<Scalar> timeList, bool sort)
{
  timeList_ = timeList;
  if (sort) {
    std::sort(timeList_.begin(), timeList_.end());
    timeList_.erase(std::unique(timeList_.begin(), timeList_.end()),
                    timeList_.end());
  }
}

template <class Scalar>
void TimeEventList<Scalar>::addTime(Scalar time)
{
  if (timeList_.size() == 0) {
    timeList_.push_back(time);
    setTimeScale();
    return;
  }

  auto it = std::upper_bound(timeList_.begin(), timeList_.end(), time);
  if (timeList_.size() == 1) {
    // Only add if a "different" time.
    if (std::abs(timeList_.front() - time) >= absTol_)
      timeList_.insert(it, time);
    setTimeScale();
    return;
  }

  // Don't add if already in list. Check both ends.
  if (it == timeList_.begin()) {
    if (std::abs(timeList_.front() - time) >= absTol_)
      timeList_.insert(it, time);
  }
  else if (it == timeList_.end()) {
    if (std::abs(timeList_.back() - time) >= absTol_)
      timeList_.insert(it, time);
  }
  else if (std::abs(*(it - 1) - time) >= absTol_ &&
           std::abs(*(it)-time) >= absTol_) {
    timeList_.insert(it, time);
  }
  setTimeScale();
}

template <class Scalar>
void TimeEventList<Scalar>::setRelTol(Scalar relTol)
{
  relTol_ = std::abs(relTol);
  setTimeScale();
}

template <class Scalar>
bool TimeEventList<Scalar>::isTime(Scalar time) const
{
  for (auto it = timeList_.begin(); it != timeList_.end(); ++it)
    if (approxEqualAbsTol(time, *it, absTol_)) return true;

  return false;
}

template <class Scalar>
Scalar TimeEventList<Scalar>::timeToNextEvent(Scalar time) const
{
  return timeOfNextEvent(time) - time;  // Neg. indicating in the past.
}

template <class Scalar>
Scalar TimeEventList<Scalar>::timeOfNextEvent(Scalar time) const
{
  if (timeList_.size() == 0) return this->getDefaultTime();

  // Check if before first event.
  if (time < timeList_.front() - absTol_) return timeList_.front();

  // Check if after or close to last event.
  if (time > timeList_.back() - absTol_)
    return std::numeric_limits<Scalar>::max();

  typename std::vector<Scalar>::const_iterator it =
      std::upper_bound(timeList_.begin(), timeList_.end(), time);
  const Scalar timeEvent = *it;

  // Check timeEvent is near time.  If so, return the next event.
  if (approxEqualAbsTol(time, timeEvent, absTol_)) return *(it + 1);

  return timeEvent;
}

template <class Scalar>
bool TimeEventList<Scalar>::eventInRange(Scalar time1, Scalar time2) const
{
  if (time1 > time2) {
    Scalar tmp = time1;
    time1      = time2;
    time2      = tmp;
  }

  if (timeList_.size() == 0) return false;

  for (auto it = timeList_.begin(); it != timeList_.end(); ++it)
    if (time1 + absTol_ < *it && *it < time2 + absTol_) return true;

  return false;
}

template <class Scalar>
void TimeEventList<Scalar>::describe(
    Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, "TimeEventList");
  l_out->setOutputToRootOnly(0);

  *l_out << "TimeEventList:"
         << "\n"
         << "  name            = " << this->getName() << "\n"
         << "  Type            = " << this->getType() << "\n"
         << "  timeScale_      = " << timeScale_ << "\n"
         << "  relTol_         = " << relTol_ << "\n"
         << "  absTol_         = " << absTol_ << "\n"
         << "  landOnExactly_  = " << landOnExactly_ << "\n"
         << "  timeList_       = ";
  if (!timeList_.empty()) {
    for (auto it = timeList_.begin(); it != timeList_.end() - 1; ++it)
      *l_out << *it << ", ";
    *l_out << *(timeList_.end() - 1) << std::endl;
  }
  else {
    *l_out << "<empty>" << std::endl;
  }
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
TimeEventList<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Event List");

  pl->setName(this->getName());
  pl->set("Name", this->getName());
  pl->set("Type", this->getType());

  pl->set("Relative Tolerance", this->getRelTol(),
          "Relative time tolerance for matching time events.");

  pl->set("Land On Exactly", this->getLandOnExactly(),
          "Should these time events be landed on exactly, i.e, adjust the "
          "timestep to hit time event, versus stepping over and keeping the "
          "time step unchanged.");

  std::vector<Scalar> times = this->getTimeList();
  std::ostringstream list;
  if (!times.empty()) {
    for (std::size_t i = 0; i < times.size() - 1; ++i) list << times[i] << ", ";
    list << times[times.size() - 1];
  }
  pl->set<std::string>("Time List", list.str(),
                       "Comma deliminated list of times");

  return pl;
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<TimeEventList<Scalar> > createTimeEventList(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto tel = Teuchos::rcp(new TimeEventList<Scalar>());
  if (pl == Teuchos::null) return tel;  // Return default TimeEventList.

  TEUCHOS_TEST_FOR_EXCEPTION(pl->get<std::string>("Type", "List") != "List",
                             std::logic_error,
                             "Error - Time Event Type != 'List'.  (='" +
                                 pl->get<std::string>("Type") + "')\n");

  pl->validateParametersAndSetDefaults(*tel->getValidParameters());

  tel->setName(pl->get("Name", "From createTimeEventList"));
  tel->setRelTol(pl->get("Relative Tolerance", tel->getRelTol()));
  tel->setLandOnExactly(pl->get("Land On Exactly", tel->getLandOnExactly()));

  std::vector<Scalar> timeList;
  std::string str = pl->get<std::string>("Time List");
  std::string delimiters(",");
  // Skip delimiters at the beginning
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find the first delimiter
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);
  while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
    // Found a token, add it to the vector
    std::string token = str.substr(lastPos, pos - lastPos);
    timeList.push_back(Scalar(std::stod(token)));
    if (pos == std::string::npos) break;

    lastPos = str.find_first_not_of(delimiters, pos);  // Skip delimiters
    pos     = str.find_first_of(delimiters, lastPos);  // Find next delimiter
  }
  tel->setTimeList(timeList);

  return tel;
}

}  // namespace Tempus
#endif  // Tempus_TimeEventList_impl_hpp
