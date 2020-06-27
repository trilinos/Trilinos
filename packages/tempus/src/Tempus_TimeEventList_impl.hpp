// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventList_impl_hpp
#define Tempus_TimeEventList_impl_hpp


namespace Tempus {

template<class Scalar>
TimeEventList<Scalar>::TimeEventList()
{
  this->setName("TimeEventList");
  setRelTol(1.0e-14);
  setTimeScale();
  setLandOnExactly(true);
}


template<class Scalar>
TimeEventList<Scalar>::TimeEventList(
  std::string name, std::vector<Scalar> timeList,
  Scalar relTol, bool landOnExactly)
{
  this->setName(name);
  setRelTol(relTol);
  setTimeScale();
  setLandOnExactly(landOnExactly);
  setTimeList(timeList);
}


template<class Scalar>
void TimeEventList<Scalar>::setTimeScale()
{
  if (timeList_.size() == 0) {
    timeScale_ = std::abs(this->getDefaultTime());
    absTol_ = relTol_*timeScale_;
    return;
  }

  timeScale_ = std::max(std::abs(timeList_.front()),
                        std::abs(timeList_.back()));
  absTol_ = relTol_*timeScale_;

  // Check if timeScale is near zero.
  if ((-absTol_ <= timeScale_ ) && (timeScale_ <= absTol_)) {
    timeScale_ = 1.0;
    absTol_ = relTol_*timeScale_;
  }
}


template<class Scalar>
void TimeEventList<Scalar>::setTimeList(std::vector<Scalar> timeList)
{
  for(auto it = std::begin(timeList); it != std::end(timeList); ++it)
    addTime(*it);
}


template<class Scalar>
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
    if (std::abs(timeList_.front()-time) >= absTol_)
      timeList_.insert(it, time);
    setTimeScale();
    return;
  }

  // Don't add if already in list. Check both ends.
  if (it == timeList_.begin()) {
    if (std::abs(timeList_.front()-time) >= absTol_)
      timeList_.insert(it, time);
  } else if (it == timeList_.end()) {
    if (std::abs(timeList_.back()-time) >= absTol_)
      timeList_.insert(it, time);
  } else if (std::abs(*(it-1) - time) >= absTol_ &&
             std::abs(*(it  ) - time) >= absTol_) {
    timeList_.insert(it, time);
  }
  setTimeScale();
}


template<class Scalar>
void TimeEventList<Scalar>::setRelTol(Scalar relTol)
{
  relTol_ = std::abs(relTol);
  setTimeScale();
}


template<class Scalar>
bool TimeEventList<Scalar>::isTime(Scalar time) const
{
  return (std::abs(timeToNextEvent(time)) <= absTol_);
}


template<class Scalar>
Scalar TimeEventList<Scalar>::timeToNextEvent(Scalar time) const
{
  return timeOfNextEvent(time) - time;  // Neg. indicating in the past.
}


template<class Scalar>
Scalar TimeEventList<Scalar>::timeOfNextEvent(Scalar time) const
{
  if (timeList_.size() == 0) return this->getDefaultTime();

  // Check if before or close to first event.
  if (timeList_.front() >= time-absTol_) return timeList_.front();

  // Check if after or close to last event.
  if (timeList_.back() <= time+absTol_) return timeList_.back();

  typename std::vector<Scalar>::const_iterator it =
    std::upper_bound(timeList_.begin(), timeList_.end(), time);

  // Check if close to left-side time event
  const Scalar timeOfLeftEvent = *(it-1);
  if (timeOfLeftEvent > time-absTol_ &&
      timeOfLeftEvent < time+absTol_) return timeOfLeftEvent;

  // Otherwise it is the next event.
  return *it;
}


template<class Scalar>
bool TimeEventList<Scalar>::eventInRange(
  Scalar time1, Scalar time2) const
{
  if (time1 > time2) {
    Scalar tmp = time1;
    time1 = time2;
    time2 = tmp;
  }

  if (timeList_.size() == 0) return false;

  // Check if range is completely outside time events.
  if (time2+absTol_ < timeList_.front() ||
       timeList_.back() < time1-absTol_) return false;

  Scalar timeEvent1 = timeOfNextEvent(time1);
  Scalar timeEvent2 = timeOfNextEvent(time2);
  // Check if the next time event is different for the two times.
  if (timeEvent1 != timeEvent2) return true;

  // Check if times bracket time event.
  if (time1-absTol_ <= timeEvent1 && timeEvent1 <= time2+absTol_) return true;

  return false;
}


template<class Scalar>
void TimeEventList<Scalar>::describe() const
{
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "TimeEventList:" << "\n"
       << "name           = " << this->getName() << "\n"
       << "timeScale_     = " << timeScale_     << "\n"
       << "relTol_        = " << relTol_        << "\n"
       << "absTol_        = " << absTol_        << "\n"
       << "landOnExactly_ = " << landOnExactly_ << "\n"
       << "timeList_      = " << std::endl;
  for (auto it = timeList_.begin(); it != timeList_.end()-1; ++it)
    *out << *it << ", ";
  *out << *(timeList_.end()-1) << "\n";
}


} // namespace Tempus
#endif // Tempus_TimeEventList_impl_hpp
