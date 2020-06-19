// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventRange_impl_hpp
#define Tempus_TimeEventRange_impl_hpp


namespace Tempus {

template<class Scalar>
TimeEventRange<Scalar>::TimeEventRange()
  : start_(this->getDefaultTime()),
    stop_ (this->getDefaultTime()),
    stride_(0.0),
    numEvents_(1),
    timeScale_(1.0),
    relTol_(1.0e-14),
    absTol_(1.0e-14),
    landOnExactly_(true)
{
  this->setName("TimeEventRange");
  setTimeRange(start_, stop_, stride_);
}


template<class Scalar>
TimeEventRange<Scalar>::TimeEventRange(
  std::string name, Scalar start, Scalar stop, Scalar stride,
  Scalar relTol, bool landOnExactly)
  : start_(start),
    stop_ (stop),
    stride_(stride),
    numEvents_(std::abs(stop-start)/stride+1),
    timeScale_(start),
    relTol_(relTol),
    absTol_(relTol*start),
    landOnExactly_(landOnExactly)
{
  this->setName(name);
  setTimeRange(start, stop, stride);
}


template<class Scalar>
TimeEventRange<Scalar>::TimeEventRange(
  std::string name, Scalar start, Scalar stop, int numEvents,
  Scalar relTol, bool landOnExactly)
  : start_(start),
    stop_ (stop),
    stride_(std::abs(stop-start)/(numEvents-1)),
    numEvents_(numEvents),
    timeScale_(start),
    relTol_(relTol),
    absTol_(relTol*start),
    landOnExactly_(landOnExactly)
{
  this->setName(name);
  setTimeRange(start, stop, numEvents);
  setRelTol(relTol);
  setLandOnExactly(landOnExactly);
}


template<class Scalar>
void TimeEventRange<Scalar>::setTimeStart(Scalar start)
{
  start_ = start;
  if (stop_ < start_) stop_ = start_;
  setTimeStride(stride_);  // Reset numEvents with the current stride.
  setTimeScale();
}


template<class Scalar>
void TimeEventRange<Scalar>::setTimeStop(Scalar stop)
{
  stop_ = stop;
  if (start_ > stop_) start_ = stop_;
  setTimeStride(stride_);  // Reset numEvents with the current stride.
  setTimeScale();
}


template<class Scalar>
void TimeEventRange<Scalar>::setTimeScale()
{
  timeScale_ = std::max(std::abs(start_), std::abs(stop_));
  absTol_ = relTol_*timeScale_;

  // Check if timeScale is near zero.
  if ((-absTol_ <= timeScale_ ) && (timeScale_ <= absTol_)) {
    timeScale_ = 1.0;
    absTol_ = relTol_*timeScale_;
  }
}


template<class Scalar>
void TimeEventRange<Scalar>::setTimeStride(Scalar stride)
{
  stride_ = stride;
  if ((start_ >= stop_-absTol_) && (start_ <= stop_+absTol_)) {
    stride_ = 0.0;
    numEvents_ = 1;
    return;
  }

  if ((stride_ > stop_ - start_) || (stride_ < 2*absTol_)) {
    stride_ = stop_ - start_;
  }

  numEvents_ = int((stop_+absTol_ - start_) / stride_) + 1;
}


template<class Scalar>
void TimeEventRange<Scalar>::setNumEvents(int numEvents)
{
  numEvents_ = numEvents;
  if (numEvents_ < 1) numEvents_ = 1;
  stride_ = (stop_ - start_)/Scalar(numEvents_-1);

  if (stride_ < 2 * absTol_) {
    setTimeStride(2*absTol_);
  }
}


template<class Scalar>
void TimeEventRange<Scalar>::setRelTol(Scalar relTol)
{
  relTol_ = std::abs(relTol);
  absTol_ = relTol_*timeScale_;
}


template<class Scalar>
bool TimeEventRange<Scalar>::isTime(Scalar time) const
{
  return (std::abs(timeToNextEvent(time)) <= absTol_);
}


template<class Scalar>
Scalar TimeEventRange<Scalar>::timeToNextEvent(Scalar time) const
{
  return timeOfNextEvent(time) - time;  // Neg. time indicates in the past.
}


template<class Scalar>
Scalar TimeEventRange<Scalar>::timeOfNextEvent(Scalar time) const
{
  // Check if before or close to first event.
  if (start_ >= time-absTol_) return start_;

  const Scalar timeOfLast = start_ + (numEvents_-1) * stride_;
  // Check if after or close to last event.
  if (timeOfLast <= time+absTol_) return timeOfLast;

  const int numStrides = (time - start_) / stride_;
  const Scalar timeOfNext = start_ + numStrides * stride_;

  // Check if close to left-side time event
  if (timeOfNext > time-absTol_ &&
      timeOfNext < time+absTol_) return timeOfNext;

  // Otherwise it is the next event.
  return timeOfNext + stride_;
}


template<class Scalar>
bool TimeEventRange<Scalar>::eventInRange(Scalar time1, Scalar time2) const
{
  if (time1 > time2) {
    Scalar tmp = time1;
    time1 = time2;
    time2 = tmp;
  }

  const Scalar timeOfLast = start_ + (numEvents_-1) * stride_;
  // Check if range is completely outside time events.
  if (time2+absTol_ < start_ || timeOfLast < time1-absTol_) return false;

  Scalar timeEvent1 = timeOfNextEvent(time1);
  Scalar timeEvent2 = timeOfNextEvent(time2);
  // Check if the next time event is different for the two times.
  if (timeEvent1 != timeEvent2) return true;

  // Check if times bracket time event.
  if (time1-absTol_ <= timeEvent1 && timeEvent1 <= time2+absTol_) return true;

  return false;
}


template<class Scalar>
void TimeEventRange<Scalar>::describe() const
{
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "TimeEventRange:" << "\n"
       << "name           = " << this->getName() << "\n"
       << "start_         = " << start_         << "\n"
       << "stop_          = " << stop_          << "\n"
       << "stride_        = " << stride_        << "\n"
       << "numEvents_     = " << numEvents_     << "\n"
       << "timeScale_     = " << timeScale_     << "\n"
       << "relTol_        = " << relTol_        << "\n"
       << "absTol_        = " << absTol_        << "\n"
       << "landOnExactly_ = " << landOnExactly_ << std::endl;
}


} // namespace Tempus
#endif // Tempus_TimeEventRange_impl_hpp
