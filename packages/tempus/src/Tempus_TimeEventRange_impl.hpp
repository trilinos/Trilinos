//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventRange_impl_hpp
#define Tempus_TimeEventRange_impl_hpp

#include "Tempus_NumericalUtils.hpp"

namespace Tempus {

template <class Scalar>
TimeEventRange<Scalar>::TimeEventRange()
  : start_(0.0),
    stop_(0.0),
    stride_(0.0),
    numEvents_(1),
    timeScale_(1.0),
    relTol_(this->getDefaultTol()),
    absTol_(this->getDefaultTol()),
    landOnExactly_(true)
{
  this->setType("Range");
  setRelTol(this->getDefaultTol()), setTimeRange(0.0, 0.0, 0.0);
  setLandOnExactly(true);
  std::ostringstream oss;
  oss << "TimeEventRange (" << start_ << "; " << stop_ << "; " << stride_
      << ")";
  this->setName(oss.str());
}

template <class Scalar>
TimeEventRange<Scalar>::TimeEventRange(Scalar start, Scalar stop, Scalar stride,
                                       std::string name, bool landOnExactly,
                                       Scalar relTol)
  : start_(start),
    stop_(stop),
    stride_(stride),
    numEvents_(0),
    timeScale_(std::max(std::abs(start_), std::abs(stop_))),
    relTol_(relTol),
    absTol_(relTol_ * timeScale_),
    landOnExactly_(landOnExactly)
{
  this->setType("Range");
  if (name == "") {
    std::ostringstream oss;
    oss << "TimeEventRange (" << start << "; " << stop << "; " << stride << ")";
    this->setName(oss.str());
  }
  else {
    this->setName(name);
  }
  setRelTol(relTol);
  setTimeRange(start, stop, stride);
  setLandOnExactly(landOnExactly);
}

template <class Scalar>
TimeEventRange<Scalar>::TimeEventRange(Scalar start, Scalar stop, int numEvents,
                                       std::string name, bool landOnExactly,
                                       Scalar relTol)
  : start_(start),
    stop_(stop),
    stride_(0.0),
    numEvents_(numEvents),
    timeScale_(std::max(std::abs(start_), std::abs(stop_))),
    relTol_(relTol),
    absTol_(relTol_ * timeScale_),
    landOnExactly_(landOnExactly)
{
  if (name == "") {
    std::ostringstream oss;
    oss << "TimeEventRange (" << start << "; " << stop << "; " << numEvents
        << ")";
    this->setName(oss.str());
  }
  else {
    this->setName(name);
  }
  setRelTol(relTol);
  setTimeRange(start, stop, numEvents);
  setLandOnExactly(landOnExactly);
}

template <class Scalar>
void TimeEventRange<Scalar>::setTimeRange(Scalar start, Scalar stop,
                                          Scalar stride)
{
  start_ = start;
  stop_  = stop;
  if (stop_ < start_) {
    Scalar tmp = start_;
    start_     = stop_;
    stop_      = tmp;
  }
  setTimeScale();
  setTimeStride(stride);
}

template <class Scalar>
void TimeEventRange<Scalar>::setTimeRange(Scalar start, Scalar stop,
                                          int numEvents)
{
  start_ = start;
  stop_  = stop;
  if (stop_ < start_) {
    Scalar tmp = start_;
    start_     = stop_;
    stop_      = tmp;
  }
  setTimeScale();
  setNumEvents(numEvents);
}

template <class Scalar>
void TimeEventRange<Scalar>::setTimeStart(Scalar start)
{
  start_ = start;
  if (stop_ < start_) stop_ = start_;
  setTimeScale();
  setTimeStride(stride_);  // Reset numEvents with the current stride.
}

template <class Scalar>
void TimeEventRange<Scalar>::setTimeStop(Scalar stop)
{
  stop_ = stop;
  if (start_ > stop_) start_ = stop_;
  setTimeScale();
  setTimeStride(stride_);  // Reset numEvents with the current stride.
}

template <class Scalar>
void TimeEventRange<Scalar>::setTimeScale()
{
  timeScale_ = std::max(std::abs(start_), std::abs(stop_));
  absTol_    = relTol_ * timeScale_;

  // Check if timeScale is near zero.
  if (approxZero(timeScale_, absTol_)) {
    timeScale_ = 1.0;
    absTol_    = relTol_ * timeScale_;
  }
}

template <class Scalar>
void TimeEventRange<Scalar>::setTimeStride(Scalar stride)
{
  stride_ = Teuchos::ScalarTraits<Scalar>::magnitude(stride);
  if (approxEqualAbsTol(start_, stop_, absTol_)) {
    stride_    = 0.0;
    numEvents_ = 1;
    return;
  }

  if ((stride_ > stop_ - start_) || (stride_ < 2 * absTol_)) {
    stride_ = stop_ - start_;
  }

  numEvents_ = int((stop_ + absTol_ - start_) / stride_) + 1;
}

template <class Scalar>
void TimeEventRange<Scalar>::setNumEvents(int numEvents)
{
  numEvents_ = numEvents;
  if (numEvents_ < 0) {  // Do not use numEvents_ to set stride!  Reset
                         // numEvents_ from stride_.
    if (stride_ < 2 * absTol_) stride_ = 2 * absTol_;
    numEvents_ = int((stop_ + absTol_ - start_) / stride_) + 1;
    return;
  }
  else if (approxEqualAbsTol(start_, stop_, absTol_)) {
    numEvents_ = 1;
    stride_    = 0.0;
    stride_    = stop_ - start_;
  }
  else {
    if (numEvents_ < 2) numEvents_ = 2;
    stride_ = (stop_ - start_) / Scalar(numEvents_ - 1);
  }

  // If stride_ is smaller than twice the absolute tolerance,
  // the time steps cannot be differentiated!
  if (stride_ <= 2 * absTol_) setTimeStride(2 * absTol_);
}

template <class Scalar>
void TimeEventRange<Scalar>::setRelTol(Scalar relTol)
{
  relTol_ = std::abs(relTol);
  setTimeScale();
}

template <class Scalar>
bool TimeEventRange<Scalar>::isTime(Scalar time) const
{
  // Check if before first event.
  if (time < start_ - absTol_) return false;

  // Check if after last event.
  const Scalar timeOfLast = start_ + (numEvents_ - 1) * stride_;
  if (time > timeOfLast + absTol_) return false;

  int numStrides = 0;
  if (!approxZero(stride_, 2 * absTol_)) numStrides = (time - start_) / stride_;

  numStrides               = std::min(std::max(0, numStrides), int(numEvents_ - 1));
  const Scalar leftBracket = start_ + numStrides * stride_;

  // Check if close to left bracket.
  if (approxEqualAbsTol(time, leftBracket, absTol_)) return true;

  // Check if close to right bracket.
  const Scalar rightBracket = leftBracket + stride_;
  if (approxEqualAbsTol(time, rightBracket, absTol_)) return true;

  return false;
}

template <class Scalar>
Scalar TimeEventRange<Scalar>::timeToNextEvent(Scalar time) const
{
  return timeOfNextEvent(time) - time;  // Neg. time indicates in the past.
}

template <class Scalar>
Scalar TimeEventRange<Scalar>::timeOfNextEvent(Scalar time) const
{
  // Check if before first event.
  if (time < start_ - absTol_) return start_;

  const Scalar timeOfLast = start_ + (numEvents_ - 1) * stride_;
  // Check if after or close to last event.
  if (time > timeOfLast - absTol_) return std::numeric_limits<Scalar>::max();

  int numStrides = 0;
  if (!approxZero(stride_, 2 * absTol_))
    numStrides = int((time - start_) / stride_) + 1;
  const Scalar timeEvent = start_ + numStrides * stride_;

  // Check timeEvent is near time.  If so, return the next event.
  if (approxEqualAbsTol(time, timeEvent, absTol_)) return timeEvent + stride_;

  return timeEvent;
}

template <class Scalar>
bool TimeEventRange<Scalar>::eventInRange(Scalar time1, Scalar time2) const
{
  if (time1 > time2) {
    Scalar tmp = time1;
    time1      = time2;
    time2      = tmp;
  }

  // Check if range is completely outside time events.
  const Scalar timeOfLast = start_ + (numEvents_ - 1) * stride_;
  if (time2 < start_ - absTol_ || timeOfLast + absTol_ < time1) return false;

  if (approxZero(stride_))
    return (time1 < start_ - absTol_ && start_ - absTol_ <= time2);

  const int strideJustBeforeTime1 = std::min(
      int(numEvents_ - 1),
      std::max(int(0), int((time1 - start_ + absTol_) / stride_ - 0.5)));

  const int strideJustAfterTime2 = std::max(
      int(0), std::min(int(numEvents_ - 1),
                       int((time2 - start_ + absTol_) / stride_ + 0.5)));

  for (int i = strideJustBeforeTime1; i <= strideJustAfterTime2; i++) {
    const Scalar timeEvent = start_ + i * stride_;
    if (time1 < timeEvent - absTol_ && timeEvent - absTol_ <= time2)
      return true;
  }

  return false;
}

template <class Scalar>
void TimeEventRange<Scalar>::describe(
    Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, "TimeEventRange");
  l_out->setOutputToRootOnly(0);

  *l_out << "TimeEventRange:"
         << "\n"
         << "  name            = " << this->getName() << "\n"
         << "  Type            = " << this->getType() << "\n"
         << "  start_          = " << start_ << "\n"
         << "  stop_           = " << stop_ << "\n"
         << "  stride_         = " << stride_ << "\n"
         << "  numEvents_      = " << numEvents_ << "\n"
         << "  timeScale_      = " << timeScale_ << "\n"
         << "  relTol_         = " << relTol_ << "\n"
         << "  absTol_         = " << absTol_ << "\n"
         << "  landOnExactly_  = " << landOnExactly_ << std::endl;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
TimeEventRange<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Event Range");

  pl->setName(this->getName());
  pl->set("Name", this->getName());
  pl->set("Type", this->getType());

  pl->set("Start Time", getTimeStart(), "Start of time range");
  pl->set("Stop Time", getTimeStop(), "Stop of time range");
  pl->set("Stride Time", getTimeStride(), "Stride of time range");

  if (getTimeStride() * Scalar(getNumEvents() - 1) -
          (getTimeStop() - getTimeStart()) <
      getAbsTol())
    pl->set("Number of Events", getNumEvents(),
            "Number of events in time range.  If specified, 'Stride Time' is "
            "reset.");

  pl->set("Relative Tolerance", getRelTol(),
          "Relative time tolerance for matching time events.");

  pl->set("Land On Exactly", getLandOnExactly(),
          "Should these time events be landed on exactly, i.e, adjust the "
          "timestep to hit time event, versus stepping over and keeping the "
          "time step unchanged.");

  return pl;
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<TimeEventRange<Scalar> > createTimeEventRange(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto ter = Teuchos::rcp(new TimeEventRange<Scalar>());
  if (pl == Teuchos::null) return ter;  // Return default TimeEventRange.

  TEUCHOS_TEST_FOR_EXCEPTION(pl->get<std::string>("Type", "Range") != "Range",
                             std::logic_error,
                             "Error - Time Event Type != 'Range'.  (='" +
                                 pl->get<std::string>("Type") + "')\n");

  auto validPL        = *ter->getValidParameters();
  bool numEventsFound = pl->isParameter("Number of Events");
  if (!numEventsFound) validPL.remove("Number of Events");

  pl->validateParametersAndSetDefaults(validPL);

  ter->setName(pl->get("Name", "From createTimeEventRange"));
  ter->setTimeStart(pl->get("Start Time", ter->getTimeStart()));
  ter->setTimeStop(pl->get("Stop Time", ter->getTimeStop()));
  ter->setTimeStride(pl->get("Stride Time", ter->getTimeStride()));
  if (numEventsFound)
    ter->setNumEvents(pl->get("Number of Events", ter->getNumEvents()));
  ter->setRelTol(pl->get("Relative Tolerance", ter->getRelTol()));
  ter->setLandOnExactly(pl->get("Land On Exactly", ter->getLandOnExactly()));

  return ter;
}

}  // namespace Tempus
#endif  // Tempus_TimeEventRange_impl_hpp
