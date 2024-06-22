//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventRangeIndex_impl_hpp
#define Tempus_TimeEventRangeIndex_impl_hpp

namespace Tempus {

template <class Scalar>
TimeEventRangeIndex<Scalar>::TimeEventRangeIndex()
  : start_(0), stop_(0), stride_(1), numEvents_(1)
{
  this->setType("Range Index");
  std::ostringstream oss;
  oss << "TimeEventRangeIndex (" << start_ << "; " << stop_ << "; " << stride_
      << ")";
  this->setName(oss.str());
  setNumEvents();
}

template <class Scalar>
TimeEventRangeIndex<Scalar>::TimeEventRangeIndex(int start, int stop,
                                                 int stride, std::string name)
  : start_(0), stop_(0), stride_(1), numEvents_(1)
{
  this->setType("Range Index");
  if (name == "") {
    std::ostringstream oss;
    oss << "TimeEventRangeIndex (" << start << "; " << stop << "; " << stride
        << ")";
    this->setName(oss.str());
  }
  else {
    this->setName(name);
  }

  this->setIndexRange(start, stop, stride);
}

template <class Scalar>
void TimeEventRangeIndex<Scalar>::setIndexStart(int start)
{
  start_ = start;
  if (stop_ < start_) {
    stop_ = start_;
  }
  setNumEvents();
}

template <class Scalar>
void TimeEventRangeIndex<Scalar>::setIndexStop(int stop)
{
  stop_ = stop;
  if (start_ > stop_) {
    start_ = stop_;
    setIndexStride(1);
  }
  setNumEvents();
}

template <class Scalar>
void TimeEventRangeIndex<Scalar>::setIndexStride(int stride)
{
  stride_ = stride;
  if (stride_ < 1) {
    stride_ = 1;
  }
  else if (stride_ > (stop_ - start_)) {
    stride_ = stop_ - start_;
  }
  setNumEvents();
}

template <class Scalar>
void TimeEventRangeIndex<Scalar>::setNumEvents()
{
  if (stride_ == 0 || start_ == stop_)
    numEvents_ = 1;
  else
    numEvents_ = int((stop_ - start_) / stride_) + 1;
}

template <class Scalar>
bool TimeEventRangeIndex<Scalar>::isIndex(int index) const
{
  const int indexOfLast = start_ + (numEvents_ - 1) * stride_;
  if (index < start_ || index > indexOfLast) return false;
  if ((index - start_) % stride_ == 0) return true;

  return false;
}

template <class Scalar>
int TimeEventRangeIndex<Scalar>::indexToNextEvent(int index) const
{
  return indexOfNextEvent(index) - index;
}

template <class Scalar>
int TimeEventRangeIndex<Scalar>::indexOfNextEvent(int index) const
{
  // Check if before the first index.
  if (index < start_) return start_;

  const int indexOfLast = start_ + (numEvents_ - 1) * stride_;
  // Check if after or equal to last index.
  if (index >= indexOfLast) return this->getDefaultIndex();

  // Check if index is an event.  If so, return next event.
  if (isIndex(index)) return index + stride_;

  const int numStrides     = (index - start_) / stride_ + 1;
  const Scalar indexOfNext = start_ + numStrides * stride_;
  return indexOfNext;
}

template <class Scalar>
bool TimeEventRangeIndex<Scalar>::eventInRangeIndex(int index1,
                                                    int index2) const
{
  if (index1 > index2) {
    int tmp = index1;
    index1  = index2;
    index2  = tmp;
  }

  // Check if range is completely outside index events.
  const Scalar indexOfLast = start_ + (numEvents_ - 1) * stride_;
  if (index2 < start_ || indexOfLast < index1) return false;

  const int strideJustBeforeIndex1 = std::min(
      int(numEvents_ - 1), std::max(int(0), int((index1 - start_) / stride_)));

  const int strideJustAfterIndex2 = std::max(
      int(0), std::min(int(numEvents_ - 1), int((index2 - start_) / stride_)));

  for (int i = strideJustBeforeIndex1; i <= strideJustAfterIndex2; i++) {
    const int indexEvent = start_ + i * stride_;
    if (index1 < indexEvent && indexEvent <= index2) return true;
  }

  return false;
}

template <class Scalar>
void TimeEventRangeIndex<Scalar>::describe(
    Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, "TimeEventRangeIndex");
  l_out->setOutputToRootOnly(0);

  *l_out << "TimeEventRangeIndex:"
         << "\n"
         << "  name       = " << this->getName() << "\n"
         << "  Type       = " << this->getType() << "\n"
         << "  start_     = " << start_ << "\n"
         << "  stop_      = " << stop_ << "\n"
         << "  stride_    = " << stride_ << "\n"
         << "  numEvents_ = " << numEvents_ << std::endl;
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
TimeEventRangeIndex<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Event Range Index");

  pl->setName(this->getName());
  pl->set("Name", this->getName());
  pl->set("Type", "Range Index");

  pl->set("Start Index", getIndexStart(), "Start of Index range");
  pl->set("Stop Index", getIndexStop(), "Stop of Index range");
  pl->set("Stride Index", getIndexStride(), "Stride of Index range");

  return pl;
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<TimeEventRangeIndex<Scalar> > createTimeEventRangeIndex(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto teri = Teuchos::rcp(new TimeEventRangeIndex<Scalar>());
  if (pl == Teuchos::null) return teri;  // Return default TimeEventRangeIndex.

  TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Type", "Range Index") != "Range Index",
      std::logic_error,
      "Error - Time Event Type != 'Range Index'.  (='" +
          pl->get<std::string>("Type") + "')\n");

  pl->validateParametersAndSetDefaults(*teri->getValidParameters());

  teri->setName(pl->get("Name", "From createTimeEventRangeIndex"));
  teri->setIndexStart(pl->get("Start Index", teri->getIndexStart()));
  teri->setIndexStop(pl->get("Stop Index", teri->getIndexStop()));
  teri->setIndexStride(pl->get("Stride Index", teri->getIndexStride()));

  return teri;
}

}  // namespace Tempus
#endif  // Tempus_TimeEventRangeIndex_impl_hpp
