// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventRangeIndex_impl_hpp
#define Tempus_TimeEventRangeIndex_impl_hpp


namespace Tempus {

template<class Scalar>
TimeEventRangeIndex<Scalar>::TimeEventRangeIndex()
  : start_(0), stop_(0), stride_(1)
{
  this->setName("TimeEventRangeIndex");
  setNumEvents();
}


template<class Scalar>
TimeEventRangeIndex<Scalar>::TimeEventRangeIndex(
  std::string name, int start, int stop, int stride)
{
  this->setName(name);
  this->setIndexRange(start, stop, stride);
}


template<class Scalar>
void TimeEventRangeIndex<Scalar>::setIndexStart(int start)
{
  start_ = start;
  if (stop_ < start_) {
    stop_ = start_;
  }
  setNumEvents();
}


template<class Scalar>
void TimeEventRangeIndex<Scalar>::setIndexStop(int stop)
{
  stop_ = stop;
  if (start_ > stop_) {
    start_ = stop_;
    setIndexStride(1);
  }
  setNumEvents();
}


template<class Scalar>
void TimeEventRangeIndex<Scalar>::setIndexStride(int stride)
{
  stride_ = stride;
  if (stride_ < 1) {
    stride_ = 1;
  } else if (stride_ > (stop_ - start_)) {
    stride_ = stop_ - start_;
  }
  setNumEvents();
}


template<class Scalar>
void TimeEventRangeIndex<Scalar>::setNumEvents()
{
  if (stride_ == 0 || start_ == stop_)
    numEvents_ = 1;
  else
    numEvents_ = int((stop_ - start_) / stride_) + 1;
}


template<class Scalar>
bool TimeEventRangeIndex<Scalar>::isIndex(int index) const
{
  return (indexToNextEvent(index) == 0);
}


template<class Scalar>
int TimeEventRangeIndex<Scalar>::indexToNextEvent(int index) const
{
  return indexOfNextEvent(index) - index;
}


template<class Scalar>
int TimeEventRangeIndex<Scalar>::indexOfNextEvent(int index) const
{
  // Check if before or equal to the first index.
  if (index <= start_) return start_;

  const int indexOfLast = start_ + (numEvents_-1) * stride_;
  // Check if after or equal to last index.
  if (indexOfLast <= index) return indexOfLast;

  // check if index is an event.
  if ((index - start_) % stride_ == 0) return index;

  const int numStrides = (index - start_) / stride_ + 1;
  const Scalar indexOfNext = start_ + numStrides * stride_;
  return indexOfNext;
}


template<class Scalar>
bool TimeEventRangeIndex<Scalar>::eventInRangeIndex(int index1, int index2) const
{
  if (index1 > index2) {
    int tmp = index1;
    index1 = index2;
    index2 = tmp;
  }

  const Scalar indexOfLast = start_ + (numEvents_-1) * stride_;
  if (index2 < start_ || indexOfLast < index1) return false;

  Scalar indexEvent1 = indexOfNextEvent(index1);
  Scalar indexEvent2 = indexOfNextEvent(index2);
  // Check if the next index event is different for the two indices.
  if (indexEvent1 != indexEvent2) return true;

  // Check if indices bracket index event.
  if (index1 <= indexEvent1 && indexEvent1 <= index2) return true;

  return false;
}


template<class Scalar>
void TimeEventRangeIndex<Scalar>::describe() const
{
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "TimeEventRange:" << "\n"
       << "name       = " << this->getName() << "\n"
       << "start_     = " << start_     << "\n"
       << "stop_      = " << stop_      << "\n"
       << "stride_    = " << stride_    << "\n"
       << "numEvents_ = " << numEvents_ << std::endl;
}



} // namespace Tempus
#endif // Tempus_TimeEventRangeIndex_impl_hpp
