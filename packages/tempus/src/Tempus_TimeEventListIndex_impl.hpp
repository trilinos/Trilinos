// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventListIndex_impl_hpp
#define Tempus_TimeEventListIndex_impl_hpp


namespace Tempus {

template<class Scalar>
TimeEventListIndex<Scalar>::TimeEventListIndex()
{
  this->setName("TimeEventListIndex");
}


template<class Scalar>
TimeEventListIndex<Scalar>::TimeEventListIndex(
  std::string name, std::vector<int> indexList)
{
  this->setName(name);
  this->setIndexList(indexList);
}


template<class Scalar>
void TimeEventListIndex<Scalar>::setIndexList(std::vector<int> indexList, bool sort)
{
  indexList_ = indexList;
  if (sort) {
    std::sort(indexList_.begin(), indexList_.end());
    indexList_.erase(std::unique(
      indexList_.begin(), indexList_.end()), indexList_.end());
  }
}


template<class Scalar>
void TimeEventListIndex<Scalar>::addIndex(int index)
{
  if (indexList_.size() == 0) {
    indexList_.push_back(index);
    return;
  }

  std::vector<int>::iterator it;
  it = std::find(indexList_.begin(), indexList_.end(), index);
  // Check if index is already in list.
  if (it != indexList_.end()) return;

  it = std::upper_bound(indexList_.begin(), indexList_.end(), index);
  indexList_.insert(it, index);
}


template<class Scalar>
bool TimeEventListIndex<Scalar>::isIndex(int index) const
{
  return (indexToNextEvent(index) == 0);
}


template<class Scalar>
int TimeEventListIndex<Scalar>::indexToNextEvent(int index) const
{
  return indexOfNextEvent(index) - index;  // Neg. indicating in the past.
}


template<class Scalar>
int TimeEventListIndex<Scalar>::indexOfNextEvent(int index) const
{
  if (indexList_.size() == 0) return this->getDefaultIndex();

  // Check if before first event.
  if (indexList_.front() >= index) return indexList_.front();

  // Check if after last event.
  if (indexList_.back() <= index) return indexList_.back();

  std::vector<int>::const_iterator it =
    std::upper_bound(indexList_.begin(), indexList_.end(), index);

  // Check if left-side index event
  const Scalar indexOfLeftEvent = *(it-1);
  if (indexOfLeftEvent == index) return indexOfLeftEvent;

  // Otherwise it is the next event.
  return *it;
}


template<class Scalar>
bool TimeEventListIndex<Scalar>::eventInRangeIndex(int index1, int index2) const
{
  if (index1 > index2) {
    int tmp = index1;
    index1 = index2;
    index2 = tmp;
  }

  if (indexList_.size() == 0) return false;

  // Check if range is completely outside index events.
  if (index2 < indexList_.front() || indexList_.back() < index1) return false;

  Scalar indexEvent1 = indexOfNextEvent(index1);
  Scalar indexEvent2 = indexOfNextEvent(index2);
  // Check if the next index event is different for the two indices.
  if (indexEvent1 != indexEvent2) return true;

  // Check if indices bracket index event.
  if (index1 <= indexEvent1 && indexEvent1 <= index2) return true;

  return false;
}


template<class Scalar>
void TimeEventListIndex<Scalar>::describe() const
{
  Teuchos::RCP<Teuchos::FancyOStream> out =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  *out << "TimeEventListIndex:" << "\n"
       << "name       = " << this->getName() << "\n"
       << "IndexList_ = " << std::endl;
  for (auto it = indexList_.begin(); it != indexList_.end()-1; ++it)
    *out << *it << ", ";
  *out << *(indexList_.end()-1) << "\n";
}


} // namespace Tempus
#endif // Tempus_TimeEventListIndex_impl_hpp
