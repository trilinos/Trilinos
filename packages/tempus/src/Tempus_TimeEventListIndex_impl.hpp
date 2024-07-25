//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventListIndex_impl_hpp
#define Tempus_TimeEventListIndex_impl_hpp

namespace Tempus {

template <class Scalar>
TimeEventListIndex<Scalar>::TimeEventListIndex()
{
  this->setType("List Index");
  this->setName("TimeEventListIndex");
}

template <class Scalar>
TimeEventListIndex<Scalar>::TimeEventListIndex(std::vector<int> indexList,
                                               std::string name)
{
  this->setType("List Index");
  if (name == "" && !indexList.empty()) {
    std::ostringstream oss;
    oss << "TimeEventListIndex (" << indexList_.front() << ", ... ,"
        << indexList_.back() << ")";
    this->setName(oss.str());
  }
  else {
    this->setName(name);
  }

  this->setIndexList(indexList);
}

template <class Scalar>
void TimeEventListIndex<Scalar>::setIndexList(std::vector<int> indexList,
                                              bool sort)
{
  indexList_ = indexList;
  if (sort) {
    std::sort(indexList_.begin(), indexList_.end());
    indexList_.erase(std::unique(indexList_.begin(), indexList_.end()),
                     indexList_.end());
  }
}

template <class Scalar>
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

template <class Scalar>
bool TimeEventListIndex<Scalar>::isIndex(int index) const
{
  return (std::find(indexList_.begin(), indexList_.end(), index) !=
          indexList_.end());
}

template <class Scalar>
int TimeEventListIndex<Scalar>::indexToNextEvent(int index) const
{
  return indexOfNextEvent(index) - index;  // Neg. indicating in the past.
}

template <class Scalar>
int TimeEventListIndex<Scalar>::indexOfNextEvent(int index) const
{
  if (indexList_.size() == 0) return this->getDefaultIndex();

  // Check if before first event.
  if (index < indexList_.front()) return indexList_.front();

  // Check if after last event.
  if (index >= indexList_.back()) return this->getDefaultIndex();

  std::vector<int>::const_iterator it =
      std::upper_bound(indexList_.begin(), indexList_.end(), index);

  return int(*it);
}

template <class Scalar>
bool TimeEventListIndex<Scalar>::eventInRangeIndex(int index1, int index2) const
{
  if (index1 > index2) {
    int tmp = index1;
    index1  = index2;
    index2  = tmp;
  }

  if (indexList_.size() == 0) return false;

  // Check if range is completely outside index events.
  if (index2 < indexList_.front() || indexList_.back() < index1) return false;

  Scalar indexEvent1 = indexOfNextEvent(index1);
  Scalar indexEvent2 = indexOfNextEvent(index2);
  // Check if the next index event is different for the two indices.
  if (indexEvent1 != indexEvent2) return true;

  // Check if indices bracket index event.
  if (index1 < indexEvent1 && indexEvent1 <= index2) return true;

  return false;
}

template <class Scalar>
void TimeEventListIndex<Scalar>::describe(
    Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, "TimeEventListIndex");
  l_out->setOutputToRootOnly(0);

  *l_out << "TimeEventListIndex:"
         << "\n"
         << "  name       = " << this->getName() << "\n"
         << "  Type       = " << this->getType() << "\n"
         << "  IndexList_ = ";
  if (!indexList_.empty()) {
    for (auto it = indexList_.begin(); it != indexList_.end() - 1; ++it)
      *l_out << *it << ", ";
    *l_out << *(indexList_.end() - 1) << std::endl;
  }
  else {
    *l_out << "<empty>" << std::endl;
  }
}

template <class Scalar>
Teuchos::RCP<const Teuchos::ParameterList>
TimeEventListIndex<Scalar>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Event List Index");

  pl->setName(this->getName());
  pl->set("Name", this->getName());
  pl->set("Type", this->getType());

  std::ostringstream list;
  if (!indexList_.empty()) {
    for (std::size_t i = 0; i < indexList_.size() - 1; ++i)
      list << indexList_[i] << ", ";
    list << indexList_[indexList_.size() - 1];
  }
  pl->set<std::string>("Index List", list.str(),
                       "Comma deliminated list of indices");

  return pl;
}

// Nonmember constructors.
// ------------------------------------------------------------------------

template <class Scalar>
Teuchos::RCP<TimeEventListIndex<Scalar> > createTimeEventListIndex(
    Teuchos::RCP<Teuchos::ParameterList> pl)
{
  auto teli = Teuchos::rcp(new TimeEventListIndex<Scalar>());
  if (pl == Teuchos::null) return teli;  // Return default TimeEventListIndex.

  TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Type", "List Index") != "List Index",
      std::logic_error,
      "Error - Time Event Type != 'List Index'.  (='" +
          pl->get<std::string>("Type") + "')\n");

  pl->validateParametersAndSetDefaults(*teli->getValidParameters());

  teli->setName(pl->get("Name", "From createTimeEventListIndex"));

  std::vector<int> indexList;
  indexList.clear();
  std::string str = pl->get<std::string>("Index List");
  std::string delimiters(",");
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
    std::string token = str.substr(lastPos, pos - lastPos);
    indexList.push_back(int(std::stoi(token)));
    if (pos == std::string::npos) break;

    lastPos = str.find_first_not_of(delimiters, pos);
    pos     = str.find_first_of(delimiters, lastPos);
  }
  teli->setIndexList(indexList);

  return teli;
}

}  // namespace Tempus
#endif  // Tempus_TimeEventListIndex_impl_hpp
