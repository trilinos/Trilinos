// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventComposite_decl_hpp
#define Tempus_TimeEventComposite_decl_hpp

// Teuchos
#include "Teuchos_Time.hpp"

// Tempus
#include "Tempus_TimeEventBase.hpp"


namespace Tempus {


/** \brief This composite TimeEvent loops over added TimeEvents.
 *
 *  Individual TimeEvents are executed in the order in which they
 *  were added.
 */
template<class Scalar>
class TimeEventComposite : virtual public TimeEventBase<Scalar>
{
public:

  /// Constructor
  TimeEventComposite()
  {
    this->setName("TimeEventComposite");
  }

  /// Destructor
  virtual ~TimeEventComposite() {}

  /// \name Basic methods
  //@{
    /// Test if time is near a TimeEvent (within tolerance).
    virtual bool isTime(Scalar time) const
    {
      bool is_time = false;
      for(auto& e : timeEvents_) {
        if (e->isTime(time)) {
          is_time = true;
          break;
        }
      }
      return is_time;
    }

    /// How much time until the next event. Negative indicating the last event is in the past.
    virtual Scalar timeToNextEvent(Scalar time) const
    {
      return timeOfNextEvent(time) - time;
    }

    /// Time of the next event. Negative indicating the last event is in the past.
    virtual Scalar timeOfNextEvent(Scalar time) const
    {
      std::vector< std::pair<Scalar, Scalar> > timeAbs;
      for(auto& e : timeEvents_)
        timeAbs.push_back(std::make_pair(e->timeOfNextEvent(time), e->getAbsTol()));
      std::sort(timeAbs.begin(), timeAbs.end());

      if (timeAbs.size() == 0) return this->getDefaultTime();

      // Check if before or close to first event.
      if (timeAbs.front().first >= time-timeAbs.front().second)
        return timeAbs.front().first;

      // Check if after or close to last event.
      if (timeAbs.back().first <= time+timeAbs.front().second)
        return timeAbs.back().first;

      typename std::vector< std::pair<Scalar, Scalar> >::const_iterator it =
        std::upper_bound(timeAbs.begin(), timeAbs.end(), std::make_pair(time, 0.0));

      // Check if close to left-side time event
      const Scalar timeOfLeftEvent   = (*(it-1)).first;
      const Scalar absTolOfLeftEvent = (*(it-1)).second;
      if (timeOfLeftEvent > time - absTolOfLeftEvent &&
          timeOfLeftEvent < time + absTolOfLeftEvent)
        return timeOfLeftEvent;

      // Otherwise it is the next event.
      return (*it).first;
    }

    /// Test if an event occurs within the time range.
    virtual bool eventInRange(Scalar time1, Scalar time2) const
    {
      bool inRange = false;
      for(auto& e : timeEvents_) {
        if (e->eventInRange(time1, time2)) {
          inRange = true;
          break;
        }
      }
      return inRange;
    }

    /// Test if index is a time event.TimeEventBase
    virtual bool isIndex(int index) const
    {
      bool is_index = false;
      for(auto& e : timeEvents_) {
        if (e->isIndex(index)) {
          is_index = true;
          break;
        }
      }
      return is_index;
    }

    /// How many indices until the next event. Negative indicating the last event is in the past.
    virtual int indexToNextEvent(int index) const
    {
      return indexOfNextEvent(index) - index;
    }

    /// Index of the next event. Negative indicating the last event is in the past.
    virtual int indexOfNextEvent(int index) const
    {
      std::vector<int> indexList;
      for(auto& e : timeEvents_)
        indexList.push_back(e->indexOfNextEvent(index));

      std::sort(indexList.begin(), indexList.end());
      indexList.erase(std::unique(
        indexList.begin(), indexList.end()), indexList.end());

      if (indexList.size() == 0) return this->getDefaultIndex();

      // Check if before first event.
      if (indexList.front() >= index) return indexList.front();

      // Check if after last event.
      if (indexList.back() <= index) return indexList.back();

      std::vector<int>::const_iterator it =
        std::upper_bound(indexList.begin(), indexList.end(), index);

      // Check if left-side index event
      const Scalar indexOfLeftEvent = *(it-1);
      if (indexOfLeftEvent == index) return indexOfLeftEvent;

      // Otherwise it is the next event.
      return *it;

    }

    /// Test if an event occurs within the time range.
    virtual bool eventInRangeIndex(int index1, int index2) const
    {
      bool inRange = false;
      for(auto& e : timeEvents_) {
        if (e->eventInRangeIndex(index1, index2)) {
          inRange = true;
          break;
        }
      }
      return inRange;
    }
  //@}


  // Add TimeEvent to the TimeEvent vector.
  void addTimeEvent(Teuchos::RCP<TimeEventBase<Scalar> > timeEvent)
  {
    timeEvents_.push_back(timeEvent);
  }

  // Clear the TimeEvent vector.
  void clearTimeEvents()
  { timeEvents_.clear(); }

  // Return the size of the TimeEvent vector.
  std::size_t getSize() const { return timeEvents_.size(); }

  /// Describe member data.
  virtual void describe() const
  {
    Teuchos::RCP<Teuchos::FancyOStream> out =
      Teuchos::VerboseObjectBase::getDefaultOStream();
    *out << "TimeEventComposite:" << "\n"
        << "name                 = " << this->getName() << "\n"
        << "Number of TimeEvents = " << timeEvents_.size() << std::endl;
    *out << "--------------------------------------------" << std::endl;
    for(auto& e : timeEvents_) {
      (*e).describe();
      *out << "--------------------------------------------" << std::endl;
    }
  }


protected:

  std::vector<Teuchos::RCP<TimeEventBase<Scalar > > > timeEvents_;

};


} // namespace Tempus

#endif // Tempus_TimeEventComposite_decl_hpp
