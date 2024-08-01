//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventComposite_decl_hpp
#define Tempus_TimeEventComposite_decl_hpp

#include "Teuchos_Time.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_config.hpp"
#include "Tempus_TimeEventBase.hpp"
#include "Tempus_TimeEventComposite.hpp"
#include "Tempus_TimeEventRange.hpp"
#include "Tempus_TimeEventRangeIndex.hpp"
#include "Tempus_TimeEventList.hpp"
#include "Tempus_TimeEventListIndex.hpp"

namespace Tempus {

/** \brief This composite TimeEvent loops over added TimeEvents.
 *
 *  Individual TimeEvents are executed in the order in which they
 *  were added.
 */
template <class Scalar>
class TimeEventComposite : virtual public TimeEventBase<Scalar> {
 public:
  /// Default Constructor
  TimeEventComposite()
  {
    this->setType("Composite");
    this->setName("TimeEventComposite");
  }

  /// Construct with full argument list of data members.
  TimeEventComposite(std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > te,
                     std::string name = "TimeEventComposite")
  {
    this->setType("Composite");
    this->setName(name);
    this->setTimeEvents(te);
  }

  /// Destructor
  virtual ~TimeEventComposite() {}

  /// \name Basic methods
  //@{
  /// Get a copy of the current set of TimeEvents.
  virtual std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > getTimeEvents()
      const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > te = timeEvents_;
    return te;
  }

  /** \brief Set the TimeEvents.
   *
   *  This completely replaces the current set of TimeEvents with
   *  the input TimeEvents
   *
   *  \param te [in] The input set of TimeEvents.
   */
  virtual void setTimeEvents(
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > te)
  {
    timeEvents_ = te;
  }

  /** \brief Test if time is near a TimeEvent (within tolerance).
   *
   *  Return true if one of the TimeEvents in the composite is within
   *  tolerance of the input time.
   *
   *  \param time [in] The input time.
   *  \return True if time is near an event (within absolute tolerance).
   */
  virtual bool isTime(Scalar time) const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents;
    return isTime(time, timeEvents);
  }

  /** \brief Test if time is near a TimeEvent (within tolerance) plus the
   * constraining TimeEvent(s).
   *
   *  Return true if one of the TimeEvents in the composite is near
   *  the input time, and the constraining TimeEvent(s) so additional
   *  details about the event can be queried.
   *
   *  \param time       [in]  The input time.
   *  \param timeEvents [out] Vector of constraining TimeEvents.
   *  \return True if time is near an event (within absolute tolerance).
   */
  virtual bool isTime(
      Scalar time,
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    timeEvents.clear();
    for (auto& e : timeEvents_) {
      if (e->isTime(time)) timeEvents.push_back(e);
    }
    return (!timeEvents.empty());
  }

  /** \brief How much time until the next event.
   *
   *  Return the amount of time to the next event (i.e., time of
   *  next event minus the input time).
   *
   *  \param time      [in] The input time.
   *  \return The time to the next event.
   */
  virtual Scalar timeToNextEvent(Scalar time) const
  {
    return timeOfNextEvent(time) - time;
  }

  /** \brief How much time until the next event plus the constraining
   * TimeEvent(s).
   *
   *  Return the amount of time to the next event (i.e., time of
   *  next event minus the input time), and the constraining TimeEvent
   *  so additional details about the event can be queried.
   *
   *  \param time       [in]  The input time.
   *  \param timeEvents [out] Vector of constraining TimeEvent.
   *  \return The time to the next event.
   */
  virtual Scalar timeToNextEvent(
      Scalar time,
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    return timeOfNextEvent(time, timeEvents) - time;
  }

  /** \brief Return the time of the next event following the input time.
   *
   *  See timeOfNextEvent(time, timeEvent).
   *
   *  \param time [in] Input time.
   *  \return Time of the next event.
   */
  virtual Scalar timeOfNextEvent(Scalar time) const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents;
    return timeOfNextEvent(time, timeEvents);
  }

  /** \brief Return the time of the next time event and constraining
   * TimeEvent(s).
   *
   *  Returns the time of the next event that follows the input time.
   *  If the input time is before all events, the time of the first
   *  event is returned.  If the input time is after all events, the
   *  default time (a time in the distant future) is returned.  If the
   *  input time is an event time, the time of the next event is returned.
   *
   *  For TimeEventComposite, find the next time event from all the
   *  TimeEvents in the composite.
   *
   *  Additionally, output the constraining TimeEvents
   *  so additional details about the events can be queried.
   *
   *  \param time       [in]  Input time.
   *  \param timeEvents [out] Constraining TimeEvent.
   *  \return Time of the next event.
   */
  virtual Scalar timeOfNextEvent(
      Scalar time,
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    timeEvents.clear();
    typedef std::pair<Scalar, Teuchos::RCP<TimeEventBase<Scalar> > > TEPAIR;
    std::vector<TEPAIR> timeEventPair;
    for (auto& e : timeEvents_)
      timeEventPair.push_back(std::make_pair(e->timeOfNextEvent(time), e));

    if (timeEventPair.empty()) return this->getDefaultTime();

    auto compare = [](TEPAIR a, TEPAIR b) { return a.first < b.first; };
    std::stable_sort(timeEventPair.begin(), timeEventPair.end(), compare);

    // The first one is the "time of next event".
    Scalar tone = timeEventPair.front().first;

    // Check if there are multiple events that match "time of next event".
    for (auto it = timeEventPair.begin(); it != timeEventPair.end(); ++it) {
      if ((*it).second->isTime(tone)) timeEvents.push_back((*it).second);
    }

    return tone;
  }

  /** \brief Test if an event occurs within the time range.
   *
   *  Find if an event is within the input range,
   *  (time1 < event-absTol and timeEvent-absTol <= time2),
   *  including the event's absolute tolerance.  For
   *  TimeEventComposite, test each TimeEvent to determine if
   *  the input time is within the range.
   *
   *  \param time1 [in] Input time of one end of the range.
   *  \param time2 [in] Input time of the other end of the range.
   *
   *  \return True if a time event is within the range.
   */
  virtual bool eventInRange(Scalar time1, Scalar time2) const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents;
    return eventInRange(time1, time2, timeEvents);
  }

  /** \brief Test if an event occurs within the time range plus the constraining
   * TimeEvent(s).
   *
   *  Find if an event is within the input range,
   *  (time1 < event-absTol and timeEvent-absTol <= time2),
   *  including the event's absolute tolerance.  For
   *  TimeEventComposite, test each TimeEvent to determine if
   *  the input time is within the range.
   *
   *  Additionally, the constraining TimeEvents are sorted by "time of
   *  next event", and returned, so additional details about the events
   *  can be queried.
   *
   *  \param time1 [in] Input time of one end of the range.
   *  \param time2 [in] Input time of the other end of the range.
   *  \param timeEvents [out] Vector of sorted constraining TimeEvent(s).
   *
   *  \return True if a time event is within the range.
   */
  virtual bool eventInRange(
      Scalar time1, Scalar time2,
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    typedef std::pair<Scalar, Teuchos::RCP<TimeEventBase<Scalar> > > TEPAIR;
    std::vector<TEPAIR> timeEventPair;
    for (auto& e : timeEvents_) {
      if (e->eventInRange(time1, time2))
        timeEventPair.push_back(std::make_pair(e->timeOfNextEvent(time1), e));
    }

    auto compare = [](TEPAIR a, TEPAIR b) { return a.first < b.first; };
    std::stable_sort(timeEventPair.begin(), timeEventPair.end(), compare);

    timeEvents.clear();
    for (auto& e : timeEventPair) timeEvents.push_back(e.second);

    return (!timeEvents.empty());
  }

  /** \brief Test if index is a time event.
   *
   *  Return true if one of the TimeEvents in the composite is
   *  the input index.
   *
   *  \param index [in] The input index.
   *  \return True if index is an event.
   */
  virtual bool isIndex(int index) const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents;
    return isIndex(index, timeEvents);
  }

  /** \brief Test if index is a time event plus the constraining TimeEvent(s).
   *
   *  Return true if one of the TimeEvents indices in the composite is near
   *  the input index, and the constraining TimeEvent so additional details
   *  about the event can be queried.
   *
   *  \param index      [in]  The input index.
   *  \param timeEvents [out] Vector of constraining TimeEvents.
   *  \return True if index is an event.
   */
  virtual bool isIndex(
      int index,
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    timeEvents.clear();
    for (auto& e : timeEvents_) {
      if (e->isIndex(index)) timeEvents.push_back(e);
    }
    return (!timeEvents.empty());
  }

  /** \brief How many indices until the next event.
   *
   *  \param index [in] The input index.
   *  \return The number of steps (indices) to the next event.
   */
  virtual int indexToNextEvent(int index) const
  {
    return indexOfNextEvent(index) - index;
  }

  /** \brief How many indices until the next event.
   *
   *  \param index      [in] The input index.
   *  \param timeEvents [out] The constraining TimeEvent.
   *  \return The number of steps (indices) to the next event.
   */
  virtual int indexToNextEvent(
      int index,
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    return indexOfNextEvent(index, timeEvents) - index;
  }

  /** \brief Return the index of the next event following the input index.
   *
   *  Returns the index of the next event that follows the input index.
   *  If the input index is before all events, the index of the first
   *  event is returned.  If the input index is after all events, the
   *  default index (an index in the distant future) is returned.  If the
   *  input index is an event index, the index of the next event is returned.
   *
   *  \param index     [in] Input index.
   *  \return Index of the next event.
   */
  virtual int indexOfNextEvent(int index) const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents;
    return indexOfNextEvent(index, timeEvents);
  }

  /** \brief Return the index of the next event following the input index plus
   * the constraining TimeEvent(s).
   *
   *  Returns the index of the next event that follows the input index.
   *  If the input index is before all events, the index of the first
   *  event is returned.  If the input index is after all events, the
   *  default index (a index in the distant future) is returned.  If the
   *  input index is an event index, the index of the next event is returned.
   *
   *  \param index      [in] Input index.
   *  \param timeEvents [out] Vector of constraining TimeEvent(s).
   *  \return Index of the next event.
   */
  virtual int indexOfNextEvent(
      int index,
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    timeEvents.clear();
    typedef std::pair<int, Teuchos::RCP<TimeEventBase<Scalar> > > TEPAIR;
    std::vector<TEPAIR> timeEventPair;
    for (auto& e : timeEvents_)
      timeEventPair.push_back(std::make_pair(e->indexOfNextEvent(index), e));

    if (timeEventPair.size() == 0) return this->getDefaultIndex();

    auto compare = [](TEPAIR a, TEPAIR b) { return a.first < b.first; };
    std::stable_sort(timeEventPair.begin(), timeEventPair.end(), compare);

    // The first one is the "index of next event".
    int ione = timeEventPair.front().first;

    // Check if there are multiple events that match "index of next event".
    for (auto it = timeEventPair.begin(); it != timeEventPair.end(); ++it) {
      if ((*it).second->isIndex(ione)) timeEvents.push_back((*it).second);
    }

    return ione;
  }

  /** \brief Test if an event occurs within the index range.
   *
   *  Find if an event is within the input range, inclusively
   *  ( index1 <= event <= index2 ).  This may require testing
   *  each event in the TimeEvent.
   *
   *  \param index1 [in] Input index of one end of the range.
   *  \param index2 [in] Input index of the other end of the range.
   *
   *  \return True if a index event is within the range.
   */
  virtual bool eventInRangeIndex(int index1, int index2) const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents;
    return eventInRangeIndex(index1, index2, timeEvents);
  }

  /** \brief Test if an event occurs within the index range plus the
   * constraining TimeEvent(s).
   *
   *  Find if an event is within the input range, inclusively
   *  ( index1 <= event <= index2 ).  This may require testing
   *  each event in the TimeEvent.
   *
   *  Additionally, the constraining TimeEvents are sorted by "index of
   *  next event", and returned, so additional details about the events
   *  can be queried.
   *
   *  \param index1 [in] Input index of one end of the range.
   *  \param index2 [in] Input index of the other end of the range.
   *  \param timeEvents [out] Vector of constraining TimeEvents.
   *
   *  \return True if a index event is within the range.
   */
  virtual bool eventInRangeIndex(
      int index1, int index2,
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    typedef std::pair<int, Teuchos::RCP<TimeEventBase<Scalar> > > TEPAIR;
    std::vector<TEPAIR> timeEventPair;
    for (auto& e : timeEvents_) {
      if (e->eventInRangeIndex(index1, index2))
        timeEventPair.push_back(std::make_pair(e->indexOfNextEvent(index1), e));
    }

    auto compare = [](TEPAIR a, TEPAIR b) { return a.first < b.first; };
    std::stable_sort(timeEventPair.begin(), timeEventPair.end(), compare);

    timeEvents.clear();
    for (auto& e : timeEventPair) timeEvents.push_back(e.second);

    return (!timeEvents.empty());
  }

  /** \brief Return the largest absolute tolerance from all the TimeEvents.
   *
   *  \return The largest absolute tolerance of all the TimeEvents.
   */
  virtual Scalar getAbsTol() const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents;
    return getAbsTol(timeEvents);
  }

  /** \brief Return the largest absolute tolerance from all the TimeEvents plus
   * the constraining TimeEvent(s).
   *
   *  All the constraining TimeEvents have the same largest absolute
   *  tolerance (within numerical tolerance).
   *
   *  \param timeEvents [out] Vector of constraining TimeEvent(s).
   *  \return The largest absolute tolerance of all the TimeEvents.
   */
  virtual Scalar getAbsTol(
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    timeEvents.clear();
    Scalar largestAbsTol = timeEvents_.front()->getAbsTol();
    timeEvents.push_back(timeEvents_.front());
    for (auto& e : timeEvents_)
      if (e->getAbsTol() > largestAbsTol) largestAbsTol = e->getAbsTol();

    for (auto& e : timeEvents_)
      if (e->getAbsTol() - largestAbsTol < largestAbsTol * 1.0e-14)
        timeEvents.push_back(e);

    return largestAbsTol;
  }

  /** \brief Return if the time events need to be landed on exactly.
   *
   *  Will return true if any of the events requires to be landed on exactly.
   *
   *  \param LOE [in] Flag indicating if TimeEvent should land on the event
   * exactly.
   */
  virtual bool getLandOnExactly() const
  {
    std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents;
    return getLandOnExactly(timeEvents);
  }

  /** \brief Return if the time events need to be landed on exactly plus the
   * constraining TimeEvent(s).
   *
   *  Will return true if any of the events requires to be landed on exactly.
   *  All the constraining TimeEvents that require to be landed on exactly
   *  will be returned through the input vector of TimeEvents.
   *
   *  \param timeEvents [out] Vector of constraining TimeEvent(s).
   *  \return LOE Flag indicating if TimeEvent should land on the event exactly.
   */
  virtual bool getLandOnExactly(
      std::vector<Teuchos::RCP<TimeEventBase<Scalar> > >& timeEvents) const
  {
    timeEvents.clear();
    for (auto& e : timeEvents_)
      if (e->getLandOnExactly()) timeEvents.push_back(e);

    return (!timeEvents.empty());
  }
  //@}

  /** \brief Add TimeEvent to the TimeEvent vector.
   *
   *  Add the TimeEvent to the composite vector.  If the
   *  TimeEvent is already in the composite (based on the
   *  TimeEvent's name), the input TimeEvent will replace
   *  the one in the composite.
   *
   *  \param timeEvent [in] The input TimeEvent.
   */
  void add(Teuchos::RCP<TimeEventBase<Scalar> > timeEvent)
  {
    std::string name                        = timeEvent->getName();
    Teuchos::RCP<TimeEventBase<Scalar> > te = find(name);
    if (!te.is_null()) {
      // auto l_out = Teuchos::fancyOStream( Teuchos::rcpFromRef(std::cout) );
      // Teuchos::OSTab ostab(*l_out, 2, "TimeEventComposite::add");
      // l_out->setOutputToRootOnly(0);

      //*l_out << "TimeEventComposite::add: Replacing Time Event, "
      //       << name << "." << "\n";
      remove(name);
    }
    timeEvents_.push_back(timeEvent);
  }

  /** \brief Remove TimeEvent based on name.
   *
   *  If the TimeEvent is not in the composite based on the
   *  TimeEvent's name, nothing is done.
   *
   *  \param name [in] The name of the TimeEvent to remove.
   */
  void remove(std::string name)
  {
    for (std::size_t i = 0; i < timeEvents_.size(); ++i) {
      if (timeEvents_[i]->getName() == name) {
        timeEvents_.erase(timeEvents_.begin() + i);
        break;
      }
    }
    // Did not find 'name', so did nothing.
  }

  /** \brief Find TimeEvent based on name.
   *
   *  If the TimeEvent is not found, Teuchos::null is returned.
   *
   *  \param name [in] The name of the TimeEvent to find.
   *  \return RCP of TimeEvent with matching input name.
   */
  Teuchos::RCP<TimeEventBase<Scalar> > find(std::string name)
  {
    for (std::size_t i = 0; i < timeEvents_.size(); ++i)
      if (timeEvents_[i]->getName() == name) return timeEvents_[i];

    return Teuchos::null;
  }

  /// Clear the TimeEvent vector.
  void clear() { timeEvents_.clear(); }

  /// Return the size of the TimeEvent vector.
  std::size_t getSize() const { return timeEvents_.size(); }

  /// Return a string of the names of Time Events (comma separated).
  std::string getTimeEventNames() const
  {
    std::stringstream tecList;
    for (std::size_t i = 0; i < timeEvents_.size(); ++i) {
      tecList << timeEvents_[i]->getName();
      if (i < timeEvents_.size() - 1) tecList << ", ";
    }
    return tecList.str();
  }

  /// Describe member data.
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const
  {
    auto l_out = Teuchos::fancyOStream(out.getOStream());
    Teuchos::OSTab ostab(*l_out, 2, "TimeEventComposite");
    l_out->setOutputToRootOnly(0);

    *l_out << "TimeEventComposite:"
           << "\n"
           << "  name                 = " << this->getName() << "\n"
           << "  Type                 = " << this->getType() << "\n"
           << "  Number of TimeEvents = " << this->getSize() << "\n"
           << "  Time Events          = " << this->getTimeEventNames()
           << std::endl;
    *l_out << "--------------------------------------------" << std::endl;
    for (auto& e : timeEvents_) {
      (*e).describe(*l_out, verbLevel);
      *l_out << "--------------------------------------------" << std::endl;
    }
  }

  /** \brief Return a valid ParameterList with current settings.
   *
   *  The returned ParameterList will contain the current parameters
   *  and can be used to reconstruct the TimeEventComposite using
   *  createTimeEventComposite(...).  The ParameterList will have
   *  the TimeEventComposite parameters along with all the parameters
   *  for the TimeEvents contained in the composite.
   *
   * \return Teuchos::ParameterList of TimeEventComposite.
   */
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl =
        Teuchos::parameterList("Time Event Composite");

    pl->setName(this->getName());
    pl->set("Name", this->getName());
    pl->set("Type", this->getType());
    pl->set<std::string>("Time Events", this->getTimeEventNames());

    for (auto& s : timeEvents_) pl->set(s->getName(), *s->getValidParameters());

    return pl;
  }

 protected:
  std::vector<Teuchos::RCP<TimeEventBase<Scalar> > > timeEvents_;
};

// Nonmember constructor - ParameterList
// ------------------------------------------------------------------------

/** \brief TimeEventComposite nonmember constructor via ParameterList.
 *
 *  If the input ParameterList is Teuchos::null, return a default
 *  TimeEventComposite, which has no TimeEvents but TimeEvents can be added.
 *  A valid ParameterList for TimeEventComposite can be obtained
 *  from TimeEventComposite::getValidParameters().
 *
 *  Limitation: Although possible, nesting TimeEventComposite within a
 *  TimeEventComposite is not a good idea and is not supported in this
 *  constructor.
 *
 *  \param pList [in] The input ParameterList to construct from.
 *  \return Constructed TimeEventComposite.
 */
template <class Scalar>
Teuchos::RCP<TimeEventComposite<Scalar> > createTimeEventComposite(
    Teuchos::RCP<Teuchos::ParameterList> const& pList)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;

  auto tec = Teuchos::rcp(new TimeEventComposite<Scalar>());
  if (pList == Teuchos::null || pList->numParams() == 0) return tec;

  TEUCHOS_TEST_FOR_EXCEPTION(
      pList->get<std::string>("Type", "Composite") != "Composite",
      std::logic_error,
      "Error - Time Event Type != 'Composite'.  (='" +
          pList->get<std::string>("Type") + "')\n");

  tec->setName(pList->get("Name", "From createTimeEventComposite"));

  // string tokenizer
  std::vector<std::string> teList;
  teList.clear();
  std::string str = pList->get<std::string>("Time Events");
  std::string delimiters(",");
  const char* WhiteSpace = " \t\v\r\n";
  // Skip delimiters at the beginning
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find the first delimiter
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);
  while ((pos != std::string::npos) || (lastPos != std::string::npos)) {
    // Found a token, add it to the vector
    std::string token = str.substr(lastPos, pos - lastPos);

    std::size_t start = token.find_first_not_of(WhiteSpace);
    std::size_t end   = token.find_last_not_of(WhiteSpace);
    token =
        (start == end ? std::string() : token.substr(start, end - start + 1));

    teList.push_back(token);
    if (pos == std::string::npos) break;

    lastPos = str.find_first_not_of(delimiters, pos);  // Skip delimiters
    pos     = str.find_first_of(delimiters, lastPos);  // Find next delimiter
  }

  // For each sublist name tokenized, add the TimeEvent.
  for (auto teName : teList) {
    RCP<ParameterList> pl =
        Teuchos::rcp(new ParameterList(pList->sublist(teName)));

    auto timeEventType = pl->get<std::string>("Type", "Unknown");
    if (timeEventType == "Range") {
      tec->add(createTimeEventRange<Scalar>(pl));
    }
    else if (timeEventType == "Range Index") {
      tec->add(createTimeEventRangeIndex<Scalar>(pl));
    }
    else if (timeEventType == "List") {
      tec->add(createTimeEventList<Scalar>(pl));
    }
    else if (timeEventType == "List Index") {
      tec->add(createTimeEventListIndex<Scalar>(pl));
    }
    else {
      RCP<Teuchos::FancyOStream> out =
          Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out->setOutputToRootOnly(0);
      Teuchos::OSTab ostab(out, 1, "createTimeEventComposite()");
      *out << "Warning -- createTimeEventComposite() - Unknown Time Event "
              "Type!\n"
           << "'Type' = '" << timeEventType << "'\n"
           << "Should call add() with this "
           << "(app-specific?) Time Event.\n"
           << std::endl;
    }
  }

  if (tec->getSize() == 0) {
    RCP<Teuchos::FancyOStream> out =
        Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    Teuchos::OSTab ostab(out, 1, "createTimeEventComposite()");
    *out << "Warning -- createTimeEventComposite() - Did not\n"
         << "           find/recognize any TimeEvents to create!\n"
         << "           If there is a app-specific TimeEvent,\n"
         << "           explicitly add it to this TimeEventComposite.\n"
         << std::endl;
  }

  return tec;
}

}  // namespace Tempus

#endif  // Tempus_TimeEventComposite_decl_hpp
