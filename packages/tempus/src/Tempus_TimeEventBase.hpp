//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventBase_decl_hpp
#define Tempus_TimeEventBase_decl_hpp

// Teuchos
#include "Teuchos_Time.hpp"
#include "Teuchos_VerboseObject.hpp"

#include "Tempus_config.hpp"

namespace Tempus {

/** \brief This class defines time events which can be used to "trigger" an
 * action.
 *
 *  Time events are points in time and/or timestep index where an
 *  an action should occur, such as, solution output (mesh and solution)
 *  diagnostic output, restart, screen dump, in-situ visualization,
 *  user-specified, or any other action.
 *
 *  This class will store a collection time events, so that an object
 *  may query it and take appropriate action.  Time events (time and
 *  timestep index) can be specified via
 *    - start, stop and stride
 *    - list of events
 */
template <class Scalar>
class TimeEventBase {
 public:
  /// Constructor
  TimeEventBase()
    : timeEventType_("Base"),
      name_("TimeEventBase"),
      defaultTime_(std::numeric_limits<Scalar>::max()),
      defaultTol_(std::numeric_limits<Scalar>::epsilon() * Scalar(100.0)),
      defaultIndex_(std::numeric_limits<int>::max())
  {
  }

  /// Destructor
  virtual ~TimeEventBase() {}

  /// \name Basic methods
  //@{
  /** \brief Test if time is near an event (within tolerance).
   *
   *  For TimeEventBase, always return false since it has no events.
   *  For TimeEvents that are not time based
   *  (e.g., index-based events - TimeEventRangeIndex), they do not
   *  do not have a time, so return false.
   *
   *  \param time [in] The input time.
   *  \return True if time is near an event (within absolute tolerance).
   */
  virtual bool isTime(Scalar time) const { return false; }

  /** \brief How much time until the next event.
   *
   *  For TimeEventBase, the time to the next event is the default time,
   *  since TimeEventBase has no events.
   *  For TimeEvents that are not time based
   *  (e.g., index-based events - TimeEventRangeIndex), they do not
   *  do not have a time to the next event, so return the default time.
   *
   *  \return The time to the next event.
   */
  virtual Scalar timeToNextEvent(Scalar time) const { return defaultTime_; }

  /** \brief Return the time of the next event following the input time.
   *
   *  Returns the time of the next event that follows the input time.
   *  If the input time is before all events, the time of the first
   *  event is returned.  If the input time is after all events, the
   *  default time (a time in the distant future) is returned.  If the
   *  input time is an event time, the time of the next event is returned.
   *
   *  For TimeEventBase, always return the default time.
   *  For TimeEvents that are not time based
   *  (e.g., index-based events - TimeEventRangeIndex), they do not
   *  do not have a time of the next event, so return the default time.
   *
   *  \param time [in] Input time.
   *  \return Time of the next event.
   */
  virtual Scalar timeOfNextEvent(Scalar time) const { return defaultTime_; }

  /** \brief Test if an event occurs within the time range.
   *
   *  Find if an event is within the input range,
   *  (time1 < timeEvent-absTol and timeEvent-absTol <= time2),
   *  including the event's absolute tolerance.  Note, this
   *  does not include time1, but does include time2.
   *
   *  For TimeEventBase, always return false since it has no events.
   *  For TimeEvents that are not time based
   *  (e.g., index-based events - TimeEventRangeIndex), there is not
   *  not event in range, so return false.
   *
   *  \param time1 [in] Input time of one end of the range.
   *  \param time2 [in] Input time of the other end of the range.
   *
   *  \return True if a time event is within the range.
   */
  virtual bool eventInRange(Scalar time1, Scalar time2) const { return false; }

  /** \brief Test if index is an event.
   *
   *  For TimeEventBase, always return false since it has no events.
   *  For TimeEvents that are not index based
   *  (e.g., time-based events - TimeEventRange), they do not
   *  do not have an index, so return false.
   *
   *  \param index [in] The input index.
   *  \return True if index is an event.
   */
  virtual bool isIndex(int index) const { return false; }

  /** \brief How many indices until the next event.
   *
   *  For TimeEventBase, the index to the next event is the default index.
   *  For TimeEvents that are not index based
   *  (e.g., time-based events - TimeEventRange), they do not
   *  do not have an index to the next event, so return the default index.
   *
   *  \return The index to the next event.
   */
  virtual int indexToNextEvent(int index) const { return defaultIndex_; }

  /** \brief Return the index of the next event following the input index.
   *
   *  Returns the index of the next event that follows the input index.
   *  If the input index is before all events, the index of the first
   *  event is returned.  If the input index is after all events, the
   *  default index (an index in the distant future) is returned.  If the
   *  input index is an event index, the index of the next event is returned.
   *
   *  For TimeEventBase, always return the default index.
   *  For TimeEvents that are not index based
   *  (e.g., time-based events - TimeEventRange), they do not
   *  do not have an index of the next event, so return default index.
   *
   *  \param index [in] Input index.
   *  \return Index of the next event.
   */
  virtual int indexOfNextEvent(int index) const { return defaultIndex_; }

  /** \brief Test if an event occurs within the index range.
   *
   *  Find if an event is within the input range,
   *  ( index1 < event <= index2 ).
   *
   *  For TimeEventBase, always return false since it has no events.
   *  For TimeEvents that are not index based
   *  (e.g., time-based events - TimeEventRange), they do not
   *  do not have an index in range, so return false.
   *
   *  \param index1 [in] Input index of one end of the range.
   *  \param index2 [in] Input index of the other end of the range.
   *
   *  \return True if a index event is within the range.
   */
  virtual bool eventInRangeIndex(int index1, int index2) const { return false; }

  /// Describe member data.
  virtual void describe(Teuchos::FancyOStream &out,
                        const Teuchos::EVerbosityLevel verbLevel) const
  {
    auto l_out = Teuchos::fancyOStream(out.getOStream());
    Teuchos::OSTab ostab(*l_out, 2, "TimeEventBase");
    l_out->setOutputToRootOnly(0);

    *l_out << "TimeEventBase name = " << getName() << std::endl;
  }

  /** \brief Return the absolute tolerance.
   *
   *  The absolute tolerance is primarily used to determine
   *  if two times are equal (i.e., t1 is equal to t2, if
   *  t2-absTol < t1 < t2+absTol).
   *
   *  For TimeEventBase and TimeEvents that are not time based
   *  (e.g., index-based events - TimeEventRangeIndex), the absolute
   *  tolerance is the default tolerance.
   *
   *  \return The absolute tolerance.
   */
  virtual Scalar getAbsTol() const { return defaultTol_; }

  /** \brief Return if the time events need to be landed on exactly.
   *
   *  If true, this returns whether the time events need to be landed
   *  on exactly, e.g., the time step needs to be adjusted so the
   *  solution is determined at the time event.
   *
   *  If false, this indicates that time event will still occur but
   *  can be stepped over without changing the time step.
   *
   *  This is for TimeEvents that are time based, e.g., TimeEventRange
   *  and TimeEventList.  For TimeEvents that are not time based (e.g.,
   *  index-based events - TimeEventRangeIndex and TimeEventListIndex),
   *  there is no need to flag this, so return false.
   *
   *  \return LOE Flag indicating if TimeEvent should land on the time event
   * exactly.
   */
  virtual bool getLandOnExactly() const { return false; }

  //@}

  /// \name Accessor methods
  //@{
  /** \brief Return the name of the TimeEvent.
   *
   *  The name of the TimeEvent can be used to identify
   *  specific TimeEvents in order to take action related
   *  to that TimeEvent (e.g., a name of "Output Special"
   *  may indicate to output some special information).
   *
   *  \return The name of the TimeEvent.
   */
  virtual std::string getName() const { return name_; }
  /// Set the name of the TimeEvent.
  virtual void setName(std::string name) { name_ = name; }
  /** \brief Return the type of TimeEvent.
   *
   *  Each derived class of TimeEventBase has a type that
   *  can be used to identify the type of the TimeEvent.
   */
  virtual std::string getType() const { return timeEventType_; }
  /// Return the default time used for TimeEvents.
  virtual Scalar getDefaultTime() const { return defaultTime_; }
  /// Return the default tolerance used by TimeEvents.
  virtual Scalar getDefaultTol() const { return defaultTol_; }
  /// Return the default index used by TimeEvents.
  virtual int getDefaultIndex() const { return defaultIndex_; }
  //@}

  /// Return ParameterList with current values.
  virtual Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl =
        Teuchos::parameterList("Time Event Base");

    pl->setName(this->getName());
    pl->set("Name", this->getName());
    pl->set("Type", this->getType());

    return pl;
  }

 protected:
  virtual void setType(std::string s) { timeEventType_ = s; }

 private:
  std::string timeEventType_;  ///< Time Event type
  std::string name_;           ///< Name to identify the TimeEvent
  const Scalar defaultTime_;   ///< Default time
  const Scalar defaultTol_;    ///< Default tolerance
  const int defaultIndex_;     ///< Default index
};

}  // namespace Tempus

#endif  // Tempus_TimeEventBase_decl_hpp
