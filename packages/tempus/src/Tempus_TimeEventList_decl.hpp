//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventList_decl_hpp
#define Tempus_TimeEventList_decl_hpp

#include <vector>

#include "Teuchos_Time.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_config.hpp"
#include "Tempus_TimeEventBase.hpp"

namespace Tempus {

/** \brief TimeEventList specifies a list of time events.
 *
 *
 */
template <class Scalar>
class TimeEventList : virtual public TimeEventBase<Scalar> {
 public:
  /// Default constructor.
  TimeEventList();

  /// Construct with full argument list of data members.
  TimeEventList(std::vector<Scalar> timeList,
                std::string name = "TimeEventList", bool landOnExactly = true,
                Scalar relTol = std::numeric_limits<Scalar>::epsilon() *
                                Scalar(100.0));

  /// Destructor
  virtual ~TimeEventList() {}

  /// \name Basic methods
  //@{
  /** \brief Test if time is near an event (within tolerance).
   *
   *  Test if any of the events in the list is near the input time.
   *
   *  \param[in] time The input time to check if it is near an event.
   *
   *  \return Return true if time is near a TimeEvent
   *          ( timeEvent-absTol < time < timeEvent+absTol ), otherwise
   *          return false.
   */
  virtual bool isTime(Scalar time) const;

  /** \brief How much time until the next event.
   *
   *  Determine the amount of time until the next timeEvent in the
   *  list (i.e., time of next event minus the input time).
   *  If the input time is after all events, the default time
   *  (a time in the distant future) minus the input time is returned.
   *
   *  \param time      [in] The input time.
   *  \return The time to the next event.
   */
  virtual Scalar timeToNextEvent(Scalar time) const;

  /** \brief Return the time of the next event following the input time.
   *
   *  Returns the time of the next event that follows the input time.
   *  If the input time is before all events, the time of the first
   *  event is returned.  If the input time is after all events, the
   *  default time (a time in the distant future) is returned.  If the
   *  input time is an event time, the time of the next event is returned.
   *
   *  \param time [in] Input time.
   *  \return Time of the next event.
   */
  virtual Scalar timeOfNextEvent(Scalar time) const;

  /** \brief Test if an event occurs within the time range.
   *
   *  Find if an event is within the input list,
   *  (time1 < timeEvent-absTol and timeEvent-absTol <= time2),
   *  including the event's absolute tolerance.  Note, this
   *  does not include time1, but does include time2.
   *
   *  \param time1 [in] Input time of one end of the range.
   *  \param time2 [in] Input time of the other end of the range.
   *
   *  \return True if a time event is within the range.
   */
  virtual bool eventInRange(Scalar time1, Scalar time2) const;

  /// Describe member data.
  virtual void describe(Teuchos::FancyOStream &out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

  /// \name Accessor methods
  //@{
  /// Return the list of time events.
  virtual std::vector<Scalar> getTimeList() const { return timeList_; }

  /** \brief Set the list of time events.
   *
   *  This will completely replace the vector of time events.
   *  If sort=true, the vector will be sorted and duplicate
   *  entries removed (i.e., only has unique entries).
   *
   *  \param timeList [in] Vector of time event.
   *  \param sort     [in] Sort vector into ascending order, if true.
   */
  virtual void setTimeList(std::vector<Scalar> timeList, bool sort = true);

  /** \brief Add the time to event vector.
   *
   *  The input time will be inserted into the vector of
   *  events in ascending order.  If the time is already
   *  present (within tolerance), it is not added to keep
   *  the vector unique.
   *
   *  \param time [in] Time to insert to vector of events.
   */
  virtual void addTime(Scalar time);

  /// Clear the vector of all events.
  virtual void clearTimeList() { timeList_.clear(); }

  /// Return the relative tolerance.
  virtual Scalar getRelTol() const { return relTol_; }

  /** \brief Set the relative tolerance.
   *
   *  The relative tolerance is used to set the absolute
   *  tolerance along with the TimeEvent time scale i.e.,
   *  absTol_ = timeScale_ * relTol_.
   *  Also see getAbsTol() and setTimeScale().
   *
   *  \param relTol [in] The input relative tolerance.
   */
  virtual void setRelTol(Scalar relTol);

  /** \brief Return the absolute tolerance.
   *
   *  The absolute tolerance is primarily used to determine
   *  if two times are equal (i.e., t1 is equal to t2, if
   *  t2-absTol_ < t1 < t2+absTol_).
   *
   *  For TimeEventList,
   *
   *  \return The absolute tolerance.
   */
  virtual Scalar getAbsTol() const { return absTol_; }

  /** \brief Return if the time events need to be landed on exactly.
   *
   *  If true, this sets whether the time events need to be landed
   *  on exactly, e.g., the time step needs to be adjusted so the
   *  solution is determined at the time event.
   *
   *  If false, this indicates that time event will still occur but
   *  can be stepped over without changing the time step.
   *
   *  \return LOE Flag indicating if TimeEvent should land on the time event
   * exactly.
   */
  virtual bool getLandOnExactly() const { return landOnExactly_; }

  /// Set if the time event should be landed on exactly.
  virtual void setLandOnExactly(bool LOE) { landOnExactly_ = LOE; }
  //@}

  /** \brief Return a valid ParameterList with current settings.
   *
   *  The returned ParameterList will contain the current parameters
   *  and can be used to reconstruct the TimeEventList using
   *  createTimeEventList(...).  The ParameterList will have
   *  the TimeEventList parameters along with all the parameters
   *  for the TimeEvents contained in the composite.
   *
   * \return Teuchos::ParameterList of TimeEventList.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

 protected:
  /** \brief Set the time scale for the time events
   *
   *  This sets the reference time scale, which is solely
   *  determined from the time events in timeList_ (this is
   *  why it is a protected function).  It tries to find
   *  an appropriate scale for the time events, e.g.,
   *  max(t0, ..., tn).
   */
  virtual void setTimeScale();

  std::vector<Scalar> timeList_;  ///< Sorted and unique list of time events.

  Scalar timeScale_;    ///< A reference time scale.
  Scalar relTol_;       ///< Relative time tolerance for matching time events.
  Scalar absTol_;       ///< Absolute time tolerance, relTol_*timeScale_.
  bool landOnExactly_;  ///< Should these time events be landed on exactly, i.e,
                        ///< adjust the timestep to hit time event, versus
                        ///< stepping over and keeping the time step unchanged.
};

// Nonmember Contructors
// ------------------------------------------------------------------------

/** \brief Nonmember Constructor via ParameterList.
 *
 *  If the input ParameterList is Teuchos::null, return a default
 *  TimeEventList.  A valid ParameterList can be obtained
 *  from getValidParameters().
 *
 *  \param pList [in] The input ParameterList to construct from.
 *  \return Constructed TimeEventList.
 */
template <class Scalar>
Teuchos::RCP<TimeEventList<Scalar> > createTimeEventList(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_TimeEventList_decl_hpp
