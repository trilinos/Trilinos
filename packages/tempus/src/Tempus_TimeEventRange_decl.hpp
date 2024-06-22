//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_TimeEventRange_decl_hpp
#define Tempus_TimeEventRange_decl_hpp

#include <tuple>

#include "Teuchos_Time.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tempus_config.hpp"
#include "Tempus_TimeEventBase.hpp"

namespace Tempus {

/** \brief TimeEventRange specifies a start, stop and stride time.
 *
 *
 */
template <class Scalar>
class TimeEventRange : virtual public TimeEventBase<Scalar> {
 public:
  /// Default constructor.
  TimeEventRange();

  /// Construct from start, stop and stride.
  TimeEventRange(Scalar start, Scalar stop, Scalar stride,
                 std::string name = "", bool landOnExactly = true,
                 Scalar relTol = std::numeric_limits<Scalar>::epsilon() *
                                 Scalar(100.0));

  /// Construct from start, stop and number of events.
  TimeEventRange(Scalar start, Scalar stop, int numEvents,
                 std::string name = "", bool landOnExactly = true,
                 Scalar relTol = std::numeric_limits<Scalar>::epsilon() *
                                 Scalar(100.0));

  /// Destructor
  virtual ~TimeEventRange() {}

  /// \name Basic methods
  //@{
  /** \brief Test if time is near an event (within tolerance).
   *
   *  Test if any of the events in the range is near the input time.
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
   *  range (i.e., time of next event minus the input time).
   *  If the input time is after all events, the default time
   *  (a time in the distant future) minus the input time is returned.
   *
   *  \param[in] time The input time to measure to the next TimeEvent.
   *
   *  \return The amount of time to next TimeEvent ( nextTimeEvent-time ).
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
   *  Find if an event is within the input range,
   *  (time1 < timeEvent-absTol and timeEvent-absTol <= time2),
   *  including the event's absolute tolerance.  Note, this
   *  does not include time1, but does include time2.
   *
   *  \param time1 [in] Input time of one end of the range.
   *  \param time2 [in] Input time of the other end of the range.
   *  \return True if an time event is within the range.
   */
  virtual bool eventInRange(Scalar time1, Scalar time2) const;
  //@}

  /// \name Accessor methods
  //@{
  /** \brief Set the range of time events from start, stop and stride.
   *
   *  This will completely replace the time range.
   *
   *  \param start  [in] The start of the time range.
   *  \param stop   [in] The stop of the time range.
   *  \param stride [in] The stride of the time range.
   */
  virtual void setTimeRange(Scalar start, Scalar stop, Scalar stride);

  /** \brief Set the range of time events from start, stop and number of events.
   *
   *  This will completely replace the time range
   *
   *  \param start     [in] The start of the time range.
   *  \param stop      [in] The stop of the time range.
   *  \param numEvents [in] The number of events in time range.
   */
  virtual void setTimeRange(Scalar start, Scalar stop, int numEvents);

  /// Return the start of the time range.
  virtual Scalar getTimeStart() const { return start_; }
  /// Set the start of the time range.
  virtual void setTimeStart(Scalar start);

  /// Return the stop of the time range.
  virtual Scalar getTimeStop() const { return stop_; }
  /// Set the stop of the time range.
  virtual void setTimeStop(Scalar stop);

  /// Return the stride of the time range.
  virtual Scalar getTimeStride() const { return stride_; }

  /** \brief Set the stride of the time range.
   *
   *  If start_ = stop_, then stride_ = 0.0, set numEvents_ = 1, and return.
   *  If stride_ > stop_-start_ or stride_ < 2*absTol_, then stride =
   * stop_-start_.
   *
   *  Reset numEvents_ = (stop_-start_)/stride_ + 1.
   *
   *  \param stride [in] The time stride for the time range.
   */
  virtual void setTimeStride(Scalar stride);

  /// Return the number of time events in the time range.
  virtual int getNumEvents() const { return numEvents_; }

  /** \brief Set the number of time events.
   *
   *  - If numEvents_ < 0, then reset numEvents from start_, stop_ and stride_.
   *  - ElseIf start_ = stop_, numEvents = 1, and reset stride = 0.
   *  - ElseIf numEvents_ < 2, numEvents = 2, and stride = stop_ - start_.
   *  - Else stride = (stop_ - start_)/(numEvents_-1).
   *
   *  If the resulting stride is less than twice the absolute tolerance,
   *  the stride is set to 2*absTol_.
   *
   *  \param numEvents [in] The number of events in time range.
   */
  virtual void setNumEvents(int numEvents);

  /// Return the relative tolerance.
  virtual Scalar getRelTol() const { return relTol_; }

  /** \brief Set the relative tolerance.
   *
   *  The relative tolerance is used to set the absolute
   *  tolerance along with the TimeEvent time scale, i.e.,
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
   *  \return The absolute tolerance.
   */
  virtual Scalar getAbsTol() const { return absTol_; }

  /** \brief Set if the time events need to be landed on exactly.
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

  /// Describe member data.
  virtual void describe(Teuchos::FancyOStream &out,
                        const Teuchos::EVerbosityLevel verbLevel) const;

  /** \brief Return a valid ParameterList with current settings.
   *
   *  The returned ParameterList will contain the current parameters
   *  and can be used to reconstruct the TimeEventRange using
   *  createTimeEventRange(...).  The ParameterList will have
   *  the TimeEventRange parameters along with all the parameters
   *  for the TimeEvents contained in the composite.
   *
   * \return Teuchos::ParameterList of TimeEventRange.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

 protected:
  /** \brief Set the time scale for the time events
   *
   *  This sets the reference time scale, which is solely
   *  determined from the time events in the range (this is
   *  why it is a protected function).  It tries to find
   *  an appropriate scale for the time events, e.g.,
   *  max(t0, ... , tn).
   */
  virtual void setTimeScale();

  Scalar start_;        ///< Start of time range.
  Scalar stop_;         ///< Stop of time range.
  Scalar stride_;       ///< Stride of time range.
  unsigned numEvents_;  ///< Number of events in time range.

  Scalar timeScale_;    ///< A reference time scale, max(abs(start_,stop_)).
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
 *  TimeEventRange.  A valid ParameterList can be obtained
 *  from getValidParameters().
 *
 *  \param pList [in] The input ParameterList to construct from.
 *  \return Constructed TimeEventRange.
 */
template <class Scalar>
Teuchos::RCP<TimeEventRange<Scalar> > createTimeEventRange(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_TimeEventRange_decl_hpp
