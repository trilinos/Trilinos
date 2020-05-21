// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventRange_decl_hpp
#define Tempus_TimeEventRange_decl_hpp

#include <tuple>

// Teuchos
#include "Teuchos_Time.hpp"

// Tempus
#include "Tempus_TimeEventBase.hpp"


namespace Tempus {


/** \brief TimeEventRange specifies a start, stop and stride time.
 *
 *
 */
template<class Scalar>
class TimeEventRange : virtual public TimeEventBase<Scalar>
{
public:

  /// Default constructor.
  TimeEventRange();

  /// Construct with full argument list of data members.
  TimeEventRange(std::string name, Scalar start, Scalar stop, Scalar stride,
                 Scalar relTol, bool landOnExactly);

  /// Construct with full argument list of data members.
  TimeEventRange(std::string name, Scalar start, Scalar stop, int numEvents,
                 Scalar relTol, bool landOnExactly);

  /// Destructor
  virtual ~TimeEventRange() {}

  /// \name Basic methods
  //@{
    /// Test if time is near a TimeEvent (within tolerance).
    virtual bool isTime(Scalar time) const;

    /// How much time until the next event. Negative indicating the last event is in the past.
    virtual Scalar timeToNextEvent(Scalar time) const;

    /// Time of the next event. Negative indicating the last event is in the past.
    virtual Scalar timeOfNextEvent(Scalar time) const;

    /// Test if an event occurs within the time range.
    virtual bool eventInRange(Scalar time1, Scalar time2) const;
  //@}

  /// \name Accessor methods
  //@{
    virtual void setTimeRange(Scalar start, Scalar stop, Scalar stride)
    { setTimeStart(start); setTimeStop(stop); setTimeStride(stride); }
    virtual void setTimeRange(Scalar start, Scalar stop, int numEvents)
    { setTimeStart(start); setTimeStop(stop); setNumEvents(numEvents); }

    virtual Scalar getTimeStart() const { return start_; }
    virtual void setTimeStart(Scalar start);

    virtual Scalar getTimeStop() const { return stop_; }
    virtual void setTimeStop(Scalar stop);

    virtual Scalar getTimeStride() const { return stride_; }
    virtual void setTimeStride(Scalar stride);

    virtual int getNumEvents() const { return numEvents_; }
    virtual void setNumEvents(int numEvents);

    virtual Scalar getRelTol() const { return relTol_; }
    virtual void setRelTol(Scalar relTol);

    virtual Scalar getAbsTol() const { return absTol_; }

    virtual bool getLandOnExactly() const { return landOnExactly_; }
    virtual void setLandOnExactly(bool LOE) { landOnExactly_ = LOE; }

    /// Describe member data.
    virtual void describe() const;
  //@}


protected:

  virtual void setTimeScale();

  Scalar start_;         ///< Start of time range.
  Scalar stop_;          ///< Stop of time range.
  Scalar stride_;        ///< Stride of time range.
  unsigned numEvents_;   ///< Number of events in time range.

  Scalar timeScale_;     ///< A reference time scale, max(abs(start_,stop_)).
  Scalar relTol_;        ///< Relative time tolerance for matching time events.
  Scalar absTol_;        ///< Absolute time tolerance, relTol_*timeScale_.
  bool   landOnExactly_; ///< Should these time events be landed on exactly, i.e, adjust the timestep to hit time event, versus stepping over and keeping the time step unchanged.

};


} // namespace Tempus

#endif // Tempus_TimeEventRange_decl_hpp
