// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventList_decl_hpp
#define Tempus_TimeEventList_decl_hpp

#include <vector>

// Teuchos
#include "Teuchos_Time.hpp"

// Tempus
#include "Tempus_TimeEventBase.hpp"


namespace Tempus {


/** \brief TimeEventList specifies a list of time events.
 *
 *
 */
template<class Scalar>
class TimeEventList : virtual public TimeEventBase<Scalar>
{
public:

  /// Default constructor.
  TimeEventList();

  /// Construct with full argument list of data members.
  TimeEventList(std::string name, std::vector<Scalar> timeList,
                 Scalar relTol, bool landOnExactly);

  /// Destructor
  virtual ~TimeEventList() {}

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

    /// Describe member data.
    virtual void describe() const;
  //@}

  /// \name Accessor methods
  //@{
    virtual std::vector<Scalar> getTimeList() const { return timeList_; }
    virtual void setTimeList(std::vector<Scalar> timeList);
    virtual void addTime(Scalar time);
    virtual void clearTimeList() { timeList_.clear(); }

    virtual Scalar getRelTol() const { return relTol_; }
    virtual void setRelTol(Scalar relTol);

    virtual Scalar getAbsTol() const { return absTol_; }

    virtual bool getLandOnExactly() const { return landOnExactly_; }
    virtual void setLandOnExactly(bool LOE) { landOnExactly_ = LOE; }
  //@}


protected:

  virtual void setTimeScale();

  std::vector<Scalar> timeList_; ///< Sorted and unique list of time events.

  Scalar timeScale_;     ///< A reference time scale, max(abs(start_,stop_)).
  Scalar relTol_;        ///< Relative time tolerance for matching time events.
  Scalar absTol_;        ///< Absolute time tolerance, relTol_*timeScale_.
  bool   landOnExactly_; ///< Should these time events be landed on exactly, i.e, adjust the timestep to hit time event, versus stepping over and keeping the time step unchanged.
};


} // namespace Tempus

#endif // Tempus_TimeEventList_decl_hpp
