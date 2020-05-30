// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeEventBase_decl_hpp
#define Tempus_TimeEventBase_decl_hpp

// Teuchos
#include "Teuchos_Time.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace Tempus {


/** \brief This class defines time events which can be used to "trigger" an action.
 *  Time events are points in time and/or timestep index where an
 *  an action should occur, such as, solution output (mesh and solution)
 *  diagnostic output, restart, screen dump, in-situ visualization,
 *  user-specified, or any other action.
 *
 *  This class will store a collection time events, so that an object
 *  may quiry it and take appropriate action.  Time events (time and
 *  timestep index) can be specified via
 *    - start, stop and stride
 *    - list of events
 *
 */
template<class Scalar>
class TimeEventBase
{
public:

  /// Constructor
  TimeEventBase()
    : name_("TimeEventBase"),
      defaultTime_ (-std::numeric_limits<Scalar>::max()*1.0e-16),
      defaultTol_  ( std::numeric_limits<Scalar>::min()),
      defaultIndex_( std::numeric_limits<int>::min())
  {}

  /// Destructor
  virtual ~TimeEventBase() {}

  /// \name Basic methods
  //@{
    /// Test if time is near a TimeEvent (within tolerance).
    virtual bool isTime(Scalar time) const
    { return false; }

    virtual Scalar getAbsTol() const
    { return defaultTol_; }

    /// How much time until the next event. Negative indicating the last event is in the past.
    virtual Scalar timeToNextEvent(Scalar time) const
    { return defaultTime_; }

    /// Time of the next event. Negative indicating the last event is in the past.
    virtual Scalar timeOfNextEvent(Scalar time) const
    { return defaultTime_; }

    /// Test if an event occurs within the time range.
    virtual bool eventInRange(Scalar time1, Scalar time2) const
    { return false; }

    /// Test if index is a time event.
    virtual bool isIndex(int index) const
    { return false; }

    /// How many indices until the next event. Negative indicating the last event is in the past.
    virtual int indexToNextEvent(int index) const
    { return defaultIndex_; }

    /// Index of the next event. Negative indicating the last event is in the past.
    virtual int indexOfNextEvent(int index) const
    { return defaultIndex_; }

    /// Test if an event occurs within the time range.
    virtual bool eventInRangeIndex(int index1, int index2) const
    { return false; }

    /// Describe member data.
    virtual void describe() const
    {
      Teuchos::RCP<Teuchos::FancyOStream> out =
        Teuchos::VerboseObjectBase::getDefaultOStream();
      *out << "TimeEventBase name = " << getName() << std::endl;
    }
  //@}

  /// \name Accessor methods
  //@{
    virtual std::string getName() const { return name_; }
    virtual void setName(std::string name) { name_ = name; }

    virtual Scalar getDefaultTime () const { return defaultTime_; }
    virtual Scalar getDefaultTol  () const { return defaultTol_; }
    virtual int    getDefaultIndex() const { return defaultIndex_; }
  //@}


private:

  std::string  name_;
  const Scalar defaultTime_;
  const Scalar defaultTol_;
  const int    defaultIndex_;

};


} // namespace Tempus

#endif // Tempus_TimeEventBase_decl_hpp
