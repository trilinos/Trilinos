/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef MLAPI_TIMEOBJECT_H
#define MLAPI_TIMEOBJECT_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "MLAPI_Error.h"
#include "MLAPI_Workspace.h"
#include "Epetra_Time.h"

namespace MLAPI {

/*!
\class TimeObject

\brief Class to track time spent in an object.

\author Marzio Sala, SNL 9214

\date Last updated on Feb-05.

*/

class TimeObject {

public:

  //! Constructor, set counter to 0.0.
  TimeObject() :
    Time_(GetEpetra_Comm())
  {
    Time_.ResetStartTime();
    TotalTime_ = 0.0;
  }

  //! Destructor.
  ~TimeObject() {};

  //! Resets the internal timer.
  inline void ResetTimer() const
  {
    Time_.ResetStartTime();
  }

  //! Updates the internal timer with the time spent since the last call to ResetTimer().
  inline void UpdateTime() const
  {
    TotalTime_ += Time_.ElapsedTime();
  }

  //! Updates the internal timer with input value \c t.
  inline void UpdateTime(double t) const
  {
    TotalTime_ += t;
  }

  //! Returns the internally stored counter.
  inline double GetTime() const
  {
    return(TotalTime_);
  }

protected:

  //! Object used to track time.
  mutable Epetra_Time Time_;
  //! Internal counter.
  mutable double TotalTime_;

}; // class TimeObject

} // namespace MLPI

#endif // MLAPI_TIMEOBJECT_H
