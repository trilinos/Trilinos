#ifndef _PETRA_TIME_H_
#define _PETRA_TIME_H_

//! Petra_Time:  The Petra Timing Class.
/*! The Petra_Time class is a wrapper that encapsulates the general
  information needed getting timing information.  Currently it return
  the elapsed time for each calling processor..
  A Petra_Comm object is required for building all Petra_Time objects.
  
  Petra_Time support both serial execution and (via MPI) parallel 
  distributed memory execution.  It is meant to insulate the user from
  the specifics of timing across a variety of platforms.
*/

#include "Petra_Petra.h"
#include "Petra_Comm.h"
#ifdef PETRA_MPI
#include "mpi.h"
#elif ICL
#include <time.h>
#else
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif

class Petra_Time {
    
  public:
  //! Petra_Time Constructor.
  /*! Creates a Petra_Time instance. This instance can be queried for
      elapsed time on the calling processor.  StartTime is also set
      for use with the ElapsedTime function.
  */
  Petra_Time(const Petra_Comm & Comm);

  //! Petra_Time Copy Constructor.
  /*! Makes an exact copy of an existing Petra_Time instance.
  */
  Petra_Time(const Petra_Time& Time);

  //! Petra_Time wall-clock time function.
  /*! Returns the wall-clock time in seconds.  A code section can be 
      timed by putting it between two calls to WallTime and taking the
      difference of the times.
  */
  double WallTime(void) const;

  //! Petra_Time function to reset the start time for a timer object.
  /*! Resets the start time for the timer object to the current time
      A code section can be 
      timed by putting it between a call to ResetStartTime and ElapsedTime.
  */
  void ResetStartTime(void);

  //! Petra_Time elapsed time function.
  /*! Returns the elapsed time in seconds since the timer object was
      constructed, or since the ResetStartTime function was called. 
      A code section can be 
      timed by putting it between the Petra_Time constructor and a call to 
      ElapsedTime, or between a call to ResetStartTime and ElapsedTime.
  */
  double ElapsedTime(void) const;

  //! Petra_Time Destructor.
  /*! Completely deletes a Petra_Time object.  
  */
  virtual ~Petra_Time(void);

 private:

  double StartTime_;
  const Petra_Comm * Comm_;
  
};

#endif /* _PETRA_TIME_H_ */
