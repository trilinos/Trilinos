//@HEADER
// ************************************************************************
// 
//                 Kokkos: A Fast Kernel Package
//              Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_TIME_H
#define KOKKOS_TIME_H

#include "Kokkos_ConfigDefs.hpp"

namespace Kokkos {

//! Kokkos_Time:  The Kokkos Timing Class.
/*! The Kokkos_Time class is a wrapper that encapsulates the general
  information needed getting timing information.  Currently it return
  the elapsed time for each calling processor..
  A Kokkos_Comm object is required for building all Kokkos_Time objects.
  
  Kokkos_Time support both serial execution and (via MPI) parallel 
  distributed memory execution.  It is meant to insulate the user from
  the specifics of timing across a variety of platforms.
*/
  class Time {
    
  public:
    //! Time Constructor.
    /*! Creates a Time instance. This instance can be queried for
      elapsed time on the calling processor.  StartTime is also set
      for use with the ElapsedTime function.
    */
    Time() {startTime_ = wallTime();};

    //! Time Copy Constructor.
    /*! Makes an exact copy of an existing Time instance.
     */
    Time(const Time& time)
      : startTime_(time.startTime_) {};

    //! Time wall-clock time function.
    /*! Returns the wall-clock time in seconds.  A code section can be 
      timed by putting it between two calls to wallTime and taking the
      difference of the times.
    */
    double wallTime() const {
#ifdef ICL

      clock_t start;
      //double duration;

      start = clock();
      return (double)( start ) / CLOCKS_PER_SEC;

#else

#ifndef MINGW
      struct timeval tp;
      static long start=0, startu;
      if (!start) {
	gettimeofday(&tp, NULL);
	start = tp.tv_sec;
	startu = tp.tv_usec;
	return(0.0);
      }
      gettimeofday(&tp, NULL);
      return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );
#else
      return( (double) clock() / CLOCKS_PER_SEC );
#endif

#endif

    };

    //! Kokkos_Time function to reset the start time for a timer object.
    /*! Resets the start time for the timer object to the current time
      A code section can be 
      timed by putting it between a call to ResetStartTime and ElapsedTime.
    */
     void resetStartTime(){
      startTime_ = wallTime();
      return;
    };

    //! Kokkos_Time elapsed time function.
    /*! Returns the elapsed time in seconds since the timer object was
      constructed, or since the ResetStartTime function was called. 
      A code section can be 
      timed by putting it between the Kokkos_Time constructor and a call to 
      ElapsedTime, or between a call to ResetStartTime and ElapsedTime.
    */
    double elapsedTime() const {
      return(wallTime()-startTime_);
    };

    //! Kokkos_Time Destructor.
    /*! Completely deletes a Kokkos_Time object.  
     */
    virtual ~Time(){};

  private:

    double startTime_;  
  };
} // namespace Kokkos
#endif /* KOKKOS_TIME_H */
