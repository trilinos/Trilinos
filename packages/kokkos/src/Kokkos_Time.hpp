
//@HEADER
// ************************************************************************
// 
//          Kokkos: A Fast Kernel Package
//              Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef KOKKOS_TIME_H
#define KOKKOS_TIME_H

//! Kokkos_Time:  The Kokkos Timing Class.
/*! The Kokkos_Time class is a wrapper that encapsulates the general
  information needed getting timing information.  Currently it return
  the elapsed time for each calling processor..
  A Kokkos_Comm object is required for building all Kokkos_Time objects.
  
  Kokkos_Time support both serial execution and (via MPI) parallel 
  distributed memory execution.  It is meant to insulate the user from
  the specifics of timing across a variety of platforms.
*/

#include "Kokkos_Object.hpp"
#include "Kokkos_Comm.hpp"

#ifdef KOKKOS_MPI
#include "mpi.hpp"
#elif ICL
#include <time.hpp>
#else
#include <sys/time.hpp>
#ifndef MINGW
#include <sys/resource.hpp>
#endif
#endif

class Kokkos_Time: public Kokkos_Object {
    
  public:
  //! Kokkos_Time Constructor.
  /*! Creates a Kokkos_Time instance. This instance can be queried for
      elapsed time on the calling processor.  StartTime is also set
      for use with the ElapsedTime function.
  */
  Kokkos_Time(const Kokkos_Comm & Comm);

  //! Kokkos_Time Copy Constructor.
  /*! Makes an exact copy of an existing Kokkos_Time instance.
  */
  Kokkos_Time(const Kokkos_Time& Time);

  //! Kokkos_Time wall-clock time function.
  /*! Returns the wall-clock time in seconds.  A code section can be 
      timed by putting it between two calls to WallTime and taking the
      difference of the times.
  */
  double WallTime(void) const;

  //! Kokkos_Time function to reset the start time for a timer object.
  /*! Resets the start time for the timer object to the current time
      A code section can be 
      timed by putting it between a call to ResetStartTime and ElapsedTime.
  */
  void ResetStartTime(void);

  //! Kokkos_Time elapsed time function.
  /*! Returns the elapsed time in seconds since the timer object was
      constructed, or since the ResetStartTime function was called. 
      A code section can be 
      timed by putting it between the Kokkos_Time constructor and a call to 
      ElapsedTime, or between a call to ResetStartTime and ElapsedTime.
  */
  double ElapsedTime(void) const;

  //! Kokkos_Time Destructor.
  /*! Completely deletes a Kokkos_Time object.  
  */
  virtual ~Kokkos_Time(void);

 private:

  double StartTime_;
  const Kokkos_Comm * Comm_;
  
};

#endif /* KOKKOS_TIME_H */
