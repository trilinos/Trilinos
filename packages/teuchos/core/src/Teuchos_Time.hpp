// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_TIME_HPP_
#define _TEUCHOS_TIME_HPP_

/*! \file Teuchos_Time.hpp
    \brief Basic wall-clock timer class
*/

#include "Teuchos_ConfigDefs.hpp"

#include <ctime>
#ifdef HAVE_MPI
#include "mpi.h"
#else
#if ICL || defined(_WIN32)
#include <time.h>
#else
#include <sys/time.h>
#ifndef MINGW
#include <sys/resource.h>
#endif
#endif
#endif


namespace Teuchos {


/// \class Time
/// \brief Wall-clock timer.
///
/// To time a section of code, place it in between calls to start()
/// and stop().  It is better to access this class through the
/// TimeMonitor class (which see) for exception safety and correct
/// behavior in reentrant code.
class TEUCHOSCORE_LIB_DLL_EXPORT Time {
public:
  /// \brief Constructor
  ///
  /// \param name [in] Name of the timer.
  /// \param start [in] If true, start the timer upon creating it.  By
  ///   default, the timer only starts running when you call start().
  Time (const std::string& name, bool start = false);

  /// \brief Current wall-clock time in seconds.
  ///
  /// This is only useful for measuring time intervals.  The absolute
  /// value returned is measured relative to some arbitrary time in
  /// the past.
  static double wallTime ();

  /// \brief Start the timer, if the timer is enabled (see disable()).
  ///
  /// \param reset [in] If true, reset the timer's total elapsed time
  ///   to zero before starting the timer.  By default, the timer
  ///   accumulates the total elapsed time for all start() ... stop()
  ///   sequences.
  void start (bool reset = false);

  //! Stop the timer, if the timer is enabled (see disable()).
  double stop ();

  //! "Disable" this timer, so that it ignores calls to start() and stop().
  void disable ();

  //! "Enable" this timer, so that it (again) respects calls to start() and stop().
  void enable ();

  //! Whether the timer is enabled (see disable()).
  bool isEnabled () const {
    return enabled_;
  }

  /// \brief The total time in seconds accumulated by this timer.
  ///
  /// \param readCurrentTime [in] If true, return the current elapsed
  ///   time since the first call to start() when the timer was
  ///   enabled, whether or not the timer is running or enabled.  If
  ///   false, return the total elapsed time as of the last call to
  ///   stop() when the timer was enabled.
  ///
  /// \note If start() has never been called when the timer was
  ///   enabled, and if readCurrentTime is true, this method will
  ///   return wallTime(), regardless of the actual start time.
  double totalElapsedTime (bool readCurrentTime = false) const;

  //! Reset the cummulative time and call count.
  void reset ();

  /// \brief Whether the timer is currently running.
  ///
  /// "Currently running" means either that start() has been called
  /// without an intervening stop() call, or that the timer was
  /// created already running and stop() has not since been called.
  bool isRunning() const {
    return isRunning_;
  }

  //! The name of this timer.
  const std::string& name() const {
    return name_;
  }

  /// \brief Increment the number of times this timer has been called,
  ///   if the timer is enabled (see disable()).
  void incrementNumCalls();

  //! The number of times this timer has been called while enabled.
  int numCalls() const {return numCalls_;}

private:
  double startTime_;
  double totalTime_;
  bool isRunning_;
  bool enabled_;
  std::string name_;
  int numCalls_;
};


} // namespace Teuchos


#endif // TEUCHOS_TIME_HPP_
