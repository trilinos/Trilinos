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


/** \brief Basic wall-clock timer class. 
 *
 *  To time a section of code, place it in between calls to start() and stop(). 
 *
 *  \note For std::exception safety and correct behavior in reentrant code, this class should
 * generally be used only through the Teuchos::TimeMonitor mechanism. 
 *
 */

/* ADDITIONAL COMMENTS:
 * Initial version by Mike Heroux and Kris Kampshoff. 
 * Modified as follows by Kevin Long, 9/29/03:
 * <ul>
 * <li> There is no need to define explicit copy ctor and dtor for this class.
 * <li> The wallTime() method returns the same value for every instance of this class, so
 * it would be best to make it a static method.
 * <li> Removed the communicator data member. Cross-processor comparisons of timings
 * can be done by the TimeMonitor.
 * </ul>
 */ 

class TEUCHOS_LIB_DLL_EXPORT Time
{

public:
  /** \brief Construct with a descriptive name */
  Time(const std::string& name, bool start = false);
  
  /** \brief Returns current wall-clock time in seconds.*/
  static double wallTime();
  
  /** \brief Starts the timer */
  void start(bool reset = false);

  /** \brief Stops the timer */
  double stop();

  /** \brief Returns the total time (in seconds) accumulated by this timer.
   *
   * <b>This should be called only when the clock is stopped.</b> */
  double totalElapsedTime(bool readCurrentTime = false) const;

  /** \brief Resets the cummulative time and number of times this timer
   * has been called.  Does not affect any other state. */
  void reset() {totalTime_ = 0; numCalls_ = 0;}

  /** \brief Indicates if this timer is currently running, i.e., if it has been started but
   * not yet stopped. 
   *
   *	It is necessary to know if a timer is running to avoid incorrectly starting or 
   *  stopping in reentrant code. */
  bool isRunning() const {return isRunning_;}

  /** \brief Return the name of this timer */
  const std::string& name() const {return name_;}

  /** \brief Increment the number of times this timer has been called */
  void incrementNumCalls() {numCalls_++;}

  /** \brief Return the number of times this timer has been called */
  int numCalls() const {return numCalls_;}

private:

  double startTime_;
  double totalTime_;
  bool isRunning_;
  std::string name_;
  int numCalls_;

};


} // namespace Teuchos


#endif // TEUCHOS_TIME_HPP_
