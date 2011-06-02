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

#ifndef TEUCHOS_TIMEMONITOR_HPP
#define TEUCHOS_TIMEMONITOR_HPP


/*! \file Teuchos_TimeMonitor.hpp
 *
 * \brief Timer class that starts when constructed and stops when the
 * destructor is called
 */

/** \example TimeMonitor/cxx_main.cpp
 *
 * This is an example of how to use the Teuchos::TimeMonitor class.
 */


#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_PerformanceMonitorBase.hpp"
#include "Teuchos_Time.hpp"


/** \brief Defines a static non-member function that returns a time monitor.
 */ 
#define TEUCHOS_TIMER(funcName, strName) \
  static Teuchos::Time& funcName() \
  {static Teuchos::RCP<Time> rtn = \
      Teuchos::TimeMonitor::getNewCounter(strName); return *rtn;}


/** \brief Defines a timer for a specific function (with differentiator).
 *
 * Same as TEUCHOS_FUNC_TIME_MONITOR(...) except required when used more than
 * once in the same function (like a block of code).
 */
#define TEUCHOS_FUNC_TIME_MONITOR_DIFF( FUNCNAME, DIFF ) \
  static Teuchos::RCP<Teuchos::Time> DIFF ## blabla_localTimer; \
  if(!DIFF ## blabla_localTimer.get()) { \
    std::ostringstream oss; \
    oss << FUNCNAME; \
    DIFF ## blabla_localTimer = Teuchos::TimeMonitor::getNewCounter(oss.str()); \
  } \
  Teuchos::TimeMonitor DIFF ## blabla_localTimeMonitor(*DIFF ## blabla_localTimer)


/** \brief Defines a timer for a specific function.
 *
 * Note that the name of the timer can be formated with stream inserts.
 * For example, we can define a time monitor for a function as follows:
 
 \code

 template<typename Scalar>
 void foo()
 {
 TEUCHOS_FUNC_TIME_MONITOR(
 "foo<"<<Teuchos::ScalarTraits<Scalar>::name()<<">()"
 );
 ...
 }

 \endcode

 * The timer can then be printed at the end of the program using

 \code

 Teuchos::TimeMonitor::summarize(std::cout);

 \endcode
 
*/
#define TEUCHOS_FUNC_TIME_MONITOR( FUNCNAME ) \
  TEUCHOS_FUNC_TIME_MONITOR_DIFF( FUNCNAME, main )


namespace Teuchos {


/** \brief A timer class that starts when constructed and stops when the
 * destructor is called.
 *
 * Termination upon destruction lets this timer behave
 * correctly even if scope is exited because of an std::exception. 
 *
 * NOTE: It is critical that this class only be used to time functions that
 * are called only within the main program and not at pre-program setup or
 * post-program teardown!
 *
 * \note Teuchos::TimeMonitor uses the Teuchos::Time class internally.
 */
class TEUCHOS_LIB_DLL_EXPORT TimeMonitor : public PerformanceMonitorBase<Time>
{
public:

  /** \name Constructor/Destructor */
  //@{
 
  /** \brief Constructor starts timer */
  TimeMonitor(Time& timer, bool reset=false)
    : PerformanceMonitorBase<Time>(timer, reset)
    {
      if (!isRecursiveCall()) counter().start(reset);
    }
 
  /** \brief Destructor causes timer to stop */
  ~TimeMonitor()
    {
      if (!isRecursiveCall()) counter().stop();
    }

  //@}

  /** \name Static functions */
  //@{
 
  /** \brief Wrapping of getNewCounter() for backwards compatibiity with old
   * code.
   */
  static Teuchos::RCP<Time> getNewTimer(const std::string& name)
    {return getNewCounter(name);}

  /** \brief Reset the global timers to zero.
   *
   * <b>Preconditions:</b><ul>
   * <li>None of the timers must currently be running!
   * </ul>
   */
  static void zeroOutTimers();
 
  /** \brief Print summary statistics for a group of timers. 
   *
   * Timings are gathered from all processors 
   *
   * \note This method <b>must</b> be called by all processors */
  static void summarize(
    std::ostream &out=std::cout, 
    const bool alwaysWriteLocal=false,
    const bool writeGlobalStats=true,
    const bool writeZeroTimers=true
    );

  //@}

};


} // namespace Teuchos


#endif // TEUCHOS_TIMEMONITOR_H
