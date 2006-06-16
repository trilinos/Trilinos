// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_TIMEMONITOR_H
#define TEUCHOS_TIMEMONITOR_H

/*! \file Teuchos_TimeMonitor.hpp
    \brief Timer class that starts when constructed and stops when the destructor is called
*/
/** \example TimeMonitor/cxx_main.cpp
    This is an example of how to use the Teuchos::TimeMonitor class.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_PerformanceMonitorBase.hpp"
#include "Teuchos_Time.hpp"

/** \brief Defines a static non-member function that returns a time monitor.
 */ 
#define TEUCHOS_TIMER(funcName, strName) \
static Teuchos::Time& funcName()       \
{static Teuchos::RefCountPtr<Time> rtn =        \
Teuchos::TimeMonitor::getNewCounter(strName); return *rtn;}

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
static Teuchos::RefCountPtr<Teuchos::Time>  blabla_localTimer; \
if(!blabla_localTimer.get()) { \
  std::ostringstream oss; \
  oss << FUNCNAME; \
  blabla_localTimer = Teuchos::TimeMonitor::getNewCounter(oss.str()); \
} \
Teuchos::TimeMonitor blabla_localTimeMonitor(*blabla_localTimer)


namespace Teuchos
{
  /** 
   * \brief A timer class that starts when constructed and stops when the 
   * destructor is called. 
   *
   * Termination upon destruction lets this timer behave
   * correctly even if scope is exited because of an exception. 
   *
   * \note Teuchos::TimeMonitor uses the Teuchos::Time class internally.
   */
  class TimeMonitor : public PerformanceMonitorBase<Time>
    {
    public:
 
      /** \brief Constructor starts timer */
      TimeMonitor(Time& timer, bool reset=false)
        : PerformanceMonitorBase<Time>(timer, reset)
        {
          if (!isRecursiveCall()) counter().start(reset);
        }

      /** \brief Destructor causes timer to stop */
      inline ~TimeMonitor()
        {
          if (!isRecursiveCall()) counter().stop();
        }

      /** Wrapping of getNewCounter() for backwards compatibiity with old code*/
      static Teuchos::RefCountPtr<Time> getNewTimer(const string& name)
      {return getNewCounter(name);}

      /** \brief Print summary statistics for a group of timers. 
       * Timings are gathered from all processors 
       *
       * \note This method <b>must</b> be called by all processors */
      static void summarize(ostream &out=std::cout, 
                            bool alwaysWriteLocal=false,
                            bool writeGlobalStats=true);

    private:
    };

}
#endif
