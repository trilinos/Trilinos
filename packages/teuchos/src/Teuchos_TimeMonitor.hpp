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

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos
{

  using std::string;

  /** 
   * \brief A timer class that starts when constructed and stops when the 
   * destructor is called. 
   *
   * Termination upon destruction lets this timer behave
   * correctly even if scope is exited because of an exception. 
   *
   * \note Teuchos::TimeMonitor uses the Teuchos::Time class internally.
   */
  class TimeMonitor
    {
    public:
 
      /** \brief Constructor starts timer */
      TimeMonitor(Time& timer)
        : timer_(timer), isRoot_(!timer.isRunning())
        {
          if (isRoot_) timer_.start();
        }

      /** \brief Destructor causes timer to stop */
      inline ~TimeMonitor()
        {
          if (isRoot_) timer_.stop();
        }

      /** \brief Print summary statistics for a group of timers. Timings are gathered
       * from all processors */
      static void summarize();

      /** \brief Create a new timer with the given name, and append it to the list of
       * timers to be used */
      static RefCountPtr<Time> getNewTimer(const string& name);
    private:
      
      Time& timer_;
      bool isRoot_;

      /** collect summary timings from all processors */
      static void gatherTimings(const Array<double>& timings,
                                Array<double>& minTime,
                                Array<double>& avgTime,
                                Array<double>& maxTime);
      
      static Array<RefCountPtr<Time> > timers_;
    };

}
#endif
