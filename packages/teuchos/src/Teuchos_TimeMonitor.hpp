#ifndef TEUCHOS_TIMEMONITOR_H
#define TEUCHOS_TIMEMONITOR_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Array.hpp"

namespace Teuchos
{

  using std::string;

  /** 
   * Teuchos::TimeMonitor is a timer that starts when constructed and stops when the 
   * destructor is called. Termination upon destruction lets this timer behave
   * correctly even if scope is exited because of an exception. 
   */
  class TimeMonitor
    {
    public:
      /** Constructor starts timer */
      TimeMonitor(Time& timer)
        : timer_(timer), isRoot_(!timer.isRunning())
        {
          if (isRoot_) timer_.start();
        }

      /** Destructor causes timer to stop */
      inline ~TimeMonitor()
        {
          if (isRoot_) timer_.stop();
        }

      /** Print summary statistics for a group of timers. Timings are gathered
       * from all processors */
      static void summarize();

      /** Create a new timer with the given name, and append it to the list of
       * timers to be  */
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
