// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_TIME_HPP_
#define _TEUCHOS_TIME_HPP_

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos 
{

  /** Teuchos::Time class is a wall-clock timer. For exception safety and correct
   * behavior in reentrant code, this class should
   * generally be used only through the Teuchos::TimeMonitor mechanism. 
   *
   * To time a section of code, place it in between calls to start() 
   * and stop(). 
   *
   * Initial version by Mike Heroux and Kris Campshaw. 
   * Modified as follows by Kevin Long, 9/29/03:
   * <ul>
   * <li> There is no need to define explicit copy ctor and dtor for this class.
   * <li> The wallTime() method returns the same value for every instance of this class, so
   * it would be best to make it a static method.
   * <li> Removed the communicator data member. Cross-processor comparisons of timings
   * can be done by the TimeMonitor.
   * </ul>
   */ 


  class Time
  {

  public:
    /** Construct with a descriptive name */
    Time(const string& name, bool start = false);
  
    /** returns current wall-clock time in seconds.*/
    static double wallTime();
  
    /** starts the timer */
    void start();

    /** stop the timer */
    double stop();

    /** returns the total time accumulated by this timer. Note that this should be called
     * only when the clock is stopped. */
    double totalElapsedTime() const {return totalTime_;}

    /** indicates if this timer is currently running, i.e., if it has been started but
     * not yet stopped. It is necessary to know if a timer is running to avoid 
     * incorrectly starting or stopping in reentrant code. */
    bool isRunning() const {return isRunning_;}

    /** return the name of this timer */
    const string& name() const {return name_;}
    
  private:
    double startTime_;

    double totalTime_;

    bool isRunning_;

    string name_;
  };

} // namespace Teuchos

// #include "Teuchos_Time.cpp"

#endif // TEUCHOS_TIME_HPP_
