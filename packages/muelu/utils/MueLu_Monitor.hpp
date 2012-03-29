#ifndef MUELU_MONITOR_HPP
#define MUELU_MONITOR_HPP

#include <string>
#include <stack>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_VerboseObject.hpp"

#include "MueLu_Level.hpp"

namespace MueLu {

  // Helper function. Similar to Teuchos::TimeMonitor::summarize().
  ArrayRCP<double> ReduceMaxMinAvg(double localValue, Teuchos::Comm<int> const &comm, int rootNode = 0);

  //! Integration of Teuchos::TimeMonitor with MueLu verbosity system
  class TimerMonitor : public BaseClass {

  public:

    TimerMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel = Timings0)
    {

      // Inherit verbosity from 'object'
      SetVerbLevel(object.GetVerbLevel());
      setOStream(object.getOStream());

      if (IsPrint(timerLevel) && 
          /* disable timer if never printed: */ (IsPrint(RuntimeTimings) || (!IsPrint(NoTimeReport)))) {

        if (!IsPrint(NoTimeReport)) {
          // Start the timer
          timer_ = Teuchos::TimeMonitor::getNewTimer("MueLu: " + msg); //TODO: Is it useful to create the timer with this static function? 
	                                                               //      I think that the timer will be register on the TimerMonitor list in the next line anyway.
          timeMonitor_ = rcp(new Teuchos::TimeMonitor(*timer_));
        } else {
          timer_ = rcp(new Teuchos::Time("MueLu: " + msg));
          timer_->start();
        }
      }
    }

    ~TimerMonitor() {
      if (timer_ != Teuchos::null) {

        // Stop the timer
        timeMonitor_ = Teuchos::null;
        timer_->stop(); // needed is timeMonitor_ not used (NoTimeReport option)
        
        if (IsPrint(RuntimeTimings)) {
          //FIXME: creates lot of barriers. An option to report time of proc0 only instead would be nice
          //FIXME: MPI_COMM_WORLD only... BTW, it is also the case in Teuchos::TimeMonitor...
          ArrayRCP<double> stats = ReduceMaxMinAvg(timer_->totalElapsedTime(), *Teuchos::DefaultComm<int>::getComm ());
          
          //FIXME: Not very important for now, but timer will be printed even if verboseLevel of Monitor/Object changed
          //       between Monitor constructor and destructor.
          GetOStream(RuntimeTimings, 0) << "Timer: " << " max=" << stats[0] << " min=" << stats[1] << " avg=" << stats[2] << std::endl;
        }
      }
    }

  protected:
    TimerMonitor() { }
    
  private:
    RCP<Teuchos::Time> timer_; // keep a reference on the timer to print stats if RuntimeTimings=ON
    RCP<Teuchos::TimeMonitor> timeMonitor_;
  };
  
  //! Manage indentation of output using OSTab and verbosity level
  class PrintMonitor : public BaseClass {

  public:

    PrintMonitor(const BaseClass& object, const std::string& msg, MsgType msgLevel = Runtime0) {

      // Inherit verbosity from 'object'
      SetVerbLevel(object.GetVerbLevel());
      setOStream(object.getOStream());
      
      // Print description and new indent
      if (IsPrint(msgLevel)) {
        GetOStream(msgLevel, 0) << msg << std::endl;
        tab_ = rcp(new Teuchos::OSTab(getOStream()));
      }

    }

    ~PrintMonitor() { }

  protected:
    PrintMonitor() { }

  private:
    RCP<Teuchos::OSTab> tab_;
  };
  
  // Main Monitor
  // A timer is created only if 'timerLevel' (Timings0 by default) is enable
  class Monitor: public BaseClass {
  public:
    Monitor(const BaseClass& object, const std::string & msg, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0) 
      : printMonitor_(object, msg + " (" + object.description() + ")", msgLevel),
        timerMonitor_(object, object.ShortClassName() + ": " + msg,    timerLevel)
    { }

  private:
    PrintMonitor printMonitor_;
    TimerMonitor timerMonitor_;
  };

  // Another version of the Monitor, that does not repeat the object description
  // A timer is created only if 'timerLevel' (Timings1 by default) is enable
  class SubMonitor: public BaseClass {
  public:
    SubMonitor(const BaseClass& object, const std::string & msg, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1) 
      : printMonitor_(object, msg, msgLevel),
        timerMonitor_(object, object.ShortClassName() + ": " + msg, timerLevel)
    { }

  private:
    PrintMonitor printMonitor_;
    TimerMonitor timerMonitor_;
  };

  // Factory monitor
  // Similar to Monitor but add a timer level by level
  class FactoryMonitor: public Monitor {
  public:
    FactoryMonitor(const BaseClass& object, const std::string & msg, int levelID, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0) 
      : Monitor(object, msg, msgLevel, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimerMonitor_ = rcp(new TimerMonitor(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    FactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0) 
      : Monitor(object, msg, msgLevel, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimerMonitor_ = rcp(new TimerMonitor(object, object.ShortClassName() + ": " +  msg + " (level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
      }
    }

  private:
    RCP<TimerMonitor> levelTimerMonitor_;
  };

  // Factory monitor
  // Similar to SubMonitor but add a timer level by level
  class SubFactoryMonitor: public SubMonitor {
  public:
    SubFactoryMonitor(const BaseClass& object, const std::string & msg, int levelID, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1) 
      : SubMonitor(object, msg, msgLevel, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimerMonitor_ = rcp(new TimerMonitor(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    SubFactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1) 
      : SubMonitor(object, msg, msgLevel, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimerMonitor_ = rcp(new TimerMonitor(object, object.ShortClassName() + ": " +  msg + " (level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
      }
    }
  private:
    RCP<TimerMonitor> levelTimerMonitor_;
  };

  //! This class wraps a Teuchos::Time and maintains a mutually exclusive property between wrapped timers.
  //! When a MutuallyExclusiveTimer is running, other timers are not running.
  //! Timers have three state: running, stopped or paused to enforce the mutually exclusive property.
  //! When the running timer is stopped, the last active timer is restarted. A stack of timers is used internally to support this functionality. 
  //! This class is useful to exclude from a timer the execution time of a subroutine.
  //!
  //! Example:
  //!
  //! Note: Only one timer can be active at a time but all timers can be inactive at the same time. Timers cannot be destroyed when they are in 'paused'.

  template<class TagName> //! The template parameter of this class can be used to define several set of mutually exclusive timer.
  class MutuallyExclusiveTimer : public BaseClass {
    
  public:
    
    MutuallyExclusiveTimer(const std::string &name, bool start=false)
      : timer_(name, false),  // second argument is false in any case, because if start==true, 
	                      // timer has to be started by MutuallyExclusiveTimer::start() instead of Teuchos::Time::start().
	isPaused_(false)
    { 
      if (start == true) timer_.start();
    }
    
    ~MutuallyExclusiveTimer() {
      // This timer can only be destroyed if it is not in the stack
      if (isPaused()) {
	// error message because cannot throw an exception in destructor
	GetOStream(Errors, 0) << "MutuallyExclusiveTimer::~MutuallyExclusiveTimer(): Error: destructor called on a paused timer." << std::endl;
	//TODO: Even if timing results will be wrong, the timer can be removed from the stack to avoid a segmentation fault.
      }
      
      stop(); // if isRunning(), remove from the stack, resume previous timer
    }
    
    //! Starts the timer. If a MutuallyExclusiveTimer is running, it will be stopped.
    //! Precondition: timer is not already paused
    //! Postcondition: timer is running. Other MutuallyExclusiveTimer are paused or stop.
    void start(bool reset=false) {
      TEUCHOS_TEST_FOR_EXCEPTION(isPaused(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTimer::start(): timer is paused. Use resume().");
      
      if (isRunning()) { return; } // If timer is already running, do not pause/push-in-the-stack/start the timer. 
                                   // Otherwise, something bad will happen when this.stop() will be called
      
      // pause currently running timer
      if (!timerStack_.empty()) {
	timerStack_.top().pause();
      }
      
      // start this timer
      timer_.start(reset);
      timerStack_.push(this);
    }
    
    // {@ Functions that can only be called on the most recent timer (= running timer or last paused timer)
    
    //!	Stops the timer. The previous MutuallyExclusiveTimer that has been paused when this timer was started will be resumed.
    // stop() can be called on an already stopped timer or on the currently running timer
    double stop() {
      TEUCHOS_TEST_FOR_EXCEPTION(isPaused(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTimer::start(): timer is paused. Use resume().");
      if (!isRunning()) { return; } // stop() can be called on stopped timer
      
      // Here, timer is running, so it is the head of the stack
      TopOfTheStack();
      
      timerStack_.pop();
      if (!timerStack_.empty()) {
	timerStack_.top().resume();
      }
      
    }
    
    //! Pause running timer. Used internally by start().
    void pause() {
      TopOfTheStack();
      
      timer_.stop();
      isPaused_ = true;
    }
    
    //! Resume paused timer. Used internally by stop()
    //! Precondition: timer is at the top of the stack
    //! Timer is not reset
    void resume() {
      TopOfTheStack();
      
      // no 'shortcut' test necessary: 
      // - if timer is stop, it is in pause (cannot be stop and not in pause because this timer is the head of the stack).
      // - if timer is running, nothing is changed by this function.

      timer_.start(false);
      isPaused_ = false;
    }

    // @}    

    //@{

    bool isRunning() {
      if (timer_.isRunning()) {
	TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.top() == this, Exceptions::RuntimeError, 
				   "MueLu::MutuallyExclusiveTimer::isRunning(): this timer is active so it is supposed to be the head of the stack");
      }
      return timer_.isRunning();
    }
    
    bool isPaused() {
      return isPaused_;
    }

    //@}

    //     //!
    //     static RCP<MutuallyExclusiveTimer> getNewTimer() {
    //       return rcp(new )
    //     }

  private:
    
    // This constructor is not public. Use MutuallyExclusiveTimer::getNewTimer() instead.
    // It is not public because I'm concerned that users will used Teuchos::Time::start()/stop() 
    //  instead of MutuallyExclusiveTimer::start()/stop() if they have access to the underlying Teuchos::Time object.
    // MutuallyExclusiveTimer(RCP<Teuchos::Time> timer)
    //       : timer_(timer) 
    //     { }
    
    MutuallyExclusiveTimer() { }
    
    Teuchos::Time timer_;
    bool isPaused_;

    //! Stack of started timers (active or paused timers).
    // - empty when no active timer
    // - head is the active timer
    // - other timers are timers paused to enforce the mutually exclusive property of the timer set.
    static std::stack<MutuallyExclusiveTimer*> timerStack_;

    void TopOfTheStack() {
      TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.empty(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTimer::TopOfTheStack(): timer is not the head of the stack (stack is empty).");
      TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.top() == this, Exceptions::RuntimeError,  "MueLu::MutuallyExclusiveTimer::TopOfTheStack(): timer is not the head of the stack.");
      TEUCHOS_TEST_FOR_EXCEPTION(!(isRunning() || isPaused()), Exceptions::RuntimeError,  "MueLu::MutuallyExclusiveTimer::TopOfTheStack(): head of the stack timer is neither active nor paused.");
    }
  
  };

  //TODO: test integrity of the stack:
  // Head = running or paused
  // Other timers of the stack = paused

} // namespace MueLu

#endif // MUELU_MONITOR_HPP
