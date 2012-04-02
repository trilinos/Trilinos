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

  //! This class wraps a Teuchos::Time and maintains a mutually exclusive property between wrapped timers.
  //! When a MutuallyExclusiveTime is running, other timers are not running.
  //! Timers have three state: running, stopped or paused to enforce the mutually exclusive property.
  //! When the running timer is stopped, the last active timer is restarted. A stack of timers is used internally to support this functionality. 
  //! This class is useful to exclude from a timer the execution time of a subroutine.
  //!
  //! Example:
  //!
  //! Note: Only one timer can be active at a time but all timers can be inactive at the same time. Timers cannot be destroyed when they are in 'paused'.

  //TODO: inheritence from PerformanceMonitorBase<Time> ?

  template<class TagName> //! The template parameter of this class can be used to define several set of mutually exclusive timer.
  class MutuallyExclusiveTime : public BaseClass {
    
  public:
    
    MutuallyExclusiveTime(const std::string &name, bool start=false)
      : timer_(rcp(new Teuchos::Time(name, false))),  // second argument is false in any case, because if start==true, 
                                                      // timer has to be started by MutuallyExclusiveTime::start() instead of Teuchos::Time::start().
	isPaused_(false)
    { 
      if (start == true) timer_->start();
    }
    
    ~MutuallyExclusiveTime() {
      // This timer can only be destroyed if it is not in the stack
      if (isPaused()) {
	// error message because cannot throw an exception in destructor
	GetOStream(Errors, 0) << "MutuallyExclusiveTime::~MutuallyExclusiveTime(): Error: destructor called on a paused timer." << std::endl;
	//TODO: Even if timing results will be wrong, the timer can be removed from the stack to avoid a segmentation fault.
      }
      
      stop(); // if isRunning(), remove from the stack, resume previous timer
    }
    
    //! Starts the timer. If a MutuallyExclusiveTime is running, it will be stopped.
    //! Precondition: timer is not already paused
    //! Postcondition: timer is running. Other MutuallyExclusiveTime are paused or stop.
    void start(bool reset=false) {
      TEUCHOS_TEST_FOR_EXCEPTION(isPaused(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::start(): timer is paused. Use resume().");
      
      if (isRunning()) { return; } // If timer is already running, do not pause/push-in-the-stack/start the timer. 
                                   // Otherwise, something bad will happen when this.stop() will be called
      
      // pause currently running timer
      if (!timerStack_.empty()) {
	timerStack_.top()->pause();
      }
      
      // start this timer
      timer_->start(reset);
      timerStack_.push(this);
    }
    
    // {@ Functions that can only be called on the most recent timer (= running timer or last paused timer)
    
    //!	Stops the timer. The previous MutuallyExclusiveTime that has been paused when this timer was started will be resumed.
    // stop() can be called on an already stopped timer or on the currently running timer
    double stop() {
      TEUCHOS_TEST_FOR_EXCEPTION(isPaused(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::start(): timer is paused. Use resume().");
      if (!isRunning()) { return timer_->stop(); } // stop() can be called on stopped timer
      
      // Here, timer is running, so it is the head of the stack
      TopOfTheStack();
      
      timerStack_.pop();
      double r = timer_->stop();

      if (!timerStack_.empty()) {
	timerStack_.top()->resume();
      }
      
      return r;
    }
    
    //! Pause running timer. Used internally by start().
    void pause() {
      if (isPaused()) // calling twice pause() is allowed
        return;

      TopOfTheStack();
      
      timer_->stop();
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

      timer_->start(false);
      isPaused_ = false;
    }

    // @}    

    //@{

    bool isRunning() {
      if (timer_->isRunning()) {
	TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.top() != this, Exceptions::RuntimeError, 
				   "MueLu::MutuallyExclusiveTime::isRunning(): this timer is active so it is supposed to be the head of the stack");
      }
      return timer_->isRunning();
    }
    
    bool isPaused() {
      TEUCHOS_TEST_FOR_EXCEPTION(isPaused_ && timer_->isRunning(), Exceptions::RuntimeError, "");
      return isPaused_;
    }

    //@}

    //! Return a new MutuallyExclusiveTime that is register with the Teuchos::TimeMonitor (for timer summary)
    // Note: this function is provide by the Timer, not the monitor (!= Teuchos)
    static RCP<MutuallyExclusiveTime<TagName> > getNewTimer(const std::string& name) {
      RCP<MutuallyExclusiveTime<TagName> > timer = rcp(new MutuallyExclusiveTime<TagName>(Teuchos::TimeMonitor::getNewTimer(name)));
      return timer;
    }

    //! Increment the number of times this timer has been called. 
    void incrementNumCalls() { timer_->incrementNumCalls(); }

  private:
    
    // This constructor is not public because I'm concerned that users will used Teuchos::Time::start()/stop() 
    // instead of MutuallyExclusiveTime::start()/stop() if they have access to the underlying Teuchos::Time object.
    MutuallyExclusiveTime(RCP<Teuchos::Time> timer)
      : timer_(timer), isPaused_(false)
    { }
    
    // MutuallyExclusiveTime() { }
    
    RCP<Teuchos::Time> timer_; // using an RCP allows to use Teuchos::TimeMonitor to keep track of the timer
    bool isPaused_;

    //! Stack of started timers (active or paused timers).
    // - empty when no active timer
    // - head is the active timer
    // - other timers are timers paused to enforce the mutually exclusive property of the timer set.
    static std::stack<MutuallyExclusiveTime<TagName>*> timerStack_;

    void TopOfTheStack() {
      TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.empty(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::TopOfTheStack(): timer is not the head of the stack (stack is empty).");
      TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.top() != this, Exceptions::RuntimeError,  "MueLu::MutuallyExclusiveTime::TopOfTheStack(): timer is not the head of the stack.");
      TEUCHOS_TEST_FOR_EXCEPTION(!(isRunning() || isPaused()), Exceptions::RuntimeError,  "MueLu::MutuallyExclusiveTime::TopOfTheStack(): head of the stack timer is neither active nor paused.");
    }

  //TODO: test integrity of the stack:
  // Head = running or paused
  // Other timers of the stack = paused
  
  };

  // Helper function. Similar to Teuchos::TimeMonitor::summarize().
  ArrayRCP<double> ReduceMaxMinAvg(double localValue, Teuchos::Comm<int> const &comm, int rootNode = 0);

  //! Integration of Teuchos::TimeMonitor with MueLu verbosity system
  class TimeMonitor : public BaseClass {

  public:

    TimeMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel = Timings0)
    {

      // Inherit verbosity from 'object'
      SetVerbLevel(object.GetVerbLevel());
      setOStream(object.getOStream());

      if (IsPrint(timerLevel) && 
          /* disable timer if never printed: */ (IsPrint(RuntimeTimings) || (!IsPrint(NoTimeReport)))) {

        if (!IsPrint(NoTimeReport)) {
          // TODO: there is no function to register a timer in Teuchos::TimeMonitor after the creation of the timer. But would be useful...
          timer_ = Teuchos::TimeMonitor::getNewTimer("MueLu: " + msg);
        } else {
          timer_ = rcp(new Teuchos::Time("MueLu: " + msg));
        }

        // Start the timer (this is what is done by Teuchos::TimeMonitor)
        timer_->start();
        timer_->incrementNumCalls();
      }
    }

    ~TimeMonitor() {
      if (timer_ != Teuchos::null) {

        // Stop the timer
        timer_->stop();
        
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
    TimeMonitor() { }
    
  private:
    RCP<Teuchos::Time> timer_;
  };

  //TODO: code duplication MutuallyExclusiveTimeMonitor / TimeMonitor

  template <class TagName>
  class MutuallyExclusiveTimeMonitor : public BaseClass {

  public:

    MutuallyExclusiveTimeMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel = Timings0)
    {
      // Inherit verbosity from 'object'
      SetVerbLevel(object.GetVerbLevel());
      setOStream(object.getOStream());

      if (IsPrint(timerLevel) && 
          /* disable timer if never printed: */ (IsPrint(RuntimeTimings) || (!IsPrint(NoTimeReport)))) {

        if (!IsPrint(NoTimeReport)) {
          timer_ = MutuallyExclusiveTime<TagName>::getNewTimer("MueLu: " + msg + " (ME)");
        } else {
          timer_ = rcp(new MutuallyExclusiveTime<TagName>("MueLu: " + msg + " (ME)"));
        }

        timer_->start();
        timer_->incrementNumCalls();
      }
    }

    ~MutuallyExclusiveTimeMonitor() {
      if (timer_ != Teuchos::null) {

        // Stop the timer
        timer_->stop();
        
        if (IsPrint(RuntimeTimings)) {
          //FIXME: creates lot of barriers. An option to report time of proc0 only instead would be nice
          //FIXME: MPI_COMM_WORLD only... BTW, it is also the case in Teuchos::TimeMonitor...
	  //TODO          ArrayRCP<double> stats = ReduceMaxMinAvg(timer_->totalElapsedTime(), *Teuchos::DefaultComm<int>::getComm ());
          
          //FIXME: Not very important for now, but timer will be printed even if verboseLevel of Monitor/Object changed
          //       between Monitor constructor and destructor.
          //TODO GetOStream(RuntimeTimings, 0) << "Timer: " << " max=" << stats[0] << " min=" << stats[1] << " avg=" << stats[2] << std::endl;
        }
      }
    }

  protected:
    MutuallyExclusiveTimeMonitor() { }
    
  private:
    RCP<MutuallyExclusiveTime<TagName> > timer_; // keep a reference on the timer to print stats if RuntimeTimings=ON //TODO:use base class instead
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
    TimeMonitor timerMonitor_;
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
    TimeMonitor timerMonitor_;
  };

  // Factory monitor
  // Similar to Monitor but add a timer level by level
  class FactoryMonitor: public Monitor {
  public:
    FactoryMonitor(const BaseClass& object, const std::string & msg, int levelID, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0) 
      : Monitor(object, msg, msgLevel, timerLevel),
        timerMonitorExclusive_(object, object.ShortClassName() + ": " + msg, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
	levelTimeMonitorExclusive_ = rcp(new MutuallyExclusiveTimeMonitor<Level>(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    //TODO: code factorization
    FactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0) 
      : Monitor(object, msg, msgLevel, timerLevel),
        timerMonitorExclusive_(object, object.ShortClassName() + ": " + msg, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " +  msg + " (level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
	levelTimeMonitorExclusive_ = rcp(new MutuallyExclusiveTimeMonitor<Level>(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
      }
    }

  private:
    RCP<TimeMonitor> levelTimeMonitor_;

    MutuallyExclusiveTimeMonitor<FactoryBase>  timerMonitorExclusive_;
    RCP<MutuallyExclusiveTimeMonitor<Level> >  levelTimeMonitorExclusive_;
  };

  // Factory monitor
  // Similar to SubMonitor but add a timer level by level
  class SubFactoryMonitor: public SubMonitor {
  public:
    SubFactoryMonitor(const BaseClass& object, const std::string & msg, int levelID, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1) 
      : SubMonitor(object, msg, msgLevel, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    SubFactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1) 
      : SubMonitor(object, msg, msgLevel, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " +  msg + " (level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
      }
    }
  private:
    RCP<TimeMonitor> levelTimeMonitor_;
  };

} // namespace MueLu

#endif // MUELU_MONITOR_HPP


