#ifndef MUELU_MONITOR_HPP
#define MUELU_MONITOR_HPP

#include <string>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_VerboseObject.hpp"

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
          timer_ = Teuchos::TimeMonitor::getNewTimer("MueLu: " + msg);
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
  // - A timer is created only if Timings0 == true
  // - Timer is printed on the go if RuntimeTimings)) == true
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
  // A timer is created only if Timings1 == true
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
  // Add a timer level by level
  class FactoryMonitor: public Monitor {
  public:
    FactoryMonitor(const BaseClass& object, int levelID, const std::string & msg, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0) 
      : Monitor(object, msg, msgLevel, Timings0)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimerMonitor_ = rcp(new TimerMonitor(object, object.ShortClassName() + "(" + Teuchos::Utils::toString(levelID) + "): " + msg, timerLevel));
      }
    }

  private:
    RCP<TimerMonitor> levelTimerMonitor_;
  };

  // Factory monitor
  // Add a timer level by level
  class SubFactoryMonitor: public SubMonitor {
  public:
    SubFactoryMonitor(const BaseClass& object, int levelID, const std::string & msg, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1) 
      : SubMonitor(object, msg, msgLevel, Timings0)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimerMonitor_ = rcp(new TimerMonitor(object, object.ShortClassName() + "(" + Teuchos::Utils::toString(levelID) + "): " + msg, timerLevel));
      }
    }

  private:
    RCP<TimerMonitor> levelTimerMonitor_;
  };

} // namespace MueLu

#endif // MUELU_MONITOR_HPP
