#ifndef MUELU_MONITOR_HPP
#define MUELU_MONITOR_HPP

#include <string>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_VerboseObject.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_MutuallyExclusiveTime.hpp"
#include "MueLu_TimeMonitor.hpp"

namespace MueLu {

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
  // A timer is created only if 'timerLevel' (Timings0 by default) is enabled
  class Monitor: public BaseClass {
  public:
    Monitor(const BaseClass& object, const std::string & msg, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0) 
      : printMonitor_(object, msg + " (" + object.description() + ")", msgLevel),
        timerMonitor_(object, object.ShortClassName() + ": " + msg + " (total)",    timerLevel)
    { }

  private:
    PrintMonitor printMonitor_;
    TimeMonitor timerMonitor_;
  };

  // Another version of the Monitor, that does not repeat the object description
  // A timer is created only if 'timerLevel' (Timings1 by default) is enabled
  class SubMonitor: public BaseClass {
  public:
    SubMonitor(const BaseClass& object, const std::string & msg, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1) 
      : printMonitor_(object, msg, msgLevel),
        timerMonitor_(object, object.ShortClassName() + ": " + msg + " (sub, total)", timerLevel)
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
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (total, level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
	levelTimeMonitorExclusive_ = rcp(new MutuallyExclusiveTimeMonitor<Level>(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    //TODO: code factorization
    FactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime0, MsgType timerLevel = Timings0) 
      : Monitor(object, msg, msgLevel, timerLevel),
        timerMonitorExclusive_(object, object.ShortClassName() + ": " + msg, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " +  msg + " (total, level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
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
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (sub, total, level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
      }
    }

    SubFactoryMonitor(const BaseClass& object, const std::string & msg, const Level & level, MsgType msgLevel = Runtime1, MsgType timerLevel = Timings1) 
      : SubMonitor(object, msg, msgLevel, timerLevel)
    { 
      if (IsPrint(TimingsByLevel)) {
        levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " +  msg + " (sub, total, level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
      }
    }
  private:
    RCP<TimeMonitor> levelTimeMonitor_;
  };

} // namespace MueLu

#endif // MUELU_MONITOR_HPP
