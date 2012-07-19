#ifndef MUELU_TIMEMONITOR_HPP
#define MUELU_TIMEMONITOR_HPP

#include <string>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"
#include "MueLu_VerboseObject.hpp"
#include "MueLu_MutuallyExclusiveTime.hpp"

namespace MueLu {

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
          timer_ = MutuallyExclusiveTime<TagName>::getNewTimer("MueLu: " + msg /*+ " (MutuallyExclusive)" */);
        } else {
          timer_ = rcp(new MutuallyExclusiveTime<TagName>     ("MueLu: " + msg /*+ " (MutuallyExclusive)" */));
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
  
} // namespace MueLu

#endif // MUELU_TIMEMONITOR_HPP


