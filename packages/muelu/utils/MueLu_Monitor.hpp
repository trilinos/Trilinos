#ifndef MUELU_MONITOR_HPP
#define MUELU_MONITOR_HPP

#include <string>
#include <Teuchos_Time.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VerboseObject.hpp"

//TODO: change timing default value to true

namespace MueLu {

  class Monitor : public BaseClass {

  public:

    Monitor(const BaseClass& object, const std::string & descr, bool timing = false) {
      init(object, descr + " (" + object.description() + ")", timing);
    }

    ~Monitor() {
      if (timer_ != Teuchos::null) {
        timer_->stop();
        //MemUtils::ReportTimeAndMemory(*timer_, comm_);
      }
    }

  protected:
    Monitor() {};

    void init (const VerboseObject& verbObject, const std::string & descr, bool timing = false) {
      // Set Monitor description
      description_ = descr;

      // Inherits the verbLevel and ostream from verbObject
      setVerbLevel(verbObject.getVerbLevel()); // Teuchos (unused)
      SetVerbLevel(verbObject.GetVerbLevel()); // MueLu
      setOStream(verbObject.getOStream());
      
      // Print description
      GetOStream(Runtime0, 0) << description_ << std::endl;
      // verbObject.describe(getOStream(), getVerbLevel()); //TODO

      tab_ = rcp(new Teuchos::OSTab(getOStream()));

      // Start timer
      if (timing) {
        timer_ = rcp(new Teuchos::Time(description_));
        timer_->start(true);
      }
    }
    
  private:
    std::string description_;
    RCP<Teuchos::Time> timer_; //TODO use TimeMonitor?
    RCP<Teuchos::OSTab> tab_;
  };

  class SubMonitor: public Monitor {
  public:
    SubMonitor(const BaseClass& object, const std::string & descr, bool timing = false) {
      init(object, descr, timing);
    }
  };

}

#endif // MUELU_MONITOR_HPP

// TODO       describe(GetOStream(Parameters0), getVerbLevel());
