// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
#include "Teuchos_StackedTimer.hpp"
#include <sstream>
#endif

namespace MueLu {

/*! @class TimeMonitor

    @brief Integrates Teuchos::TimeMonitor with MueLu verbosity system.
*/
class TimeMonitor : public BaseClass {
 public:
  TimeMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel = Timings0);

  ~TimeMonitor();

 protected:
  TimeMonitor();

 private:
  RCP<Teuchos::Time> timer_;
};  // class TimeMonitor

// TODO: code duplication MutuallyExclusiveTimeMonitor / TimeMonitor

/*! @class MutuallyExclusiveTimeMonitor

    @brief Similar to TimeMonitor, but uses MutuallyExclusiveTime objects.

*/
template <class TagName>
class MutuallyExclusiveTimeMonitor : public BaseClass {
 public:
  /*! @brief Constructor

     @param[in] object      Reference to the class instance that is creating this MutuallyExclusiveTimeMonitor.
     @param[in] msg         String that indicates what the Monitor is monitoring, e.g., "Build"
     @param[in] timerLevel  Governs whether timing information should be *gathered*.  Setting this to NoTimeReport prevents the creation of timers.
 */
  MutuallyExclusiveTimeMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel = Timings0) {
    // Inherit props from 'object'
    SetVerbLevel(object.GetVerbLevel());
    SetProcRankVerbose(object.GetProcRankVerbose());

    if (IsPrint(timerLevel) &&
        /* disable timer if never printed: */ (IsPrint(RuntimeTimings) || (!IsPrint(NoTimeReport)))) {
      if (!IsPrint(NoTimeReport)) {
        timer_ = MutuallyExclusiveTime<TagName>::getNewTimer("MueLu: " + msg /*+ " (MutuallyExclusive)" */);
      } else {
        timer_ = rcp(new MutuallyExclusiveTime<TagName>("MueLu: " + msg /*+ " (MutuallyExclusive)" */));
      }

      timer_->start();
      timer_->incrementNumCalls();
    }
  }

  ~MutuallyExclusiveTimeMonitor() {
    // Stop the timer if present
    if (timer_ != Teuchos::null)
      timer_->stop();
  }

 protected:
  MutuallyExclusiveTimeMonitor() {}

 private:
  RCP<MutuallyExclusiveTime<TagName> > timer_;  // keep a reference on the timer to print stats if RuntimeTimings=ON //TODO:use base class instead
};

extern template class MutuallyExclusiveTimeMonitor<FactoryBase>;
extern template class MutuallyExclusiveTimeMonitor<Level>;

}  // namespace MueLu

#endif  // MUELU_TIMEMONITOR_HPP
