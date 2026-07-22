// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TimeMonitor.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

namespace {
//! Certain output formats must be observed
std::string cleanupLabel(const std::string& label) {
  // Our main interest is to remove double quotes from the label
  // since tracing formats can be written in JSON/YAML
  auto clean_label                     = label;
  std::vector<std::string> bad_values  = {"\""};
  std::vector<std::string> good_values = {"'"};
  for (size_t i = 0; i < bad_values.size(); ++i) {
    const auto& bad_value  = bad_values[i];
    const auto& good_value = good_values[i];
    size_t pos             = 0;
    while ((pos = clean_label.find(bad_value, pos)) != std::string::npos) {
      clean_label.replace(pos, bad_value.size(), good_value);
      pos += good_value.size();
    }
  }
  return clean_label;
}
}  // namespace

TimeMonitor::TimeMonitor(const BaseClass& object, const std::string& msg, MsgType timerLevel) {
  // Inherit props from 'object'
  SetVerbLevel(object.GetVerbLevel());
  SetProcRankVerbose(object.GetProcRankVerbose());

#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
  useStackedTimer_ = !Teuchos::TimeMonitor::stackedTimerNameIsDefault();
#else
  useStackedTimer_ = false;
#endif

  if (IsPrint(timerLevel) &&
      /* disable timer if never printed: */ (IsPrint(RuntimeTimings) || (!IsPrint(NoTimeReport)))) {
    label_ = cleanupLabel("MueLu: " + msg);

    if (!IsPrint(NoTimeReport)) {
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
      if (useStackedTimer_) {
        const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
        stackedTimer->start(label_);
      } else
#endif
      {
        // TODO: there is no function to register a timer in Teuchos::TimeMonitor after the creation of the timer. But would be useful...
        static std::map<std::string, RCP<Teuchos::Time>> timers;
        auto it = timers.find(label_);
        if (it == timers.end()) {
          timer_         = Teuchos::TimeMonitor::getNewTimer(label_);
          timers[label_] = timer_;
        } else {
          timer_ = it->second;
          // Start the timer (this is what is done by Teuchos::TimeMonitor)
          timer_->incrementNumCalls();
          timer_->start();
        }
      }
    } else {
      timer_ = rcp(new Teuchos::Time(label_));
      // Start the timer (this is what is done by Teuchos::TimeMonitor)
      timer_->incrementNumCalls();
      timer_->start();
    }
  }
}  // TimeMonitor::TimeMonitor()

TimeMonitor::TimeMonitor() {}

TimeMonitor::~TimeMonitor() {
  // Stop the timer if present

#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
  if (useStackedTimer_ && (!label_.empty())) {
    try {
      const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
      stackedTimer->stop(label_);
    } catch (std::runtime_error&) {
      std::ostringstream warning;
      warning << "\n*********************************************************************\n"
                 "WARNING: Overlapping timers detected!\n"
                 "A TimeMonitor timer was stopped before a nested subtimer was\n"
                 "stopped. This is not allowed by the StackedTimer. This corner case\n"
                 "typically occurs if the TimeMonitor is stored in an RCP and the RCP is\n"
                 "assigned to a new timer. To disable this warning, either fix the\n"
                 "ordering of timer creation and destuction or disable the StackedTimer\n";
      std::cout << warning.str() << std::endl;
      Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);
    }
  } else
#endif
  {
    if (timer_ != Teuchos::null) {
      timer_->stop();
    }
  }
}  // TimeMonitor::~TimeMonitor()

}  // namespace MueLu
