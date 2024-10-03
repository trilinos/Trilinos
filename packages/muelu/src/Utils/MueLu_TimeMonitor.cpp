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

  if (IsPrint(timerLevel) &&
      /* disable timer if never printed: */ (IsPrint(RuntimeTimings) || (!IsPrint(NoTimeReport)))) {
    if (!IsPrint(NoTimeReport)) {
      // TODO: there is no function to register a timer in Teuchos::TimeMonitor after the creation of the timer. But would be useful...
      timer_ = Teuchos::TimeMonitor::getNewTimer(cleanupLabel("MueLu: " + msg));
    } else {
      timer_ = rcp(new Teuchos::Time(cleanupLabel("MueLu: " + msg)));
    }

    // Start the timer (this is what is done by Teuchos::TimeMonitor)
    timer_->start();
    timer_->incrementNumCalls();
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
    const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
    if (nonnull(stackedTimer))
      stackedTimer->start(timer_->name());
#endif
  }
}  // TimeMonitor::TimeMonitor()

TimeMonitor::TimeMonitor() {}

TimeMonitor::~TimeMonitor() {
  // Stop the timer if present
  if (timer_ != Teuchos::null) {
    timer_->stop();
#ifdef HAVE_TEUCHOS_ADD_TIME_MONITOR_TO_STACKED_TIMER
    try {
      const auto stackedTimer = Teuchos::TimeMonitor::getStackedTimer();
      if (nonnull(stackedTimer))
        stackedTimer->stop(timer_->name());
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
#endif
  }
}  // TimeMonitor::~TimeMonitor()

template class MutuallyExclusiveTimeMonitor<FactoryBase>;
template class MutuallyExclusiveTimeMonitor<Level>;

}  // namespace MueLu
