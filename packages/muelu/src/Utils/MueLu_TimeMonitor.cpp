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
    label_ = cleanupLabel("MueLu: " + msg);

    if (!IsPrint(NoTimeReport)) {
      timer_ = rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(label_)));
    }
  }
}  // TimeMonitor::TimeMonitor()

TimeMonitor::TimeMonitor() {}

TimeMonitor::~TimeMonitor() {}

}  // namespace MueLu
