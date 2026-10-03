// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_Monitor.hpp"
int MueLu::FactoryMonitor::timerIdentifier_ = 0;

namespace MueLu {
PrintMonitor::PrintMonitor(const BaseClass& object, const std::string& msg, MsgType msgLevel)
  : object_(object) {
  tabbed = false;
  if (object_.IsPrint(msgLevel)) {
    // Print description and new indent
    object_.GetOStream(msgLevel, 0) << msg << std::endl;
    object_.getOStream()->pushTab();
    tabbed = true;
  }
}

PrintMonitor::~PrintMonitor() {
  if (tabbed) object_.getOStream()->popTab();
}

Monitor::Monitor(const BaseClass& object, const std::string& msg, MsgType msgLevel, MsgType timerLevel)
  : printMonitor_(object, msg + " (" + object.description() + ")", msgLevel)
  , timerMonitor_(object, object.ShortClassName() + ": " + msg, timerLevel) {}

Monitor::Monitor(const BaseClass& object, const std::string& msg, const std::string& label, MsgType msgLevel, MsgType timerLevel)
  : printMonitor_(object, label + msg + " (" + object.description() + ")", msgLevel)
  , timerMonitor_(object, label + object.ShortClassName() + ": " + msg, timerLevel) {}

Monitor::~Monitor() = default;

SubMonitor::SubMonitor(const BaseClass& object, const std::string& msg, MsgType msgLevel, MsgType timerLevel)
  : printMonitor_(object, msg, msgLevel)
  , timerMonitor_(object, object.ShortClassName() + ": " + msg, timerLevel) {}

SubMonitor::SubMonitor(const BaseClass& object, const std::string& msg, const std::string& label, MsgType msgLevel, MsgType timerLevel)
  : printMonitor_(object, label + msg, msgLevel)
  , timerMonitor_(object, label + object.ShortClassName() + ": " + msg, timerLevel) {}

SubMonitor::~SubMonitor() = default;

FactoryMonitor::FactoryMonitor(const BaseClass& object, const std::string& msg, int levelID, MsgType msgLevel, MsgType timerLevel)
  : Monitor(object, msg, msgLevel, timerLevel) {
}

FactoryMonitor::FactoryMonitor(const BaseClass& object, const std::string& msg, const Level& level, MsgType msgLevel, MsgType timerLevel)
  : Monitor(object, msg, FormattingHelper::getColonLabel(level.getObjectLabel()), msgLevel, timerLevel) {
}

FactoryMonitor::~FactoryMonitor() = default;

SubFactoryMonitor::SubFactoryMonitor(const BaseClass& object, const std::string& msg, int levelID, MsgType msgLevel, MsgType timerLevel)
  : SubMonitor(object, msg, msgLevel, timerLevel) {
}

SubFactoryMonitor::SubFactoryMonitor(const BaseClass& object, const std::string& msg, const Level& level, MsgType msgLevel, MsgType timerLevel)
  : SubMonitor(object, msg, FormattingHelper::getColonLabel(level.getObjectLabel()), msgLevel, timerLevel) {
}

SubFactoryMonitor::~SubFactoryMonitor() = default;

}  // namespace MueLu
