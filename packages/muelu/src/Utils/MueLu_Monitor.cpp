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
  , timerMonitor_(object, object.ShortClassName() + ": " + msg + " (total)", timerLevel) {}

Monitor::Monitor(const BaseClass& object, const std::string& msg, const std::string& label, MsgType msgLevel, MsgType timerLevel)
  : printMonitor_(object, label + msg + " (" + object.description() + ")", msgLevel)
  , timerMonitor_(object, label + object.ShortClassName() + ": " + msg + " (total)", timerLevel) {}

Monitor::~Monitor() = default;

SubMonitor::SubMonitor(const BaseClass& object, const std::string& msg, MsgType msgLevel, MsgType timerLevel)
  : printMonitor_(object, msg, msgLevel)
  , timerMonitor_(object, object.ShortClassName() + ": " + msg + " (sub, total)", timerLevel) {}

SubMonitor::SubMonitor(const BaseClass& object, const std::string& msg, const std::string& label, MsgType msgLevel, MsgType timerLevel)
  : printMonitor_(object, label + msg, msgLevel)
  , timerMonitor_(object, label + object.ShortClassName() + ": " + msg + " (sub, total)", timerLevel) {}

SubMonitor::~SubMonitor() = default;

FactoryMonitor::FactoryMonitor(const BaseClass& object, const std::string& msg, int levelID, MsgType msgLevel, MsgType timerLevel)
  : Monitor(object, msg, msgLevel, timerLevel)
  , timerMonitorExclusive_(object, object.ShortClassName() + ": " + msg, timerLevel) {
  if (object.IsPrint(TimingsByLevel)) {
    if (Teuchos::TimeMonitor::getStackedTimer().is_null())
      levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (total, level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
    levelTimeMonitorExclusive_ = rcp(new MutuallyExclusiveTimeMonitor<Level>(object, object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
  }
}

FactoryMonitor::FactoryMonitor(const BaseClass& object, const std::string& msg, const Level& level, MsgType msgLevel, MsgType timerLevel)
  : Monitor(object, msg, FormattingHelper::getColonLabel(level.getObjectLabel()), msgLevel, timerLevel)
  , timerMonitorExclusive_(object, FormattingHelper::getColonLabel(level.getObjectLabel()) + object.ShortClassName() + ": " + msg, timerLevel) {
  if (object.IsPrint(TimingsByLevel)) {
    std::string label = FormattingHelper::getColonLabel(level.getObjectLabel());
    if (Teuchos::TimeMonitor::getStackedTimer().is_null())
      levelTimeMonitor_ = rcp(new TimeMonitor(object, label + object.ShortClassName() + ": " + msg + " (total, level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
    levelTimeMonitorExclusive_ = rcp(new MutuallyExclusiveTimeMonitor<Level>(object, label + object.ShortClassName() + ": " + msg + " (level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
  }
}

FactoryMonitor::~FactoryMonitor() = default;

SubFactoryMonitor::SubFactoryMonitor(const BaseClass& object, const std::string& msg, int levelID, MsgType msgLevel, MsgType timerLevel)
  : SubMonitor(object, msg, msgLevel, timerLevel) {
  if (object.IsPrint(TimingsByLevel) && Teuchos::TimeMonitor::getStackedTimer().is_null())
    levelTimeMonitor_ = rcp(new TimeMonitor(object, object.ShortClassName() + ": " + msg + " (sub, total, level=" + Teuchos::Utils::toString(levelID) + ")", timerLevel));
}

SubFactoryMonitor::SubFactoryMonitor(const BaseClass& object, const std::string& msg, const Level& level, MsgType msgLevel, MsgType timerLevel)
  : SubMonitor(object, msg, FormattingHelper::getColonLabel(level.getObjectLabel()), msgLevel, timerLevel) {
  if (object.IsPrint(TimingsByLevel) && Teuchos::TimeMonitor::getStackedTimer().is_null()) {
    std::string label = FormattingHelper::getColonLabel(level.getObjectLabel());
    levelTimeMonitor_ = rcp(new TimeMonitor(object, label + object.ShortClassName() + ": " + msg + " (sub, total, level=" + Teuchos::Utils::toString(level.GetLevelID()) + ")", timerLevel));
  }
}

SubFactoryMonitor::~SubFactoryMonitor() = default;

}  // namespace MueLu
