// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <sstream>
#include <string>

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_BaseClass.hpp"
#include "MueLu_Monitor.hpp"

namespace {

class TimeMonitorTestObject : public MueLu::BaseClass {
 public:
  std::string ShortClassName() const override { return "TimeMonitorUnitTest"; }
};

class ScopedStackedTimer {
 public:
  explicit ScopedStackedTimer(const Teuchos::RCP<Teuchos::StackedTimer>& replacement)
    : original_(Teuchos::TimeMonitor::getStackedTimer()) {
    Teuchos::TimeMonitor::setStackedTimer(replacement);
  }

  ~ScopedStackedTimer() {
    Teuchos::TimeMonitor::setStackedTimer(original_);
  }

 private:
  Teuchos::RCP<Teuchos::StackedTimer> original_;
};

TEUCHOS_UNIT_TEST(TimeMonitor, FlatAndLevelTimers) {
  // Create timers for 5 different levels and a total timer and check that
  // they are all created and called the correct number of times
  constexpr int numLevels      = 5;
  const std::string totalLabel = "MueLu: TimeMonitorUnitTest: Computing Ac (total)";

  ScopedStackedTimer stackedTimerGuard(Teuchos::null);

  TimeMonitorTestObject object;
  object.SetVerbLevel(MueLu::Timings);

  for (int levelID = 0; levelID < numLevels; ++levelID) {
    MueLu::FactoryMonitor monitor(object, "Computing Ac", levelID);
  }

  const auto totalTimer = Teuchos::TimeMonitor::lookupCounter(totalLabel);
  // Check if timer is registered
  TEST_ASSERT(!totalTimer.is_null());
  if (!totalTimer.is_null())
    // Check if timer is called the correct number of times
    TEST_EQUALITY(totalTimer->numCalls(), numLevels);

  const auto comm = Teuchos::DefaultComm<int>::getComm();
  std::ostringstream report;
  Teuchos::TimeMonitor::report(comm.ptr(), report, "MueLu: TimeMonitorUnitTest: Computing Ac");
  if (comm->getRank() == 0)
    // Confirm that the total timer is present in the report
    TEST_ASSERT(report.str().find(totalLabel) != std::string::npos);

  for (int levelID = 0; levelID < numLevels; ++levelID) {
    // Construct expected name of level timer
    const std::string levelLabel = "MueLu: TimeMonitorUnitTest: Computing Ac (total, level=" +
                                   std::to_string(levelID) + ")";
    const auto levelTimer = Teuchos::TimeMonitor::lookupCounter(levelLabel);
    TEST_ASSERT(!levelTimer.is_null());
    if (!levelTimer.is_null())
      TEST_EQUALITY(levelTimer->numCalls(), 1);
    if (comm->getRank() == 0)
      // Confirm that each level timer is present in the report
      TEST_ASSERT(report.str().find(levelLabel) != std::string::npos);
  }
}

}  // namespace
