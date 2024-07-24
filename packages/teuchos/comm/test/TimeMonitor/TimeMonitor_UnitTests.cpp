// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"


namespace {


int someFunction()
{
  TEUCHOS_FUNC_TIME_MONITOR("someFunction1");
  TEUCHOS_FUNC_TIME_MONITOR_DIFF("someFunction2", diff);
  return 5;
}


TEUCHOS_UNIT_TEST( TimeMonitor, someFunction_timed )
{
  const int rtn = someFunction();
  TEST_EQUALITY(rtn, 5);
}

TEUCHOS_UNIT_TEST( TimeMonitor, formatting)
{
    using namespace Teuchos;

    // Zero out previous timers
    TimeMonitor::zeroOutTimers();

    std::ostringstream out1, out2;
    const std::string filter = "";

    // Check std::fixed formatting
    out1 << std::fixed;
    TimeMonitor::report(out1, filter);
    bool test1 = (out1.str().find("someFunction1    0.0000 (0)") != std::string::npos);
    TEST_EQUALITY(test1, true);

    // Check std::scientific formatting
    out2 << std::scientific;
    TimeMonitor::report(out2, filter);
    bool test2 = (out2.str().find("someFunction1    0.0000e+00 (0)") != std::string::npos);
    TEST_EQUALITY(test2, true);
}


} // namespace
