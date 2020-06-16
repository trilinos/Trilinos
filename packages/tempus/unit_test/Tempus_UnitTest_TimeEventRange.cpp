// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_DefaultComm.hpp"

#include "Tempus_TimeEventRange.hpp"

#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <cmath>
static double PI = M_PI;

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventRange");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  TEST_FLOATING_EQUALITY(te->getTimeStart (), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 1);

  TEST_FLOATING_EQUALITY(te->getRelTol(), 1.0e-14, 1.0e-14);
  TEST_COMPARE(te->getLandOnExactly(), ==, true);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Construction_Stride)
{
  auto te = rcp(new Tempus::TimeEventRange<double>(
                "TestName", 0.0, PI, 1.0, 1.0e-14, true));

  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 0.0 , 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 1.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 4);

  auto teRange2 = rcp(new Tempus::TimeEventRange<double>(
      "teRange2", -PI/2.0, PI/2.0, PI/4.0, 1.0e-14, true));

  TEST_FLOATING_EQUALITY(teRange2->timeToNextEvent(0.1), PI/4.0-0.1, 1.0e-14);

}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Construction_NumEvents)
{
  auto te = rcp(new Tempus::TimeEventRange<double>(
                "TestName", 0.0, PI, 5, 1.0e-14, true));

  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 0.0,      1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), PI,     1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), PI/4.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 5);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Basic_Accessors)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());

  te->setRelTol(0.1);
  TEST_FLOATING_EQUALITY(te->getRelTol(), 0.1, 1.0e-14);
  te->setRelTol(1.0e-14);
  te->setLandOnExactly(false);
  TEST_COMPARE(te->getLandOnExactly(), ==, false);
  te->setLandOnExactly(true);

  // Reset start after stop.
  te->setTimeStart(1.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 1);

  // Reset stop.
  te->setTimeStop(4.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 2);

  // Reset stride.
  te->setTimeStride(0.5);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.5, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 7);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Stride)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeStart(1.0);
  te->setTimeStop(4.0);
  te->setTimeStride(0.5);

  // Negative stride should be reset to stop_-start_.
  te->setTimeStride(-0.5);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 2);

  // Large stride should be reset to stop_-start_.
  te->setTimeStride(5.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 2);

  // Stride smaller than relative tolerance should be reset to stop_-start_.
  te->setTimeStride(1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 2);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, setTimeRange)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());

  // Set with time range.
  te->setTimeRange(0.0, PI, 1.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart (), 0.0,  1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop  (), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 1.0,  1.0e-14);
  TEST_COMPARE(te->getNumEvents (), ==, 4);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, isTime)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeRange(0.0, PI, 1.0);

  // Test isTime.
  //   Around first event.
  TEST_COMPARE(te->isTime(-10.0e-14), ==, false);   // Just outside tolerance.
  TEST_COMPARE(te->isTime( -1.0e-14), ==, true );   // Just inside tolerance.
  TEST_COMPARE(te->isTime(  0.0    ), ==, true );   // Right on timeEvent.
  TEST_COMPARE(te->isTime(  1.0e-14), ==, true );   // Just inside tolerance.
  TEST_COMPARE(te->isTime( 10.0e-14), ==, false);   // Just outside tolerance.

  //   Around mid event.
  TEST_COMPARE(te->isTime(1.0 + -10.0e-14), ==, false); // Just outside tolerance.
  TEST_COMPARE(te->isTime(1.0 +  -1.0e-14), ==, true ); // Just inside tolerance.
  TEST_COMPARE(te->isTime(1.0 +   0.0    ), ==, true ); // Right on timeEvent.
  TEST_COMPARE(te->isTime(1.0 +   1.0e-14), ==, true ); // Just inside tolerance.
  TEST_COMPARE(te->isTime(1.0 +  10.0e-14), ==, false); // Just outside tolerance.

  //   Around last event.
  TEST_COMPARE(te->isTime(3.0 + -10.0e-14), ==, false); // Just outside tolerance.
  TEST_COMPARE(te->isTime(3.0 +  -1.0e-14), ==, true ); // Just inside tolerance.
  TEST_COMPARE(te->isTime(3.0 +   0.0    ), ==, true ); // Right on timeEvent.
  TEST_COMPARE(te->isTime(3.0 +   1.0e-14), ==, true ); // Just inside tolerance.
  TEST_COMPARE(te->isTime(3.0 +  10.0e-14), ==, false); // Just outside tolerance.
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, timeToNextEvent)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeRange(0.0, PI, 1.0);

  // Test timeToNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-10.0e-14),  1.0e-13, 1.0e-14);     // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( -1.0e-14),  1.0e-14, 1.0e-14);     // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(  0.0    ),  0.0    , 1.0e-14);     // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(  1.0e-14), -1.0e-14, 1.0e-14);     // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( 10.0e-14), 1.0-1.0e-13, 1.0e-14);  // Just outside tolerance.

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+ -10.0e-14),  1.0e-13, 1.0e-02);    // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+  -1.0e-14),  1.0e-14, 1.0e-01);    // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+   0.0    ),  0.0    , 1.0e-02);    // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+   1.0e-14), -1.0e-14, 1.0e-01);    // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0+  10.0e-14), 1.0-1.0e-13, 1.0e-14); // Just outside tolerance.

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+ -10.0e-14),  1.0e-13, 1.0e-02);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+  -1.0e-14),  1.0e-14, 1.0e-01);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+   0.0    ),  0.0    , 1.0e-02);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+   1.0e-14), -1.0e-14, 1.0e-01);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0+  10.0e-14), -1.0e-13, 1.0e-02);  // Just outside tolerance.
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, timeOfNextEvent)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeRange(0.0, PI, 1.0);

  // Test timeOfNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-10.0e-14), 0.0, 1.0e-14);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( -1.0e-14), 0.0, 1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(  0.0    ), 0.0, 1.0e-14);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(  1.0e-14), 0.0, 1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( 10.0e-14), 1.0, 1.0e-14);  // Just outside tolerance.

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+ -10.0e-14), 1.0, 1.0e-14);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+  -1.0e-14), 1.0, 1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+   0.0    ), 1.0, 1.0e-14);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+   1.0e-14), 1.0, 1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0+  10.0e-14), 2.0, 1.0e-14);  // Just outside tolerance.

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+ -10.0e-14), 3.0, 1.0e-14);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+  -1.0e-14), 3.0, 1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+   0.0    ), 3.0, 1.0e-14);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+   1.0e-14), 3.0, 1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0+  10.0e-14), 3.0, 1.0e-14);  // Just outside tolerance.
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, eventInRange)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeRange(0.0, PI, 1.0);

  // Test eventInRange.
  //   Right end.
  //   Around first event.
  TEST_COMPARE(te->eventInRange(-1.0, -10.0e-14), ==, false);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(-1.0,  -1.0e-14), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(-1.0,   0.0    ), ==, true );  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(-1.0,   1.0e-14), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(-1.0,  10.0e-14), ==, true );  // Just outside tolerance.

  //   Around mid event.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 + -10.0e-14), ==, false);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 +  -1.0e-14), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 +   0.0    ), ==, true );  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 +   1.0e-14), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 +  10.0e-14), ==, true );  // Just outside tolerance.

  //   Around last event.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 + -10.0e-14), ==, false);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 +  -1.0e-14), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 +   0.0    ), ==, true );  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 +   1.0e-14), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 +  10.0e-14), ==, true );  // Just outside tolerance.

  //   Left end.
  //   Around first event.
  TEST_COMPARE(te->eventInRange(-10.0e-14, 0.5), ==, true );  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange( -1.0e-14, 0.5), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(  0.0    , 0.5), ==, true );  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(  1.0e-14, 0.5), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange( 10.0e-14, 0.5), ==, false);  // Just outside tolerance.

  //   Around mid event.
  TEST_COMPARE(te->eventInRange(1.0 + -10.0e-14, 1.5), ==, true );  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(1.0 +  -1.0e-14, 1.5), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(1.0 +   0.0    , 1.5), ==, true );  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(1.0 +   1.0e-14, 1.5), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(1.0 +  10.0e-14, 1.5), ==, false);  // Just outside tolerance.

  //   Around last event.
  TEST_COMPARE(te->eventInRange(3.0 + -10.0e-14, 4.0), ==, true );  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(3.0 +  -1.0e-14, 4.0), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(3.0 +   0.0    , 4.0), ==, true );  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(3.0 +   1.0e-14, 4.0), ==, true );  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(3.0 +  10.0e-14, 4.0), ==, false);  // Just outside tolerance.
}


} // namespace Tempus_Test
