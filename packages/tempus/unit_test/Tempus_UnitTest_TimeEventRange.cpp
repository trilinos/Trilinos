//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"
#include "Tempus_TimeEventRange.hpp"

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

  TEST_COMPARE(te->getName(), ==, "TimeEventRange (0; 0; 0)");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  TEST_FLOATING_EQUALITY(te->getTimeStart(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 1);

  TEST_FLOATING_EQUALITY(
      te->getRelTol(), std::numeric_limits<double>::epsilon() * 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(
      te->getAbsTol(), std::numeric_limits<double>::epsilon() * 100.0, 1.0e-14);

  TEST_COMPARE(te->getLandOnExactly(), ==, true);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Full_Construction_Stride)
{
  auto te = rcp(new Tempus::TimeEventRange<double>(0.0, PI, 1.0, "TestName",
                                                   true, 1.0e-14));

  // te->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 1.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 4);

  auto teRange2 = rcp(new Tempus::TimeEventRange<double>(
      -PI / 2.0, PI / 2.0, PI / 4.0, "teRange2", true, 1.0e-14));

  TEST_FLOATING_EQUALITY(teRange2->timeToNextEvent(0.1), PI / 4.0 - 0.1,
                         1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, Full_Construction_NumEvents)
{
  auto te = rcp(new Tempus::TimeEventRange<double>(0.0, PI, 5, "TestName", true,
                                                   1.0e-14));

  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), PI / 4.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 5);
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
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 1);

  // Reset stop.
  te->setTimeStop(4.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 2);

  // Reset stride.
  te->setTimeStride(0.5);
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.5, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 7);
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
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 0.5, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 7);

  // Large stride should be reset to stop_-start_.
  te->setTimeStride(5.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 2);

  // Stride smaller than relative tolerance should be reset to stop_-start_.
  te->setTimeStride(1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 3.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, setTimeRange)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());

  // Set with time range.
  te->setTimeRange(0.0, PI, 1.0);
  TEST_FLOATING_EQUALITY(te->getTimeStart(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStop(), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->getTimeStride(), 1.0, 1.0e-14);
  TEST_COMPARE(te->getNumEvents(), ==, 4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, isTime)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeRange(0.0, PI, 1.0);
  te->setRelTol(1.0e-14);

  //   Before first event.  (This is exactly one stride before range.)
  TEST_COMPARE(te->isTime(-1.0 + -10.0e-14), ==,
               false);  // Just outside tolerance.
  TEST_COMPARE(te->isTime(-1.0 + -0.1e-14), ==,
               false);                              // Just inside tolerance.
  TEST_COMPARE(te->isTime(-1.0 + 0.0), ==, false);  // Right on timeEvent.
  TEST_COMPARE(te->isTime(-1.0 + 0.1e-14), ==,
               false);  // Just inside tolerance.
  TEST_COMPARE(te->isTime(-1.0 + 10.0e-14), ==,
               false);  // Just outside tolerance.

  //   Around first event.
  TEST_COMPARE(te->isTime(-10.0e-14), ==, false);  // Just outside tolerance.
  TEST_COMPARE(te->isTime(-0.1e-14), ==, true);    // Just inside tolerance.
  TEST_COMPARE(te->isTime(0.0), ==, true);         // Right on timeEvent.
  TEST_COMPARE(te->isTime(0.1e-14), ==, true);     // Just inside tolerance.
  TEST_COMPARE(te->isTime(10.0e-14), ==, false);   // Just outside tolerance.

  //   Around mid event.
  TEST_COMPARE(te->isTime(1.0 + -10.0e-14), ==,
               false);                                 // Just outside tolerance.
  TEST_COMPARE(te->isTime(1.0 + -0.1e-14), ==, true);  // Just inside tolerance.
  TEST_COMPARE(te->isTime(1.0 + 0.0), ==, true);       // Right on timeEvent.
  TEST_COMPARE(te->isTime(1.0 + 0.1e-14), ==, true);   // Just inside tolerance.
  TEST_COMPARE(te->isTime(1.0 + 10.0e-14), ==,
               false);  // Just outside tolerance.

  //   Around last event.
  TEST_COMPARE(te->isTime(3.0 + -10.0e-14), ==,
               false);                                 // Just outside tolerance.
  TEST_COMPARE(te->isTime(3.0 + -0.1e-14), ==, true);  // Just inside tolerance.
  TEST_COMPARE(te->isTime(3.0 + 0.0), ==, true);       // Right on timeEvent.
  TEST_COMPARE(te->isTime(3.0 + 0.1e-14), ==, true);   // Just inside tolerance.
  TEST_COMPARE(te->isTime(3.0 + 10.0e-14), ==,
               false);  // Just outside tolerance.

  //   After last event.  (This is exactly one stride after range.)
  TEST_COMPARE(te->isTime(4.0 + -10.0e-14), ==,
               false);  // Just outside tolerance.
  TEST_COMPARE(te->isTime(4.0 + -0.1e-14), ==,
               false);                                 // Just inside tolerance.
  TEST_COMPARE(te->isTime(4.0 + 0.0), ==, false);      // Right on timeEvent.
  TEST_COMPARE(te->isTime(4.0 + 0.1e-14), ==, false);  // Just inside tolerance.
  TEST_COMPARE(te->isTime(4.0 + 10.0e-14), ==,
               false);  // Just outside tolerance.
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, timeToNextEvent)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeRange(0.0, PI, 1.0);
  te->setRelTol(1.0e-14);

  // Test timeToNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-10.0e-14), 10.0e-14,
                         1.0e-14);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-0.1e-14), 1.0 + 0.1e-14,
                         1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), 1.0 + 0.0,
                         1.0e-14);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.1e-14), 1.0 - 0.1e-14,
                         1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(10.0e-14), 1.0 - 10.0e-14,
                         1.0e-14);  // Just outside tolerance.

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0 + -10.0e-14), 10.0e-14,
                         1.0e-02);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0 + -0.1e-14), 1.0 + 0.1e-14,
                         1.0e-01);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0 + 0.0), 1.0 + 0.0,
                         1.0e-02);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0 + 0.1e-14), 1.0 - 0.1e-14,
                         1.0e-01);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(1.0 + 10.0e-14), 1.0 - 10.0e-14,
                         1.0e-14);  // Just outside tolerance.

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0 + -10.0e-14), 10.0e-14,
                         1.0e-02);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0 + -0.1e-14),
                         te->getDefaultTime(),
                         1.0e-01);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0 + 0.0), te->getDefaultTime(),
                         1.0e-02);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0 + 0.1e-14),
                         te->getDefaultTime(),
                         1.0e-01);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.0 + 10.0e-14),
                         te->getDefaultTime(),
                         1.0e-02);  // Just outside tolerance.
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, timeOfNextEvent)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeRange(0.0, PI, 1.0);
  te->setRelTol(1.0e-14);

  // Test timeOfNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-10.0e-14), 0.0,
                         1.0e-14);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-0.1e-14), 1.0,
                         1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(0.0), 1.0,
                         1.0e-14);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(0.1e-14), 1.0,
                         1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(10.0e-14), 1.0,
                         1.0e-14);  // Just outside tolerance.

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0 + -10.0e-14), 1.0,
                         1.0e-14);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0 + -0.1e-14), 2.0,
                         1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0 + 0.0), 2.0,
                         1.0e-14);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0 + 0.1e-14), 2.0,
                         1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.0 + 10.0e-14), 2.0,
                         1.0e-14);  // Just outside tolerance.

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0 + -10.0e-14), 3.0,
                         1.0e-14);  // Just outside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0 + -0.1e-14),
                         te->getDefaultTime(),
                         1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0 + 0.0), te->getDefaultTime(),
                         1.0e-14);  // Right on timeEvent.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0 + 0.1e-14),
                         te->getDefaultTime(),
                         1.0e-14);  // Just inside tolerance.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0 + 10.0e-14),
                         te->getDefaultTime(),
                         1.0e-14);  // Just outside tolerance.
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, eventInRange)
{
  auto te = rcp(new Tempus::TimeEventRange<double>());
  te->setTimeRange(0.0, PI, 1.0);
  te->setRelTol(1.0e-14);

  // Test eventInRange.
  //   Right end of input range.
  //   Around first event.
  TEST_COMPARE(te->eventInRange(-1.0, -10.0e-14), ==,
               false);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(-1.0, -0.1e-14), ==,
               true);                                   // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(-1.0, 0.0), ==, true);  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(-1.0, 0.1e-14), ==,
               true);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(-1.0, 10.0e-14), ==,
               true);  // Just outside tolerance.

  //   Around mid event.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 + -10.0e-14), ==,
               false);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 + -0.1e-14), ==,
               true);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 + 0.0), ==,
               true);  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 + 0.1e-14), ==,
               true);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(0.5, 1.0 + 10.0e-14), ==,
               true);  // Just outside tolerance.

  //   Around last event.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 + -10.0e-14), ==,
               false);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 + -0.1e-14), ==,
               true);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 + 0.0), ==,
               true);  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 + 0.1e-14), ==,
               true);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(2.5, 3.0 + 10.0e-14), ==,
               true);  // Just outside tolerance.

  //   Left end of input range.
  //   Around first event.
  TEST_COMPARE(te->eventInRange(-10.0e-14, 0.5), ==,
               true);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(-0.1e-14, 0.5), ==,
               false);                                  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(0.0, 0.5), ==, false);  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(0.1e-14, 0.5), ==,
               false);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(10.0e-14, 0.5), ==,
               false);  // Just outside tolerance.

  //   Around mid event.
  TEST_COMPARE(te->eventInRange(1.0 + -10.0e-14, 1.5), ==,
               true);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(1.0 + -0.1e-14, 1.5), ==,
               false);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(1.0 + 0.0, 1.5), ==,
               false);  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(1.0 + 0.1e-14, 1.5), ==,
               false);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(1.0 + 10.0e-14, 1.5), ==,
               false);  // Just outside tolerance.

  //   Around last event.
  TEST_COMPARE(te->eventInRange(3.0 + -10.0e-14, 4.0), ==,
               true);  // Just outside tolerance.
  TEST_COMPARE(te->eventInRange(3.0 + -0.1e-14, 4.0), ==,
               false);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(3.0 + 0.0, 4.0), ==,
               false);  // Right on timeEvent.
  TEST_COMPARE(te->eventInRange(3.0 + 0.1e-14, 4.0), ==,
               false);  // Just inside tolerance.
  TEST_COMPARE(te->eventInRange(3.0 + 10.0e-14, 4.0), ==,
               false);  // Just outside tolerance.
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, getValidParameters)
{
  auto ter = rcp(new Tempus::TimeEventRange<double>());

  auto pl = ter->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Type"), ==, "Range");
  TEST_COMPARE(pl->get<std::string>("Name"), ==, "TimeEventRange (0; 0; 0)");
  TEST_FLOATING_EQUALITY(pl->get<double>("Start Time"), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Stop Time"), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Stride Time"), 0.0, 1.0e-14);
  TEST_COMPARE(pl->get<int>("Number of Events"), ==, 1);
  TEST_FLOATING_EQUALITY(pl->get<double>("Relative Tolerance"),
                         std::numeric_limits<double>::epsilon() * 100.0,
                         1.0e-14);
  TEST_COMPARE(pl->get<bool>("Land On Exactly"), ==, true);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, createTimeEventRange)
{
  // Construct parameterList similar to getValidParameters().
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Event Range");

  pl->set("Name", "Unit Test Time Event Range");
  pl->set("Type", "Range");
  pl->set("Start Time", -0.1);
  pl->set("Stop Time", 1.1);
  pl->set("Stride Time", 0.1);
  pl->set("Relative Tolerance", 1.0e-10);
  pl->set("Land On Exactly", false);

  // Construct TimeEventRange from ParameterList.
  auto ter = Tempus::createTimeEventRange<double>(pl);

  // ter->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(ter->getName(), ==, "Unit Test Time Event Range");
  TEST_COMPARE(ter->getType(), ==, "Range");
  TEST_FLOATING_EQUALITY(ter->getTimeStart(), -0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(ter->getTimeStop(), 1.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(ter->getTimeStride(), 0.1, 1.0e-14);
  TEST_COMPARE(ter->getNumEvents(), ==, 13);
  TEST_FLOATING_EQUALITY(ter->getRelTol(), 1.0e-10, 1.0e-14);
  TEST_COMPARE(ter->getLandOnExactly(), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRange, SingleEventAtZero)
{
  auto ter = rcp(new Tempus::TimeEventRange<double>(0.0, 0.0, 0.0,
                                                    "SingleEventAtZero", true));
  ter->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(ter->getNumEvents(), ==, 1);

  TEST_COMPARE(ter->isTime(0.0), ==, true);
  TEST_FLOATING_EQUALITY(ter->timeToNextEvent(-1.0), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(ter->timeOfNextEvent(-1.0), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(ter->timeToNextEvent(0.0), ter->getDefaultTime(),
                         1.0e-14);
  TEST_FLOATING_EQUALITY(ter->timeOfNextEvent(0.0), ter->getDefaultTime(),
                         1.0e-14);
  TEST_COMPARE(ter->eventInRange(-1.0, 1.0), ==, true);
  TEST_COMPARE(ter->eventInRange(0.0, 1.0), ==, false);
  TEST_COMPARE(ter->eventInRange(0.0, 0.0), ==, false);
}

}  // namespace Tempus_Unit_Test
