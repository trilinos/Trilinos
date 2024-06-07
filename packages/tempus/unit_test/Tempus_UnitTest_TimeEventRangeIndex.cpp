//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_TimeEventRangeIndex.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventRangeIndex (0; 0; 1)");

  // Test when everything is zero.
  TEST_COMPARE(te->getIndexStart(), ==, 0);
  TEST_COMPARE(te->getIndexStop(), ==, 0);
  TEST_COMPARE(te->getIndexStride(), ==, 1);
  TEST_COMPARE(te->getNumEvents(), ==, 1);

  // Check base class defaults (functions not implemented in
  // TimeEventRangeIndex).
  TEST_COMPARE(te->isTime(1.0), ==, false);
  TEST_COMPARE(te->timeToNextEvent(1.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->timeOfNextEvent(1.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->eventInRange(1.0, 4.0), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, Construction)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>(-1, 10, 2, "TestName"));

  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test when everything is zero.
  TEST_COMPARE(te->getIndexStart(), ==, -1);
  TEST_COMPARE(te->getIndexStop(), ==, 10);
  TEST_COMPARE(te->getIndexStride(), ==, 2);
  TEST_COMPARE(te->getNumEvents(), ==, 6);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, Basic_Accesors)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>());

  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  // Reset start after stop.
  te->setIndexStart(1);
  TEST_COMPARE(te->getIndexStart(), ==, 1);
  TEST_COMPARE(te->getIndexStop(), ==, 1);
  TEST_COMPARE(te->getIndexStride(), ==, 1);
  TEST_COMPARE(te->getNumEvents(), ==, 1);

  // Reset stop.
  te->setIndexStop(5);
  TEST_COMPARE(te->getIndexStart(), ==, 1);
  TEST_COMPARE(te->getIndexStop(), ==, 5);
  TEST_COMPARE(te->getIndexStride(), ==, 1);
  TEST_COMPARE(te->getNumEvents(), ==, 5);

  // Reset stride.
  te->setIndexStride(2);
  TEST_COMPARE(te->getIndexStart(), ==, 1);
  TEST_COMPARE(te->getIndexStop(), ==, 5);
  TEST_COMPARE(te->getIndexStride(), ==, 2);
  TEST_COMPARE(te->getNumEvents(), ==, 3);

  // Set with index range.
  te->setIndexRange(-5, 5, 3);
  TEST_COMPARE(te->getIndexStart(), ==, -5);
  TEST_COMPARE(te->getIndexStop(), ==, 5);
  TEST_COMPARE(te->getIndexStride(), ==, 3);
  TEST_COMPARE(te->getNumEvents(), ==, 4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, Stride)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>());
  te->setIndexStart(1);
  te->setIndexStop(5);
  TEST_COMPARE(te->getIndexStart(), ==, 1);
  TEST_COMPARE(te->getIndexStop(), ==, 5);

  // Negative stride should be reset to stop_-start_.
  te->setIndexStride(-1);
  TEST_COMPARE(te->getIndexStride(), ==, 1);
  TEST_COMPARE(te->getNumEvents(), ==, 5);

  // Large stride should be reset to stop_-start_.
  te->setIndexStride(5);
  TEST_COMPARE(te->getIndexStride(), ==, 4);
  TEST_COMPARE(te->getNumEvents(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, isIndex)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>());
  te->setIndexRange(-5, 5, 3);

  // Test isIndex.
  TEST_COMPARE(te->isIndex(-6), ==, false);  // Around first event.
  TEST_COMPARE(te->isIndex(-5), ==, true);
  TEST_COMPARE(te->isIndex(-4), ==, false);

  TEST_COMPARE(te->isIndex(0), ==, false);  // Around mid event.
  TEST_COMPARE(te->isIndex(1), ==, true);
  TEST_COMPARE(te->isIndex(2), ==, false);

  TEST_COMPARE(te->isIndex(3), ==, false);  // Around last event.
  TEST_COMPARE(te->isIndex(4), ==, true);
  TEST_COMPARE(te->isIndex(5), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, indexToNextEvent)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>());
  te->setIndexRange(-5, 5, 3);

  // Test indexToNextEvent.
  TEST_COMPARE(te->indexToNextEvent(-9), ==, 4);  //   Around first event.
  TEST_COMPARE(te->indexToNextEvent(-5), ==, 3);
  TEST_COMPARE(te->indexToNextEvent(-4), ==, 2);

  TEST_COMPARE(te->indexToNextEvent(-1), ==, 2);  //   Around mid event.
  TEST_COMPARE(te->indexToNextEvent(1), ==, 3);
  TEST_COMPARE(te->indexToNextEvent(3), ==, 1);

  TEST_COMPARE(te->indexToNextEvent(2), ==, 2);  //   Around last event.
  TEST_COMPARE(te->indexToNextEvent(4), ==, te->getDefaultIndex() - 4);
  TEST_COMPARE(te->indexToNextEvent(8), ==, te->getDefaultIndex() - 8);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, indexOfNextEvent)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>());
  te->setIndexRange(-5, 5, 3);

  // Test indexOfNextEvent.
  TEST_COMPARE(te->indexOfNextEvent(-9), ==, -5);  //   Around first event.
  TEST_COMPARE(te->indexOfNextEvent(-5), ==, -2);
  TEST_COMPARE(te->indexOfNextEvent(-4), ==, -2);

  TEST_COMPARE(te->indexOfNextEvent(-1), ==, 1);  //   Around mid event.
  TEST_COMPARE(te->indexOfNextEvent(1), ==, 4);
  TEST_COMPARE(te->indexOfNextEvent(3), ==, 4);

  TEST_COMPARE(te->indexOfNextEvent(2), ==, 4);  //   Around last event.
  TEST_COMPARE(te->indexOfNextEvent(4), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(8), ==, te->getDefaultIndex());
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, eventInRangeIndex)
{
  auto te = rcp(new Tempus::TimeEventRangeIndex<double>());
  te->setIndexRange(-5, 5, 3);

  // Test eventInRangeIndex.
  //   Right end.
  TEST_COMPARE(te->eventInRangeIndex(-9, -6), ==,
               false);  // Around first event.
  TEST_COMPARE(te->eventInRangeIndex(-9, -5), ==, true);
  TEST_COMPARE(te->eventInRangeIndex(-9, -4), ==, true);

  TEST_COMPARE(te->eventInRangeIndex(-1, 0), ==, false);  // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-1, 1), ==, true);
  TEST_COMPARE(te->eventInRangeIndex(-1, 2), ==, true);

  TEST_COMPARE(te->eventInRangeIndex(2, 3), ==, false);  // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(2, 4), ==, true);
  TEST_COMPARE(te->eventInRangeIndex(2, 5), ==, true);

  //   Left end.
  TEST_COMPARE(te->eventInRangeIndex(-6.0, -3), ==,
               true);  // Around first event.
  TEST_COMPARE(te->eventInRangeIndex(-5.0, -3), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(-4.0, -3), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(-3, 0), ==, true);  // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-2, 0), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(-1, 0), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(3, 8), ==, true);  // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(4, 8), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(5, 8), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, getValidParameters)
{
  auto teri = rcp(new Tempus::TimeEventRangeIndex<double>());

  auto pl = teri->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Type"), ==, "Range Index");
  TEST_COMPARE(pl->get<std::string>("Name"), ==,
               "TimeEventRangeIndex (0; 0; 1)");
  TEST_COMPARE(pl->get<int>("Start Index"), ==, 0);
  TEST_COMPARE(pl->get<int>("Stop Index"), ==, 0);
  TEST_COMPARE(pl->get<int>("Stride Index"), ==, 1);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventRangeIndex, createTimeEventRange)
{
  // Construct parameterList similar to getValidParameters().
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Event Range Index");

  pl->set("Name", "Unit Test Time Event Range Index");
  pl->set("Type", "Range Index");
  pl->set("Start Index", -1);
  pl->set("Stop Index", 11);
  pl->set("Stride Index", 2);

  // Construct TimeEventRangeIndex from ParameterList.
  auto teri = Tempus::createTimeEventRangeIndex<double>(pl);

  teri->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(teri->getName(), ==, "Unit Test Time Event Range Index");
  TEST_COMPARE(teri->getType(), ==, "Range Index");
  TEST_COMPARE(teri->getIndexStart(), ==, -1);
  TEST_COMPARE(teri->getIndexStop(), ==, 11);
  TEST_COMPARE(teri->getIndexStride(), ==, 2);
  TEST_COMPARE(teri->getNumEvents(), ==, 7);
}

}  // namespace Tempus_Unit_Test
