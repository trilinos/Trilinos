//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_TimeEventComposite.hpp"
#include "Tempus_TimeEventRange.hpp"
#include "Tempus_TimeEventRangeIndex.hpp"
#include "Tempus_TimeEventList.hpp"
#include "Tempus_TimeEventListIndex.hpp"

static double PI = M_PI;

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

// TimeEventRanges for testing.
// ------------------------------------------------------------
Teuchos::RCP<Tempus::TimeEventRange<double> > getTestRange1()
{
  auto teRange1 = rcp(new Tempus::TimeEventRange<double>(
      0.0, PI, 1.0, "teRange1", true, 1.0e-14));
  return teRange1;
}

Teuchos::RCP<Tempus::TimeEventRange<double> > getTestRange2()
{
  auto teRange2 = rcp(new Tempus::TimeEventRange<double>(
      -PI / 2.0, PI / 2.0, PI / 4.0, "teRange2", true, 1.0e-14));
  return teRange2;
}

Teuchos::RCP<Tempus::TimeEventRange<double> > getTestRange3()
{
  auto teRange3 = rcp(new Tempus::TimeEventRange<double>(
      4.0, 10.0, 4.0, "teRange3", true, 1.0e-14));
  return teRange3;
}

// TimeEventLists for testing.
// ------------------------------------------------------------
Teuchos::RCP<Tempus::TimeEventList<double> > getTestList1()
{
  std::vector<double> testList1;
  testList1.push_back(-1.0);
  testList1.push_back(0.0);
  testList1.push_back(5.0);
  testList1.push_back(2.0);
  testList1.push_back(PI);

  auto teList1 = rcp(
      new Tempus::TimeEventList<double>(testList1, "teList1", true, 1.0e-14));
  return teList1;
}

Teuchos::RCP<Tempus::TimeEventList<double> > getTestList2()
{
  std::vector<double> testList2;
  testList2.push_back(-0.5);
  testList2.push_back(1.25);
  testList2.push_back(4.95);
  testList2.push_back(12.34);

  auto teList2 = rcp(
      new Tempus::TimeEventList<double>(testList2, "teList2", true, 1.0e-14));
  return teList2;
}

Teuchos::RCP<Tempus::TimeEventList<double> > getTestList3()
{
  std::vector<double> testList3;
  testList3.push_back(-5.0);
  testList3.push_back(-PI);

  auto teList3 = rcp(
      new Tempus::TimeEventList<double>(testList3, "teList3", true, 1.0e-14));
  return teList3;
}

// TimeEventRangeIndices for testing.
// ------------------------------------------------------------
Teuchos::RCP<Tempus::TimeEventRangeIndex<double> > getTestRangeIndex1()
{
  auto teRangeIndex1 =
      rcp(new Tempus::TimeEventRangeIndex<double>(-1, 10, 3, "teRangeIndex1"));
  return teRangeIndex1;
}

Teuchos::RCP<Tempus::TimeEventRangeIndex<double> > getTestRangeIndex2()
{
  auto teRangeIndex2 =
      rcp(new Tempus::TimeEventRangeIndex<double>(-5, 8, 4, "teRangeIndex2"));
  return teRangeIndex2;
}

Teuchos::RCP<Tempus::TimeEventRangeIndex<double> > getTestRangeIndex3()
{
  auto teRangeIndex3 =
      rcp(new Tempus::TimeEventRangeIndex<double>(10, 17, 5, "teRangeIndex3"));
  return teRangeIndex3;
}

// TimeEventRangeIndices for testing.
// ------------------------------------------------------------
Teuchos::RCP<Tempus::TimeEventListIndex<double> > getTestListIndex1()
{
  std::vector<int> testListIndex1;
  testListIndex1.push_back(-2);
  testListIndex1.push_back(0);
  testListIndex1.push_back(7);
  testListIndex1.push_back(3);
  testListIndex1.push_back(-5);

  auto teListIndex1 = rcp(
      new Tempus::TimeEventListIndex<double>(testListIndex1, "teListIndex1"));
  return teListIndex1;
}

Teuchos::RCP<Tempus::TimeEventListIndex<double> > getTestListIndex2()
{
  std::vector<int> testListIndex2;
  testListIndex2.push_back(2);

  auto teListIndex2 = rcp(
      new Tempus::TimeEventListIndex<double>(testListIndex2, "teListIndex2"));
  return teListIndex2;
}

Teuchos::RCP<Tempus::TimeEventListIndex<double> > getTestListIndex3()
{
  std::vector<int> testListIndex3;
  testListIndex3.push_back(14);
  testListIndex3.push_back(9);

  auto teListIndex3 = rcp(
      new Tempus::TimeEventListIndex<double>(testListIndex3, "teListIndex3"));

  return teListIndex3;
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, Default_Construction)
{
  auto tec = rcp(new Tempus::TimeEventComposite<double>());

  TEST_COMPARE(tec->getType(), ==, "Composite");
  TEST_COMPARE(tec->getName(), ==, "TimeEventComposite");
  TEST_COMPARE(tec->isTime(0.0), ==, false);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(0.0), tec->getDefaultTime(),
                         1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(0.0), tec->getDefaultTime(),
                         1.0e-14);
  TEST_COMPARE(tec->eventInRange(0.0, 1.0), ==, false);

  TEST_COMPARE(tec->isIndex(0), ==, false);
  TEST_COMPARE(tec->indexToNextEvent(0), ==, tec->getDefaultIndex());
  TEST_COMPARE(tec->indexOfNextEvent(0), ==, tec->getDefaultIndex());
  TEST_COMPARE(tec->eventInRange(0, 10), ==, false);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 0);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, Full_Construction)
{
  std::vector<Teuchos::RCP<Tempus::TimeEventBase<double> > > timeEvents;
  timeEvents.push_back(getTestRange1());
  timeEvents.push_back(getTestList1());
  timeEvents.push_back(getTestRangeIndex1());
  timeEvents.push_back(getTestListIndex1());

  auto tec = rcp(new Tempus::TimeEventComposite<double>(
      timeEvents, "Test TimeEventComposite"));
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->getType(), ==, "Composite");
  TEST_COMPARE(tec->getName(), ==, "Test TimeEventComposite");

  TEST_COMPARE(tec->isTime(3.0), ==, true);
  TEST_COMPARE(tec->isTime(2.0), ==, true);
  TEST_COMPARE(tec->isIndex(2), ==, true);
  TEST_COMPARE(tec->isIndex(3), ==, true);

  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(-2.5), 1.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(0.5), 0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(4.5), 0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(7.5), tec->getDefaultTime(),
                         1.0e-14);

  TEST_COMPARE(tec->indexToNextEvent(-6), ==, 1);
  TEST_COMPARE(tec->indexToNextEvent(1), ==, 1);
  TEST_COMPARE(tec->indexToNextEvent(7), ==, 1);
  TEST_COMPARE(tec->indexToNextEvent(9), ==, tec->getDefaultIndex() - 9);

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-PI), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-0.5), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(2.5), 3.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(7.5), tec->getDefaultTime(),
                         1.0e-14);

  TEST_COMPARE(tec->indexOfNextEvent(-6), ==, -5);
  TEST_COMPARE(tec->indexOfNextEvent(1), ==, 2);
  TEST_COMPARE(tec->indexOfNextEvent(7), ==, 8);
  TEST_COMPARE(tec->indexOfNextEvent(9), ==, tec->getDefaultIndex());

  TEST_COMPARE(tec->eventInRange(-5.0, -2.0), ==, false);
  TEST_COMPARE(tec->eventInRange(-2.0, -0.5), ==, true);
  TEST_COMPARE(tec->eventInRange(1.2, 1.8), ==, false);
  TEST_COMPARE(tec->eventInRange(3.1, 4.0), ==, true);
  TEST_COMPARE(tec->eventInRange(4.5, 6.0), ==, true);

  TEST_COMPARE(tec->eventInRangeIndex(-8, -6), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(1, 1), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(5, 7), ==, true);
  TEST_COMPARE(tec->eventInRangeIndex(8, 10), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(12, 14), ==, false);

  timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, setTimeEvents)
{
  std::vector<Teuchos::RCP<Tempus::TimeEventBase<double> > > timeEvents;
  timeEvents.push_back(getTestRange1());
  timeEvents.push_back(getTestList1());
  timeEvents.push_back(getTestRangeIndex1());
  timeEvents.push_back(getTestListIndex1());

  auto tec = rcp(new Tempus::TimeEventComposite<double>());

  tec->setTimeEvents(timeEvents);

  auto timeEventsSet = tec->getTimeEvents();
  TEST_COMPARE(timeEventsSet.size(), ==, 4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, One_TimeEventRange)
{
  auto tec      = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange1 = getTestRange1();

  tec->add(teRange1);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isTime(0.0), ==, true);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(0.0), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(1.1), 2.0, 1.0e-14);
  TEST_COMPARE(tec->eventInRange(0.0, 1.0), ==, true);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoOverlapping_TimeEventRanges)
{
  auto tec      = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange1 = getTestRange1();
  auto teRange2 = getTestRange2();

  tec->add(teRange1);
  tec->add(teRange2);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isTime(0.0), ==, true);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(0.1), PI / 4.0 - 0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(1.1), PI / 2.0, 1.0e-14);
  TEST_COMPARE(tec->eventInRange(0.0, 1.0), ==, true);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoSeparated_TimeEventRanges)
{
  auto tec      = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange3 = getTestRange3();
  auto teRange2 = getTestRange2();

  tec->add(teRange3);
  tec->add(teRange2);
  tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isTime(4.0), ==, true);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(1.1), PI / 2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(3.0), 4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(5.0), 8.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(10.0), tec->getDefaultTime(),
                         1.0e-14);
  TEST_COMPARE(tec->eventInRange(-3.0, -2.0), ==, false);
  TEST_COMPARE(tec->eventInRange(-1.0, 0.0), ==, true);
  TEST_COMPARE(tec->eventInRange(1.0, 2.0), ==, true);
  TEST_COMPARE(tec->eventInRange(2.0, 3.0), ==, false);
  TEST_COMPARE(tec->eventInRange(5.0, 7.0), ==, false);
  TEST_COMPARE(tec->eventInRange(7.0, 9.0), ==, true);
  TEST_COMPARE(tec->eventInRange(9.0, 11.0), ==, false);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, One_TimeEventList)
{
  auto tec     = rcp(new Tempus::TimeEventComposite<double>());
  auto teList1 = getTestList1();

  tec->add(teList1);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isTime(2.0), ==, true);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(3.5), 1.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-2.0), -1.0, 1.0e-14);
  TEST_COMPARE(tec->eventInRange(4.99, 10.0), ==, true);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoOverlapping_TimeEventLists)
{
  auto tec     = rcp(new Tempus::TimeEventComposite<double>());
  auto teList1 = getTestList1();
  auto teList2 = getTestList2();

  tec->add(teList1);
  tec->add(teList2);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isTime(1.25), ==, true);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(4.0), 0.95, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(6.5), 12.34, 1.0e-14);
  TEST_COMPARE(tec->eventInRange(0.1, 1.0), ==, false);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoSeparated_TimeEventLists)
{
  auto tec     = rcp(new Tempus::TimeEventComposite<double>());
  auto teList3 = getTestList3();
  auto teList2 = getTestList2();

  tec->add(teList3);
  tec->add(teList2);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isTime(4.0), ==, false);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-8.9), -5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-0.3), 1.25, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-4.0), -PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(20.0), tec->getDefaultTime(),
                         1.0e-14);
  TEST_COMPARE(tec->eventInRange(-6.0, -4.0), ==, true);
  TEST_COMPARE(tec->eventInRange(-3.0, 0.0), ==, true);
  TEST_COMPARE(tec->eventInRange(2.0, 3.0), ==, false);
  TEST_COMPARE(tec->eventInRange(4.9, 5.1), ==, true);
  TEST_COMPARE(tec->eventInRange(12.0, 12.4), ==, true);
  TEST_COMPARE(tec->eventInRange(14.0, 15.0), ==, false);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, One_TimeEventRangeIndex)
{
  auto tec           = rcp(new Tempus::TimeEventComposite<double>());
  auto teRangeIndex1 = getTestRangeIndex1();

  tec->add(teRangeIndex1);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isIndex(5), ==, true);
  TEST_COMPARE(tec->indexToNextEvent(3), ==, 2);
  TEST_COMPARE(tec->indexOfNextEvent(3), ==, 5);
  TEST_COMPARE(tec->eventInRangeIndex(3, 9), ==, true);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoOverlapping_TimeEventRangeIndex)
{
  auto tec           = rcp(new Tempus::TimeEventComposite<double>());
  auto teRangeIndex1 = getTestRangeIndex1();
  auto teRangeIndex2 = getTestRangeIndex2();

  tec->add(teRangeIndex1);
  tec->add(teRangeIndex2);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isIndex(-1), ==, true);
  TEST_COMPARE(tec->indexToNextEvent(-2), ==, 1);
  TEST_COMPARE(tec->indexOfNextEvent(2), ==, 3);
  TEST_COMPARE(tec->eventInRangeIndex(0, 1), ==, false);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoSeparate_TimeEventRangeIndex)
{
  auto tec           = rcp(new Tempus::TimeEventComposite<double>());
  auto teRangeIndex3 = getTestRangeIndex3();
  auto teRangeIndex2 = getTestRangeIndex2();

  tec->add(teRangeIndex3);
  tec->add(teRangeIndex2);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isIndex(15), ==, true);
  TEST_COMPARE(tec->indexOfNextEvent(9), ==, 10);
  TEST_COMPARE(tec->indexOfNextEvent(-6), ==, -5);
  TEST_COMPARE(tec->indexOfNextEvent(6), ==, 7);
  TEST_COMPARE(tec->indexOfNextEvent(16), ==, tec->getDefaultIndex());
  TEST_COMPARE(tec->eventInRangeIndex(-3, -2), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(-2, 0), ==, true);
  TEST_COMPARE(tec->eventInRangeIndex(1, 2), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(6, 7), ==, true);
  TEST_COMPARE(tec->eventInRangeIndex(7, 8), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(10, 13), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(14, 20), ==, true);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, One_TimeEventListIndex)
{
  auto tec          = rcp(new Tempus::TimeEventComposite<double>());
  auto teListIndex1 = getTestListIndex1();

  tec->add(teListIndex1);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isIndex(3), ==, true);
  TEST_COMPARE(tec->indexToNextEvent(1), ==, 2);
  TEST_COMPARE(tec->indexOfNextEvent(4), ==, 7);
  TEST_COMPARE(tec->eventInRangeIndex(1, 3), ==, true);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 1);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoOverlapping_TimeEventListIndex)
{
  auto tec          = rcp(new Tempus::TimeEventComposite<double>());
  auto teListIndex1 = getTestListIndex1();
  auto teListIndex2 = getTestListIndex2();

  tec->add(teListIndex1);
  tec->add(teListIndex2);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isIndex(2), ==, true);
  TEST_COMPARE(tec->indexToNextEvent(0), ==, 2);
  TEST_COMPARE(tec->indexOfNextEvent(1), ==, 2);
  TEST_COMPARE(tec->eventInRangeIndex(-1, 3), ==, true);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoSeparate_TimeEventListIndex)
{
  auto tec          = rcp(new Tempus::TimeEventComposite<double>());
  auto teListIndex3 = getTestListIndex3();
  auto teListIndex2 = getTestListIndex2();

  tec->add(teListIndex3);
  tec->add(teListIndex2);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isIndex(14), ==, true);
  TEST_COMPARE(tec->indexOfNextEvent(2), ==, 9);
  TEST_COMPARE(tec->indexOfNextEvent(5), ==, 9);
  TEST_COMPARE(tec->indexOfNextEvent(19), ==, tec->getDefaultIndex());
  TEST_COMPARE(tec->eventInRangeIndex(0, 1), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(3, 10), ==, true);
  TEST_COMPARE(tec->eventInRangeIndex(15, 20), ==, false);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 2);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, OneOfEach_TimeEvent)
{
  auto tec           = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange1      = getTestRange1();
  auto teList1       = getTestList1();
  auto teRangeIndex1 = getTestRangeIndex1();
  auto teListIndex1  = getTestListIndex1();

  tec->add(teRange1);
  tec->add(teList1);
  tec->add(teRangeIndex1);
  tec->add(teListIndex1);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->isTime(3.0), ==, true);
  TEST_COMPARE(tec->isTime(2.0), ==, true);
  TEST_COMPARE(tec->isIndex(2), ==, true);
  TEST_COMPARE(tec->isIndex(3), ==, true);

  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(-2.5), 1.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(0.5), 0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(4.5), 0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(7.5), tec->getDefaultTime(),
                         1.0e-14);

  TEST_COMPARE(tec->indexToNextEvent(-6), ==, 1);
  TEST_COMPARE(tec->indexToNextEvent(1), ==, 1);
  TEST_COMPARE(tec->indexToNextEvent(7), ==, 1);
  TEST_COMPARE(tec->indexToNextEvent(9), ==, tec->getDefaultIndex() - 9);

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-PI), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-0.5), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(2.5), 3.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(7.5), tec->getDefaultTime(),
                         1.0e-14);

  TEST_COMPARE(tec->indexOfNextEvent(-6), ==, -5);
  TEST_COMPARE(tec->indexOfNextEvent(1), ==, 2);
  TEST_COMPARE(tec->indexOfNextEvent(7), ==, 8);
  TEST_COMPARE(tec->indexOfNextEvent(9), ==, tec->getDefaultIndex());

  TEST_COMPARE(tec->eventInRange(-5.0, -2.0), ==, false);
  TEST_COMPARE(tec->eventInRange(-2.0, -0.5), ==, true);
  TEST_COMPARE(tec->eventInRange(1.2, 1.8), ==, false);
  TEST_COMPARE(tec->eventInRange(3.1, 4.0), ==, true);
  TEST_COMPARE(tec->eventInRange(4.5, 6.0), ==, true);

  TEST_COMPARE(tec->eventInRangeIndex(-8, -6), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(1, 1), ==, false);
  TEST_COMPARE(tec->eventInRangeIndex(5, 7), ==, true);
  TEST_COMPARE(tec->eventInRangeIndex(7, 10), ==, true);
  TEST_COMPARE(tec->eventInRangeIndex(12, 14), ==, false);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, OneOfEach_Plus_TimeEvent)
{
  auto tec           = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange1      = getTestRange1();
  auto teList1       = getTestList1();
  auto teRangeIndex1 = getTestRangeIndex1();
  auto teListIndex1  = getTestListIndex1();

  tec->add(teRange1);
  tec->add(teList1);
  tec->add(teRangeIndex1);
  tec->add(teListIndex1);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  // Constraining TimeEvent(s)
  std::vector<Teuchos::RCP<Tempus::TimeEventBase<double> > > teCons;

  TEST_COMPARE(tec->isTime(3.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRange1");

  TEST_COMPARE(tec->isTime(-1.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teList1");

  TEST_COMPARE(tec->isIndex(2, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRangeIndex1");

  TEST_COMPARE(tec->isIndex(-2, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teListIndex1");

  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(-2.5, teCons), 1.5, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teList1");

  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(0.5, teCons), 0.5, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRange1");

  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(4.5, teCons), 0.5, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teList1");

  TEST_FLOATING_EQUALITY(tec->timeToNextEvent(7.5, teCons),
                         tec->getDefaultTime(), 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 0);

  TEST_COMPARE(tec->indexToNextEvent(-6, teCons), ==, 1);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teListIndex1");

  TEST_COMPARE(tec->indexToNextEvent(1, teCons), ==, 1);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRangeIndex1");

  TEST_COMPARE(tec->indexToNextEvent(7, teCons), ==, 1);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRangeIndex1");

  TEST_COMPARE(tec->indexToNextEvent(9, teCons), ==,
               tec->getDefaultIndex() - 9);
  TEST_COMPARE(teCons.size(), ==, 0);

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-PI, teCons), -1.0, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teList1");

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-0.5, teCons), 0.0, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 2);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRange1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teList1");

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(2.5, teCons), 3.0, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRange1");

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(7.5, teCons),
                         tec->getDefaultTime(), 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 0);

  TEST_COMPARE(tec->indexOfNextEvent(-6, teCons), ==, -5);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teListIndex1");

  TEST_COMPARE(tec->indexOfNextEvent(1, teCons), ==, 2);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRangeIndex1");

  TEST_COMPARE(tec->indexOfNextEvent(7, teCons), ==, 8);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(teCons[0]->getName(), ==, "teRangeIndex1");

  TEST_COMPARE(tec->indexOfNextEvent(9, teCons), ==, tec->getDefaultIndex());
  TEST_COMPARE(teCons.size(), ==, 0);

  TEST_COMPARE(tec->eventInRange(-5.0, -2.0, teCons), ==, false);
  TEST_COMPARE(teCons.size(), ==, 0);
  TEST_COMPARE(tec->eventInRange(-2.0, -0.5, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(tec->eventInRange(1.2, 1.8, teCons), ==, false);
  TEST_COMPARE(teCons.size(), ==, 0);
  TEST_COMPARE(tec->eventInRange(3.1, 4.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(tec->eventInRange(4.5, 6.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 1);

  TEST_COMPARE(tec->eventInRangeIndex(-8, -6, teCons), ==, false);
  TEST_COMPARE(teCons.size(), ==, 0);
  TEST_COMPARE(tec->eventInRangeIndex(1, 1, teCons), ==, false);
  TEST_COMPARE(teCons.size(), ==, 0);
  TEST_COMPARE(tec->eventInRangeIndex(5, 7, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 1);
  TEST_COMPARE(tec->eventInRangeIndex(8, 10, teCons), ==, false);
  TEST_COMPARE(teCons.size(), ==, 0);
  TEST_COMPARE(tec->eventInRangeIndex(12, 14, teCons), ==, false);
  TEST_COMPARE(teCons.size(), ==, 0);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, Multiple_Simultaneous_Events)
{
  // Define Events
  auto ter1 = rcp(new Tempus::TimeEventRange<double>(0.0, 20.0, 5.0, "ter1",
                                                     true, 1.0e-14));
  auto ter2 = rcp(new Tempus::TimeEventRange<double>(0.0, 20.0, 2.0, "ter2",
                                                     false, 1.0e-14));

  std::vector<double> testList1;
  testList1.push_back(0.0);
  testList1.push_back(4.0);
  testList1.push_back(5.0);
  testList1.push_back(9.0);
  testList1.push_back(20.0);
  auto tel1 =
      rcp(new Tempus::TimeEventList<double>(testList1, "tel1", true, 1.0e-14));

  std::vector<double> testList2;
  testList2.push_back(0.0);
  testList2.push_back(3.0);
  testList2.push_back(5.0);
  testList2.push_back(13.0);
  testList2.push_back(20.0);
  auto tel2 =
      rcp(new Tempus::TimeEventList<double>(testList2, "tel2", false, 1.0e-14));

  auto teri1 =
      rcp(new Tempus::TimeEventRangeIndex<double>(0, 200, 50, "teri1"));
  auto teri2 =
      rcp(new Tempus::TimeEventRangeIndex<double>(0, 200, 20, "teri2"));

  std::vector<int> testListIndex1;
  testListIndex1.push_back(0);
  testListIndex1.push_back(40);
  testListIndex1.push_back(50);
  testListIndex1.push_back(90);
  testListIndex1.push_back(200);
  auto teli1 =
      rcp(new Tempus::TimeEventListIndex<double>(testListIndex1, "teli1"));

  std::vector<int> testListIndex2;
  testListIndex2.push_back(0);
  testListIndex2.push_back(30);
  testListIndex2.push_back(50);
  testListIndex2.push_back(130);
  testListIndex2.push_back(200);
  auto teli2 =
      rcp(new Tempus::TimeEventListIndex<double>(testListIndex2, "teli2"));

  auto tec = rcp(new Tempus::TimeEventComposite<double>());
  tec->add(ter1);
  tec->add(ter2);
  tec->add(tel1);
  tec->add(tel2);
  tec->add(teri1);
  tec->add(teri2);
  tec->add(teli1);
  tec->add(teli2);

  // tec->describe(out, Teuchos::VERB_EXTREME);

  // Constraining TimeEvent(s)
  std::vector<Teuchos::RCP<Tempus::TimeEventBase<double> > > teCons;

  TEST_COMPARE(tec->isTime(0.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "ter2");
  TEST_COMPARE(teCons[2]->getName(), ==, "tel1");
  TEST_COMPARE(teCons[3]->getName(), ==, "tel2");

  TEST_COMPARE(tec->isTime(5.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 3);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "tel1");
  TEST_COMPARE(teCons[2]->getName(), ==, "tel2");

  TEST_COMPARE(tec->isTime(10.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 2);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "ter2");

  TEST_COMPARE(tec->isTime(20.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "ter2");
  TEST_COMPARE(teCons[2]->getName(), ==, "tel1");
  TEST_COMPARE(teCons[3]->getName(), ==, "tel2");

  TEST_COMPARE(tec->isIndex(0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teri2");
  TEST_COMPARE(teCons[2]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[3]->getName(), ==, "teli2");

  TEST_COMPARE(tec->isIndex(50, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 3);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[2]->getName(), ==, "teli2");

  TEST_COMPARE(tec->isIndex(100, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 2);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teri2");

  TEST_COMPARE(tec->isIndex(200, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teri2");
  TEST_COMPARE(teCons[2]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[3]->getName(), ==, "teli2");

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(-1.0, teCons), 0.0, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "ter2");
  TEST_COMPARE(teCons[2]->getName(), ==, "tel1");
  TEST_COMPARE(teCons[3]->getName(), ==, "tel2");

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(4.0, teCons), 5.0, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 3);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "tel1");
  TEST_COMPARE(teCons[2]->getName(), ==, "tel2");

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(9.0, teCons), 10.0, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 2);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "ter2");

  TEST_FLOATING_EQUALITY(tec->timeOfNextEvent(19.0, teCons), 20.0, 1.0e-14);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "ter2");
  TEST_COMPARE(teCons[2]->getName(), ==, "tel1");
  TEST_COMPARE(teCons[3]->getName(), ==, "tel2");

  TEST_COMPARE(tec->indexOfNextEvent(-1, teCons), ==, 0);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teri2");
  TEST_COMPARE(teCons[2]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[3]->getName(), ==, "teli2");

  TEST_COMPARE(tec->indexOfNextEvent(40, teCons), ==, 50);
  TEST_COMPARE(teCons.size(), ==, 3);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[2]->getName(), ==, "teli2");

  TEST_COMPARE(tec->indexOfNextEvent(90, teCons), ==, 100);
  TEST_COMPARE(teCons.size(), ==, 2);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teri2");

  TEST_COMPARE(tec->indexOfNextEvent(190, teCons), ==, 200);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teri2");
  TEST_COMPARE(teCons[2]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[3]->getName(), ==, "teli2");

  // Sorted order is still input order.
  TEST_COMPARE(tec->eventInRange(-1.0, 21.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_COMPARE(teCons[1]->getName(), ==, "ter2");
  TEST_COMPARE(teCons[2]->getName(), ==, "tel1");
  TEST_COMPARE(teCons[3]->getName(), ==, "tel2");

  // Sorted order is based on "time of next event".
  TEST_COMPARE(tec->eventInRange(0.0, 21.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter2");
  TEST_FLOATING_EQUALITY(teCons[0]->timeOfNextEvent(0.0), 2.0, 1.0e-14);
  TEST_COMPARE(teCons[1]->getName(), ==, "tel2");
  TEST_FLOATING_EQUALITY(teCons[1]->timeOfNextEvent(0.0), 3.0, 1.0e-14);
  TEST_COMPARE(teCons[2]->getName(), ==, "tel1");
  TEST_FLOATING_EQUALITY(teCons[2]->timeOfNextEvent(0.0), 4.0, 1.0e-14);
  TEST_COMPARE(teCons[3]->getName(), ==, "ter1");
  TEST_FLOATING_EQUALITY(teCons[3]->timeOfNextEvent(0.0), 5.0, 1.0e-14);

  TEST_COMPARE(tec->eventInRange(7.0, 21.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter2");
  TEST_FLOATING_EQUALITY(teCons[0]->timeOfNextEvent(7.0), 8.0, 1.0e-14);
  TEST_COMPARE(teCons[1]->getName(), ==, "tel1");
  TEST_FLOATING_EQUALITY(teCons[1]->timeOfNextEvent(7.0), 9.0, 1.0e-14);
  TEST_COMPARE(teCons[2]->getName(), ==, "ter1");
  TEST_FLOATING_EQUALITY(teCons[2]->timeOfNextEvent(7.0), 10.0, 1.0e-14);
  TEST_COMPARE(teCons[3]->getName(), ==, "tel2");
  TEST_FLOATING_EQUALITY(teCons[3]->timeOfNextEvent(7.0), 13.0, 1.0e-14);

  TEST_COMPARE(tec->eventInRange(19.0, 21.0, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "ter1");
  TEST_FLOATING_EQUALITY(teCons[0]->timeOfNextEvent(19.0), 20.0, 1.0e-14);
  TEST_COMPARE(teCons[1]->getName(), ==, "ter2");
  TEST_FLOATING_EQUALITY(teCons[1]->timeOfNextEvent(19.0), 20.0, 1.0e-14);
  TEST_COMPARE(teCons[2]->getName(), ==, "tel1");
  TEST_FLOATING_EQUALITY(teCons[2]->timeOfNextEvent(19.0), 20.0, 1.0e-14);
  TEST_COMPARE(teCons[3]->getName(), ==, "tel2");
  TEST_FLOATING_EQUALITY(teCons[3]->timeOfNextEvent(19.0), 20.0, 1.0e-14);

  // Sorted order is still input order.
  TEST_COMPARE(tec->eventInRangeIndex(-10, 210, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[1]->getName(), ==, "teri2");
  TEST_COMPARE(teCons[2]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[3]->getName(), ==, "teli2");

  // Sorted order is based on "time of next event".
  TEST_COMPARE(tec->eventInRangeIndex(0, 210, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri2");
  TEST_COMPARE(teCons[0]->indexOfNextEvent(0), ==, 20);
  TEST_COMPARE(teCons[1]->getName(), ==, "teli2");
  TEST_COMPARE(teCons[1]->indexOfNextEvent(0), ==, 30);
  TEST_COMPARE(teCons[2]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[2]->indexOfNextEvent(0), ==, 40);
  TEST_COMPARE(teCons[3]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[3]->indexOfNextEvent(0), ==, 50);

  TEST_COMPARE(tec->eventInRangeIndex(70, 210, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri2");
  TEST_COMPARE(teCons[0]->indexOfNextEvent(70), ==, 80);
  TEST_COMPARE(teCons[1]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[1]->indexOfNextEvent(70), ==, 90);
  TEST_COMPARE(teCons[2]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[2]->indexOfNextEvent(70), ==, 100);
  TEST_COMPARE(teCons[3]->getName(), ==, "teli2");
  TEST_COMPARE(teCons[3]->indexOfNextEvent(70), ==, 130);

  TEST_COMPARE(tec->eventInRangeIndex(190, 210, teCons), ==, true);
  TEST_COMPARE(teCons.size(), ==, 4);
  TEST_COMPARE(teCons[0]->getName(), ==, "teri1");
  TEST_COMPARE(teCons[0]->indexOfNextEvent(190), ==, 200);
  TEST_COMPARE(teCons[1]->getName(), ==, "teri2");
  TEST_COMPARE(teCons[1]->indexOfNextEvent(190), ==, 200);
  TEST_COMPARE(teCons[2]->getName(), ==, "teli1");
  TEST_COMPARE(teCons[2]->indexOfNextEvent(190), ==, 200);
  TEST_COMPARE(teCons[3]->getName(), ==, "teli2");
  TEST_COMPARE(teCons[3]->indexOfNextEvent(190), ==, 200);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, remove)
{
  // Construct ParmeterList for testing.
  auto tec  = rcp(new Tempus::TimeEventComposite<double>());
  auto ter  = rcp(new Tempus::TimeEventRange<double>());
  auto teri = rcp(new Tempus::TimeEventRangeIndex<double>());
  auto tel  = rcp(new Tempus::TimeEventList<double>());
  auto teli = rcp(new Tempus::TimeEventListIndex<double>());

  ter->setName("Test Range");
  teri->setName("Test Range Index");
  tel->setName("Test List");
  teli->setName("Test List Index");

  tec->add(ter);
  tec->add(teri);
  tec->add(tel);
  tec->add(teli);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  tec->remove("Blah Blah");
  TEST_COMPARE(tec->getSize(), ==, 4);

  tec->remove("Test Range Index");
  TEST_COMPARE(tec->getSize(), ==, 3);
  auto names = tec->getTimeEventNames();
  TEST_COMPARE(names, ==, "Test Range, Test List, Test List Index");

  tec->remove("Test List Index");
  TEST_COMPARE(tec->getSize(), ==, 2);
  names = tec->getTimeEventNames();
  TEST_COMPARE(names, ==, "Test Range, Test List");

  tec->remove("Test Range");
  TEST_COMPARE(tec->getSize(), ==, 1);
  names = tec->getTimeEventNames();
  TEST_COMPARE(names, ==, "Test List");

  tec->remove("Test List");
  TEST_COMPARE(tec->getSize(), ==, 0);
  names = tec->getTimeEventNames();
  TEST_COMPARE(names, ==, "");
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, find)
{
  // Construct ParmeterList for testing.
  auto tec  = rcp(new Tempus::TimeEventComposite<double>());
  auto ter  = rcp(new Tempus::TimeEventRange<double>());
  auto teri = rcp(new Tempus::TimeEventRangeIndex<double>());
  auto tel  = rcp(new Tempus::TimeEventList<double>());
  auto teli = rcp(new Tempus::TimeEventListIndex<double>());

  ter->setName("Test Range");
  teri->setName("Test Range Index");
  tel->setName("Test List");
  teli->setName("Test List Index");

  tec->add(ter);
  tec->add(teri);
  tec->add(tel);
  tec->add(teli);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  auto teTest = tec->find("Test Range");
  TEST_COMPARE(teTest->getName(), ==, "Test Range");

  teTest = tec->find("Test Range Index");
  TEST_COMPARE(teTest->getName(), ==, "Test Range Index");

  teTest = tec->find("Test List");
  TEST_COMPARE(teTest->getName(), ==, "Test List");

  teTest = tec->find("Test List Index");
  TEST_COMPARE(teTest->getName(), ==, "Test List Index");

  teTest = tec->find("Blah Blah");
  TEST_COMPARE(teTest, ==, Teuchos::null);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, clear)
{
  // Construct ParmeterList for testing.
  auto tec  = rcp(new Tempus::TimeEventComposite<double>());
  auto ter  = rcp(new Tempus::TimeEventRange<double>());
  auto teri = rcp(new Tempus::TimeEventRangeIndex<double>());
  auto tel  = rcp(new Tempus::TimeEventList<double>());
  auto teli = rcp(new Tempus::TimeEventListIndex<double>());

  tec->add(ter);
  tec->add(teri);
  tec->add(tel);
  tec->add(teli);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tec->getSize(), ==, 4);
  tec->clear();
  TEST_COMPARE(tec->getSize(), ==, 0);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, getValidParameters)
{
  // Construct ParmeterList for testing.
  auto tec  = rcp(new Tempus::TimeEventComposite<double>());
  auto ter  = rcp(new Tempus::TimeEventRange<double>());
  auto teri = rcp(new Tempus::TimeEventRangeIndex<double>());
  auto tel  = rcp(new Tempus::TimeEventList<double>());
  auto teli = rcp(new Tempus::TimeEventListIndex<double>());

  ter->setName("Test Range");
  teri->setName("Test Range Index");
  tel->setName("Test List");
  teli->setName("Test List Index");

  tec->add(ter);
  tec->add(teri);
  tec->add(tel);
  tec->add(teli);
  // tec->describe(out, Teuchos::VERB_EXTREME);

  auto timeEvents = tec->getTimeEvents();
  TEST_COMPARE(timeEvents.size(), ==, 4);

  auto pl = tec->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Name"), ==, "TimeEventComposite");
  TEST_COMPARE(pl->get<std::string>("Type"), ==, "Composite");
  TEST_COMPARE(pl->get<std::string>("Time Events"), ==,
               "Test Range, Test Range Index, Test List, Test List Index");
  TEST_COMPARE(pl->isSublist("Test Range"), ==, true);
  TEST_COMPARE(pl->isSublist("Test Range Index"), ==, true);
  TEST_COMPARE(pl->isSublist("Test List"), ==, true);
  TEST_COMPARE(pl->isSublist("Test List Index"), ==, true);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(
        unusedParameters.str(), ==,
        "WARNING: Parameter \"Test Range\"    [unused] is unused\n"
        "WARNING: Parameter \"Test Range Index\"    [unused] is unused\n"
        "WARNING: Parameter \"Test List\"    [unused] is unused\n"
        "WARNING: Parameter \"Test List Index\"    [unused] is unused\n");
  }

  auto terPL = pl->sublist("Test Range");
  TEST_COMPARE(terPL.get<std::string>("Type"), ==, "Range");
  TEST_COMPARE(terPL.get<std::string>("Name"), ==, "Test Range");
  TEST_FLOATING_EQUALITY(terPL.get<double>("Start Time"), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(terPL.get<double>("Stop Time"), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(terPL.get<double>("Stride Time"), 0.0, 1.0e-14);
  TEST_COMPARE(terPL.get<int>("Number of Events"), ==, 1);
  TEST_FLOATING_EQUALITY(terPL.get<double>("Relative Tolerance"),
                         std::numeric_limits<double>::epsilon() * 100.0,
                         1.0e-14);

  TEST_COMPARE(terPL.get<bool>("Land On Exactly"), ==, true);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    terPL.unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }

  auto teriPL = pl->sublist("Test Range Index");
  TEST_COMPARE(teriPL.get<std::string>("Type"), ==, "Range Index");
  TEST_COMPARE(teriPL.get<std::string>("Name"), ==, "Test Range Index");
  TEST_COMPARE(teriPL.get<int>("Start Index"), ==, 0);
  TEST_COMPARE(teriPL.get<int>("Stop Index"), ==, 0);
  TEST_COMPARE(teriPL.get<int>("Stride Index"), ==, 1);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    teriPL.unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }

  auto telPL = pl->sublist("Test List");
  TEST_COMPARE(telPL.get<std::string>("Type"), ==, "List");
  TEST_COMPARE(telPL.get<std::string>("Name"), ==, "Test List");
  TEST_FLOATING_EQUALITY(telPL.get<double>("Relative Tolerance"),
                         std::numeric_limits<double>::epsilon() * 100.0,
                         1.0e-14);
  TEST_COMPARE(telPL.get<bool>("Land On Exactly"), ==, true);
  TEST_COMPARE(telPL.get<std::string>("Time List"), ==, "");

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    telPL.unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }

  auto teliPL = pl->sublist("Test List Index");
  TEST_COMPARE(teliPL.get<std::string>("Type"), ==, "List Index");
  TEST_COMPARE(teliPL.get<std::string>("Name"), ==, "Test List Index");
  TEST_COMPARE(teliPL.get<std::string>("Index List"), ==, "");

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    teliPL.unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, createTimeEventComposite)
{
  using Teuchos::ParameterList;
  using Teuchos::rcp_const_cast;

  {  // Construct from default ParameterList.
    auto tecTmp = rcp(new Tempus::TimeEventComposite<double>());
    auto pl     = rcp_const_cast<ParameterList>(tecTmp->getValidParameters());

    auto tec = Tempus::createTimeEventComposite<double>(pl);
    // tec->describe(out, Teuchos::VERB_EXTREME);

    TEST_COMPARE(tec->getName(), ==, "TimeEventComposite");
    TEST_COMPARE(tec->getType(), ==, "Composite");
    TEST_COMPARE(tec->getSize(), ==, 0);
    TEST_COMPARE(tec->getTimeEventNames(), ==, "");
  }

  {  // Construct with TimeEvents.
    auto tecTmp = rcp(new Tempus::TimeEventComposite<double>());
    auto ter    = rcp(new Tempus::TimeEventRange<double>());
    auto teri   = rcp(new Tempus::TimeEventRangeIndex<double>());
    auto tel    = rcp(new Tempus::TimeEventList<double>());
    auto teli   = rcp(new Tempus::TimeEventListIndex<double>());

    ter->setName("Test Range");
    teri->setName("Test Range Index");
    tel->setName("Test List");
    teli->setName("Test List Index");

    tecTmp->add(ter);
    tecTmp->add(teri);
    tecTmp->add(tel);
    tecTmp->add(teli);

    auto pl = rcp_const_cast<ParameterList>(tecTmp->getValidParameters());

    auto tec = Tempus::createTimeEventComposite<double>(pl);
    // tec->describe(out, Teuchos::VERB_EXTREME);

    TEST_COMPARE(tec->getName(), ==, "TimeEventComposite");
    TEST_COMPARE(tec->getType(), ==, "Composite");
    TEST_COMPARE(tec->getSize(), ==, 4);
    TEST_COMPARE(tec->getTimeEventNames(), ==,
                 "Test Range, Test Range Index, Test List, Test List Index");
  }

  {  // Construct with non-Tempus TimeEvent.
    auto tecTmp = rcp(new Tempus::TimeEventComposite<double>());
    auto pl     = rcp_const_cast<ParameterList>(tecTmp->getValidParameters());

    auto nonTempusTE = Teuchos::parameterList("Application Time Event");
    nonTempusTE->set<std::string>("Name", "Application Time Event");
    nonTempusTE->set<std::string>("Type", "Application Time Event Type");
    nonTempusTE->set<double>("Secret Sauce", 1.2345);
    pl->set("Application Time Event", *nonTempusTE);
    pl->set("Time Events", "Application Time Event");

    auto tec = Tempus::createTimeEventComposite<double>(pl);
    // tec->describe(out, Teuchos::VERB_EXTREME);

    TEST_COMPARE(tec->getName(), ==, "TimeEventComposite");
    TEST_COMPARE(tec->getType(), ==, "Composite");
    TEST_COMPARE(tec->getSize(), ==, 0);
    TEST_COMPARE(tec->getTimeEventNames(), ==, "");
  }
}

}  // namespace Tempus_Unit_Test
