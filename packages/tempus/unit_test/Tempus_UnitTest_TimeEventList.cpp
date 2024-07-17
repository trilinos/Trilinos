//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"
#include "Tempus_TimeEventList.hpp"

static double PI = M_PI;

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventList<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventList");

  TEST_COMPARE(te->getTimeList().size(), ==, 0);
  TEST_FLOATING_EQUALITY(
      te->getRelTol(), std::numeric_limits<double>::epsilon() * 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(
      te->getAbsTol(), std::numeric_limits<double>::epsilon() * 100.0, 1.0e-14);

  TEST_COMPARE(te->getLandOnExactly(), ==, true);

  // Check base class defaults (functions not implemented in TimeEventList).
  TEST_COMPARE(te->isIndex(1), ==, false);
  TEST_COMPARE(te->indexToNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(1, 4), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, Full_Construction)
{
  std::vector<double> testVector;
  testVector.push_back(-1.0);
  testVector.push_back(0.0);
  testVector.push_back(5.0);
  testVector.push_back(2.0);
  testVector.push_back(PI);

  auto te = rcp(
      new Tempus::TimeEventList<double>(testVector, "TestName", true, 1.0e-14));

  TEST_COMPARE(te->getName(), ==, "TestName");
  TEST_FLOATING_EQUALITY(te->getRelTol(), 1.0e-14, 1.0e-14);
  TEST_COMPARE(te->getLandOnExactly(), ==, true);

  auto testList = te->getTimeList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_FLOATING_EQUALITY(testList[0], -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[1], 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[2], 2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[3], PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[4], 5.0, 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, Basic_Accessors)
{
  auto te = rcp(new Tempus::TimeEventList<double>());

  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");
  te->setRelTol(0.1);
  TEST_FLOATING_EQUALITY(te->getRelTol(), 0.1, 1.0e-14);
  te->setRelTol(1.0e-14);
  te->setLandOnExactly(false);
  TEST_COMPARE(te->getLandOnExactly(), ==, false);

  // Test addTime.
  te->addTime(0.0);
  te->addTime(PI);
  te->addTime(-1.0);
  te->addTime(2.0);
  te->addTime(5.0);

  // Add times that should not be duplicated.
  te->addTime(0.0);
  te->addTime(PI);

  auto testList = te->getTimeList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_FLOATING_EQUALITY(testList[0], -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[1], 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[2], 2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[3], PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[4], 5.0, 1.0e-14);

  // Test that two events within relative tolerance are added or not.
  te->addTime(2.0 + 1.0e-14);
  TEST_COMPARE(te->getTimeList().size(), ==, 5);
  te->addTime(2.0 + 1.0e-13);
  TEST_COMPARE(te->getTimeList().size(), ==, 6);

  // Test setTimeList()
  te->clearTimeList();
  te->setTimeList(testList);
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_FLOATING_EQUALITY(testList[0], -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[1], 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[2], 2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[3], PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[4], 5.0, 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, isTime)
{
  auto te = rcp(new Tempus::TimeEventList<double>());
  te->addTime(0.0);
  te->addTime(PI);
  te->addTime(-1.0);
  te->addTime(2.0);
  te->addTime(5.0);
  te->setRelTol(1.0e-14);

  // Test isTime.
  //   Around first event.
  TEST_COMPARE(te->isTime(-10.0e-14), ==, false);  // Just outside tolerance.
  TEST_COMPARE(te->isTime(-0.1e-14), ==, true);    // Just inside tolerance.
  TEST_COMPARE(te->isTime(0.0), ==, true);         // Right on timeEvent.
  TEST_COMPARE(te->isTime(0.1e-14), ==, true);     // Just inside tolerance.
  TEST_COMPARE(te->isTime(10.0e-14), ==, false);   // Just outside tolerance.

  //   Around mid event.
  TEST_COMPARE(te->isTime(PI + -10.0e-14), ==,
               false);                                // Just outside tolerance.
  TEST_COMPARE(te->isTime(PI + -0.1e-14), ==, true);  // Just inside tolerance.
  TEST_COMPARE(te->isTime(PI + 0.0), ==, true);       // Right on timeEvent.
  TEST_COMPARE(te->isTime(PI + 0.1e-14), ==, true);   // Just inside tolerance.
  TEST_COMPARE(te->isTime(PI + 10.0e-14), ==,
               false);  // Just outside tolerance.

  //   Around last event.
  TEST_COMPARE(te->isTime(5.0 + -10.0e-14), ==,
               false);                                 // Just outside tolerance.
  TEST_COMPARE(te->isTime(5.0 + -0.1e-14), ==, true);  // Just inside tolerance.
  TEST_COMPARE(te->isTime(5.0 + 0.0), ==, true);       // Right on timeEvent.
  TEST_COMPARE(te->isTime(5.0 + 0.1e-14), ==, true);   // Just inside tolerance.
  TEST_COMPARE(te->isTime(5.0 + 10.0e-14), ==,
               false);  // Just outside tolerance.
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, timeToNextEvent)
{
  std::vector<double> testList;
  testList.push_back(0.0);
  testList.push_back(PI);
  testList.push_back(-1.0);
  testList.push_back(2.0);
  testList.push_back(5.0);

  auto te = rcp(
      new Tempus::TimeEventList<double>(testList, "testList", true, 1.0e-14));

  // Test timeToNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 + -10.0e-14), 1.0e-13,
                         1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 + -0.1e-14), 1.0 + 0.1e-14,
                         1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 + 0.0), 1.0 + 0.0, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 + 0.1e-14), 1.0 - 0.1e-14,
                         1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 + 10.0e-14), 1.0 - 1.0e-13,
                         1.0e-14);

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI + -10.0e-14), 1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI + -0.1e-14), 5.0 - PI + 0.1e-14,
                         1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI + 0.0), 5.0 - PI + 0.0,
                         1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI + 0.1e-14), 5.0 - PI - 0.1e-14,
                         1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI + 10.0e-14), 5.0 - PI - 1.0e-13,
                         1.0e-14);

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0 + -10.0e-14), 1.0e-13,
                         1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0 + -0.1e-14),
                         te->getDefaultTime(), 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0 + 0.0), te->getDefaultTime(),
                         1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0 + 0.1e-14),
                         te->getDefaultTime(), 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0 + 10.0e-14),
                         te->getDefaultTime(), 1.0e-02);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, timeOfNextEvent)
{
  std::vector<double> testList;
  testList.push_back(0.0);
  testList.push_back(PI);
  testList.push_back(-1.0);
  testList.push_back(2.0);
  testList.push_back(5.0);

  auto te = rcp(
      new Tempus::TimeEventList<double>(testList, "testList", true, 1.0e-14));

  // Test timeOfNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 + -10.0e-14), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 + -0.1e-14), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 + 0.0), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 + 0.1e-14), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 + 10.0e-14), 0.0, 1.0e-14);

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0 + -10.0e-14), 2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0 + -0.1e-14), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0 + 0.0), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0 + 0.1e-14), PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0 + 10.0e-14), PI, 1.0e-14);

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0 + -10.0e-14), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0 + -0.1e-14),
                         te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0 + 0.0), te->getDefaultTime(),
                         1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0 + 0.1e-14),
                         te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0 + 10.0e-14),
                         te->getDefaultTime(), 1.0e-14);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, eventInRange)
{
  std::vector<double> testList;
  testList.push_back(0.0);
  testList.push_back(PI);
  testList.push_back(-1.0);
  testList.push_back(2.0);
  testList.push_back(5.0);

  auto te = rcp(
      new Tempus::TimeEventList<double>(testList, "testList", true, 1.0e-14));

  // Test eventInRange.
  //   Right end.
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 + -10.0e-14), ==,
               false);  // Around first event.
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 + -0.1e-14), ==, true);
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 + 0.0), ==, true);
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 + 0.1e-14), ==, true);
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 + 10.0e-14), ==, true);

  TEST_COMPARE(te->eventInRange(3.0, PI + -10.0e-14), ==,
               false);  // Around mid event.
  TEST_COMPARE(te->eventInRange(3.0, PI + -0.1e-14), ==, true);
  TEST_COMPARE(te->eventInRange(3.0, PI + 0.0), ==, true);
  TEST_COMPARE(te->eventInRange(3.0, PI + 0.1e-14), ==, true);
  TEST_COMPARE(te->eventInRange(3.0, PI + 10.0e-14), ==, true);

  TEST_COMPARE(te->eventInRange(4.5, 5.0 + -10.0e-14), ==,
               false);  // Around last event.
  TEST_COMPARE(te->eventInRange(4.5, 5.0 + -0.1e-14), ==, true);
  TEST_COMPARE(te->eventInRange(4.5, 5.0 + 0.0), ==, true);
  TEST_COMPARE(te->eventInRange(4.5, 5.0 + 0.1e-14), ==, true);
  TEST_COMPARE(te->eventInRange(4.5, 5.0 + 10.0e-14), ==, true);

  //   Left end.
  TEST_COMPARE(te->eventInRange(-1.0 + -10.0e-14, -0.5), ==,
               true);  // Around first event.
  TEST_COMPARE(te->eventInRange(-1.0 + -0.1e-14, -0.5), ==, false);
  TEST_COMPARE(te->eventInRange(-1.0 + 0.0, -0.5), ==, false);
  TEST_COMPARE(te->eventInRange(-1.0 + 0.1e-14, -0.5), ==, false);
  TEST_COMPARE(te->eventInRange(-1.0 + 10.0e-14, -0.5), ==, false);

  TEST_COMPARE(te->eventInRange(PI + -10.0e-14, 3.5), ==,
               true);  // Around mid event.
  TEST_COMPARE(te->eventInRange(PI + -0.1e-14, 3.5), ==, false);
  TEST_COMPARE(te->eventInRange(PI + 0.0, 3.5), ==, false);
  TEST_COMPARE(te->eventInRange(PI + 0.1e-14, 3.5), ==, false);
  TEST_COMPARE(te->eventInRange(PI + 10.0e-14, 3.5), ==, false);

  TEST_COMPARE(te->eventInRange(5.0 + -10.0e-14, 6.0), ==,
               true);  // Around last event.
  TEST_COMPARE(te->eventInRange(5.0 + -0.1e-14, 6.0), ==, false);
  TEST_COMPARE(te->eventInRange(5.0 + 0.0, 6.0), ==, false);
  TEST_COMPARE(te->eventInRange(5.0 + 0.1e-14, 6.0), ==, false);
  TEST_COMPARE(te->eventInRange(5.0 + 10.0e-14, 6.0), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, getValidParameters)
{
  auto tel = rcp(new Tempus::TimeEventList<double>());

  auto pl = tel->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Type"), ==, "List");
  TEST_COMPARE(pl->get<std::string>("Name"), ==, "TimeEventList");
  TEST_FLOATING_EQUALITY(pl->get<double>("Relative Tolerance"),
                         std::numeric_limits<double>::epsilon() * 100.0,
                         1.0e-14);
  TEST_COMPARE(pl->get<bool>("Land On Exactly"), ==, true);
  TEST_COMPARE(pl->get<std::string>("Time List"), ==, "");

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, createTimeEventList)
{
  // Construct parameterList similar to getValidParameters().
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Event List");

  pl->set("Name", "Unit Test Time Event List");
  pl->set("Type", "List");
  pl->set("Relative Tolerance", 1.0e-10);
  pl->set("Land On Exactly", false);

  std::vector<double> times;
  times.push_back(-0.1);
  times.push_back(0.1);
  times.push_back(0.5);
  times.push_back(1.1);
  std::ostringstream list;
  for (std::size_t i = 0; i < times.size() - 1; ++i) list << times[i] << ", ";
  list << times[times.size() - 1];
  pl->set<std::string>("Time List", list.str());

  // Construct TimeEventList from ParameterList.
  auto tel = Tempus::createTimeEventList<double>(pl);

  tel->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tel->getName(), ==, "Unit Test Time Event List");
  TEST_COMPARE(tel->getType(), ==, "List");
  TEST_FLOATING_EQUALITY(tel->getRelTol(), 1.0e-10, 1.0e-14);
  TEST_COMPARE(tel->getLandOnExactly(), ==, false);
  auto teList = tel->getTimeList();
  TEST_FLOATING_EQUALITY(teList[0], -0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(teList[1], 0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(teList[2], 0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(teList[3], 1.1, 1.0e-14);
}

}  // namespace Tempus_Unit_Test
