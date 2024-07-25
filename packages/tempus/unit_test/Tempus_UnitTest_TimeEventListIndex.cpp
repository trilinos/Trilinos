//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_TimeEventListIndex.hpp"

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventListIndex<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventListIndex");

  TEST_COMPARE(te->getIndexList().size(), ==, 0);

  TEST_COMPARE(te->isIndex(-1), ==, false);
  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->isIndex(1), ==, false);

  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(2, -1), ==, false);

  // Check base class defaults (functions not implemented in
  // TimeEventListIndex).
  TEST_COMPARE(te->isTime(1.0), ==, false);
  TEST_COMPARE(te->timeToNextEvent(1.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->timeOfNextEvent(1.0), ==, te->getDefaultTime());
  TEST_COMPARE(te->eventInRange(1.0, 4.0), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, Construction)
{
  std::vector<int> testVector;
  testVector.push_back(-2);
  testVector.push_back(0);
  testVector.push_back(7);
  testVector.push_back(3);
  testVector.push_back(-5);

  auto te = rcp(new Tempus::TimeEventListIndex<double>(testVector, "TestName"));

  TEST_COMPARE(te->getName(), ==, "TestName");

  auto testList = te->getIndexList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_COMPARE(testList[0], ==, -5);
  TEST_COMPARE(testList[1], ==, -2);
  TEST_COMPARE(testList[2], ==, 0);
  TEST_COMPARE(testList[3], ==, 3);
  TEST_COMPARE(testList[4], ==, 7);

  // Test adding a duplicate event index.
  te->addIndex(3);
  TEST_COMPARE(te->getIndexList().size(), ==, 5);
  te->addIndex(1);
  TEST_COMPARE(te->getIndexList().size(), ==, 6);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, isIndex)
{
  auto te = rcp(new Tempus::TimeEventListIndex<double>());

  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  // Test isIndex with one element.
  te->addIndex(-5);
  TEST_COMPARE(te->isIndex(-6), ==, false);
  TEST_COMPARE(te->isIndex(-5), ==, true);
  TEST_COMPARE(te->isIndex(-4), ==, false);

  // Test isIndex with two elements.
  te->addIndex(1);
  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->isIndex(1), ==, true);
  TEST_COMPARE(te->isIndex(2), ==, false);

  // Test addIndex.
  te->addIndex(-2);
  te->addIndex(4);
  te->addIndex(-9);
  auto testList = te->getIndexList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_COMPARE(testList[0], ==, -9);
  TEST_COMPARE(testList[1], ==, -5);
  TEST_COMPARE(testList[2], ==, -2);
  TEST_COMPARE(testList[3], ==, 1);
  TEST_COMPARE(testList[4], ==, 4);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, indexToNextEvent)
{
  std::vector<int> testListIndex;
  testListIndex.push_back(-5);
  testListIndex.push_back(1);
  testListIndex.push_back(-2);
  testListIndex.push_back(4);
  testListIndex.push_back(-9);

  auto te =
      rcp(new Tempus::TimeEventListIndex<double>(testListIndex, "teListIndex"));

  // Test indexToNextEvent.
  //   Around first event.
  TEST_COMPARE(te->indexToNextEvent(-12), ==, 3);
  TEST_COMPARE(te->indexToNextEvent(-9), ==, 4);
  TEST_COMPARE(te->indexToNextEvent(-8), ==, 3);

  //   Around mid event.
  TEST_COMPARE(te->indexToNextEvent(-4), ==, 2);
  TEST_COMPARE(te->indexToNextEvent(-2), ==, 3);
  TEST_COMPARE(te->indexToNextEvent(0), ==, 1);

  //   Around last event.
  TEST_COMPARE(te->indexToNextEvent(2), ==, 2);
  TEST_COMPARE(te->indexToNextEvent(4), ==, te->getDefaultIndex() - 4);
  TEST_COMPARE(te->indexToNextEvent(9), ==, te->getDefaultIndex() - 9);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, indexOfNextEvent)
{
  std::vector<int> testListIndex;
  testListIndex.push_back(-5);
  testListIndex.push_back(1);
  testListIndex.push_back(-2);
  testListIndex.push_back(4);
  testListIndex.push_back(-9);

  auto te =
      rcp(new Tempus::TimeEventListIndex<double>(testListIndex, "teListIndex"));

  // Test indexOfNextEvent.
  //   Around first event.
  TEST_COMPARE(te->indexOfNextEvent(-12), ==, -9);
  TEST_COMPARE(te->indexOfNextEvent(-9), ==, -5);
  TEST_COMPARE(te->indexOfNextEvent(-8), ==, -5);

  //   Around mid event.
  TEST_COMPARE(te->indexOfNextEvent(-4), ==, -2);
  TEST_COMPARE(te->indexOfNextEvent(-2), ==, 1);
  TEST_COMPARE(te->indexOfNextEvent(0), ==, 1);

  //   Around last event.
  TEST_COMPARE(te->indexOfNextEvent(2), ==, 4);
  TEST_COMPARE(te->indexOfNextEvent(4), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(9), ==, te->getDefaultIndex());
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, eventInRangeIndex)
{
  std::vector<int> testListIndex;
  testListIndex.push_back(-5);
  testListIndex.push_back(1);
  testListIndex.push_back(-2);
  testListIndex.push_back(4);
  testListIndex.push_back(-9);

  auto te =
      rcp(new Tempus::TimeEventListIndex<double>(testListIndex, "teListIndex"));

  // Test eventInRangeIndex.
  //   Right end.
  TEST_COMPARE(te->eventInRangeIndex(-12.0, -10), ==,
               false);  // Around first event.
  TEST_COMPARE(te->eventInRangeIndex(-12.0, -9), ==, true);
  TEST_COMPARE(te->eventInRangeIndex(-12.0, -8), ==, true);

  TEST_COMPARE(te->eventInRangeIndex(-4, -3), ==, false);  // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-4, -2), ==, true);
  TEST_COMPARE(te->eventInRangeIndex(-4, -1), ==, true);

  TEST_COMPARE(te->eventInRangeIndex(3, 3), ==, false);  // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(3, 4), ==, true);
  TEST_COMPARE(te->eventInRangeIndex(3, 6), ==, true);

  //   Left end.
  TEST_COMPARE(te->eventInRangeIndex(-12, -7), ==,
               true);  // Around first event.
  TEST_COMPARE(te->eventInRangeIndex(-9, -7), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(-8, -7), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(-3, 0), ==, true);  // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-2, 0), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(-1, 0), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(3, 8), ==, true);  // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(4, 8), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(5, 8), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, getValidParameters)
{
  auto teli = rcp(new Tempus::TimeEventListIndex<double>());

  auto pl = teli->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Type"), ==, "List Index");
  TEST_COMPARE(pl->get<std::string>("Name"), ==, "TimeEventListIndex");
  TEST_COMPARE(pl->get<std::string>("Index List"), ==, "");

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, createTimeEventListIndex)
{
  // Construct parameterList similar to getValidParameters().
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Teuchos::parameterList("Time Event List Index");

  pl->set("Name", "Unit Test Time Event List Index");
  pl->set("Type", "List Index");

  std::vector<int> indices;
  indices.push_back(-99);
  indices.push_back(13);
  indices.push_back(97);
  indices.push_back(101);
  std::ostringstream list;
  for (std::size_t i = 0; i < indices.size() - 1; ++i)
    list << indices[i] << ", ";
  list << indices[indices.size() - 1];
  pl->set<std::string>("Index List", list.str());

  // Construct TimeEventListIndex from ParameterList.
  auto teli = Tempus::createTimeEventListIndex<double>(pl);

  teli->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(teli->getName(), ==, "Unit Test Time Event List Index");
  TEST_COMPARE(teli->getType(), ==, "List Index");
  auto teList = teli->getIndexList();
  TEST_COMPARE(teList[0], ==, -99);
  TEST_COMPARE(teList[1], ==, 13);
  TEST_COMPARE(teList[2], ==, 97);
  TEST_COMPARE(teList[3], ==, 101);
}

}  // namespace Tempus_Unit_Test
