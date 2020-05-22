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

#include "Tempus_TimeEventListIndex.hpp"

#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"


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
  TEST_COMPARE(te->isIndex(te->getDefaultIndex()), ==, true );
  TEST_COMPARE(te->isIndex( 1), ==, false);

  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(2, -1), ==, false);

  // Check base class defaults (functions not implemented in TimeEventListIndex).
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
  testVector.push_back( 0);
  testVector.push_back( 7);
  testVector.push_back( 3);
  testVector.push_back(-5);

  auto te = rcp(new Tempus::TimeEventListIndex<double>(
                "TestName", testVector));

  TEST_COMPARE(te->getName(), ==, "TestName");

  auto testList = te->getIndexList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_COMPARE(testList[0], ==, -5);
  TEST_COMPARE(testList[1], ==, -2);
  TEST_COMPARE(testList[2], ==,  0);
  TEST_COMPARE(testList[3], ==,  3);
  TEST_COMPARE(testList[4], ==,  7);

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
  TEST_COMPARE(te->isIndex(-5), ==, true );
  TEST_COMPARE(te->isIndex(-4), ==, false);

  // Test isIndex with two elements.
  te->addIndex(1);
  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->isIndex(1), ==, true );
  TEST_COMPARE(te->isIndex(2), ==, false);

  // Test addIndex.
  te->addIndex(-2);
  te->addIndex( 4);
  te->addIndex(-9);
  auto testList = te->getIndexList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_COMPARE(testList[0], ==, -9);
  TEST_COMPARE(testList[1], ==, -5);
  TEST_COMPARE(testList[2], ==, -2);
  TEST_COMPARE(testList[3], ==,  1);
  TEST_COMPARE(testList[4], ==,  4);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, indexToNextEvent)
{
  std::vector<int> testListIndex;
  testListIndex.push_back(-5);
  testListIndex.push_back( 1);
  testListIndex.push_back(-2);
  testListIndex.push_back( 4);
  testListIndex.push_back(-9);

  auto te = rcp(new Tempus::TimeEventListIndex<double>(
                "teListIndex", testListIndex));

  // Test indexToNextEvent.
  //   Around first event.
  TEST_COMPARE(te->indexToNextEvent(-12), ==,  3);
  TEST_COMPARE(te->indexToNextEvent( -9), ==,  0);
  TEST_COMPARE(te->indexToNextEvent( -8), ==,  3);

  //   Around mid event.
  TEST_COMPARE(te->indexToNextEvent(-4), ==,  2);
  TEST_COMPARE(te->indexToNextEvent(-2), ==,  0);
  TEST_COMPARE(te->indexToNextEvent( 0), ==,  1);

  //   Around last event.
  TEST_COMPARE(te->indexToNextEvent(2), ==,  2);
  TEST_COMPARE(te->indexToNextEvent(4), ==,  0);
  TEST_COMPARE(te->indexToNextEvent(9), ==, -5);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, indexOfNextEvent)
{
  std::vector<int> testListIndex;
  testListIndex.push_back(-5);
  testListIndex.push_back( 1);
  testListIndex.push_back(-2);
  testListIndex.push_back( 4);
  testListIndex.push_back(-9);

  auto te = rcp(new Tempus::TimeEventListIndex<double>(
                "teListIndex", testListIndex));

  // Test indexOfNextEvent.
  //   Around first event.
  TEST_COMPARE(te->indexOfNextEvent(-12), ==, -9);
  TEST_COMPARE(te->indexOfNextEvent( -9), ==, -9);
  TEST_COMPARE(te->indexOfNextEvent( -8), ==, -5);

  //   Around mid event.
  TEST_COMPARE(te->indexOfNextEvent(-4), ==, -2);
  TEST_COMPARE(te->indexOfNextEvent(-2), ==, -2);
  TEST_COMPARE(te->indexOfNextEvent( 0), ==,  1);

  //   Around last event.
  TEST_COMPARE(te->indexOfNextEvent(2), ==,  4);
  TEST_COMPARE(te->indexOfNextEvent(4), ==,  4);
  TEST_COMPARE(te->indexOfNextEvent(9), ==,  4);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventListIndex, eventInRangeIndex)
{
  std::vector<int> testListIndex;
  testListIndex.push_back(-5);
  testListIndex.push_back( 1);
  testListIndex.push_back(-2);
  testListIndex.push_back( 4);
  testListIndex.push_back(-9);

  auto te = rcp(new Tempus::TimeEventListIndex<double>(
                "teListIndex", testListIndex));

  // Test eventInRangeIndex.
  //   Right end.
  TEST_COMPARE(te->eventInRangeIndex(-12.0, -10), ==, false);   // Around first event.
  TEST_COMPARE(te->eventInRangeIndex(-12.0,  -9), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-12.0,  -8), ==, true );

  TEST_COMPARE(te->eventInRangeIndex(-4, -3), ==, false);   // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-4, -2), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-4, -1), ==, true );

  TEST_COMPARE(te->eventInRangeIndex(3, 3), ==, false);   // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(3, 4), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(3, 6), ==, true );

  //   Left end.
  TEST_COMPARE(te->eventInRangeIndex(-12, -7), ==, true );   // Around first event.
  TEST_COMPARE(te->eventInRangeIndex( -9, -7), ==, true );
  TEST_COMPARE(te->eventInRangeIndex( -8, -7), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(-3, 0), ==, true );   // Around mid event.
  TEST_COMPARE(te->eventInRangeIndex(-2, 0), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(-1, 0), ==, false);

  TEST_COMPARE(te->eventInRangeIndex(3, 8), ==, true );   // Around last event.
  TEST_COMPARE(te->eventInRangeIndex(4, 8), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(5, 8), ==, false);
}


} // namespace Tempus_Test
