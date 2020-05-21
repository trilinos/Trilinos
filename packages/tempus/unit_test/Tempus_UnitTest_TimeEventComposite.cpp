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

#include "Tempus_TimeEventComposite.hpp"
#include "Tempus_TimeEventRange.hpp"
#include "Tempus_TimeEventRangeIndex.hpp"
#include "Tempus_TimeEventList.hpp"
#include "Tempus_TimeEventListIndex.hpp"


#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <cmath>
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
    "teRange1", 0.0, PI, 1.0, 1.0e-14, true));
  return teRange1;
}

Teuchos::RCP<Tempus::TimeEventRange<double> > getTestRange2()
{
  auto teRange2 = rcp(new Tempus::TimeEventRange<double>(
    "teRange2", -PI/2.0, PI/2.0, PI/4.0, 1.0e-14, true));
  return teRange2;
}

Teuchos::RCP<Tempus::TimeEventRange<double> > getTestRange3()
{
  auto teRange3 = rcp(new Tempus::TimeEventRange<double>(
    "teRange3", 4.0, 10.0, 4.0, 1.0e-14, true));
  return teRange3;
}


// TimeEventLists for testing.
// ------------------------------------------------------------
Teuchos::RCP<Tempus::TimeEventList<double> > getTestList1()
{
  std::vector<double> testList1;
  testList1.push_back(-1.0);
  testList1.push_back( 0.0);
  testList1.push_back( 5.0);
  testList1.push_back( 2.0);
  testList1.push_back(  PI);

  auto teList1 = rcp(new Tempus::TimeEventList<double>(
                     "teList1", testList1, 1.0e-14, true));
  return teList1;
}

Teuchos::RCP<Tempus::TimeEventList<double> > getTestList2()
{
  std::vector<double> testList2;
  testList2.push_back(-0.5);
  testList2.push_back( 1.25);
  testList2.push_back( 4.95);
  testList2.push_back(12.34);

  auto teList2 = rcp(new Tempus::TimeEventList<double>(
                     "teList2", testList2, 1.0e-14, true));
  return teList2;
}

Teuchos::RCP<Tempus::TimeEventList<double> > getTestList3()
{
  std::vector<double> testList3;
  testList3.push_back(-5.0);
  testList3.push_back(-PI);

  auto teList3 = rcp(new Tempus::TimeEventList<double>(
                     "teList3", testList3, 1.0e-14, true));
  return teList3;
}


// TimeEventRangeIndices for testing.
// ------------------------------------------------------------
Teuchos::RCP<Tempus::TimeEventRangeIndex<double> > getTestRangeIndex1()
{
  auto teRangeIndex1 = rcp(new Tempus::TimeEventRangeIndex<double>(
                           "teRangeIndex1", -1, 10, 3));
  return teRangeIndex1;
}

Teuchos::RCP<Tempus::TimeEventRangeIndex<double> > getTestRangeIndex2()
{
  auto teRangeIndex2 = rcp(new Tempus::TimeEventRangeIndex<double>(
                           "teRangeIndex2", -5, 8, 4));
  return teRangeIndex2;
}

Teuchos::RCP<Tempus::TimeEventRangeIndex<double> > getTestRangeIndex3()
{
  auto teRangeIndex3 = rcp(new Tempus::TimeEventRangeIndex<double>(
                           "teRangeIndex3", 10, 17, 5));
  return teRangeIndex3;
}


// TimeEventRangeIndices for testing.
// ------------------------------------------------------------
Teuchos::RCP<Tempus::TimeEventListIndex<double> > getTestListIndex1()
{
  std::vector<int> testListIndex1;
  testListIndex1.push_back(-2);
  testListIndex1.push_back( 0);
  testListIndex1.push_back( 7);
  testListIndex1.push_back( 3);
  testListIndex1.push_back(-5);

  auto teListIndex1 = rcp(new Tempus::TimeEventListIndex<double>(
                          "teListIndex1", testListIndex1));
  return teListIndex1;
}

Teuchos::RCP<Tempus::TimeEventListIndex<double> > getTestListIndex2()
{
  std::vector<int> testListIndex2;
  testListIndex2.push_back( 2);

  auto teListIndex2 = rcp(new Tempus::TimeEventListIndex<double>(
                          "teListIndex2", testListIndex2));
  return teListIndex2;
}

Teuchos::RCP<Tempus::TimeEventListIndex<double> > getTestListIndex3()
{
  std::vector<int> testListIndex3;
  testListIndex3.push_back(14);
  testListIndex3.push_back( 9);

  auto teListIndex3 = rcp(new Tempus::TimeEventListIndex<double>(
                          "teListIndex3", testListIndex3));

  return teListIndex3;
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());

  TEST_COMPARE(te->isTime(0.0), ==, false);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, false);

  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRange(0, 10), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, One_TimeEventRange)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange1 = getTestRange1();

  te->addTimeEvent(teRange1);

  TEST_COMPARE(te->isTime(0.0), ==, true );
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.1), 2.0, 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, true );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoOverlapping_TimeEventRanges)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange1 = getTestRange1();
  auto teRange2 = getTestRange2();

  te->addTimeEvent(teRange1);
  te->addTimeEvent(teRange2);
  //te->describe();

  TEST_COMPARE(te->isTime(0.0), ==, true );
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.1), PI/4.0-0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.1), PI/2.0, 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, true );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoSeparated_TimeEventRanges)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange3 = getTestRange3();
  auto teRange2 = getTestRange2();

  te->addTimeEvent(teRange3);
  te->addTimeEvent(teRange2);
  //te->describe();

  TEST_COMPARE(te->isTime(4.0), ==, true );
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(1.1), PI/2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(3.0),    4.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0),    8.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(10.0),   8.0, 1.0e-14);
  TEST_COMPARE(te->eventInRange(-3.0, -2.0), ==, false);
  TEST_COMPARE(te->eventInRange(-1.0,  0.0), ==, true );
  TEST_COMPARE(te->eventInRange( 1.0,  2.0), ==, true );
  TEST_COMPARE(te->eventInRange( 2.0,  3.0), ==, false);
  TEST_COMPARE(te->eventInRange( 5.0,  7.0), ==, false);
  TEST_COMPARE(te->eventInRange( 7.0,  9.0), ==, true );
  TEST_COMPARE(te->eventInRange( 9.0, 11.0), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, One_TimeEventList)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teList1 = getTestList1();

  te->addTimeEvent(teList1);

  TEST_COMPARE(te->isTime(2.0), ==, true );
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(3.5), 1.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-2.0), -1.0, 1.0e-14);
  TEST_COMPARE(te->eventInRange(4.99, 10.0), ==, true );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoOverlapping_TimeEventLists)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teList1 = getTestList1();
  auto teList2 = getTestList2();

  te->addTimeEvent(teList1);
  te->addTimeEvent(teList2);
  //te->describe();

  TEST_COMPARE(te->isTime(1.25), ==, true );
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(4.0), 0.95, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(6.5), 12.34, 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.1, 1.0), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoSeparated_TimeEventLists)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teList3 = getTestList3();
  auto teList2 = getTestList2();

  te->addTimeEvent(teList3);
  te->addTimeEvent(teList2);
  //te->describe();

  TEST_COMPARE(te->isTime(4.0), ==, false);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-8.9),  -5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-0.3),  1.25, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-4.0),   -PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(20.0), 12.34, 1.0e-14);
  TEST_COMPARE(te->eventInRange(-6.0, -4.0), ==, true );
  TEST_COMPARE(te->eventInRange(-3.0,  0.0), ==, true );
  TEST_COMPARE(te->eventInRange( 2.0,  3.0), ==, false);
  TEST_COMPARE(te->eventInRange( 4.9,  5.1), ==, true );
  TEST_COMPARE(te->eventInRange(12.0, 12.4), ==, true );
  TEST_COMPARE(te->eventInRange(14.0, 15.0), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, One_TimeEventRangeIndex)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teRangeIndex1 = getTestRangeIndex1();

  te->addTimeEvent(teRangeIndex1);

  TEST_COMPARE(te->isIndex(5), ==, true );
  TEST_COMPARE(te->indexToNextEvent(3), ==, 2);
  TEST_COMPARE(te->indexOfNextEvent(3), ==, 5);
  TEST_COMPARE(te->eventInRangeIndex(3, 9), ==, true );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoOverlapping_TimeEventRangeIndex)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teRangeIndex1 = getTestRangeIndex1();
  auto teRangeIndex2 = getTestRangeIndex2();

  te->addTimeEvent(teRangeIndex1);
  te->addTimeEvent(teRangeIndex2);
  //te->describe();

  TEST_COMPARE(te->isIndex(-1), ==, true );
  TEST_COMPARE(te->indexToNextEvent(-2), ==, 1);
  TEST_COMPARE(te->indexOfNextEvent( 2), ==, 2);
  TEST_COMPARE(te->eventInRangeIndex(0, 1), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoSeparate_TimeEventRangeIndex)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teRangeIndex3 = getTestRangeIndex3();
  auto teRangeIndex2 = getTestRangeIndex2();

  te->addTimeEvent(teRangeIndex3);
  te->addTimeEvent(teRangeIndex2);
  //te->describe();

  TEST_COMPARE(te->isIndex(15), ==, true );
  TEST_COMPARE(te->indexOfNextEvent( 9), ==, 10);
  TEST_COMPARE(te->indexOfNextEvent(-6), ==, -5);
  TEST_COMPARE(te->indexOfNextEvent( 6), ==,  7);
  TEST_COMPARE(te->indexOfNextEvent(16), ==, 15);
  TEST_COMPARE(te->eventInRangeIndex(-3, -2), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(-1,  0), ==, true );
  TEST_COMPARE(te->eventInRangeIndex( 1,  2), ==, false);
  TEST_COMPARE(te->eventInRangeIndex( 7,  8), ==, true );
  TEST_COMPARE(te->eventInRangeIndex( 8,  9), ==, false);
  TEST_COMPARE(te->eventInRangeIndex(10, 13), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(14, 20), ==, true );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, One_TimeEventListIndex)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teListIndex1 = getTestListIndex1();

  te->addTimeEvent(teListIndex1);

  TEST_COMPARE(te->isIndex(3), ==, true );
  TEST_COMPARE(te->indexToNextEvent(1), ==, 2);
  TEST_COMPARE(te->indexOfNextEvent(4), ==, 7);
  TEST_COMPARE(te->eventInRangeIndex(1, 3), ==, true );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoOverlapping_TimeEventListIndex)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teListIndex1 = getTestListIndex1();
  auto teListIndex2 = getTestListIndex2();

  te->addTimeEvent(teListIndex1);
  te->addTimeEvent(teListIndex2);
  //te->describe();

  TEST_COMPARE(te->isIndex(2), ==, true );
  TEST_COMPARE(te->indexToNextEvent(0), ==, 0);
  TEST_COMPARE(te->indexOfNextEvent(1), ==, 2);
  TEST_COMPARE(te->eventInRangeIndex(-1, 3), ==, true );
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, TwoSeparate_TimeEventListIndex)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teListIndex3 = getTestListIndex3();
  auto teListIndex2 = getTestListIndex2();

  te->addTimeEvent(teListIndex3);
  te->addTimeEvent(teListIndex2);
  //te->describe();

  TEST_COMPARE(te->isIndex(14), ==, true );
  TEST_COMPARE(te->indexOfNextEvent(2), ==, 2);
  TEST_COMPARE(te->indexOfNextEvent(5), ==, 9);
  TEST_COMPARE(te->indexOfNextEvent(19), ==, 14);
  TEST_COMPARE(te->eventInRangeIndex( 0,  1), ==, false);
  TEST_COMPARE(te->eventInRangeIndex( 3, 10), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(15, 20), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventComposite, OneOfEach_TimeEvent)
{
  auto te = rcp(new Tempus::TimeEventComposite<double>());
  auto teRange1      = getTestRange1();
  auto teList1       = getTestList1();
  auto teRangeIndex1 = getTestRangeIndex1();
  auto teListIndex1  = getTestListIndex1();

  te->addTimeEvent(teRange1);
  te->addTimeEvent(teList1);
  te->addTimeEvent(teRangeIndex1);
  te->addTimeEvent(teListIndex1);

  TEST_COMPARE(te->isTime (3.0), ==, true );
  TEST_COMPARE(te->isTime (2.0), ==, true );
  TEST_COMPARE(te->isIndex(  2), ==, true );
  TEST_COMPARE(te->isIndex(  3), ==, true );

  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-2.5),  1.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( 0.5),  0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( 4.5),  0.5, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent( 7.5), -2.5, 1.0e-14);

  TEST_COMPARE(te->indexToNextEvent(-6), ==,  1);
  TEST_COMPARE(te->indexToNextEvent( 1), ==,  1);
  TEST_COMPARE(te->indexToNextEvent( 7), ==,  0);
  TEST_COMPARE(te->indexToNextEvent( 9), ==, -1);

  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( -PI), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-0.5),  0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( 2.5),  3.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent( 7.5),  5.0, 1.0e-14);

  TEST_COMPARE(te->indexOfNextEvent(-6), ==, -5);
  TEST_COMPARE(te->indexOfNextEvent( 1), ==,  2);
  TEST_COMPARE(te->indexOfNextEvent( 7), ==,  7);
  TEST_COMPARE(te->indexOfNextEvent( 9), ==,  8);

  TEST_COMPARE(te->eventInRange(-5.0, -2.0), ==, false);
  TEST_COMPARE(te->eventInRange(-2.0, -0.5), ==, true );
  TEST_COMPARE(te->eventInRange( 1.2,  1.8), ==, false);
  TEST_COMPARE(te->eventInRange( 3.1,  4.0), ==, true );
  TEST_COMPARE(te->eventInRange( 4.5,  6.0), ==, true );

  TEST_COMPARE(te->eventInRangeIndex(-8, -6), ==, false);
  TEST_COMPARE(te->eventInRangeIndex( 1,  1), ==, false);
  TEST_COMPARE(te->eventInRangeIndex( 5,  7), ==, true );
  TEST_COMPARE(te->eventInRangeIndex( 8, 10), ==, true );
  TEST_COMPARE(te->eventInRangeIndex(12, 14), ==, false);

}


} // namespace Tempus_Test
