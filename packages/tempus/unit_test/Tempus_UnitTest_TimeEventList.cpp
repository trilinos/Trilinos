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

#include "Tempus_TimeEventList.hpp"

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
TEUCHOS_UNIT_TEST(TimeEventList, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventList<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventList");

  TEST_COMPARE(te->getTimeList().size(), ==, 0);
  TEST_FLOATING_EQUALITY(te->getRelTol(), 1.0e-14, 1.0e-14);
  TEST_COMPARE(te->getLandOnExactly(), ==, true);

  // Check base class defaults (functions not implemented in TimeEventList).
  TEST_COMPARE(te->isIndex(1), ==, false);
  TEST_COMPARE(te->indexToNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(1,4), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, Construction)
{
  std::vector<double> testVector;
  testVector.push_back(-1.0);
  testVector.push_back( 0.0);
  testVector.push_back( 5.0);
  testVector.push_back( 2.0);
  testVector.push_back(  PI);

  auto te = rcp(new Tempus::TimeEventList<double>(
                "TestName", testVector, 1.0e-14, true));

  TEST_COMPARE(te->getName(), ==, "TestName");
  TEST_FLOATING_EQUALITY(te->getRelTol(), 1.0e-14, 1.0e-14);
  TEST_COMPARE(te->getLandOnExactly(), ==, true);

  auto testList = te->getTimeList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_FLOATING_EQUALITY(testList[0], -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[1],  0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[2],  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[3],   PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[4],  5.0, 1.0e-14);
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
  te->addTime( 2.0);
  te->addTime( 5.0);

  // Add times that should not be duplicated.
  te->addTime(0.0);
  te->addTime(PI);

  auto testList = te->getTimeList();
  TEST_COMPARE(testList.size(), ==, 5);
  TEST_FLOATING_EQUALITY(testList[0], -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[1],  0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[2],  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[3],   PI, 1.0e-14);
  TEST_FLOATING_EQUALITY(testList[4],  5.0, 1.0e-14);

  // Test that two events within relative tolerance are added or not.
  te->addTime( 2.0 + 1.0e-14);
  TEST_COMPARE(te->getTimeList().size(), ==, 5);
  te->addTime( 2.0 + 1.0e-13);
  TEST_COMPARE(te->getTimeList().size(), ==, 6);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, isTime)
{
  auto te = rcp(new Tempus::TimeEventList<double>());
  te->addTime( 0.0);
  te->addTime(  PI);
  te->addTime(-1.0);
  te->addTime( 2.0);
  te->addTime( 5.0);

  // Test isTime.
  //   Around first event.
  TEST_COMPARE(te->isTime(-10.0e-14), ==, false);   // Just outside tolerance.
  TEST_COMPARE(te->isTime( -1.0e-14), ==, true );   // Just inside tolerance.
  TEST_COMPARE(te->isTime(  0.0    ), ==, true );   // Right on timeEvent.
  TEST_COMPARE(te->isTime(  1.0e-14), ==, true );   // Just inside tolerance.
  TEST_COMPARE(te->isTime( 10.0e-14), ==, false);   // Just outside tolerance.

  //   Around mid event.
  TEST_COMPARE(te->isTime(PI + -10.0e-14), ==, false); // Just outside tolerance.
  TEST_COMPARE(te->isTime(PI +  -1.0e-14), ==, true ); // Just inside tolerance.
  TEST_COMPARE(te->isTime(PI +   0.0    ), ==, true ); // Right on timeEvent.
  TEST_COMPARE(te->isTime(PI +   1.0e-14), ==, true ); // Just inside tolerance.
  TEST_COMPARE(te->isTime(PI +  10.0e-14), ==, false); // Just outside tolerance.

  //   Around last event.
  TEST_COMPARE(te->isTime(5.0 + -10.0e-14), ==, false); // Just outside tolerance.
  TEST_COMPARE(te->isTime(5.0 +  -1.0e-14), ==, true ); // Just inside tolerance.
  TEST_COMPARE(te->isTime(5.0 +   0.0    ), ==, true ); // Right on timeEvent.
  TEST_COMPARE(te->isTime(5.0 +   1.0e-14), ==, true ); // Just inside tolerance.
  TEST_COMPARE(te->isTime(5.0 +  10.0e-14), ==, false); // Just outside tolerance.
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, timeToNextEvent)
{
  std::vector<double> testList;
  testList.push_back( 0.0);
  testList.push_back(  PI);
  testList.push_back(-1.0);
  testList.push_back( 2.0);
  testList.push_back( 5.0);

  auto te = rcp(new Tempus::TimeEventList<double>(
                "testList", testList, 1.0e-14, true));

  // Test timeToNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 + -10.0e-14),     1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 +  -1.0e-14),     1.0e-14, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 +   0.0    ),     0.0    , 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 +   1.0e-14),    -1.0e-14, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(-1.0 +  10.0e-14), 1.0-1.0e-13, 1.0e-14);

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI + -10.0e-14),        1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI +  -1.0e-14),        1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI +   0.0    ),        0.0    , 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI +   1.0e-14),       -1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(PI +  10.0e-14), 5.0-PI-1.0e-13, 1.0e-14);

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+ -10.0e-14),  1.0e-13, 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+  -1.0e-14),  1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+   0.0    ),  0.0    , 1.0e-02);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+   1.0e-14), -1.0e-14, 1.0e-01);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(5.0+  10.0e-14), -1.0e-13, 1.0e-02);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, timeOfNextEvent)
{
  std::vector<double> testList;
  testList.push_back( 0.0);
  testList.push_back(  PI);
  testList.push_back(-1.0);
  testList.push_back( 2.0);
  testList.push_back( 5.0);

  auto te = rcp(new Tempus::TimeEventList<double>(
                "testList", testList, 1.0e-14, true));

  // Test timeOfNextEvent.
  //   Around first event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 + -10.0e-14), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 +  -1.0e-14), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 +   0.0    ), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 +   1.0e-14), -1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(-1.0 +  10.0e-14),  0.0, 1.0e-14);

  //   Around mid event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+ -10.0e-14),  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+  -1.0e-14),  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+   0.0    ),  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+   1.0e-14),  2.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(2.0+  10.0e-14), PI, 1.0e-14);

  //   Around last event.
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+ -10.0e-14), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+  -1.0e-14), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+   0.0    ), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+   1.0e-14), 5.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(5.0+  10.0e-14), 5.0, 1.0e-14);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventList, eventInRange)
{
  std::vector<double> testList;
  testList.push_back( 0.0);
  testList.push_back(  PI);
  testList.push_back(-1.0);
  testList.push_back( 2.0);
  testList.push_back( 5.0);

  auto te = rcp(new Tempus::TimeEventList<double>(
                "testList", testList, 1.0e-14, true));

  // Test eventInRange.
  //   Right end.
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 + -10.0e-14), ==, false);   // Around first event.
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 +  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 +   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 +   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(-2.0, -1.0 +  10.0e-14), ==, true );

  TEST_COMPARE(te->eventInRange(3.0, PI + -10.0e-14), ==, false);   // Around mid event.
  TEST_COMPARE(te->eventInRange(3.0, PI +  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(3.0, PI +   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(3.0, PI +   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(3.0, PI +  10.0e-14), ==, true );

  TEST_COMPARE(te->eventInRange(4.5, 5.0 + -10.0e-14), ==, false);   // Around last event.
  TEST_COMPARE(te->eventInRange(4.5, 5.0 +  -1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(4.5, 5.0 +   0.0    ), ==, true );
  TEST_COMPARE(te->eventInRange(4.5, 5.0 +   1.0e-14), ==, true );
  TEST_COMPARE(te->eventInRange(4.5, 5.0 +  10.0e-14), ==, true );

  //   Left end.
  TEST_COMPARE(te->eventInRange(-1.0 + -10.0e-14, -0.5), ==, true );   // Around first event.
  TEST_COMPARE(te->eventInRange(-1.0 +  -1.0e-14, -0.5), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0 +   0.0    , -0.5), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0 +   1.0e-14, -0.5), ==, true );
  TEST_COMPARE(te->eventInRange(-1.0 +  10.0e-14, -0.5), ==, false );

  TEST_COMPARE(te->eventInRange(PI + -10.0e-14, 3.5), ==, true );   // Around mid event.
  TEST_COMPARE(te->eventInRange(PI +  -1.0e-14, 3.5), ==, true );
  TEST_COMPARE(te->eventInRange(PI +   0.0    , 3.5), ==, true );
  TEST_COMPARE(te->eventInRange(PI +   1.0e-14, 3.5), ==, true );
  TEST_COMPARE(te->eventInRange(PI +  10.0e-14, 3.5), ==, false);

  TEST_COMPARE(te->eventInRange(5.0 + -10.0e-14, 6.0), ==, true );   // Around last event.
  TEST_COMPARE(te->eventInRange(5.0 +  -1.0e-14, 6.0), ==, true );
  TEST_COMPARE(te->eventInRange(5.0 +   0.0    , 6.0), ==, true );
  TEST_COMPARE(te->eventInRange(5.0 +   1.0e-14, 6.0), ==, true );
  TEST_COMPARE(te->eventInRange(5.0 +  10.0e-14, 6.0), ==, false);
}


} // namespace Tempus_Test
