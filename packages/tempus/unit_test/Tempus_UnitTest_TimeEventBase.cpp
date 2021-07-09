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

#include "Tempus_TimeEventBase.hpp"

#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"


namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventBase, Default_Construction)
{
  auto te = rcp(new Tempus::TimeEventBase<double>());

  TEST_COMPARE(te->getName(), ==, "TimeEventBase");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  TEST_COMPARE(te->isTime(0.0), ==, false);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(0.0), te->getDefaultTime(), 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, false);

  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRange(0, 10), ==, false);

  // Check base class defaults (functions not implemented in TimeEventRange).
  TEST_COMPARE(te->isIndex(1), ==, false);
  TEST_COMPARE(te->indexToNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(1,4), ==, false);
}


} // namespace Tempus_Test
