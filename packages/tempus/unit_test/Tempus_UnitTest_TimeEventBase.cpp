//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_TimeEventBase.hpp"

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

  TEST_COMPARE(te->getType(), ==, "Base");

  TEST_COMPARE(te->getName(), ==, "TimeEventBase");
  te->setName("TestName");
  TEST_COMPARE(te->getName(), ==, "TestName");

  TEST_COMPARE(te->isTime(0.0), ==, false);
  TEST_FLOATING_EQUALITY(
      te->getAbsTol(), std::numeric_limits<double>::epsilon() * 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeToNextEvent(0.0), te->getDefaultTime(),
                         1.0e-14);
  TEST_FLOATING_EQUALITY(te->timeOfNextEvent(0.0), te->getDefaultTime(),
                         1.0e-14);
  TEST_FLOATING_EQUALITY(te->getDefaultTol(), te->getAbsTol(), 1.0e-14);
  TEST_COMPARE(te->eventInRange(0.0, 1.0), ==, false);

  TEST_COMPARE(te->isIndex(0), ==, false);
  TEST_COMPARE(te->indexToNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(0), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRange(0, 10), ==, false);

  // Check base class defaults.
  TEST_COMPARE(te->isIndex(1), ==, false);
  TEST_COMPARE(te->indexToNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->indexOfNextEvent(1), ==, te->getDefaultIndex());
  TEST_COMPARE(te->eventInRangeIndex(1, 4), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeEventBase, getValidParameters)
{
  auto teb = rcp(new Tempus::TimeEventBase<double>());

  auto pl = teb->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Type"), ==, "Base");
  TEST_COMPARE(pl->get<std::string>("Name"), ==, "TimeEventBase");

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

}  // namespace Tempus_Unit_Test
