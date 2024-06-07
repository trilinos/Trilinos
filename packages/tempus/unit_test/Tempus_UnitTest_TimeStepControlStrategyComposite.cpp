//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#include "Tempus_UnitTest_Utils.hpp"

#include "Tempus_TimeStepControl.hpp"
#include "Tempus_TimeStepControlStrategyBasicVS.hpp"
#include "Tempus_TimeStepControlStrategyConstant.hpp"
#include "Tempus_TimeStepControlStrategyIntegralController.hpp"
#include "Tempus_TimeStepControlStrategyComposite.hpp"

namespace Tempus_Unit_Test {

using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::sublist;

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyComposite, Default_Construction_1)
{
  auto tscsc = rcp(new Tempus::TimeStepControlStrategyComposite<double>());
  TEUCHOS_TEST_FOR_EXCEPT(
      tscsc->isInitialized());  // Should NOT be initialized.

  auto tscsConstant =
      rcp(new Tempus::TimeStepControlStrategyConstant<double>());
  tscsc->addStrategy(tscsConstant);
  tscsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscsc->isInitialized());

  TEST_COMPARE(tscsc->size(), ==, 1);

  std::vector<Teuchos::RCP<Tempus::TimeStepControlStrategy<double>>>
      strategies = tscsc->getStrategies();

  auto strategyConstant =
      rcp_dynamic_cast<Tempus::TimeStepControlStrategyConstant<double>>(
          strategies[0]);

  TEUCHOS_TEST_FOR_EXCEPT(strategyConstant->getStepType() != "Constant");
  TEUCHOS_TEST_FOR_EXCEPT(strategyConstant->getConstantTimeStep() != 0.0);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyComposite, Default_Construction_2)
{
  auto tscsc = rcp(new Tempus::TimeStepControlStrategyComposite<double>());
  TEUCHOS_TEST_FOR_EXCEPT(
      tscsc->isInitialized());  // Should NOT be initialized.

  auto tscsBasicVS = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
  tscsc->addStrategy(tscsBasicVS);
  auto tscsIntCtrl =
      rcp(new Tempus::TimeStepControlStrategyIntegralController<double>());
  tscsc->addStrategy(tscsIntCtrl);
  tscsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscsc->isInitialized());

  TEST_COMPARE(tscsc->size(), ==, 2);

  std::vector<Teuchos::RCP<Tempus::TimeStepControlStrategy<double>>>
      strategies = tscsc->getStrategies();

  auto strategyBasicVS =
      rcp_dynamic_cast<Tempus::TimeStepControlStrategyBasicVS<double>>(
          strategies[0]);
  TEUCHOS_TEST_FOR_EXCEPT(strategyBasicVS->getStrategyType() != "Basic VS");

  auto strategyIntCtrl = rcp_dynamic_cast<
      Tempus::TimeStepControlStrategyIntegralController<double>>(strategies[1]);
  TEUCHOS_TEST_FOR_EXCEPT(strategyIntCtrl->getStrategyType() !=
                          "Integral Controller");
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyComposite, Create_Construction)
{
  // Construct ParmeterList for testing.
  auto tscsc_temp   = rcp(new Tempus::TimeStepControlStrategyComposite<double>());
  auto tscs_BasicVS = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
  auto tscs_IC =
      rcp(new Tempus::TimeStepControlStrategyIntegralController<double>());

  tscsc_temp->addStrategy(tscs_BasicVS);
  tscsc_temp->addStrategy(tscs_IC);
  tscsc_temp->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscsc_temp->isInitialized());

  auto pl =
      rcp_const_cast<Teuchos::ParameterList>(tscsc_temp->getValidParameters());

  auto tscsc = Tempus::createTimeStepControlStrategyComposite<double>(pl);

  TEST_COMPARE(tscsc->size(), ==, 2);

  std::vector<Teuchos::RCP<Tempus::TimeStepControlStrategy<double>>>
      strategies = tscsc->getStrategies();

  auto strategyBasicVS =
      Teuchos::rcp_dynamic_cast<Tempus::TimeStepControlStrategyBasicVS<double>>(
          strategies[0]);

  TEUCHOS_TEST_FOR_EXCEPT(strategyBasicVS->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(strategyBasicVS->getAmplFactor() != 1.75);
  TEUCHOS_TEST_FOR_EXCEPT(strategyBasicVS->getReductFactor() != 0.5);
  TEUCHOS_TEST_FOR_EXCEPT(strategyBasicVS->getMinEta() != 0.0);
  TEUCHOS_TEST_FOR_EXCEPT(strategyBasicVS->getMaxEta() != 1.0e+16);

  auto strategyIC = Teuchos::rcp_dynamic_cast<
      Tempus::TimeStepControlStrategyIntegralController<double>>(strategies[1]);

  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getStepType() != "Variable");
  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getController() != "PID");
  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getKI() != 0.58);
  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getKP() != 0.21);
  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getKD() != 0.10);
  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getSafetyFactor() != 0.90);
  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getSafetyFactorAfterReject() != 0.90);
  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getFacMax() != 5.0);
  TEUCHOS_TEST_FOR_EXCEPT(strategyIC->getFacMin() != 0.5);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyComposite, getValidParameters)
{
  // Construct ParmeterList for testing.
  auto tscsc_temp   = rcp(new Tempus::TimeStepControlStrategyComposite<double>());
  auto tscs_BasicVS = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
  auto tscs_IC =
      rcp(new Tempus::TimeStepControlStrategyIntegralController<double>());

  tscsc_temp->addStrategy(tscs_BasicVS);
  tscsc_temp->addStrategy(tscs_IC);
  tscsc_temp->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tscsc_temp->isInitialized());

  auto pl = tscsc_temp->getValidParameters();

  TEST_COMPARE(pl->get<std::string>("Strategy Type"), ==, "Composite");
  TEST_COMPARE(pl->get<std::string>("Strategy List"), ==,
               "Basic VS, Integral Controller");
  TEST_COMPARE(pl->isSublist("Basic VS"), ==, true);
  TEST_COMPARE(pl->isSublist("Integral Controller"), ==, true);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(
        unusedParameters.str(), ==,
        "WARNING: Parameter \"Basic VS\"    [unused] is unused\n"
        "WARNING: Parameter \"Integral Controller\"    [unused] is unused\n");
  }

  auto BasicVS_PL = pl->sublist("Basic VS");
  TEST_COMPARE(BasicVS_PL.get<std::string>("Strategy Type"), ==, "Basic VS");
  TEST_FLOATING_EQUALITY(BasicVS_PL.get<double>("Amplification Factor"), 1.75,
                         1.0e-14);
  TEST_FLOATING_EQUALITY(BasicVS_PL.get<double>("Reduction Factor"), 0.5,
                         1.0e-14);
  TEST_FLOATING_EQUALITY(
      BasicVS_PL.get<double>("Minimum Value Monitoring Function"), 0.0,
      1.0e-14);
  TEST_FLOATING_EQUALITY(
      BasicVS_PL.get<double>("Maximum Value Monitoring Function"), 1.0e+16,
      1.0e-14);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    BasicVS_PL.unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }

  auto IntCtrl_PL = pl->sublist("Integral Controller");
  TEST_COMPARE(IntCtrl_PL.get<std::string>("Strategy Type"), ==,
               "Integral Controller");
  TEST_COMPARE(IntCtrl_PL.get<std::string>("Controller Type"), ==, "PID");
  TEST_FLOATING_EQUALITY(IntCtrl_PL.get<double>("KI"), 0.58, 1.0e-14);
  TEST_FLOATING_EQUALITY(IntCtrl_PL.get<double>("KP"), 0.21, 1.0e-14);
  TEST_FLOATING_EQUALITY(IntCtrl_PL.get<double>("KD"), 0.1, 1.0e-14);
  TEST_FLOATING_EQUALITY(IntCtrl_PL.get<double>("Safety Factor"), 0.9, 1.0e-14);
  TEST_FLOATING_EQUALITY(
      IntCtrl_PL.get<double>("Safety Factor After Step Rejection"), 0.9,
      1.0e-14);
  TEST_FLOATING_EQUALITY(IntCtrl_PL.get<double>("Maximum Safety Factor"), 5.0,
                         1.0e-14);
  TEST_FLOATING_EQUALITY(IntCtrl_PL.get<double>("Minimum Safety Factor"), 0.5,
                         1.0e-14);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    IntCtrl_PL.unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControlStrategyComposite, Setting_Strategies_PLs)
{
  {  // Test with default ParameterList
    auto pl = Tempus::getTimeStepControlStrategyCompositePL<double>();

    auto tscsc = Tempus::createTimeStepControlStrategyComposite<double>(pl);
    TEUCHOS_TEST_FOR_EXCEPT(!tscsc->isInitialized());

    std::vector<Teuchos::RCP<Tempus::TimeStepControlStrategy<double>>>
        strategies = tscsc->getStrategies();

    // Default strategy is "Constant"
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStepType() == "Constant"));
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStrategyType() == "Composite"));
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStrategies().size() == 1));
    TEUCHOS_TEST_FOR_EXCEPT(!(strategies[0]->getStepType() == "Constant"));
    TEUCHOS_TEST_FOR_EXCEPT(!(strategies[0]->getStrategyType() == "Constant"));
  }

  {  // Test with empty "Strategy List"
    auto pl = Tempus::getTimeStepControlStrategyCompositePL<double>();
    pl->set("Strategy List", "");
    pl->remove("Constant");

    auto tscsc = Tempus::createTimeStepControlStrategyComposite<double>(pl);
    TEUCHOS_TEST_FOR_EXCEPT(
        tscsc->isInitialized());  // Should NOT be initialized!

    // Default strategy is "Variable"
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStepType() == "Variable"));
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStrategyType() == "Composite"));
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStrategies().size() == 0));
  }

  {  // Test with non-Tempus strategy
    auto pl = Tempus::getTimeStepControlStrategyCompositePL<double>();
    pl->remove("Constant");
    pl->set("Strategy List", "Application Strategy");

    auto nonTempusStrategyPL = Teuchos::parameterList("Application Strategy");
    nonTempusStrategyPL->set<std::string>("Strategy Type",
                                          "Application Strategy");
    nonTempusStrategyPL->set<double>("Secret Sauce", 1.2345);
    pl->set("Application Strategy", *nonTempusStrategyPL);

    auto tscsc = Tempus::createTimeStepControlStrategyComposite<double>(pl);
    TEUCHOS_TEST_FOR_EXCEPT(
        tscsc->isInitialized());  // Should NOT be initialized!

    // Default strategy is "Variable"
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStepType() == "Variable"));
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStrategyType() == "Composite"));
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStrategies().size() == 0));
  }

  {  // Test with two Tempus strategies and a non-Tempus strategy

    // Setup ParameterList for this test.
    auto temp = rcp(new Tempus::TimeStepControlStrategyComposite<double>());
    auto tscsBasicVS =
        rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
    temp->addStrategy(tscsBasicVS);
    auto tscsIntCtrl =
        rcp(new Tempus::TimeStepControlStrategyIntegralController<double>());
    temp->addStrategy(tscsIntCtrl);
    temp->initialize();

    auto pl =
        rcp_const_cast<Teuchos::ParameterList>(temp->getValidParameters());
    auto sList = pl->get<std::string>("Strategy List");
    pl->set("Strategy List", sList + ", Application Strategy");

    auto nonTempusStrategyPL = Teuchos::parameterList("Application Strategy");
    nonTempusStrategyPL->set<std::string>("Strategy Type",
                                          "Application Strategy");
    nonTempusStrategyPL->set<double>("Secret Sauce", 1.2345);
    pl->set("Application Strategy", *nonTempusStrategyPL);

    auto tscsc = Tempus::createTimeStepControlStrategyComposite<double>(pl);
    TEUCHOS_TEST_FOR_EXCEPT(!tscsc->isInitialized());

    // Default strategy is "Constant"
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStepType() == "Variable"));
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStrategyType() == "Composite"));
    TEUCHOS_TEST_FOR_EXCEPT(!(tscsc->getStrategies().size() == 2));
  }
}

}  // namespace Tempus_Unit_Test
