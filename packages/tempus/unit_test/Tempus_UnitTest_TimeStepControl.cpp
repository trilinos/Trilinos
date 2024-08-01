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
#include "Tempus_TimeStepControlStrategyConstant.hpp"
#include "Tempus_TimeStepControlStrategyBasicVS.hpp"
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
TEUCHOS_UNIT_TEST(TimeStepControl, Default_Construction)
{
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  // Test the get functions (i.e., defaults).
  TEST_COMPARE(tsc->getStepType(), ==, "Constant");
  TEST_FLOATING_EQUALITY(tsc->getInitTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getFinalTime(), 1.0e+99, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMinTimeStep(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getInitTimeStep(), 1.0e+99, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxTimeStep(), 1.0e+99, 1.0e-14);
  TEST_COMPARE(tsc->getInitIndex(), ==, 0);
  TEST_COMPARE(tsc->getFinalIndex(), ==, 1000000);
  TEST_FLOATING_EQUALITY(tsc->getMaxAbsError(), 1.0e-08, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxRelError(), 1.0e-08, 1.0e-14);
  TEST_COMPARE(tsc->getMaxFailures(), ==, 10);
  TEST_COMPARE(tsc->getMaxConsecFailures(), ==, 5);
  TEST_COMPARE(tsc->getNumTimeSteps(), ==, -1);
  TEST_COMPARE(tsc->getPrintDtChanges(), ==, true);
  TEST_COMPARE(tsc->getOutputExactly(), ==, true);
  TEST_COMPARE(tsc->getOutputIndexInterval(), ==, 1000000);
  TEST_FLOATING_EQUALITY(tsc->getOutputTimeInterval(), 1.0e+99, 1.0e-14);
  auto tec = tsc->getTimeEvents();
  TEST_COMPARE(tec->getSize(), ==, 2);
  TEST_COMPARE(tec->getTimeEventNames(), ==,
               "Output Index Interval, Output Time Interval");

  // Test the set functions.
  tsc->setInitTime(1.0);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setFinalTime(100.0);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setMinTimeStep(0.01);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setInitTimeStep(0.02);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setMaxTimeStep(0.05);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setInitIndex(-100);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setFinalIndex(100);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setMaxAbsError(1.0e-06);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setMaxRelError(1.0e-06);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setMaxFailures(8);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setMaxConsecFailures(4);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setNumTimeSteps(-1);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setPrintDtChanges(false);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setOutputExactly(false);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setOutputIndexInterval(9);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());
  tsc->setOutputTimeInterval(0.1);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  TEST_COMPARE(tsc->getStepType(), ==, "Constant");
  TEST_FLOATING_EQUALITY(tsc->getInitTime(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getFinalTime(), 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMinTimeStep(), 0.01, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getInitTimeStep(), 0.02, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxTimeStep(), 0.05, 1.0e-14);
  TEST_COMPARE(tsc->getInitIndex(), ==, -100);
  TEST_COMPARE(tsc->getFinalIndex(), ==, 100);
  TEST_FLOATING_EQUALITY(tsc->getMaxAbsError(), 1.0e-06, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxRelError(), 1.0e-06, 1.0e-14);
  TEST_COMPARE(tsc->getMaxFailures(), ==, 8);
  TEST_COMPARE(tsc->getMaxConsecFailures(), ==, 4);
  TEST_COMPARE(tsc->getNumTimeSteps(), ==, -1);
  TEST_COMPARE(tsc->getPrintDtChanges(), ==, false);
  TEST_COMPARE(tsc->getOutputExactly(), ==, false);
  TEST_COMPARE(tsc->getOutputIndexInterval(), ==, 9);
  TEST_FLOATING_EQUALITY(tsc->getOutputTimeInterval(), 0.1, 1.0e-14);

  auto tecTmp = rcp(new Tempus::TimeEventComposite<double>());
  auto ter    = rcp(new Tempus::TimeEventRange<double>());
  auto tel    = rcp(new Tempus::TimeEventList<double>());
  ter->setName("Test Range");
  tel->setName("Test List");
  tecTmp->add(ter);
  tecTmp->add(tel);
  tsc->setTimeEvents(tecTmp);
  tec = tsc->getTimeEvents();
  TEST_COMPARE(tec->getSize(), ==, 2);
  TEST_COMPARE(tec->getTimeEventNames(), ==, "Test Range, Test List");
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, Full_Construction)
{
  std::vector<int> outputIndices;
  outputIndices.push_back(7);
  outputIndices.push_back(11);
  outputIndices.push_back(13);

  std::vector<double> outputTimes;
  outputTimes.push_back(0.3);
  outputTimes.push_back(0.7);
  outputTimes.push_back(1.3);
  outputTimes.push_back(1.7);

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

  auto tscsc = rcp(new Tempus::TimeStepControlStrategyConstant<double>());

  auto tsc = rcp(new Tempus::TimeStepControl<double>(
      1.0,           /* initTime_ */
      100.0,         /* finalTime_ */
      0.01,          /* minTimeStep_ */
      0.02,          /* initTimeStep_ */
      0.05,          /* maxTimeStep_ */
      -100,          /* initIndex_ */
      100,           /* finalIndex_ */
      1.0e-06,       /* maxAbsError_ */
      1.0e-06,       /* maxRelError_ */
      8,             /* maxFailures_ */
      4,             /* maxConsecFailures_ */
      -1,            /* numTimeSteps_ */
      false,         /* printDtChanges_ */
      false,         /* outputExactly_ */
      outputIndices, /* outputIndices_ */
      outputTimes,   /* outputTimes_ */
      9,             /* outputIndexInterval_ */
      0.011,         /* outputTimeInterval_ */
      tecTmp,        /* timeEvent_ */
      tscsc          /* stepControlStrategy_ */
      ));

  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  TEST_COMPARE(tsc->getStepType(), ==, "Constant");
  TEST_FLOATING_EQUALITY(tsc->getInitTime(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getFinalTime(), 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMinTimeStep(), 0.01, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getInitTimeStep(), 0.02, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxTimeStep(), 0.05, 1.0e-14);
  TEST_COMPARE(tsc->getInitIndex(), ==, -100);
  TEST_COMPARE(tsc->getFinalIndex(), ==, 100);
  TEST_FLOATING_EQUALITY(tsc->getMaxAbsError(), 1.0e-06, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxRelError(), 1.0e-06, 1.0e-14);
  TEST_COMPARE(tsc->getMaxFailures(), ==, 8);
  TEST_COMPARE(tsc->getMaxConsecFailures(), ==, 4);
  TEST_COMPARE(tsc->getNumTimeSteps(), ==, -1);
  TEST_COMPARE(tsc->getPrintDtChanges(), ==, false);
  TEST_COMPARE(tsc->getOutputExactly(), ==, false);
  TEST_COMPARE(tsc->getOutputIndexInterval(), ==, 9);
  TEST_FLOATING_EQUALITY(tsc->getOutputTimeInterval(), 0.011, 1.0e-14);
  auto tec = tsc->getTimeEvents();
  TEST_COMPARE(tec->getSize(), ==, 8);
  TEST_COMPARE(
      tec->getTimeEventNames(), ==,
      "Test Range, Test Range Index, Test List, Test List Index, Output Time "
      "Interval, Output Time List, Output Index Interval, Output Index List");
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, createTimeStepControl)
{
  Teuchos::RCP<Teuchos::ParameterList> pl =
      Tempus::getTimeStepControlPL<double>();

  pl->set<double>("Initial Time", 1.0);
  pl->set<double>("Final Time", 100.0);
  pl->set<double>("Minimum Time Step", 0.01);
  pl->set<double>("Initial Time Step", 0.02);
  pl->set<double>("Maximum Time Step", 0.05);
  pl->set<int>("Initial Time Index", -100);
  pl->set<int>("Final Time Index", 100);
  pl->set<double>("Maximum Absolute Error", 1.0e-06);
  pl->set<double>("Maximum Relative Error", 1.0e-06);
  pl->set<int>("Maximum Number of Stepper Failures", 8);
  pl->set<int>("Maximum Number of Consecutive Stepper Failures", 4);
  pl->set<int>("Number of Time Steps", -1);
  pl->set<bool>("Print Time Step Changes", false);
  pl->set<bool>("Output Exactly On Output Times", false);
  pl->set<std::string>("Output Index List", "7, 11, 13");
  pl->set<std::string>("Output Time List", "0.3, 0.7, 1.3, 1.7");
  pl->set<int>("Output Index Interval", 9);
  pl->set<double>("Output Time Interval", 0.011);

  auto tscs   = rcp(new Tempus::TimeStepControlStrategyConstant<double>());
  auto tscsPL = tscs->getValidParameters();
  pl->set("Time Step Control Strategy", *tscsPL);

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
  auto tecPL = rcp_const_cast<ParameterList>(tec->getValidParameters());
  pl->set("Time Step Control Events", *tecPL);

  auto tsc = Tempus::createTimeStepControl<double>(pl);
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  tsc->describe(out, Teuchos::VERB_EXTREME);

  TEST_COMPARE(tsc->getStepType(), ==, "Constant");
  TEST_FLOATING_EQUALITY(tsc->getInitTime(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getFinalTime(), 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMinTimeStep(), 0.01, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getInitTimeStep(), 0.02, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxTimeStep(), 0.05, 1.0e-14);
  TEST_COMPARE(tsc->getInitIndex(), ==, -100);
  TEST_COMPARE(tsc->getFinalIndex(), ==, 100);
  TEST_FLOATING_EQUALITY(tsc->getMaxAbsError(), 1.0e-06, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxRelError(), 1.0e-06, 1.0e-14);
  TEST_COMPARE(tsc->getMaxFailures(), ==, 8);
  TEST_COMPARE(tsc->getMaxConsecFailures(), ==, 4);
  TEST_COMPARE(tsc->getNumTimeSteps(), ==, -1);
  TEST_COMPARE(tsc->getPrintDtChanges(), ==, false);
  TEST_COMPARE(tsc->getOutputExactly(), ==, false);
  TEST_COMPARE(tsc->getOutputIndices()[0], ==, 7);
  TEST_COMPARE(tsc->getOutputIndices()[1], ==, 11);
  TEST_COMPARE(tsc->getOutputIndices()[2], ==, 13);
  TEST_FLOATING_EQUALITY(tsc->getOutputTimes()[0], 0.3, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getOutputTimes()[1], 0.7, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getOutputTimes()[2], 1.3, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getOutputTimes()[3], 1.7, 1.0e-14);
  TEST_COMPARE(tsc->getOutputIndexInterval(), ==, 9);
  TEST_FLOATING_EQUALITY(tsc->getOutputTimeInterval(), 0.011, 1.0e-14);

  tec = tsc->getTimeEvents();
  TEST_COMPARE(tec->getSize(), ==, 8);
  TEST_COMPARE(
      tec->getTimeEventNames(), ==,
      "Output Index List, Output Index Interval, Output Time List, Output Time "
      "Interval, Test Range, Test Range Index, Test List, Test List Index");
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, Accessors)
{
  auto tsc      = rcp(new Tempus::TimeStepControl<double>());
  int iFirst    = 0;
  int iLast     = 101;
  int iStep     = 17;
  double dFirst = 0.0;
  double dLast  = 0.989;
  double dStep  = 0.01;

  tsc->setInitTime(dFirst);
  TEST_COMPARE(tsc->getInitTime(), ==, dFirst);
  tsc->setFinalTime(dLast);
  TEST_COMPARE(tsc->getFinalTime(), ==, dLast);
  tsc->setMinTimeStep(dStep);
  TEST_COMPARE(tsc->getMinTimeStep(), ==, dStep);
  tsc->setInitTimeStep(dStep);
  TEST_COMPARE(tsc->getInitTimeStep(), ==, dStep);
  tsc->setMaxTimeStep(dLast);
  TEST_COMPARE(tsc->getMaxTimeStep(), ==, dLast);
  tsc->setInitIndex(iFirst);
  TEST_COMPARE(tsc->getInitIndex(), ==, iFirst);
  tsc->setFinalIndex(iLast);
  TEST_COMPARE(tsc->getFinalIndex(), ==, iLast);
  tsc->setMaxAbsError(dStep);
  TEST_COMPARE(tsc->getMaxAbsError(), ==, dStep);
  tsc->setMaxRelError(dStep);
  TEST_COMPARE(tsc->getMaxRelError(), ==, dStep);
  tsc->setOutputExactly(false);
  TEST_COMPARE(tsc->getOutputExactly(), ==, false);
  tsc->setOutputExactly(true);
  TEST_COMPARE(tsc->getOutputExactly(), ==, true);

  std::vector<int> iVSet{0, 1, 2, 3, 5, 8, 13, 21, 34};
  tsc->setOutputIndices(iVSet);
  TEUCHOS_TEST_FOR_EXCEPT(tsc->getOutputIndices() != iVSet);

  tsc->setOutputIndexInterval(iStep);
  TEST_COMPARE(tsc->getOutputIndexInterval(), ==, iStep);
  tsc->setOutputTimeInterval(dStep);
  TEST_COMPARE(tsc->getOutputTimeInterval(), ==, dStep);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, setOutputTimes)
{
  auto tsc = rcp(new Tempus::TimeStepControl<double>());

  std::vector<double> times_in;
  times_in.push_back(0.0000000000000000e-11);
  times_in.push_back(0.1001384570000000e-11);
  times_in.push_back(0.2002769140000000e-11);
  times_in.push_back(0.3004153710000000e-11);
  times_in.push_back(0.4005538280000000e-11);
  times_in.push_back(0.5006922850000000e-11);
  times_in.push_back(0.6008307420000000e-11);
  times_in.push_back(0.7009691990000000e-11);
  times_in.push_back(0.8011076560000000e-11);
  times_in.push_back(0.9012461130000000e-11);
  times_in.push_back(1.0013845700000000e-11);

  tsc->setOutputTimes(times_in);
  tsc->initialize();
  auto times_out = tsc->getOutputTimes();
  double maxDiff = 0.0;

  // out << "\n  times_in, times_out = " << std::endl;
  for (size_t i = 0; i < times_in.size(); ++i) {
    // out << std::setw(25) << std::setprecision(16) << times_in[i] << ","
    //     << std::setw(25) << std::setprecision(16) << times_out[i]
    //     << std::endl;
    maxDiff = std::max(std::fabs(times_in[i] - times_out[i]), maxDiff);
  }
  // out << "  maxDiff = " << maxDiff << std::endl;

  TEST_COMPARE(maxDiff, <, 1.0e-25);

  times_in.clear();
  times_in.push_back(0.00000000000000000000000000000000);
  times_in.push_back(0.00000000000100138457000000009381);
  times_in.push_back(0.00000000000200276914000000018762);
  times_in.push_back(0.00000000000300415371000000007949);
  times_in.push_back(0.00000000000400553828000000037525);
  times_in.push_back(0.00000000000500692284999999986321);
  times_in.push_back(0.00000000000600830742000000015898);
  times_in.push_back(0.00000000000700969198999999964694);
  times_in.push_back(0.00000000000801107656000000075050);
  times_in.push_back(0.00000000000901246112999999943067);
  times_in.push_back(0.00000000001001384569999999972643);

  tsc->setOutputTimes(times_in);
  tsc->initialize();
  times_out = tsc->getOutputTimes();
  maxDiff   = 0.0;

  // out << "\n  times_in, times_out = " << std::endl;
  for (size_t i = 0; i < times_in.size(); ++i) {
    // out << std::setw(30) << std::setprecision(20) << times_in[i] << ","
    //     << std::setw(30) << std::setprecision(20) << times_out[i]
    //     << std::endl;
    maxDiff = std::max(std::fabs(times_in[i] - times_out[i]), maxDiff);
  }
  // out << "  maxDiff = " << maxDiff << std::endl;

  TEST_COMPARE(maxDiff, <, 1.0e-25);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, getOutputIndicesandIntervals)
{
  auto tsc                     = rcp(new Tempus::TimeStepControl<double>());
  int setOutputTimeIndex       = 17;
  double setOutputTimeInterval = 1.101001000100001e-7;

  tsc->setFinalTime(1.0);
  tsc->setOutputIndexInterval(setOutputTimeIndex);
  tsc->setOutputTimeInterval(setOutputTimeInterval);

  int getOutputTimeIndex       = tsc->getOutputIndexInterval();
  double getOutputTimeInterval = tsc->getOutputTimeInterval();
  TEST_COMPARE(getOutputTimeInterval, ==, setOutputTimeInterval);
  TEST_COMPARE(getOutputTimeIndex, ==, setOutputTimeIndex);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, timeInRange)
{
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  // Testing lambda function
  auto testTimeInRange = [=](double initTime, double finalTime) {
    tsc->setInitTime(initTime);
    tsc->setFinalTime(finalTime);
    tsc->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

    const int i                = (initTime == 0)
                                     ? 0
                                     : 1 + (int)std::floor(std::log10(std::fabs(initTime)));
    const double absTolInit10  = std::pow(10, i - 10);
    const double absTolInit15  = std::pow(10, i - 15);
    const int j                = (finalTime == 0)
                                     ? 0
                                     : 1 + (int)std::floor(std::log10(std::fabs(finalTime)));
    const double absTolFinal10 = std::pow(10, j - 10);
    const double absTolFinal15 = std::pow(10, j - 15);

    // Classic unit testing of critical values.
    if (initTime == 0.0) {
      TEUCHOS_TEST_FOR_EXCEPT(tsc->timeInRange(initTime - 0.1));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(
          tsc->timeInRange(initTime - 0.1 * std::fabs(initTime)));
    }
    TEUCHOS_TEST_FOR_EXCEPT(tsc->timeInRange(initTime - absTolInit10));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->timeInRange(initTime - absTolInit15));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->timeInRange(initTime));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->timeInRange(initTime + absTolInit15));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->timeInRange(initTime + absTolInit10));
    TEUCHOS_TEST_FOR_EXCEPT(
        !tsc->timeInRange(initTime + 0.3 * (std::fabs(finalTime - initTime))));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->timeInRange(finalTime - absTolFinal10));
    TEUCHOS_TEST_FOR_EXCEPT(tsc->timeInRange(finalTime - absTolFinal15));
    TEUCHOS_TEST_FOR_EXCEPT(tsc->timeInRange(finalTime));
    TEUCHOS_TEST_FOR_EXCEPT(tsc->timeInRange(finalTime + absTolFinal15));
    TEUCHOS_TEST_FOR_EXCEPT(tsc->timeInRange(finalTime + absTolFinal10));
    if (finalTime == 0.0) {
      TEUCHOS_TEST_FOR_EXCEPT(tsc->timeInRange(finalTime + 0.1));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(
          tsc->timeInRange(finalTime + 0.1 * std::fabs(finalTime)));
    }
  };

  // Test with initTime = 0.0
  testTimeInRange(0.0, 1.0);

  // Test with finalTime = 0.0
  testTimeInRange(-1.0, 0.0);

  // Test large and small times
  testTimeInRange(9.9e-20, 3.3e+20);
  testTimeInRange(-1.9e+20, 2.3e-20);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, indexInRange)
{
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  // Testing lambda function
  auto testIndexInRange = [=](double initIndex, double finalIndex) {
    tsc->setInitIndex(initIndex);
    tsc->setFinalIndex(finalIndex);
    tsc->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

    // Classic unit testing of critical values.
    TEUCHOS_TEST_FOR_EXCEPT(tsc->indexInRange(initIndex - 7));
    TEUCHOS_TEST_FOR_EXCEPT(tsc->indexInRange(initIndex - 1));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->indexInRange(initIndex));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->indexInRange(initIndex + 1));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->indexInRange(
        initIndex + (int)0.3 * (std::fabs(finalIndex - initIndex))));
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->indexInRange(finalIndex - 1));
    TEUCHOS_TEST_FOR_EXCEPT(tsc->indexInRange(finalIndex));
    TEUCHOS_TEST_FOR_EXCEPT(tsc->indexInRange(finalIndex + 1));
    TEUCHOS_TEST_FOR_EXCEPT(tsc->indexInRange(finalIndex + 7));
  };

  // Test with initIndex = 0.0
  testIndexInRange(0, 10);

  // Test with finalIndex = 0.0
  testIndexInRange(-10, 0);

  // Test large and small indices
  testIndexInRange(-190000, 20);
  testIndexInRange(-19, 200000);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, getValidParameters)
{
  std::vector<int> outputIndices;
  outputIndices.push_back(7);
  outputIndices.push_back(11);
  outputIndices.push_back(13);

  std::vector<double> outputTimes;
  outputTimes.push_back(0.3);
  outputTimes.push_back(0.7);
  outputTimes.push_back(1.3);
  outputTimes.push_back(1.7);

  auto tec = rcp(new Tempus::TimeEventComposite<double>());

  auto tscsc = rcp(new Tempus::TimeStepControlStrategyConstant<double>());

  auto tsc = rcp(new Tempus::TimeStepControl<double>(
      1.0, 100.0, 0.01, 0.02, 0.05, -100, 100, 1.0e-06, 1.0e-06, 8, 4, -1,
      false, false, outputIndices, outputTimes, 9, 0.011, tec, tscsc));
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  auto pl = tsc->getValidParameters();

  TEST_FLOATING_EQUALITY(pl->get<double>("Initial Time"), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Final Time"), 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Minimum Time Step"), 0.01, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Initial Time Step"), 0.02, 1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Maximum Time Step"), 0.05, 1.0e-14);
  TEST_COMPARE(pl->get<int>("Initial Time Index"), ==, -100);
  TEST_COMPARE(pl->get<int>("Final Time Index"), ==, 100);
  TEST_FLOATING_EQUALITY(pl->get<double>("Maximum Absolute Error"), 1.0e-06,
                         1.0e-14);
  TEST_FLOATING_EQUALITY(pl->get<double>("Maximum Relative Error"), 1.0e-06,
                         1.0e-14);
  TEST_COMPARE(pl->get<int>("Maximum Number of Stepper Failures"), ==, 8);
  TEST_COMPARE(pl->get<int>("Maximum Number of Consecutive Stepper Failures"),
               ==, 4);
  TEST_COMPARE(pl->get<int>("Number of Time Steps"), ==, -1);
  TEST_COMPARE(pl->get<bool>("Print Time Step Changes"), ==, false);
  TEST_COMPARE(pl->get<bool>("Output Exactly On Output Times"), ==, false);
  TEST_COMPARE(pl->get<std::string>("Output Index List"), ==, "7, 11, 13");
  TEST_COMPARE(pl->get<std::string>("Output Time List"), ==,
               "0.3, 0.7, 1.3, 1.7");
  TEST_COMPARE(pl->get<int>("Output Index Interval"), ==, 9);
  TEST_FLOATING_EQUALITY(pl->get<double>("Output Time Interval"), 0.011,
                         1.0e-14);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    pl->unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==,
                 "WARNING: Parameter \"Time Step Control Strategy\"    "
                 "[unused] is unused\n");
  }

  auto tscs_PL = pl->sublist("Time Step Control Strategy");
  TEST_COMPARE(tscs_PL.get<std::string>("Strategy Type"), ==, "Constant");
  TEST_FLOATING_EQUALITY(tscs_PL.get<double>("Time Step"), 0.0, 1.0e-14);

  {  // Ensure that parameters are "used", excluding sublists.
    std::ostringstream unusedParameters;
    tscs_PL.unused(unusedParameters);
    TEST_COMPARE(unusedParameters.str(), ==, "");
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, setNumTimeSteps)
{
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  tsc->setInitTime(0.0);
  tsc->setFinalTime(100.0);
  tsc->setMinTimeStep(0.01);
  tsc->setInitTimeStep(0.02);
  tsc->setMaxTimeStep(0.05);
  tsc->setNumTimeSteps(-1);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  TEST_FLOATING_EQUALITY(tsc->getInitTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getFinalTime(), 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMinTimeStep(), 0.01, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getInitTimeStep(), 0.02, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxTimeStep(), 0.05, 1.0e-14);
  TEST_COMPARE(tsc->getNumTimeSteps(), ==, -1);

  tsc->setNumTimeSteps(100);
  tsc->initialize();

  TEST_FLOATING_EQUALITY(tsc->getInitTime(), 0.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getFinalTime(), 100.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMinTimeStep(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getInitTimeStep(), 1.0, 1.0e-14);
  TEST_FLOATING_EQUALITY(tsc->getMaxTimeStep(), 1.0, 1.0e-14);
  TEST_COMPARE(tsc->getNumTimeSteps(), ==, 100);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, SetDtAfterOutput_Variable)
{
  // Setup the SolutionHistory --------------------------------
  auto model           = rcp(new Tempus_Test::SinCosModel<double>());
  auto inArgsIC        = model->getNominalValues();
  auto icSolution      = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX<double>(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);
  double dt = 0.9;
  solutionHistory->getCurrentState()->setTimeStep(dt);

  // Setup the TimeStepControl --------------------------------
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  std::vector<double> outputTimes;
  double outputTime = 0.8;
  outputTimes.push_back(outputTime);
  tsc->setOutputTimes(outputTimes);
  tsc->setOutputExactly(true);
  auto tscs = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
  tsc->setTimeStepControlStrategy(tscs);
  tsc->setMinTimeStep(dt / 2.0);
  tsc->setInitTimeStep(dt);
  tsc->setMaxTimeStep(2.0 * dt);
  tsc->setPrintDtChanges(true);
  tsc->initialize();
  TEST_COMPARE(tsc->getOutputExactly(), ==, true);
  Tempus::Status status = Tempus::Status::WORKING;

  // ** First Timestep ** //
  // Set dt to hit outputTime.
  solutionHistory->initWorkingState();
  auto currentState = solutionHistory->getCurrentState();
  auto workingState = solutionHistory->getWorkingState();

  tsc->setNextTimeStep(solutionHistory, status);

  TEST_FLOATING_EQUALITY(workingState->getTimeStep(), outputTime, 1.0e-14);
  TEST_FLOATING_EQUALITY(workingState->getTime(), outputTime, 1.0e-14);
  TEST_COMPARE(workingState->getOutput(), ==, true);

  // ** Successful timestep !! ** //
  workingState->setSolutionStatus(Tempus::Status::PASSED);

  solutionHistory->promoteWorkingState();

  // ** Second Timestep ** //
  // Set dt to timestep before output.
  solutionHistory->initWorkingState();
  currentState = solutionHistory->getCurrentState();
  workingState = solutionHistory->getWorkingState();

  tsc->setNextTimeStep(solutionHistory, status);

  TEST_FLOATING_EQUALITY(workingState->getTimeStep(), dt, 1.0e-14);
  TEST_FLOATING_EQUALITY(currentState->getTime() + workingState->getTimeStep(),
                         workingState->getTime(), 1.0e-14);

  TEST_COMPARE(workingState->getOutput(), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, SetDtAfterOutput_Constant)
{
  // Setup the SolutionHistory --------------------------------
  auto model           = rcp(new Tempus_Test::SinCosModel<double>());
  auto inArgsIC        = model->getNominalValues();
  auto icSolution      = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX<double>(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);
  double dt = 0.9;
  solutionHistory->getCurrentState()->setTimeStep(dt);

  // Setup the TimeStepControl --------------------------------
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  std::vector<double> outputTimes;
  double outputTime = 0.8;
  outputTimes.push_back(outputTime);
  tsc->setOutputTimes(outputTimes);
  tsc->setMinTimeStep(dt / 2.0);
  tsc->setInitTimeStep(dt);
  tsc->setMaxTimeStep(2.0 * dt);
  tsc->setOutputExactly(false);
  tsc->initialize();
  TEST_COMPARE(tsc->getOutputExactly(), ==, false);
  Tempus::Status status = Tempus::Status::WORKING;

  // Set dt to hit outputTime for first timestep.
  solutionHistory->initWorkingState();
  auto currentState = solutionHistory->getCurrentState();
  auto workingState = solutionHistory->getWorkingState();

  tsc->setNextTimeStep(solutionHistory, status);
  double timeN = workingState->getTime();
  TEST_COMPARE(timeN, ==, dt);
  // TEST_COMPARE( std::fabs(timeN-dt)/dt, <, 1.0e-15);
  TEST_COMPARE(workingState->getOutput(), ==, true);

  // ** Successful timestep !! ** //
  workingState->setSolutionStatus(Tempus::Status::PASSED);

  solutionHistory->promoteWorkingState();

  // Set dt to timestep before output for second timestep.
  solutionHistory->initWorkingState();
  // Set local RCPs for WS and CS after initialize.
  currentState = solutionHistory->getCurrentState();
  workingState = solutionHistory->getWorkingState();

  tsc->setNextTimeStep(solutionHistory, status);
  timeN = workingState->getTime();
  TEST_COMPARE((timeN), ==, 2 * dt);

  double dtN = workingState->getTimeStep();
  TEST_COMPARE(dt, ==, dtN);

  TEST_COMPARE(workingState->getOutput(), ==, false);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, ConstantTimeStep_Roundoff)
{
  // Setup the SolutionHistory --------------------------------
  auto model           = rcp(new Tempus_Test::SinCosModel<double>());
  auto inArgsIC        = model->getNominalValues();
  auto icSolution      = rcp_const_cast<Thyra::VectorBase<double>>(inArgsIC.get_x());
  auto icState         = Tempus::createSolutionStateX<double>(icSolution);
  auto solutionHistory = rcp(new Tempus::SolutionHistory<double>());
  solutionHistory->addState(icState);
  double dt = 1.0e-04;
  solutionHistory->getCurrentState()->setTimeStep(dt);

  // Setup the TimeStepControl --------------------------------
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  std::vector<double> outputTimes;
  double outputTime = 0.8;
  outputTimes.push_back(outputTime);
  tsc->setOutputTimes(outputTimes);
  tsc->setOutputExactly(true);
  tsc->setTimeStepControlStrategy();
  tsc->setInitTimeStep(dt);
  tsc->initialize();
  Tempus::Status status = Tempus::Status::WORKING;

  // Take 10000 timesteps.
  for (int i = 0; i < 10000; ++i) {
    solutionHistory->initWorkingState();
    tsc->setNextTimeStep(solutionHistory, status);

    // ** Successful timestep !! ** //
    solutionHistory->getWorkingState()->setSolutionStatus(
        Tempus::Status::PASSED);

    solutionHistory->promoteWorkingState();
  }

  auto currentState = solutionHistory->getCurrentState();
  double time       = currentState->getTime();
  TEST_COMPARE(std::fabs(time - 1.0), <, 1.0e-15);
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, Setting_Strategies_PLs)
{
  {  // Test without "Time Step Control Strategy"
    auto pl = Tempus::getTimeStepControlPL<double>();
    pl->remove("Time Step Control Strategy");

    auto tsc = Tempus::createTimeStepControl<double>(pl, false);
    TEUCHOS_TEST_FOR_EXCEPT(tsc->isInitialized());
    tsc->initialize();
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

    // Default strategy is "Constant"
    TEUCHOS_TEST_FOR_EXCEPT(!(tsc->getStepType() == "Constant"));
    TEUCHOS_TEST_FOR_EXCEPT(
        !(tsc->getTimeStepControlStrategy()->getStrategyType() == "Constant"));
  }

  {  // Test with "Basic VS" strategy
    auto pl = Tempus::getTimeStepControlPL<double>();
    pl->remove("Time Step Control Strategy");
    pl->set("Time Step Control Strategy",
            *(Tempus::getTimeStepControlStrategyBasicVS_PL<double>()));

    auto tsc = Tempus::createTimeStepControl<double>(pl);
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

    // Strategy should be "Basic VS"
    TEUCHOS_TEST_FOR_EXCEPT(!(tsc->getStepType() == "Variable"));
    TEUCHOS_TEST_FOR_EXCEPT(
        !(tsc->getTimeStepControlStrategy()->getStrategyType() == "Basic VS"));
  }

  {  // Test with "Integral Controller" strategy
    auto pl = Tempus::getTimeStepControlPL<double>();
    pl->remove("Time Step Control Strategy");
    pl->set(
        "Time Step Control Strategy",
        *(Tempus::getTimeStepControlStrategyIntegralControllerPL<double>()));

    auto tsc = Tempus::createTimeStepControl<double>(pl);
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

    // Strategy should be "Integral Controller"
    TEUCHOS_TEST_FOR_EXCEPT(!(tsc->getStepType() == "Variable"));
    TEUCHOS_TEST_FOR_EXCEPT(
        !(tsc->getTimeStepControlStrategy()->getStrategyType() ==
          "Integral Controller"));
  }

  {  // Test with "Composite" strategy
    auto pl = Tempus::getTimeStepControlPL<double>();
    pl->remove("Time Step Control Strategy");
    pl->set("Time Step Control Strategy",
            *(Tempus::getTimeStepControlStrategyCompositePL<double>()));

    auto tsc = Tempus::createTimeStepControl<double>(pl);
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

    // Strategy should be "Composite"
    TEUCHOS_TEST_FOR_EXCEPT(!(tsc->getStepType() == "Constant"));
    TEUCHOS_TEST_FOR_EXCEPT(
        !(tsc->getTimeStepControlStrategy()->getStrategyType() == "Composite"));
  }

  {  // Test with a non-Tempus strategy
    auto pl = Tempus::getTimeStepControlPL<double>();
    pl->remove("Time Step Control Strategy");

    auto nonTempusStrategyPL =
        Teuchos::parameterList("Time Step Control Strategy");
    nonTempusStrategyPL->set<std::string>("Strategy Type",
                                          "Application Strategy");
    nonTempusStrategyPL->set<double>("Secret Sauce", 1.2345);

    pl->set("Time Step Control Strategy", *nonTempusStrategyPL);

    auto tsc = Tempus::createTimeStepControl<double>(pl);
    TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

    // Without finding a Tempus Strategy, the strategy should be the default
    // strategy "Constant"
    TEUCHOS_TEST_FOR_EXCEPT(!(tsc->getStepType() == "Constant"));
    TEUCHOS_TEST_FOR_EXCEPT(
        !(tsc->getTimeStepControlStrategy()->getStrategyType() == "Constant"));
  }
}

// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, Cast_Composite)
{
  // Setup Composite for this test.
  auto temp        = rcp(new Tempus::TimeStepControlStrategyComposite<double>());
  auto tscsBasicVS = rcp(new Tempus::TimeStepControlStrategyBasicVS<double>());
  temp->addStrategy(tscsBasicVS);
  auto tscsIntCtrl =
      rcp(new Tempus::TimeStepControlStrategyIntegralController<double>());
  temp->addStrategy(tscsIntCtrl);
  temp->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!temp->isInitialized());

  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  tsc->setTimeStepControlStrategy(temp);
  tsc->initialize();
  TEUCHOS_TEST_FOR_EXCEPT(!tsc->isInitialized());

  Teuchos::RCP<Tempus::TimeStepControlStrategy<double>> strategy =
      tsc->getTimeStepControlStrategy();

  // Test that we can cast this to a TimeStepControlStrategyComposite.
  auto tscsc =
      rcp_dynamic_cast<Tempus::TimeStepControlStrategyComposite<double>>(
          strategy);

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

}  // namespace Tempus_Unit_Test
