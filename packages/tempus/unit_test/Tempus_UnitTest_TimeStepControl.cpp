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

#include "Thyra_VectorStdOps.hpp"

#include "Tempus_TimeStepControl.hpp"

#include "../TestUtils/Tempus_ConvergenceTestUtils.hpp"

#include <fstream>
#include <vector>

namespace Tempus_Unit_Test {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;
using Teuchos::ParameterList;
using Teuchos::sublist;
using Teuchos::getParametersFromXmlFile;

// Comment out any of the following tests to exclude from build/run.
#define SETOUTPUTTIMES
#define SETANDGETOUTPUTINDICESANDINTERVALS


#ifdef SETOUTPUTTIMES
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, setOutputTimes)
{
  auto tsc = rcp(new Tempus::TimeStepControl<double>());

  //auto tscPL = tsc->getParameterList();
  //std::cout << "tscPL = \n" << *tscPL << std::endl;

  std::vector<double> times_in;
  //std::cout << "Test 1" << std::endl;
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

  //std::cout << "\n  times_in, times_out = " << std::endl;
  for (size_t i=0; i < times_in.size(); ++i) {
    //std::cout << std::setw(25) << std::setprecision(16) << times_in[i] << ","
    //          << std::setw(25) << std::setprecision(16) << times_out[i]
    //          << std::endl;
    maxDiff = std::max(std::abs(times_in[i] - times_out[i]), maxDiff);
  }
  //std::cout << "  maxDiff = " << maxDiff << std::endl;

  TEST_COMPARE(maxDiff, <, 1.0e-25);


  //std::cout << "Test 2" << std::endl;
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
  maxDiff = 0.0;

  //std::cout << "\n  times_in, times_out = " << std::endl;
  for (size_t i=0; i < times_in.size(); ++i) {
    //std::cout << std::setw(30) << std::setprecision(20) << times_in[i] << ","
    //          << std::setw(30) << std::setprecision(20) << times_out[i]
    //          << std::endl;
    maxDiff = std::max(std::abs(times_in[i] - times_out[i]), maxDiff);
  }
  //std::cout << "  maxDiff = " << maxDiff << std::endl;

  TEST_COMPARE(maxDiff, <, 1.0e-25);

  //tscPL = tsc->getParameterList();
  //std::cout << "tscPL = \n" << *tscPL << std::endl;
}
#endif // SETOUTPUTTIMES


#ifdef SETANDGETOUTPUTINDICESANDINTERVALS
// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, getOutputIndicesandIntervals){
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  int setOutputTimeIndex = 17;
  double setOutputTimeInterval = 1.101001000100001e-7;

  tsc->setOutputIndexInterval(setOutputTimeIndex);
  tsc->setOutputTimeInterval(setOutputTimeInterval);

  int getOutputTimeIndex = tsc->getOutputIndexInterval();
  double getOutputTimeInterval = tsc->getOutputTimeInterval();
  TEST_COMPARE(getOutputTimeInterval, ==, setOutputTimeInterval);
  TEST_COMPARE(getOutputTimeIndex, ==, setOutputTimeIndex);
}
#endif // SETANDGETOUTPUTINDICESANDINTERVALS


} // namespace Tempus_Test
