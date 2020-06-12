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

#include "Tempus_Types.hpp"
#include "Tempus_TimeStepControl.hpp"

#include "../TestModels/SinCosModel.hpp"
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


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, Accessors)
{
  auto tsc = rcp(new Tempus::TimeStepControl<double>());
  int    iSet = 17;
  double dSet = 0.989;

  tsc->setInitTime(dSet); TEST_COMPARE( tsc->getInitTime(), ==, dSet);
  tsc->setFinalTime(dSet); TEST_COMPARE( tsc->getFinalTime(), ==, dSet);
  tsc->setMinTimeStep(dSet); TEST_COMPARE( tsc->getMinTimeStep(), ==, dSet);
  tsc->setInitTimeStep(dSet); TEST_COMPARE( tsc->getInitTimeStep(), ==, dSet);
  tsc->setMaxTimeStep(dSet); TEST_COMPARE( tsc->getMaxTimeStep(), ==, dSet);
  tsc->setInitIndex(iSet); TEST_COMPARE( tsc->getInitIndex(), ==, iSet);
  tsc->setFinalIndex(iSet); TEST_COMPARE( tsc->getFinalIndex(), ==, iSet);
  tsc->setMaxAbsError(dSet); TEST_COMPARE( tsc->getMaxAbsError(), ==, dSet);
  tsc->setMaxRelError(dSet); TEST_COMPARE( tsc->getMaxRelError(), ==, dSet);
  tsc->setMinOrder(iSet); TEST_COMPARE( tsc->getMinOrder(), ==, iSet);
  tsc->setInitOrder(iSet); TEST_COMPARE( tsc->getInitOrder(), ==, iSet);
  tsc->setMaxOrder(iSet); TEST_COMPARE( tsc->getMaxOrder(), ==, iSet);
  tsc->setStepType("Constant"); TEST_COMPARE( tsc->getStepType(), ==, "Constant");
  tsc->setStepType("Variable"); TEST_COMPARE( tsc->getStepType(), ==, "Variable");
  tsc->setOutputExactly(false); TEST_COMPARE( tsc->getOutputExactly(), ==, false);
  tsc->setOutputExactly(true); TEST_COMPARE( tsc->getOutputExactly(), ==, true);

  std::vector<int> iVSet{ 0, 1, 1, 2, 3, 5, 8, 13, 21, 34 };
  tsc->setOutputIndices(iVSet); TEUCHOS_TEST_FOR_EXCEPT(tsc->getOutputIndices() != iVSet);


  tsc->setOutputIndexInterval(iSet); TEST_COMPARE( tsc->getOutputIndexInterval(), ==, iSet);
  tsc->setOutputTimeInterval(dSet); TEST_COMPARE( tsc->getOutputTimeInterval(), ==, dSet);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, SetDtAfterOutput_Variable)
{
  // Setup the SolutionHistory --------------------------------
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =model->getNominalValues();
  auto icSolution =rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = rcp(new Tempus::SolutionState<double>(icSolution));
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
  tsc->setStepType("Variable");
  TEST_COMPARE(tsc->getOutputExactly(), ==, true);
  Tempus::Status status = Tempus::Status::WORKING;

  // ** First Timestep ** //
  // Set dt to hit outputTime.
  // If last step is PASSED (i.e., WS is null), then initWorkingState()
  //   * Deep Copy CS to WS (maybe new RCP; may recycle old RCP for WS).
  solutionHistory->initWorkingState();
  // Need to reset local RCPs for WS and CS after initialize.
  auto currentState = solutionHistory->getCurrentState();
  auto workingState = solutionHistory->getWorkingState();

  tsc->getNextTimeStep(solutionHistory, status);
  double timeNM1 = currentState->getTime();
  double timeN   = workingState->getTime();
  TEST_COMPARE(timeN, ==, outputTime);
  //TEST_COMPARE( std::abs(timeN-outputTime)/outputTime, <, 1.0e-15);
  TEST_COMPARE(workingState->getOutput(), ==, true);

  // ** Successful timestep !! ** //
  workingState->setSolutionStatus(Tempus::Status::PASSED);

  // If workingState PASSED, then RCP CS = RCP WS and RCP WS = null.
  solutionHistory->promoteWorkingState();

  // ** Second Timestep ** //
  // Set dt to timestep before output.
  solutionHistory->initWorkingState();
  // Set local RCPs for WS and CS after initialize.
  currentState = solutionHistory->getCurrentState();
  workingState = solutionHistory->getWorkingState();

  tsc->getNextTimeStep(solutionHistory, status);
  timeNM1 = currentState->getTime();
  timeN   = workingState->getTime();
  TEST_COMPARE( (timeNM1 + dt), ==, timeN);

  double dtN = workingState->getTimeStep();
  TEST_COMPARE( dt, ==, dtN);

  TEST_COMPARE(workingState->getOutput(), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, SetDtAfterOutput_Constant)
{
  // Setup the SolutionHistory --------------------------------
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =model->getNominalValues();
  auto icSolution =rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = rcp(new Tempus::SolutionState<double>(icSolution));
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
  tsc->setStepType("Constant");
  tsc->setOutputExactly(false);
  TEST_COMPARE(tsc->getOutputExactly(), ==, false);
  Tempus::Status status = Tempus::Status::WORKING;

  // Set dt to hit outputTime for first timestep.
  // If last step is PASSED (i.e., WS is null), then initWorkingState()
  //   * Deep Copy CS to WS (maybe new RCP; may recycle old RCP for WS).
  solutionHistory->initWorkingState();
  // Need to reset local RCPs for WS and CS after initialize.
  auto currentState = solutionHistory->getCurrentState();
  auto workingState = solutionHistory->getWorkingState();

  tsc->getNextTimeStep(solutionHistory, status);
  double timeN   = workingState->getTime();
  TEST_COMPARE(timeN, ==, dt);
  //TEST_COMPARE( std::abs(timeN-dt)/dt, <, 1.0e-15);
  TEST_COMPARE(workingState->getOutput(), ==, true);

  // ** Successful timestep !! ** //
  workingState->setSolutionStatus(Tempus::Status::PASSED);

  // If workingState PASSED, then RCP CS = RCP WS and RCP WS = null.
  solutionHistory->promoteWorkingState();

  // Set dt to timestep before output for second timestep.
  solutionHistory->initWorkingState();
  // Set local RCPs for WS and CS after initialize.
  currentState = solutionHistory->getCurrentState();
  workingState = solutionHistory->getWorkingState();

  tsc->getNextTimeStep(solutionHistory, status);
  timeN   = workingState->getTime();
  TEST_COMPARE( (timeN), ==, 2*dt);

  double dtN = workingState->getTimeStep();
  TEST_COMPARE( dt, ==, dtN);

  TEST_COMPARE(workingState->getOutput(), ==, false);
}


// ************************************************************
// ************************************************************
TEUCHOS_UNIT_TEST(TimeStepControl, ConstantTimeStep_Roundoff)
{
  // Setup the SolutionHistory --------------------------------
  auto model   = rcp(new Tempus_Test::SinCosModel<double>());
  Thyra::ModelEvaluatorBase::InArgs<double> inArgsIC =model->getNominalValues();
  auto icSolution =rcp_const_cast<Thyra::VectorBase<double> >(inArgsIC.get_x());
  auto icState = rcp(new Tempus::SolutionState<double>(icSolution));
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
  tsc->setStepType("Constant");
  tsc->setTimeStepControlStrategy();
  tsc->setMinOrder(1);
  tsc->setMaxOrder(1);
  tsc->setInitTimeStep(dt);
  Tempus::Status status = Tempus::Status::WORKING;


  // Take 10000 timesteps.
  for (int i=0; i < 10000; ++i) {
    solutionHistory->initWorkingState();
    tsc->getNextTimeStep(solutionHistory, status);

    // ** Successful timestep !! ** //
    solutionHistory->getWorkingState()->setSolutionStatus(Tempus::Status::PASSED);

    solutionHistory->promoteWorkingState();
  }

  auto currentState = solutionHistory->getCurrentState();
  double time = currentState->getTime();
  TEST_COMPARE( std::abs(time-1.0), <, 1.0e-15);
}


} // namespace Tempus_Test
