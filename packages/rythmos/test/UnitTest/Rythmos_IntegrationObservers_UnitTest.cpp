//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_MockIntegrationObserver.hpp"
#include "Rythmos_MockStepperDecorator.hpp"
#include "Rythmos_TimeRange.hpp"
#include "Rythmos_StepControlInfo.hpp"

#include "Rythmos_CompositeIntegrationObserver.hpp"
#include "Rythmos_LoggingIntegrationObserver.hpp"

namespace Rythmos {

TEUCHOS_UNIT_TEST( Rythmos_IntegrationObservers, MockIntegrationObserver) {
  
  RCP<MockIntegrationObserver<double> > 
    observer = createMockIntegrationObserver<double>();

  std::list<std::string> call_stack;
  call_stack.push_back(observer->nameCloneIntegrationObserver_);
  call_stack.push_back(observer->nameResetIntegrationObserver_);
  call_stack.push_back(observer->nameObserveStartTimeIntegration_);
  call_stack.push_back(observer->nameObserveStartTimeStep_);
  call_stack.push_back(observer->nameObserveCompletedTimeStep_);
  call_stack.push_back(observer->nameObserveStartTimeStep_);
  call_stack.push_back(observer->nameObserveFailedTimeStep_);
  call_stack.push_back(observer->nameObserveStartTimeStep_);
  call_stack.push_back(observer->nameObserveCompletedTimeStep_);
  call_stack.push_back(observer->nameObserveStartTimeStep_);
  call_stack.push_back(observer->nameObserveCompletedTimeStep_);
  call_stack.push_back(observer->nameObserveEndTimeIntegration_);
  call_stack.push_back(observer->nameResetIntegrationObserver_);

  // Objects needed for observer calls
  const TimeRange<double> timeRange;
  MockStepperDecorator<double> stepper; 
  StepControlInfo<double> stepControlInfo;

  // Test a correct series of calls
  {
    observer->setCallStack(call_stack);

    observer->cloneIntegrationObserver();
    observer->resetIntegrationObserver(timeRange);
    observer->observeStartTimeIntegration(stepper);
    observer->observeStartTimeStep(stepper, stepControlInfo, 0);
    observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
    observer->observeStartTimeStep(stepper, stepControlInfo, 0);
    observer->observeFailedTimeStep(stepper, stepControlInfo, 0);
    observer->observeStartTimeStep(stepper, stepControlInfo, 0);
    observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
    observer->observeStartTimeStep(stepper, stepControlInfo, 0);
    observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
    observer->observeEndTimeIntegration(stepper);
    observer->resetIntegrationObserver(timeRange);
    
    // Stack should be empty, check that it throws gracefully (no seg fault) 
    TEST_THROW(observer->observeCompletedTimeStep(stepper, stepControlInfo, 0),
	       std::logic_error);
  }

  // Reset call stack and test incorrect series of calls
  {
    observer->setCallStack(call_stack);

    observer->cloneIntegrationObserver();

    TEST_THROW(observer->observeCompletedTimeStep(stepper, stepControlInfo, 0),
	       std::logic_error);
  }

}


TEUCHOS_UNIT_TEST( Rythmos_IntegrationObservers, LoggingIntegrationObserver) {
  
  RCP<LoggingIntegrationObserver<double> > 
    observer = createLoggingIntegrationObserver<double>();

  std::list<std::string> call_stack;
  call_stack.push_back(observer->nameCloneIntegrationObserver_);
  call_stack.push_back(observer->nameResetIntegrationObserver_);
  call_stack.push_back(observer->nameObserveStartTimeIntegration_);
  call_stack.push_back(observer->nameObserveStartTimeStep_);
  call_stack.push_back(observer->nameObserveCompletedTimeStep_);
  call_stack.push_back(observer->nameObserveStartTimeStep_);
  call_stack.push_back(observer->nameObserveFailedTimeStep_);
  call_stack.push_back(observer->nameObserveStartTimeStep_);
  call_stack.push_back(observer->nameObserveCompletedTimeStep_);
  call_stack.push_back(observer->nameObserveStartTimeStep_);
  call_stack.push_back(observer->nameObserveCompletedTimeStep_);
  call_stack.push_back(observer->nameObserveEndTimeIntegration_);
  call_stack.push_back(observer->nameResetIntegrationObserver_);

  // Objects needed for observer calls
  const TimeRange<double> timeRange;
  MockStepperDecorator<double> stepper; 
  StepControlInfo<double> stepControlInfo;

  observer->cloneIntegrationObserver();
  observer->resetIntegrationObserver(timeRange);
  observer->observeStartTimeIntegration(stepper);
  observer->observeStartTimeStep(stepper, stepControlInfo, 0);
  observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
  observer->observeStartTimeStep(stepper, stepControlInfo, 0);
  observer->observeFailedTimeStep(stepper, stepControlInfo, 0);
  observer->observeStartTimeStep(stepper, stepControlInfo, 0);
  observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
  observer->observeStartTimeStep(stepper, stepControlInfo, 0);
  observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
  observer->observeEndTimeIntegration(stepper);
  observer->resetIntegrationObserver(timeRange);
    
  const std::map<std::string,int>& counters = *(observer->getCounters());
  const std::list<std::string>& order = *(observer->getOrder());

  TEST_EQUALITY(counters.find(observer->nameCloneIntegrationObserver_)->second, 1);
  TEST_EQUALITY(counters.find(observer->nameResetIntegrationObserver_)->second, 2);
  TEST_EQUALITY(counters.find(observer->nameObserveStartTimeIntegration_)->second, 1);
  TEST_EQUALITY(counters.find(observer->nameObserveEndTimeIntegration_)->second, 1);
  TEST_EQUALITY(counters.find(observer->nameObserveStartTimeStep_)->second, 4);
  TEST_EQUALITY(counters.find(observer->nameObserveCompletedTimeStep_)->second, 3);
  TEST_EQUALITY(counters.find(observer->nameObserveFailedTimeStep_)->second, 1);

  TEUCHOS_ASSERT(order.size() == call_stack.size());
  std::list<std::string>::const_iterator observer_order = order.begin();  
  std::list<std::string>::const_iterator gold_standard = call_stack.begin();  
  for ( ; observer_order != order.end(); ++observer_order,++gold_standard)
    TEST_EQUALITY(*observer_order, *gold_standard);

  observer->resetLogCounters();
  TEST_EQUALITY(counters.find(observer->nameCloneIntegrationObserver_)->second, 0);
  TEST_EQUALITY(counters.find(observer->nameResetIntegrationObserver_)->second, 0);
  TEST_EQUALITY(counters.find(observer->nameObserveStartTimeIntegration_)->second, 0);
  TEST_EQUALITY(counters.find(observer->nameObserveEndTimeIntegration_)->second, 0);
  TEST_EQUALITY(counters.find(observer->nameObserveStartTimeStep_)->second, 0);
  TEST_EQUALITY(counters.find(observer->nameObserveCompletedTimeStep_)->second, 0);
  TEST_EQUALITY(counters.find(observer->nameObserveFailedTimeStep_)->second, 0);
  TEST_EQUALITY(order.size(), 0);

}

TEUCHOS_UNIT_TEST( Rythmos_IntegrationObservers, CompositeIntegrationObserver) {
  
  RCP<MockIntegrationObserver<double> > 
    mockObserver = createMockIntegrationObserver<double>();

  std::list<std::string> call_stack;
  call_stack.push_back(mockObserver->nameCloneIntegrationObserver_);
  call_stack.push_back(mockObserver->nameResetIntegrationObserver_);
  call_stack.push_back(mockObserver->nameObserveStartTimeIntegration_);
  call_stack.push_back(mockObserver->nameObserveStartTimeStep_);
  call_stack.push_back(mockObserver->nameObserveCompletedTimeStep_);
  call_stack.push_back(mockObserver->nameObserveStartTimeStep_);
  call_stack.push_back(mockObserver->nameObserveFailedTimeStep_);
  call_stack.push_back(mockObserver->nameObserveStartTimeStep_);
  call_stack.push_back(mockObserver->nameObserveCompletedTimeStep_);
  call_stack.push_back(mockObserver->nameObserveStartTimeStep_);
  call_stack.push_back(mockObserver->nameObserveCompletedTimeStep_);
  call_stack.push_back(mockObserver->nameObserveEndTimeIntegration_);
  call_stack.push_back(mockObserver->nameResetIntegrationObserver_);

  mockObserver->setCallStack(call_stack);

  RCP<LoggingIntegrationObserver<double> > 
    loggingObserver = createLoggingIntegrationObserver<double>();

  RCP<CompositeIntegrationObserver<double> > 
    observer = createCompositeIntegrationObserver<double>();

  observer->addObserver(mockObserver);
  observer->addObserver(loggingObserver);

  // Objects needed for observer calls
  const TimeRange<double> timeRange;
  MockStepperDecorator<double> stepper; 
  StepControlInfo<double> stepControlInfo;

  observer->cloneIntegrationObserver();
  observer->resetIntegrationObserver(timeRange);
  observer->observeStartTimeIntegration(stepper);
  observer->observeStartTimeStep(stepper, stepControlInfo, 0);
  observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
  observer->observeStartTimeStep(stepper, stepControlInfo, 0);
  observer->observeFailedTimeStep(stepper, stepControlInfo, 0);
  observer->observeStartTimeStep(stepper, stepControlInfo, 0);
  observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
  observer->observeStartTimeStep(stepper, stepControlInfo, 0);
  observer->observeCompletedTimeStep(stepper, stepControlInfo, 0);
  observer->observeEndTimeIntegration(stepper);
  observer->resetIntegrationObserver(timeRange);

  const std::map<std::string,int>& counters = *(loggingObserver->getCounters());
  const std::list<std::string>& order = *(loggingObserver->getOrder());

  TEST_EQUALITY(counters.find(loggingObserver->nameCloneIntegrationObserver_)->second, 1);
  TEST_EQUALITY(counters.find(loggingObserver->nameResetIntegrationObserver_)->second, 2);
  TEST_EQUALITY(counters.find(loggingObserver->nameObserveStartTimeIntegration_)->second, 1);
  TEST_EQUALITY(counters.find(loggingObserver->nameObserveEndTimeIntegration_)->second, 1);
  TEST_EQUALITY(counters.find(loggingObserver->nameObserveStartTimeStep_)->second, 4);
  TEST_EQUALITY(counters.find(loggingObserver->nameObserveCompletedTimeStep_)->second, 3);
  TEST_EQUALITY(counters.find(loggingObserver->nameObserveFailedTimeStep_)->second, 1);

  TEUCHOS_ASSERT(order.size() == call_stack.size());
  std::list<std::string>::const_iterator observer_order = order.begin();  
  std::list<std::string>::const_iterator gold_standard = call_stack.begin();  
  for ( ; observer_order != order.end(); ++observer_order,++gold_standard)
    TEST_EQUALITY(*observer_order, *gold_standard);

}

} // namespace Rythmos



