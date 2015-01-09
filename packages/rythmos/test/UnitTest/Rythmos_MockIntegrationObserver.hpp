//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
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

#ifndef Rythmos_MOCK_INTEGRATION_OBSERVER_HPP
#define Rythmos_MOCK_INTEGRATION_OBSERVER_HPP

#include "Rythmos_IntegrationObserverBase.hpp"
#include "Teuchos_Assert.hpp"
#include <list>
#include <string>

namespace Rythmos {


/** \brief Mock IntegrationOberserver that enforces call stack.
 */
template<class Scalar>
class MockIntegrationObserver : virtual public IntegrationObserverBase<Scalar>
{
public:

  MockIntegrationObserver();  

  void setCallStack(const std::list<std::string> call_stack);

  /** \name Overridden from IntegrationObserverBase */
  //@{

  RCP<IntegrationObserverBase<Scalar> > cloneIntegrationObserver() const;
  
  void 
  resetIntegrationObserver(const TimeRange<Scalar> &integrationTimeDomain);

  void observeStartTimeIntegration(const StepperBase<Scalar> &stepper);

  void observeEndTimeIntegration(const StepperBase<Scalar> &stepper);

  void observeStartTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  void observeCompletedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  void observeFailedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    );

  //@}

  /** \name string names for call stack 

      Use these strings to validate a call stack with this observer
  */
  //@{
  
  const std::string nameCloneIntegrationObserver_;
  const std::string nameResetIntegrationObserver_;
  const std::string nameObserveStartTimeIntegration_;
  const std::string nameObserveEndTimeIntegration_;
  const std::string nameObserveStartTimeStep_;
  const std::string nameObserveCompletedTimeStep_;
  const std::string nameObserveFailedTimeStep_;

  //@}

private:

  /** \brief Asserts next call on the stack is correct and removes from stack

      This is a const method so that it can be called form the
      derived IntegrationObserver methods.
  */
  void verifyCallAndPop(const std::string call) const;

private:

  mutable std::list<std::string> expectedCallStack_;

};


/** \brief Nonmember constructor.
 *
 * \relates MockIntegrationObserver
 */
template<class Scalar>
Teuchos::RCP<MockIntegrationObserver<Scalar> >
createMockIntegrationObserver()
{
  const Teuchos::RCP<MockIntegrationObserver<Scalar> > mio = 
    Teuchos::rcp(new MockIntegrationObserver<Scalar>);

  return mio;
}


// //////////////////////////////////////////////////////
// Implementations

template<typename Scalar>
MockIntegrationObserver<Scalar>::MockIntegrationObserver() :
  nameCloneIntegrationObserver_("cloneIntegrationObserver"),
  nameResetIntegrationObserver_("resetIntegrationObserver"),
  nameObserveStartTimeIntegration_("observeStartTimeIntegration"),
  nameObserveEndTimeIntegration_("observeEndTimeIntegration"),
  nameObserveStartTimeStep_("observeStartTimeStep"),
  nameObserveCompletedTimeStep_("observeCompletedTimeStep"),
  nameObserveFailedTimeStep_("observeFailedTimeStep")
{ }

template<typename Scalar>
void MockIntegrationObserver<Scalar>::
setCallStack(const std::list<std::string> expectedCallStack)
{
  expectedCallStack_ = expectedCallStack;
}

template<typename Scalar>
RCP<IntegrationObserverBase<Scalar> > 
MockIntegrationObserver<Scalar>::cloneIntegrationObserver() const
{
  this->verifyCallAndPop(nameCloneIntegrationObserver_);
  Teuchos::RCP<IntegrationObserverBase<Scalar> > observer = 
    Teuchos::rcp(new MockIntegrationObserver<Scalar>(*this));
  return observer;
}
  
template<typename Scalar>
void 
MockIntegrationObserver<Scalar>::
resetIntegrationObserver(const TimeRange<Scalar> &integrationTimeDomain)
{
  this->verifyCallAndPop(nameResetIntegrationObserver_);
}

template<typename Scalar>
void MockIntegrationObserver<Scalar>::
observeStartTimeIntegration(const StepperBase<Scalar> &stepper)
{
  this->verifyCallAndPop(nameObserveStartTimeIntegration_);
}				

template<typename Scalar>
void MockIntegrationObserver<Scalar>::
observeEndTimeIntegration(const StepperBase<Scalar> &stepper)
{
  this->verifyCallAndPop(nameObserveEndTimeIntegration_);
}				

template<typename Scalar>
void MockIntegrationObserver<Scalar>::observeStartTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
  this->verifyCallAndPop(nameObserveStartTimeStep_);
}				

template<typename Scalar>
void MockIntegrationObserver<Scalar>::observeCompletedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
  this->verifyCallAndPop(nameObserveCompletedTimeStep_);
}				

template<typename Scalar>
void MockIntegrationObserver<Scalar>::observeFailedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
  this->verifyCallAndPop(nameObserveFailedTimeStep_);
}				

template<typename Scalar>
void MockIntegrationObserver<Scalar>::verifyCallAndPop(const std::string call) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(expectedCallStack_.size() == 0, std::logic_error, "Registered the call\"" << call << "\" but the expected stack is empty! The time integration logic has broken the expected behavior of this mock time integration observer!");
  TEUCHOS_TEST_FOR_EXCEPTION(call != *expectedCallStack_.begin(), std::logic_error, "The expected IntegrationObserver call, \"" << *expectedCallStack_.begin() << "\" is not equal to the actual call recieved, \"" << call << "\".  The time integration logic has broken the expected behavior of this mock time integration observer!");
  expectedCallStack_.pop_front();
}



} // namespace Rythmos


#endif //Rythmos_MOCK_INTEGRATION_OBSERVER_HPP
