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

#ifndef Rythmos_LOGGING_INTEGRATION_OBSERVER_HPP
#define Rythmos_LOGGING_INTEGRATION_OBSERVER_HPP

#include "Rythmos_IntegrationObserverBase.hpp"
#include "Teuchos_RCP.hpp"
#include <map>
#include <list>
#include <string>

namespace Rythmos {


/** \brief Logging IntegrationOberserver that counts calls to observer functions and lists their order.
 */
template<class Scalar>
class LoggingIntegrationObserver : virtual public IntegrationObserverBase<Scalar>
{
public:

  LoggingIntegrationObserver();  

  void resetLogCounters();

  Teuchos::RCP<const std::map<std::string,int> > getCounters();

  Teuchos::RCP<const std::list<std::string> > getOrder();

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

  /** \name string names logged in map 

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

      This is a const method so that it can be called from the
      derived IntegrationObserver methods that are const.
  */
  void logCall(const std::string call) const;

private:

  Teuchos::RCP< std::map<std::string,int> > counters_;
  Teuchos::RCP< std::list<std::string> > order_;

};


/** \brief Nonmember constructor.
 *
 * \relates LoggingIntegrationObserver
 */
template<class Scalar>
Teuchos::RCP<LoggingIntegrationObserver<Scalar> >
createLoggingIntegrationObserver()
{
  const Teuchos::RCP<LoggingIntegrationObserver<Scalar> > lio = 
    Teuchos::rcp(new LoggingIntegrationObserver<Scalar>);

  return lio;
}


// //////////////////////////////////////////////////////
// Implementations

template<typename Scalar>
LoggingIntegrationObserver<Scalar>::LoggingIntegrationObserver() :
  nameCloneIntegrationObserver_("cloneIntegrationObserver"),
  nameResetIntegrationObserver_("resetIntegrationObserver"),
  nameObserveStartTimeIntegration_("observeStartTimeIntegration"),
  nameObserveEndTimeIntegration_("observeEndTimeIntegration"),
  nameObserveStartTimeStep_("observeStartTimeStep"),
  nameObserveCompletedTimeStep_("observeCompletedTimeStep"),
  nameObserveFailedTimeStep_("observeFailedTimeStep")
{ 
  counters_ = Teuchos::rcp(new std::map<std::string,int>);
  order_ = Teuchos::rcp(new std::list<std::string>);
  this->resetLogCounters();
}

template<typename Scalar>
void LoggingIntegrationObserver<Scalar>::
resetLogCounters()
{
  (*counters_)[nameCloneIntegrationObserver_] = 0;
  (*counters_)[nameResetIntegrationObserver_] = 0;
  (*counters_)[nameObserveStartTimeIntegration_] = 0;
  (*counters_)[nameObserveEndTimeIntegration_] = 0;
  (*counters_)[nameObserveStartTimeStep_] = 0;
  (*counters_)[nameObserveCompletedTimeStep_] = 0;
  (*counters_)[nameObserveFailedTimeStep_] = 0;
  order_->clear();
}

template<typename Scalar>
RCP<IntegrationObserverBase<Scalar> > 
LoggingIntegrationObserver<Scalar>::cloneIntegrationObserver() const
{
  logCall(nameCloneIntegrationObserver_);
  Teuchos::RCP<IntegrationObserverBase<Scalar> > observer = 
    Teuchos::rcp(new LoggingIntegrationObserver<Scalar>(*this));
  return observer;
}
  
template<typename Scalar>
void 
LoggingIntegrationObserver<Scalar>::
resetIntegrationObserver(const TimeRange<Scalar> &integrationTimeDomain)
{
  logCall(nameResetIntegrationObserver_);
}

template<typename Scalar>
void LoggingIntegrationObserver<Scalar>::
observeStartTimeIntegration(const StepperBase<Scalar> &stepper)
{
  logCall(nameObserveStartTimeIntegration_);
}				

template<typename Scalar>
void LoggingIntegrationObserver<Scalar>::
observeEndTimeIntegration(const StepperBase<Scalar> &stepper)
{
  logCall(nameObserveEndTimeIntegration_);
}				

template<typename Scalar>
void LoggingIntegrationObserver<Scalar>::observeStartTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
  logCall(nameObserveStartTimeStep_);
}				

template<typename Scalar>
void LoggingIntegrationObserver<Scalar>::observeCompletedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
  logCall(nameObserveCompletedTimeStep_);
}				

template<typename Scalar>
void LoggingIntegrationObserver<Scalar>::observeFailedTimeStep(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfo,
    const int timeStepIter
    )
{
  logCall(nameObserveFailedTimeStep_);
}				

template<typename Scalar>
RCP<const std::map<std::string,int> > LoggingIntegrationObserver<Scalar>::getCounters()
{
  return counters_;
}

template<typename Scalar>
RCP<const std::list<std::string> > LoggingIntegrationObserver<Scalar>::getOrder()
{
  return order_;
}


template<typename Scalar>
void LoggingIntegrationObserver<Scalar>::logCall(const std::string call) const
{
  (*counters_)[call] += 1;
  order_->push_back(call);
}


} // namespace Rythmos


#endif //Rythmos_LOGGING_INTEGRATION_OBSERVER_HPP
