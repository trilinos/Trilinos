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


#ifndef RYTHMOS_RAMPING_INTEGRATION_CONTROL_STRATEGY_DEF_HPP
#define RYTHMOS_RAMPING_INTEGRATION_CONTROL_STRATEGY_DEF_HPP

#include "Rythmos_RampingIntegrationControlStrategy_decl.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Rythmos {


template<class Scalar>
RCP<RampingIntegrationControlStrategy<Scalar> >
rampingIntegrationControlStrategy()
{
  RCP<RampingIntegrationControlStrategy<Scalar> >
    integrationControl =
      Teuchos::rcp(new RampingIntegrationControlStrategy<Scalar>());
  return integrationControl;
}


template<class Scalar>
RCP<RampingIntegrationControlStrategy<Scalar> >
rampingIntegrationControlStrategy( const RCP<ParameterList> &paramList )
{
  RCP<RampingIntegrationControlStrategy<Scalar> >
    integrationControl =
      Teuchos::rcp(new RampingIntegrationControlStrategy<Scalar>());
  integrationControl->setParameterList(paramList);
  return integrationControl;
}


//
// Implementation
//


// Static members


template<class Scalar>
const std::string
RampingIntegrationControlStrategy<Scalar>::take_variable_steps_name_
= "Take Variable Steps";

template<class Scalar>
const bool
RampingIntegrationControlStrategy<Scalar>::take_variable_steps_default_
= true;


template<class Scalar>
const std::string
RampingIntegrationControlStrategy<Scalar>::num_constant_steps_name_
= "Number of Initial Constant Steps";

template<class Scalar>
const int
RampingIntegrationControlStrategy<Scalar>::num_constant_steps_default_
= 0;


template<class Scalar>
const std::string
RampingIntegrationControlStrategy<Scalar>::num_ramping_steps_name_
= "Number of Ramping Steps";

template<class Scalar>
const int
RampingIntegrationControlStrategy<Scalar>::num_ramping_steps_default_
= 6;


template<class Scalar>
const std::string
RampingIntegrationControlStrategy<Scalar>::initial_dt_name_
= "Initial dt";

template<class Scalar>
const double
RampingIntegrationControlStrategy<Scalar>::initial_dt_default_
= -1.0;


template<class Scalar>
const std::string
RampingIntegrationControlStrategy<Scalar>::min_dt_name_
= "Min dt";

template<class Scalar>
const double
RampingIntegrationControlStrategy<Scalar>::min_dt_default_
= std::numeric_limits<Scalar>::min();


template<class Scalar>
const std::string
RampingIntegrationControlStrategy<Scalar>::max_dt_name_
= "Max dt";

template<class Scalar>
const double
RampingIntegrationControlStrategy<Scalar>::max_dt_default_
= std::numeric_limits<Scalar>::max();


template<class Scalar>
const std::string
RampingIntegrationControlStrategy<Scalar>::ramping_factor_name_
= "Ramping Factor";

template<class Scalar>
const double
RampingIntegrationControlStrategy<Scalar>::ramping_factor_default_
= 1.0;


template<class Scalar>
const std::string
RampingIntegrationControlStrategy<Scalar>::max_step_failures_name_
= "Maximum Number of Step Failures";

template<class Scalar>
const int
RampingIntegrationControlStrategy<Scalar>::max_step_failures_default_
= 100;



// Constructors/Initializers


template<class Scalar>
RampingIntegrationControlStrategy<Scalar>::RampingIntegrationControlStrategy()
  : take_variable_steps_(take_variable_steps_default_),
    num_constant_steps_(num_constant_steps_default_),
    num_ramping_steps_(num_ramping_steps_default_),
    initial_dt_(initial_dt_default_),
    min_dt_(min_dt_default_),
    max_dt_(max_dt_default_),
    ramping_factor_(ramping_factor_default_),
    num_step_failures_(0),
    max_step_failures_(max_step_failures_default_),
    current_dt_(-1.0)
{}


// Overridden from ParameterListAcceptor


template<class Scalar>
void RampingIntegrationControlStrategy<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  using Teuchos::as;
  using Teuchos::get;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParametersAndSetDefaults(*getValidParameters());
  this->setMyParamList(paramList);

  take_variable_steps_ = paramList->get<bool>  (take_variable_steps_name_,
                                                take_variable_steps_default_);
  num_constant_steps_  = paramList->get<int>   (num_constant_steps_name_,
                                                num_constant_steps_default_);
  num_ramping_steps_   = paramList->get<int>   (num_ramping_steps_name_,
                                                num_ramping_steps_default_);
  initial_dt_          = paramList->get<double>(initial_dt_name_,
                                                initial_dt_default_);
  min_dt_              = paramList->get<double>(min_dt_name_, min_dt_default_);
  max_dt_              = paramList->get<double>(max_dt_name_, max_dt_default_);
  ramping_factor_      = paramList->get<double>(ramping_factor_name_,
                                                ramping_factor_default_);
  max_step_failures_   = paramList->get<int>   (max_step_failures_name_,
                                                max_step_failures_default_);

  Teuchos::readVerboseObjectSublist(&*paramList,this);
}


template<class Scalar>
RCP<const ParameterList>
RampingIntegrationControlStrategy<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();

    pl->set( take_variable_steps_name_, take_variable_steps_default_,
      "Take variable time steps after '" + num_constant_steps_name_ +
      "' plus '" + num_ramping_steps_name_ + "' steps.  Variable time "
      "stepping allows the Stepper to adjust the time step through a "
      "StepControlStrategy after fixed time steps during initial constant "
      "steps and ramping steps.  If false, fixed-time steps are taken "
      "after ramping.  Fixed time stepping requires the Stepper "
      "to take the time step set by this IntegrationControlStrategy.");
    pl->set(num_constant_steps_name_, num_constant_steps_default_,
      "Number of initial constant steps to take before starting the ramping.");
    pl->set(num_ramping_steps_name_, num_ramping_steps_default_,
      "Number of ramping steps to take before handing control to "
      "variable stepper if '" + take_variable_steps_name_ +
      "' is set to true.  Otherwise take fixed-time steps.");
    pl->set(initial_dt_name_, initial_dt_default_,
	    "Initial time step.");
    pl->set(min_dt_name_, min_dt_default_,
	    "Minimum time step.");
    pl->set(max_dt_name_, max_dt_default_,
	    "Maximum time step.");
    pl->set(ramping_factor_name_, ramping_factor_default_,
	    "Time step growth factor used during ramping phase. dt_{n+1} = "
      "(ramping factor) * dt_n");
    pl->set(max_step_failures_name_, max_step_failures_default_,
	    "The maximum number of step failures before exiting with error.");
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from IntegrationControlStrategyBase


template<class Scalar>
bool RampingIntegrationControlStrategy<Scalar>::handlesFailedTimeSteps() const
{
  return true;
}


template<class Scalar>
RCP<IntegrationControlStrategyBase<Scalar> >
RampingIntegrationControlStrategy<Scalar>::cloneIntegrationControlStrategy() const
{
  RCP<RampingIntegrationControlStrategy<Scalar> >
    integrCtrlStry = rampingIntegrationControlStrategy<Scalar>();
  const RCP<const ParameterList> paramList = this->getParameterList();
  if (!is_null(paramList))
    integrCtrlStry->setParameterList(Teuchos::parameterList(*paramList));

  integrCtrlStry->take_variable_steps_ = this->take_variable_steps_;
  integrCtrlStry->num_constant_steps_  = this->num_constant_steps_;
  integrCtrlStry->num_ramping_steps_   = this->num_ramping_steps_;
  integrCtrlStry->initial_dt_          = this->initial_dt_;
  integrCtrlStry->min_dt_              = this->min_dt_;
  integrCtrlStry->max_dt_              = this->max_dt_;
  integrCtrlStry->ramping_factor_      = this->ramping_factor_;
  integrCtrlStry->current_dt_          = this->current_dt_;

  return integrCtrlStry;
}


template<class Scalar>
void
RampingIntegrationControlStrategy<Scalar>::resetIntegrationControlStrategy(
  const TimeRange<Scalar> &integrationTimeDomain
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT(integrationTimeDomain.length() > ST::zero());
#endif
  integrationTimeDomain_ = integrationTimeDomain;
  if (max_dt_ < ST::zero()) {
    max_dt_ = integrationTimeDomain_.length();
  }

  current_dt_ = initial_dt_;
}


template<class Scalar>
StepControlInfo<Scalar>
RampingIntegrationControlStrategy<Scalar>::getNextStepControlInfo(
  const StepperBase<Scalar> &stepper,
  const StepControlInfo<Scalar> &stepCtrlInfoLast,
  const int timeStepIter
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT(integrationTimeDomain_.length() > ST::zero());
#endif

  StepControlInfo<Scalar> trialStepCtrlInfo;

  if (timeStepIter < num_constant_steps_)
    current_dt_ = initial_dt_;
  else if (timeStepIter < num_constant_steps_ + num_ramping_steps_)
    current_dt_ *= ramping_factor_;

  current_dt_ = std::min(max_dt_, current_dt_);
  current_dt_ = std::max(min_dt_, current_dt_);

  num_step_failures_ = std::min(num_step_failures_-1,0);

  trialStepCtrlInfo.stepSize = current_dt_;
  if (take_variable_steps_) {
    if (timeStepIter < num_constant_steps_ + num_ramping_steps_)
      trialStepCtrlInfo.stepType = STEP_TYPE_FIXED;
    else
      trialStepCtrlInfo.stepType = STEP_TYPE_VARIABLE;
  } else {
    trialStepCtrlInfo.stepType = STEP_TYPE_FIXED;
  }

  return trialStepCtrlInfo;
}


template<class Scalar>
bool RampingIntegrationControlStrategy<Scalar>::resetForFailedTimeStep(
  const StepperBase<Scalar> &stepper,
  const StepControlInfo<Scalar> &stepCtrlInfoLast,
  const int timeStepIter,
  const StepControlInfo<Scalar> &stepCtrlInfo
  )
{
  if (current_dt_ == min_dt_) return false;
  num_step_failures_++;
  if (num_step_failures_ > max_step_failures_) return false;
  current_dt_ = std::max(min_dt_, current_dt_/ramping_factor_);
  return true;
}


//
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_RAMPING_INTEGRATION_CONTROL_STRATEGY_INSTANT(SCALAR) \
  \
  template class RampingIntegrationControlStrategy< SCALAR >; \
  \
  template RCP<RampingIntegrationControlStrategy< SCALAR > > \
  rampingIntegrationControlStrategy(); \
  \
  template RCP<RampingIntegrationControlStrategy< SCALAR > > \
  rampingIntegrationControlStrategy( const RCP<ParameterList> &paramList );


} // namespace Rythmos


#endif // RYTHMOS_RAMPING_INTEGRATION_CONTROL_STRATEGY_DEF_HPP
