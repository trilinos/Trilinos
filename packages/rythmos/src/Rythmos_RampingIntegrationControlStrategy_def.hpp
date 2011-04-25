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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
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
    integrationControl = Teuchos::rcp(new RampingIntegrationControlStrategy<Scalar>());
  return integrationControl;
}


template<class Scalar> 
RCP<RampingIntegrationControlStrategy<Scalar> >
rampingIntegrationControlStrategy( const RCP<ParameterList> &paramList )
{
  RCP<RampingIntegrationControlStrategy<Scalar> >
    integrationControl = Teuchos::rcp(new RampingIntegrationControlStrategy<Scalar>());
  integrationControl->setParameterList(paramList);
  return integrationControl;
}


//
// Implementation
//


// Static members


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



// Constructors/Initializers


template<class Scalar>
RampingIntegrationControlStrategy<Scalar>::RampingIntegrationControlStrategy() :
  num_ramping_steps_(num_ramping_steps_default_),
  initial_dt_(initial_dt_default_),
  max_dt_(max_dt_default_),
  ramping_factor_(ramping_factor_default_),
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
  TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters());
  this->setMyParamList(paramList);

  num_ramping_steps_ = paramList->get<int>(num_ramping_steps_name_);
  initial_dt_ = paramList->get<double>(initial_dt_name_);
  max_dt_ = paramList->get<double>(max_dt_name_);
  ramping_factor_ = paramList->get<double>(ramping_factor_name_);

  Teuchos::readVerboseObjectSublist(&*paramList,this);
}


template<class Scalar> 
RCP<const ParameterList>
RampingIntegrationControlStrategy<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(num_ramping_steps_name_, num_ramping_steps_default_,
      "Number of ramping steps to take before handing control to variable stepper.");
    pl->set(initial_dt_name_, initial_dt_default_,
	    "Initial teim step.");
    pl->set(max_dt_name_, max_dt_default_,
	    "Maximum time step.");
    pl->set(ramping_factor_name_, ramping_factor_default_,
	    "Time step growth factor used during ramping phase.dt_{n+1} = (ramping factor) * dt_n");
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from IntegrationControlStrategyBase


template<class Scalar>
RCP<IntegrationControlStrategyBase<Scalar> >
RampingIntegrationControlStrategy<Scalar>::cloneIntegrationControlStrategy() const
{
  RCP<RampingIntegrationControlStrategy<Scalar> >
    integrCtrlStry = rampingIntegrationControlStrategy<Scalar>();
  const RCP<const ParameterList> paramList = this->getParameterList();
  if (!is_null(paramList))
    integrCtrlStry->setParameterList(Teuchos::parameterList(*paramList));
  integrCtrlStry->num_ramping_steps_ = this->num_ramping_steps_;
  integrCtrlStry->initial_dt_ = this->initial_dt_;
  integrCtrlStry->max_dt_ = this->max_dt_;
  integrCtrlStry->ramping_factor_ = this->ramping_factor_;
  integrCtrlStry->current_dt_ = this->current_dt_;
  return integrCtrlStry;
}


template<class Scalar>
void
RampingIntegrationControlStrategy<Scalar>::resetIntegrationControlStrategy(
  const TimeRange<Scalar> &integrationTimeDomain
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef RYTHMOS_DEBUG
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

#ifdef RYTHMOS_DEBUG
  TEUCHOS_ASSERT(integrationTimeDomain_.length() > ST::zero());
#endif
  
  StepControlInfo<Scalar> trialStepCtrlInfo;
  
  if (timeStepIter == 0)
    current_dt_ = initial_dt_;
  else
    current_dt_ *= ramping_factor_;

  if (timeStepIter < num_ramping_steps_) {
    trialStepCtrlInfo.stepType = STEP_TYPE_FIXED;
    trialStepCtrlInfo.stepSize = current_dt_;
  }
  else  {
    trialStepCtrlInfo.stepType = STEP_TYPE_VARIABLE;
    trialStepCtrlInfo.stepSize = max_dt_;
  }

  return trialStepCtrlInfo;
  
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
