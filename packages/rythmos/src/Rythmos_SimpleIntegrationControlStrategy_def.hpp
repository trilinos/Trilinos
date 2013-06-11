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


#ifndef RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_DEF_HPP
#define RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_DEF_HPP

#include "Rythmos_SimpleIntegrationControlStrategy_decl.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"


namespace Rythmos {


template<class Scalar> 
RCP<SimpleIntegrationControlStrategy<Scalar> >
simpleIntegrationControlStrategy()
{
  RCP<SimpleIntegrationControlStrategy<Scalar> >
    integrationControl = Teuchos::rcp(new SimpleIntegrationControlStrategy<Scalar>());
  return integrationControl;
}


template<class Scalar> 
RCP<SimpleIntegrationControlStrategy<Scalar> >
simpleIntegrationControlStrategy( const RCP<ParameterList> &paramList )
{
  RCP<SimpleIntegrationControlStrategy<Scalar> >
    integrationControl = Teuchos::rcp(new SimpleIntegrationControlStrategy<Scalar>());
  integrationControl->setParameterList(paramList);
  return integrationControl;
}


//
// Implementation
//


// Static members


template<class Scalar> 
const std::string
SimpleIntegrationControlStrategy<Scalar>::takeVariableSteps_name_
= "Take Variable Steps";

template<class Scalar> 
const bool
SimpleIntegrationControlStrategy<Scalar>::takeVariableSteps_default_
= true;


template<class Scalar> 
const std::string
SimpleIntegrationControlStrategy<Scalar>::max_dt_name_
= "Max dt";

template<class Scalar> 
const double
SimpleIntegrationControlStrategy<Scalar>::max_dt_default_
= std::numeric_limits<Scalar>::max();


template<class Scalar> 
const std::string
SimpleIntegrationControlStrategy<Scalar>::numTimeSteps_name_
= "Number of Time Steps";

template<class Scalar> 
const int
SimpleIntegrationControlStrategy<Scalar>::numTimeSteps_default_
= -1;


template<class Scalar> 
const std::string
SimpleIntegrationControlStrategy<Scalar>::fixed_dt_name_
= "Fixed dt";

template<class Scalar> 
const double
SimpleIntegrationControlStrategy<Scalar>::fixed_dt_default_
= -1.0;


// Constructors/Initializers


template<class Scalar>
SimpleIntegrationControlStrategy<Scalar>::SimpleIntegrationControlStrategy()
  :takeVariableSteps_(takeVariableSteps_default_),
   max_dt_(max_dt_default_),
   numTimeSteps_(numTimeSteps_default_),
   fixed_dt_(fixed_dt_default_)
{}


// Overridden from ParameterListAcceptor


template<class Scalar> 
void SimpleIntegrationControlStrategy<Scalar>::setParameterList(
  RCP<ParameterList> const& paramList
  )
{
  using Teuchos::as;
  using Teuchos::get;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  TEUCHOS_TEST_FOR_EXCEPT(is_null(paramList));
  paramList->validateParameters(*getValidParameters());
  this->setMyParamList(paramList);
  takeVariableSteps_ = paramList->get(
    takeVariableSteps_name_, takeVariableSteps_ );
  if (!takeVariableSteps_) {
    numTimeSteps_ = paramList->get(numTimeSteps_name_,numTimeSteps_);
    fixed_dt_ = paramList->get(fixed_dt_name_,fixed_dt_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      numTimeSteps_ < 0 && fixed_dt_ <= ST::zero(), std::logic_error,
      "Error, when taking fixed steps, the user must set the parameters "
      "\""<<numTimeSteps_name_<<"\" > 0 or \""<<fixed_dt_name_<<"\" > 0.0!" );
  }
  else {
    max_dt_ = paramList->get(max_dt_name_,max_dt_);
  }
  Teuchos::readVerboseObjectSublist(&*paramList,this);
}


template<class Scalar> 
RCP<const ParameterList>
SimpleIntegrationControlStrategy<Scalar>::getValidParameters() const
{
  static RCP<const ParameterList> validPL;
  if (is_null(validPL) ) {
    RCP<ParameterList> pl = Teuchos::parameterList();
    pl->set(
      takeVariableSteps_name_, takeVariableSteps_default_,
      "Take variable time steps or fixed time steps.\n"
      "If set to false, then the parameter \"" + fixed_dt_name_ + "\"\n"
      "or \"" + numTimeSteps_name_ + "\" must be set!"
      );
    pl->set(
      max_dt_name_, max_dt_default_,
      "Gives the max size of the variable time steps.  This is only read and used if\n"
      "\"" + takeVariableSteps_name_ + "\" is set to true."
      );
    pl->set(
      numTimeSteps_name_, numTimeSteps_default_,
      "Gives the number of fixed time steps.  The actual step size gets computed\n"
      "on the fly given the size of the time domain.\n"
      "This is only read and used if \"" + takeVariableSteps_name_ + "\" is set to false\n"
      "and \"" + fixed_dt_name_ + "\" is set to < 0.0."
      );
    pl->set(
      fixed_dt_name_, fixed_dt_default_,
      "Gives the size of the fixed time steps.  This is only read and used if\n"
      "\"" + takeVariableSteps_name_ + "\" is set to false."
      );
    Teuchos::setupVerboseObjectSublist(&*pl);
    validPL = pl;
  }
  return validPL;
}


// Overridden from IntegrationControlStrategyBase


template<class Scalar>
RCP<IntegrationControlStrategyBase<Scalar> >
SimpleIntegrationControlStrategy<Scalar>::cloneIntegrationControlStrategy() const
{
  RCP<SimpleIntegrationControlStrategy<Scalar> >
    integrCtrlStry = simpleIntegrationControlStrategy<Scalar>();
  const RCP<const ParameterList> paramList = this->getParameterList();
  if (!is_null(paramList))
    integrCtrlStry->setParameterList(Teuchos::parameterList(*paramList));
  integrCtrlStry->takeVariableSteps_ = takeVariableSteps_;
  integrCtrlStry->max_dt_ = max_dt_;
  integrCtrlStry->numTimeSteps_ = numTimeSteps_;
  integrCtrlStry->fixed_dt_ = fixed_dt_;
  return integrCtrlStry;
}


template<class Scalar>
void
SimpleIntegrationControlStrategy<Scalar>::resetIntegrationControlStrategy(
  const TimeRange<Scalar> &integrationTimeDomain
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT(integrationTimeDomain.length() >= ST::zero());
#endif
  integrationTimeDomain_ = integrationTimeDomain;
  if (takeVariableSteps_) {
    if (max_dt_ < ST::zero()) {
      max_dt_ = integrationTimeDomain_.length();
    }
  }
  else {
    if (fixed_dt_ < ST::zero()) {
#ifdef HAVE_RYTHMOS_DEBUG
      TEUCHOS_ASSERT(numTimeSteps_ > 0);
#endif
      fixed_dt_ = integrationTimeDomain_.length()/numTimeSteps_;
    }
  }
}


template<class Scalar>
StepControlInfo<Scalar>
SimpleIntegrationControlStrategy<Scalar>::getNextStepControlInfo(
  const StepperBase<Scalar> &stepper,
  const StepControlInfo<Scalar> &stepCtrlInfoLast,
  const int timeStepIter
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

#ifdef HAVE_RYTHMOS_DEBUG
  TEUCHOS_ASSERT(integrationTimeDomain_.length() >= ST::zero());
#endif
  
  StepControlInfo<Scalar> trialStepCtrlInfo;
  
  if (takeVariableSteps_) {
    trialStepCtrlInfo.stepType = STEP_TYPE_VARIABLE;
    trialStepCtrlInfo.stepSize = max_dt_;
  }
  else {
    trialStepCtrlInfo.stepType = STEP_TYPE_FIXED;
    trialStepCtrlInfo.stepSize = fixed_dt_;
  }
  
  return trialStepCtrlInfo;
  
}

// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_INSTANT(SCALAR) \
  \
  template class SimpleIntegrationControlStrategy< SCALAR >; \
  \
  template RCP<SimpleIntegrationControlStrategy< SCALAR > > \
  simpleIntegrationControlStrategy(); \
  \
  template RCP<SimpleIntegrationControlStrategy< SCALAR > > \
  simpleIntegrationControlStrategy( const RCP<ParameterList> &paramList );
   

} // namespace Rythmos


#endif // RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_DEF_HPP
