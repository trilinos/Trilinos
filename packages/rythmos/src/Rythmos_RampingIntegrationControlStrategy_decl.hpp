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


#ifndef RYTHMOS_RAMPING_INTEGRATION_CONTROL_STRATEGY_DECL_HPP
#define RYTHMOS_RAMPING_INTEGRATION_CONTROL_STRATEGY_DECL_HPP


#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_StepControlInfo.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Rythmos {


/** \brief Controls inital ramping at a fixed or incrementing time step size
 *
 * This integration control strategy starts a run by taking a fixed step size for a number of time steps (can optionally ramp the size).  It then turns over control to a variable time stepper.
 */
template<class Scalar>
class RampingIntegrationControlStrategy
  : virtual public IntegrationControlStrategyBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \brief Constructors/Initializers. */
  //@{

  /** \brief . */
  RampingIntegrationControlStrategy();

  //@}

  /** \name Overridden from ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);

  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \brief Overridden from IntegrationControlStrategyBase */
  //@{

  /** \brief . */
  RCP<IntegrationControlStrategyBase<Scalar> >
  cloneIntegrationControlStrategy() const;

  /** \brief . */
  void resetIntegrationControlStrategy(
    const TimeRange<Scalar> &integrationTimeDomain
    );

  /** \brief . */
  StepControlInfo<Scalar>
  getNextStepControlInfo(
    const StepperBase<Scalar> &stepper,
    const StepControlInfo<Scalar> &stepCtrlInfoLast,
    const int timeStepIter
    );

  //@}

private:

  int num_ramping_steps_;
  Scalar initial_dt_;
  Scalar max_dt_;
  Scalar ramping_factor_;

  Scalar current_dt_;

  TimeRange<Scalar> integrationTimeDomain_;

  static const std::string num_ramping_steps_name_;
  static const int num_ramping_steps_default_;

  static const std::string initial_dt_name_;
  static const double initial_dt_default_;

  static const std::string max_dt_name_;
  static const double max_dt_default_;

  static const std::string ramping_factor_name_;
  static const double ramping_factor_default_;

};


/** \brief .
 *
 * \relates RampingIntegrationControlStrategy
 */
template<class Scalar> 
RCP<RampingIntegrationControlStrategy<Scalar> >
rampingIntegrationControlStrategy();


/** \brief .
 *
 * \relates RampingIntegrationControlStrategy
 */
template<class Scalar> 
RCP<RampingIntegrationControlStrategy<Scalar> >
rampingIntegrationControlStrategy( const RCP<ParameterList> &paramList );



} // namespace Rythmos


#endif // RYTHMOS_RAMPING_INTEGRATION_CONTROL_STRATEGY_DECL_HPP
