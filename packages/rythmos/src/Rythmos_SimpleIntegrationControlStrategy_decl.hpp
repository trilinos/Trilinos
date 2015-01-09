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


#ifndef RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_DECL_HPP
#define RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_DECL_HPP


#include "Rythmos_IntegrationControlStrategyBase.hpp"
#include "Rythmos_StepperBase.hpp"
#include "Rythmos_StepControlInfo.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"


namespace Rythmos {


/** \brief Base class for strategy objects that control integration by
 * selecting step sizes for a stepper.
 *
 * ToDo: Finish Implementation!
 */
template<class Scalar>
class SimpleIntegrationControlStrategy
  : virtual public IntegrationControlStrategyBase<Scalar>,
    virtual public Teuchos::ParameterListAcceptorDefaultBase
{
public:

  /** \brief Constructors/Initializers. */
  //@{

  /** \brief . */
  SimpleIntegrationControlStrategy();

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

  bool takeVariableSteps_;
  Scalar max_dt_;
  int numTimeSteps_;
  Scalar fixed_dt_;

  TimeRange<Scalar> integrationTimeDomain_;

  static const std::string takeVariableSteps_name_;
  static const bool takeVariableSteps_default_;

  static const std::string max_dt_name_;
  static const double max_dt_default_;

  static const std::string numTimeSteps_name_;
  static const int numTimeSteps_default_;

  static const std::string fixed_dt_name_;
  static const double fixed_dt_default_;

};


/** \brief .
 *
 * \relates SimpleIntegrationControlStrategy
 */
template<class Scalar> 
RCP<SimpleIntegrationControlStrategy<Scalar> >
simpleIntegrationControlStrategy();


/** \brief .
 *
 * \relates SimpleIntegrationControlStrategy
 */
template<class Scalar> 
RCP<SimpleIntegrationControlStrategy<Scalar> >
simpleIntegrationControlStrategy( const RCP<ParameterList> &paramList );



} // namespace Rythmos


#endif // RYTHMOS_SIMPLE_INTEGRATION_CONTROL_STRATEGY_DECL_HPP
