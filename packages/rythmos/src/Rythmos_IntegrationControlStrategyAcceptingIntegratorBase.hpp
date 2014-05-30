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


#ifndef RYTHMOS_INTEGRATION_CONTROL_STRATEGY_ACCEPTING_INTEGRATOR_BASE_HPP
#define RYTHMOS_INTEGRATION_CONTROL_STRATEGY_ACCEPTING_INTEGRATOR_BASE_HPP


#include "Rythmos_Types.hpp"
#include "Rythmos_IntegratorBase.hpp"
#include "Rythmos_IntegrationControlStrategyBase.hpp"

namespace Rythmos {

/** \brief Mix-in interface for integrator objects that accept an
 * integration control strategy object to be used for evaluating steps and picking
 * step-sizes & orders.
 *
 * ToDo: Finish documentation!
 */

template<class Scalar>
class IntegrationControlStrategyAcceptingIntegratorBase : virtual public IntegratorBase<Scalar>
{
public:

  /** \brief . */
  virtual void setIntegrationControlStrategy(
      const RCP<IntegrationControlStrategyBase<Scalar> >& integrationControlStrategy
      ) = 0;

  /** \brief . */
  virtual RCP<IntegrationControlStrategyBase<Scalar> > 
    getNonconstIntegrationControlStrategy() = 0;

  /** \brief . */
  virtual RCP<const IntegrationControlStrategyBase<Scalar> > 
    getIntegrationControlStrategy() const = 0;

};


} // namespace Rythmos


#endif // RYTHMOS_INTEGRATION_CONTROL_STRATEGY_ACCEPTING_INTEGRATOR_BASE_HPP
