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


#ifndef RYTHMOS_STEP_CONTROL_STRATEGY_ACCEPTING_STEPPER_BASE_HPP
#define RYTHMOS_STEP_CONTROL_STRATEGY_ACCEPTING_STEPPER_BASE_HPP


#include "Rythmos_StepperBase.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"

namespace Rythmos {


/** \brief Mix-in interface for stepper objects that accept a
 * step control strategy object to be used for evaluating steps and picking
 * step-sizes & orders.
 *
 * ToDo: Finish documentation!
 */

template<class Scalar>
class StepControlStrategyAcceptingStepperBase : virtual public StepperBase<Scalar>
{
public:

  /** \brief . */
  virtual void setStepControlStrategy(
      const Teuchos::RCP<StepControlStrategyBase<Scalar> >& stepControlStrategy
      ) = 0;

  /** \brief . */
  virtual Teuchos::RCP<StepControlStrategyBase<Scalar> > 
    getNonconstStepControlStrategy() = 0;

  /** \brief . */
  virtual Teuchos::RCP<const StepControlStrategyBase<Scalar> > 
    getStepControlStrategy() const = 0;

};


} // namespace Rythmos


#endif // RYTHMOS_STEP_CONTROL_STRATEGY_ACCEPTING_STEPPER_BASE_HPP
