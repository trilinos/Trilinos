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


#ifndef RYTHMOS_STEP_CONTROL_STRATEGY_BASE_HPP
#define RYTHMOS_STEP_CONTROL_STRATEGY_BASE_HPP


#include "Rythmos_StepperBase.hpp"

namespace Rythmos {


/** \brief step control strategy interface to be used for evaluating steps and picking
 * step-sizes & orders.
 *
 * ToDo: Finish documentation!
 */

/** \brief . */
enum AttemptedStepStatusFlag { PREDICT_AGAIN, CONTINUE_ANYWAY, REP_ERR_FAIL, REP_CONV_FAIL };

/** \brief . */
enum StepControlStrategyState { UNINITIALIZED, BEFORE_FIRST_STEP, MID_STEP, AFTER_CORRECTION, READY_FOR_NEXT_STEP };

/** \brief Convert StepControlStrategyState to string. */
inline
const char* toString( const StepControlStrategyState stepControlStrategyState )
{
  switch(stepControlStrategyState) {
    case UNINITIALIZED:
      return "UNINITIALIZED";
    case BEFORE_FIRST_STEP:
      return "BEFORE_FIRST_STEP";
    case MID_STEP:
      return "MID_STEP";
    case AFTER_CORRECTION:
      return "AFTER_CORRECTION";
    case READY_FOR_NEXT_STEP:
      return "READY_FOR_NEXT_STEP";
#ifdef HAVE_RYTHMOS_DEBUG
    default:
      TEUCHOS_TEST_FOR_EXCEPT("Invalid enum value!");
#endif
  }
  return 0; // Should never get here!
}

/** \brief The member functions in the StepControlStrategyBase move you
 * between these states in the following fashion:
 *
 * setRequestedStepSize:  valid in UNINITIALIZED and BEFORE_FIRST_STEP and READY_FOR_NEXT_STEP
 *   changes state:  UNINITIALIZED -> BEFORE_FIRST_STEP
 *
 * nextStepSize:  valid in BEFORE_FIRST_STEP and READY_FOR_NEXT_STEP
 *   changes state:  BEFORE_FIRST_STEP -> MID_STEP
 *                   READY_FOR_NEXT_STEP -> MID_STEP
 *
 * setCorrection:  valid in MID_STEP
 *   changes state:  MID_STEP -> AFTER_CORRECTION
 *
 * acceptStep:  valid in AFTER_CORRECTION
 *
 * completeStep:  valid in AFTER_CORRECTION
 *   changes state:  AFTER_CORRECTION -> READY_FOR_NEXT_STEP
 *
 * rejectStep:  valid in AFTER_CORRECTION
 *   changes state:  AFTER_CORRECTION -> READY_FOR_NEXT_STEP
 *
 */
template<class Scalar>
class StepControlStrategyBase
  : virtual public Teuchos::Describable
  , virtual public Teuchos::ParameterListAcceptor
  , virtual public Teuchos::VerboseObject<StepControlStrategyBase<Scalar> >
{
public:

  /** \brief . */
  virtual void initialize(const StepperBase<Scalar>& stepper) =0;

  /** \brief . */
  virtual void setRequestedStepSize(
      const StepperBase<Scalar>& stepper
      , const Scalar& stepSize
      , const StepSizeType& stepSizeType
      ) = 0;

  /** \brief . */
  virtual void nextStepSize(
      const StepperBase<Scalar>& stepper
      , Scalar* stepSize
      , StepSizeType* stepSizeType
      , int* order
      ) = 0;

  /** \brief . */
  virtual void setCorrection(
      const StepperBase<Scalar>& stepper
      , const RCP<const Thyra::VectorBase<Scalar> >& soln
      , const RCP<const Thyra::VectorBase<Scalar> >& ee
      , int solveStatus
      ) = 0;

  /** \brief . */
  virtual bool acceptStep(
      const StepperBase<Scalar>& stepper
      ,Scalar* LETValue
      ) = 0;

  /** \brief . */
  virtual void completeStep(
      const StepperBase<Scalar>& stepper
      ) = 0;

  /** \brief . */
  virtual AttemptedStepStatusFlag rejectStep(
      const StepperBase<Scalar>& stepper
      ) = 0;

  /** \brief . */
  virtual StepControlStrategyState getCurrentState() = 0;

  /** \brief . */
  virtual int getMaxOrder() const = 0;

  /** \brief . */
  virtual void setStepControlData(const StepperBase<Scalar>& stepper) = 0;

  /** \brief . */
  virtual bool supportsCloning() const;

  /** \brief . */
  virtual RCP<StepControlStrategyBase<Scalar> > cloneStepControlStrategyAlgorithm() const;


};

template<class Scalar>
bool StepControlStrategyBase<Scalar>::supportsCloning() const
{
  return false;
}


template<class Scalar>
RCP<StepControlStrategyBase<Scalar> >
StepControlStrategyBase<Scalar>::cloneStepControlStrategyAlgorithm() const
{
  return Teuchos::null;
}


} // namespace Rythmos


#endif // RYTHMOS_STEP_CONTROL_STRATEGY_BASE_HPP
