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

#ifndef Rythmos_STEPPER_BASE_DEF_H
#define Rythmos_STEPPER_BASE_DEF_H

#include "Rythmos_StepperBase_decl.hpp"

namespace Rythmos {

/** \brief .
 *
 * \relates StepperBase
 */
template<class Scalar>
bool isInitialized( const StepperBase<Scalar>& stepper )
{
  return stepper.getTimeRange().isValid();
}


// ///////////////////////////////
// Implementations


template<class Scalar>
bool StepperBase<Scalar>::isImplicit() const
{
  return false;
}


template<class Scalar>
bool StepperBase<Scalar>::acceptsModel() const
{
  return true;
}


template<class Scalar>
bool StepperBase<Scalar>::supportsCloning() const
{
  return false;
}


template<class Scalar>
RCP<StepperBase<Scalar> >
StepperBase<Scalar>::cloneStepperAlgorithm() const
{
  return Teuchos::null;
}


template<class Scalar>
void StepperBase<Scalar>::setInitialCondition(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
  )
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "Error, the function setIntialCondition(...) is not implemented\n"
    "in the class \"" << Teuchos::typeName(*this) << "\"!" );
  // ToDo: Remove this default implementation and make every concrete
  // subclass implement this!
}

template<class Scalar>
void StepperBase<Scalar>::setStepControlData(const StepperBase & stepper)
{
}


// 
// Explicit Instantiation macro
//
// Must be expanded from within the Rythmos namespace!
//

#define RYTHMOS_STEPPER_BASE_INSTANT(SCALAR) \
  \
  template class StepperBase< SCALAR >; \
  \
  template bool isInitialized( const StepperBase< SCALAR >& stepper );
   

} // namespace Rythmos

#endif //Rythmos_STEPPER_BASE_DEF_H
