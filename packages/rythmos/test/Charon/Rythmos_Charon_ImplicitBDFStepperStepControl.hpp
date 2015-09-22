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

#ifndef RYTHMOS_CHARON_IMPLICIT_BDF_STEPPER_STEP_CONTROL_HPP
#define RYTHMOS_CHARON_IMPLICIT_BDF_STEPPER_STEP_CONTROL_HPP

#include "Rythmos_Charon_ImplicitBDFStepperStepControlDecl.hpp"
#include "Rythmos_Charon_ImplicitBDFStepperErrWtVecCalc.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"

namespace RythmosCharon {

template<class Scalar>
CharonImplicitBDFStepperStepControl<Scalar>::CharonImplicitBDFStepperStepControl() { }
template<class Scalar>
CharonImplicitBDFStepperStepControl<Scalar>::~CharonImplicitBDFStepperStepControl() { }
template<class Scalar>
void CharonImplicitBDFStepperStepControl<Scalar>::setErrWtVecCalc(const Teuchos::RCP<Rythmos::ErrWtVecCalcBase<Scalar> >& errWtVecCalc)
{ }
template<class Scalar>
Teuchos::RCP<const Rythmos::ErrWtVecCalcBase<Scalar> > CharonImplicitBDFStepperStepControl<Scalar>::getErrWtVecCalc() const
{
  return Teuchos::null;
}
template<class Scalar>
void CharonImplicitBDFStepperStepControl<Scalar>::initialize(const Rythmos::StepperBase<Scalar>& stepper)
{ }
template<class Scalar>
void CharonImplicitBDFStepperStepControl<Scalar>::setRequestedStepSize(
    const Rythmos::StepperBase<Scalar>& stepper
    , const Scalar& stepSize
    , const Rythmos::StepSizeType& stepSizeType
    )
{ }
template<class Scalar>
void CharonImplicitBDFStepperStepControl<Scalar>::nextStepSize(
    const Rythmos::StepperBase<Scalar>& stepper
    , Scalar* stepSize
    , Rythmos::StepSizeType* stepSizeType
    , int* order 
    )
{ }
template<class Scalar>
void CharonImplicitBDFStepperStepControl<Scalar>::setCorrection(
    const Rythmos::StepperBase<Scalar>& stepper
    , const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& soln
    , const Teuchos::RCP<const Thyra::VectorBase<Scalar> >& ee
    , int solveStatus
    )
{ }
template<class Scalar>
bool CharonImplicitBDFStepperStepControl<Scalar>::acceptStep(
    const Rythmos::StepperBase<Scalar>& stepper
    ,Scalar* LETValue
    )
{ 
  return false;
}
template<class Scalar>
void CharonImplicitBDFStepperStepControl<Scalar>::completeStep(
    const Rythmos::StepperBase<Scalar>& stepper
    )
{ }
template<class Scalar>
Rythmos::AttemptedStepStatusFlag CharonImplicitBDFStepperStepControl<Scalar>::rejectStep(
    const Rythmos::StepperBase<Scalar>& stepper
    )
{ 
  return Rythmos::REP_ERR_FAIL;
}
template<class Scalar>
Rythmos::StepControlStrategyState CharonImplicitBDFStepperStepControl<Scalar>::getCurrentState()
{
  return Rythmos::UNINITIALIZED;
}
template<class Scalar>
int CharonImplicitBDFStepperStepControl<Scalar>::getMaxOrder() const
{ 
  return 0;
}
template<class Scalar>
void CharonImplicitBDFStepperStepControl<Scalar>::setStepControlData(const Rythmos::StepperBase<Scalar>& stepper)
{ }
template<class Scalar>
bool CharonImplicitBDFStepperStepControl<Scalar>::supportsCloning() const
{ 
  return false;
}
template<class Scalar>
Teuchos::RCP<Rythmos::StepControlStrategyBase<Scalar> > CharonImplicitBDFStepperStepControl<Scalar>::cloneStepControlStrategyAlgorithm() const
{
  return Teuchos::null;
}
// Overridden from Teuchos::ParameterListAcceptor
template<class Scalar>
void CharonImplicitBDFStepperStepControl<Scalar>::setParameterList( Teuchos::RCP<Teuchos::ParameterList> const& paramList )
{ }
template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList> CharonImplicitBDFStepperStepControl<Scalar>::getNonconstParameterList()
{ 
  return Teuchos::parameterList();
}
template<class Scalar>
Teuchos::RCP<Teuchos::ParameterList> CharonImplicitBDFStepperStepControl<Scalar>::unsetParameterList()
{
  return Teuchos::parameterList();
}
template<class Scalar>
Teuchos::RCP<const Teuchos::ParameterList> CharonImplicitBDFStepperStepControl<Scalar>::getValidParameters() const
{
  return Teuchos::parameterList();
}

} // namespace RythmosCharon

#endif // RYTHMOS_CHARON_IMPLICIT_BDF_STEPPER_STEP_CONTROL_HPP

