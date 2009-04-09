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

#include "Rythmos_Charon_ImplicitBDFStepperErrWtVecCalc.hpp"
#include "Rythmos_ErrWtVecCalcBase.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

namespace RythmosCharon {

CharonImplicitBDFStepperErrWtVecCalc::

CharonImplicitBDFStepperErrWtVecCalc::CharonImplicitBDFStepperErrWtVecCalc() { }
CharonImplicitBDFStepperErrWtVecCalc::~CharonImplicitBDFStepperErrWtVecCalc() { }
void CharonImplicitBDFStepperErrWtVecCalc::errWtVecSet(
    Thyra::VectorBase<double>* weight, 
    const Thyra::VectorBase<double>& vector, 
    double relTol, 
    double absTol
    ) const
{ }
// Overridden from Teuchos::ParameterListAcceptor
void CharonImplicitBDFStepperErrWtVecCalc::setParameterList( Teuchos::RCP<Teuchos::ParameterList> const& paramList )
{ }
Teuchos::RCP<Teuchos::ParameterList> CharonImplicitBDFStepperErrWtVecCalc::getNonconstParameterList()
{ 
  return Teuchos::parameterList();
}
Teuchos::RCP<Teuchos::ParameterList> CharonImplicitBDFStepperErrWtVecCalc::unsetParameterList()
{
  return Teuchos::parameterList();
}
Teuchos::RCP<const Teuchos::ParameterList> CharonImplicitBDFStepperErrWtVecCalc::getValidParameters() const
{
  return Teuchos::parameterList();
}


} // namespace RythmosCharon


