// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_SingularJacobianSolve_ItRef.H"
#include "LOCA_ErrorCheck.H"

LOCA::SingularJacobianSolve::ItRef::ItRef(Teuchos::ParameterList& params)
{
  reset(params);
}

LOCA::SingularJacobianSolve::ItRef::ItRef(
			  const LOCA::SingularJacobianSolve::ItRef& source)
{
}

LOCA::SingularJacobianSolve::ItRef::~ItRef()
{
}

LOCA::SingularJacobianSolve::Generic*
LOCA::SingularJacobianSolve::ItRef::clone() const 
{
  return new ItRef(*this);
}

LOCA::SingularJacobianSolve::Generic&
LOCA::SingularJacobianSolve::ItRef::operator=(
			  const LOCA::SingularJacobianSolve::Generic& source)
{
  return operator=(dynamic_cast<const LOCA::SingularJacobianSolve::ItRef&>(source));
}

LOCA::SingularJacobianSolve::ItRef&
LOCA::SingularJacobianSolve::ItRef::operator=(
			  const LOCA::SingularJacobianSolve::ItRef& source)
{
  return *this;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::ItRef::reset(Teuchos::ParameterList& params) 
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::ItRef::compute(
				Teuchos::ParameterList& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& input,
			        const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result) 
{
  string callingFunction = 
    "LOCA::SingularJacobianSolve::ItRef::compute()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  finalStatus = grp.applyJacobianInverse(params, input, result);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  NOX::Abstract::Vector* remainder = input.clone(NOX::ShapeCopy);

  status = grp.applyJacobian(result, *remainder);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // r = b-Ax
  remainder->update(1.0, input, -1.0);

  NOX::Abstract::Vector* refinement = input.clone(NOX::ShapeCopy);

  // Ay=r
  status = grp.applyJacobianInverse(params, *remainder, *refinement);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // x+=y
  result.update(1.0, *refinement, 1.0);
  
  delete remainder;
  delete refinement;

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::ItRef::computeMulti(
				Teuchos::ParameterList& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector*const* inputs,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector** results,
				int nVecs) 
{
  string callingFunction = 
    "LOCA::SingularJacobianSolve::ItRef::computeMulti()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  NOX::Abstract::Vector** remainders = new NOX::Abstract::Vector*[nVecs];
  NOX::Abstract::Vector** refinements = new NOX::Abstract::Vector*[nVecs];

  finalStatus = grp.applyJacobianInverseMulti(params, inputs, results, nVecs);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  for (int i=0; i<nVecs; i++) {
    remainders[i] = inputs[i]->clone(NOX::ShapeCopy);
    refinements[i] = inputs[i]->clone(NOX::ShapeCopy);

    status = grp.applyJacobian(*(results[i]), *(remainders[i]));
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);

    // r = b-Ax
    remainders[i]->update(1.0, *(inputs[i]), -1.0);
  }

  // Ay=r
  status = grp.applyJacobianInverseMulti(params, remainders, refinements, 
					 nVecs);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // x+=y
  for (int i=0; i<nVecs; i++) {
    results[i]->update(1.0, *(refinements[i]), 1.0);
    delete remainders[i];
    delete refinements[i];
  }
  
  delete [] remainders;
  delete [] refinements;

  return finalStatus;
}
