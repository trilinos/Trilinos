// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_SingularJacobianSolve_ItRef.H"

LOCA::SingularJacobianSolve::ItRef::ItRef(NOX::Parameter::List& params)
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
LOCA::SingularJacobianSolve::ItRef::reset(NOX::Parameter::List& params) 
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::ItRef::compute(
				NOX::Parameter::List& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& input,
			        const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result) 
{
  NOX::Abstract::Group::ReturnType 
    res = grp.applyJacobianInverse(params, input, result);

  NOX::Abstract::Vector* remainder = input.clone(NOX::ShapeCopy);

  res = grp.applyJacobian(result, *remainder);

  // r = b-Ax
  remainder->update(1.0, input, -1.0);

  NOX::Abstract::Vector* refinement = input.clone(NOX::ShapeCopy);

  // Ay=r
  res = grp.applyJacobianInverse(params, *remainder, *refinement);

  // x+=y
  result.update(1.0, *refinement, 1.0);
  
  delete remainder;
  delete refinement;

  return res;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::ItRef::computeMulti(
				NOX::Parameter::List& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector*const* inputs,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector** results,
				int nVecs) 
{
  NOX::Abstract::Vector** remainders = new NOX::Abstract::Vector*[nVecs];
  NOX::Abstract::Vector** refinements = new NOX::Abstract::Vector*[nVecs];

  NOX::Abstract::Group::ReturnType 
    res = grp.applyJacobianInverseMulti(params, inputs, results, nVecs);

  for (int i=0; i<nVecs; i++) {
    remainders[i] = inputs[i]->clone(NOX::ShapeCopy);
    refinements[i] = inputs[i]->clone(NOX::ShapeCopy);

    res = grp.applyJacobian(*(results[i]), *(remainders[i]));

    // r = b-Ax
    remainders[i]->update(1.0, *(inputs[i]), -1.0);
  }

  // Ay=r
  res = grp.applyJacobianInverseMulti(params, remainders, refinements, nVecs);

  // x+=y
  for (int i=0; i<nVecs; i++) {
    results[i]->update(1.0, *(refinements[i]), 1.0);
    delete remainders[i];
    delete refinements[i];
  }
  
  delete [] remainders;
  delete [] refinements;

  return res;
}
