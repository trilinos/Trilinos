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
#include "LOCA_SingularJacobianSolve_Default.H"

LOCA::SingularJacobianSolve::Default::Default(NOX::Parameter::List& params)
{
  reset(params);
}

LOCA::SingularJacobianSolve::Default::Default(
			  const LOCA::SingularJacobianSolve::Default& source)
{
}

LOCA::SingularJacobianSolve::Default::~Default()
{
}

LOCA::SingularJacobianSolve::Generic*
LOCA::SingularJacobianSolve::Default::clone() const 
{
  return new Default(*this);
}

LOCA::SingularJacobianSolve::Generic&
LOCA::SingularJacobianSolve::Default::operator=(
			  const LOCA::SingularJacobianSolve::Generic& source)
{
  return operator=(dynamic_cast<const LOCA::SingularJacobianSolve::Default&>(source));
}

LOCA::SingularJacobianSolve::Default&
LOCA::SingularJacobianSolve::Default::operator=(
			  const LOCA::SingularJacobianSolve::Default& source)
{
  return *this;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Default::reset(NOX::Parameter::List& params) 
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Default::compute(
				NOX::Parameter::List& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& input,
			        const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result) 
{
  return grp.applyJacobianInverse(params, input, result);
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Default::computeMulti(
				NOX::Parameter::List& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector*const* inputs,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector** results,
				int nVecs) 
{
  return grp.applyJacobianInverseMulti(params, inputs, results, nVecs);
}
