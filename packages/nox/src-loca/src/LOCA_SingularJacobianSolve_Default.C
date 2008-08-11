// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_SingularJacobianSolve_Default.H"

LOCA::SingularJacobianSolve::Default::Default(Teuchos::ParameterList& params)
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
LOCA::SingularJacobianSolve::Default::reset(Teuchos::ParameterList& params) 
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Default::compute(
				Teuchos::ParameterList& params,
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
				Teuchos::ParameterList& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector*const* inputs,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector** results,
				int nVecs) 
{
  return grp.applyJacobianInverseMulti(params, inputs, results, nVecs);
}
