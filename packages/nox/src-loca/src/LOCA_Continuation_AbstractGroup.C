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

#include "NOX_Parameter_List.H"
#include "LOCA_Continuation_AbstractGroup.H"

NOX::Abstract::Group::ReturnType
LOCA::Continuation::AbstractGroup::applyJacobianInverseMulti(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** outputs, int nVecs) const
{
  NOX::Abstract::Group::ReturnType res;

  for (int i=0; i<nVecs; i++) {
    res = applyJacobianInverse(params, *inputs[i], *outputs[i]);
    if (res != NOX::Abstract::Group::Ok)
      return res;
  }

  return res;
}

double
LOCA::Continuation::AbstractGroup::computeScaledDotProduct(
					 const NOX::Abstract::Vector& a,
					 const NOX::Abstract::Vector& b) const
{
  return a.dot(b);
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::AbstractGroup::computeEigenvalues(
					        NOX::Parameter::List& params)
{
  errorCheck.throwError("LOCA::Continuation::AbstractGroup::computeEigenvalues",
			       "No eigensolver defined for group");
  return NOX::Abstract::Group::Failed;
}
