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

#include "LOCA_Bifurcation_TPBord_AbstractGroup.H"
#include "LOCA_ErrorCheck.H"

double
LOCA::Bifurcation::TPBord::AbstractGroup::innerProduct(
				              const NOX::Abstract::Vector& x,
					      const NOX::Abstract::Vector& y) 
{
  return x.dot(y);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::AbstractGroup::applyBorderedJacobianInverseMulti(
			       bool trans,
			       NOX::Parameter::List& params,
			       const NOX::Abstract::Vector& a,
			       const NOX::Abstract::Vector& b,
			       const NOX::Abstract::Vector*const* vInputs,
			       double *sInputs,
			       NOX::Abstract::Vector** vResults,
			       double *sResults,
			       int nVecs) const
{
  string callingFunction = 
    "LOCA::Bifurcation::TPBord::AbstractGroup::applyBorderedJacobianInverseMulti()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<nVecs; i++) {
    status = applyBorderedJacobianInverse(trans, params, a, b, *vInputs[i], 
					  sInputs[i], *vResults[i], 
					  sResults[i]);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  return finalStatus;
}
