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

#include "NOX_Parameter_List.H"
#include "NOX_MultiVector.H"
#include "LOCA_Continuation_AbstractGroup.H"
#include "LOCA_ErrorCheck.H"

NOX::Abstract::Group::ReturnType
LOCA::Continuation::AbstractGroup::applyJacobianMultiVector(
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const
{
  string callingFunction = 
    "LOCA::Continuation::AbstractGroup::applyJacobian()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;
  
  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobian(input[i], result[i]);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::AbstractGroup::applyJacobianInverseMultiVector(
                                    NOX::Parameter::List& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const
{
  string callingFunction = 
    "LOCA::Continuation::AbstractGroup::applyJacobianInverse()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;
  
  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobianInverse(params, input[i], result[i]);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::AbstractGroup::applyJacobianInverseMulti(
			    NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** outputs, int nVecs) const
{
  string callingFunction = 
    "LOCA::Continuation::AbstractGroup::applyJacobianInverseMulti()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  for (int i=0; i<nVecs; i++) {
    status = applyJacobianInverse(params, *inputs[i], *outputs[i]);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  return finalStatus;
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
 LOCA::ErrorCheck::throwError(
		       "LOCA::Continuation::AbstractGroup::computeEigenvalues",
		       "No eigensolver defined for group");
  return NOX::Abstract::Group::Failed;
}

NOX::Abstract::Group::ReturnType
LOCA::Continuation::AbstractGroup::applyHouseholderJacobianInverse(
				         NOX::Parameter::List& params,
					 const NOX::Abstract::Vector& f,
					 const NOX::Abstract::Vector& dfdp,
					 const NOX::Abstract::Vector& ux,
					 double up, double beta,
					 NOX::Abstract::Vector& result_x,
					 double& result_p) const
{
  LOCA::ErrorCheck::throwError(
	  "LOCA::Continuation::AbstractGroup::applyHouseholderJacobianInverse",
	  "No implementation defined for group");
  return NOX::Abstract::Group::Failed;
}

void
LOCA::Continuation::AbstractGroup::scaleVector(NOX::Abstract::Vector& x) const
{
}
