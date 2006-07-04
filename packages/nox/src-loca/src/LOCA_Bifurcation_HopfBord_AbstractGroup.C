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

#include "Teuchos_ParameterList.hpp"
#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H"
#include "LOCA_ErrorCheck.H"

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::AbstractGroup::applyComplexInverseMulti(
			     Teuchos::ParameterList& params,
			     const NOX::Abstract::Vector* const* inputs_real,
			     const NOX::Abstract::Vector* const* inputs_imag,
			     double frequency,
			     NOX::Abstract::Vector** results_real,
			     NOX::Abstract::Vector** results_imag,
			     int nVecs) const
{
  string callingFunction = 
    "LOCA::Bifurcation::HopfBord::AbstractGroup::applyJacobianInverseMulti()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;
  
  for (int i=0; i<nVecs; i++) {
    status = applyComplexInverse(params, *inputs_real[i], *inputs_imag[i],
			      frequency, *results_real[i], *results_imag[i]);
    finalStatus = 
      LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						   callingFunction);
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::AbstractGroup::applyComplex(
		                      const NOX::Abstract::Vector& input_real,
				      const NOX::Abstract::Vector& input_imag,
				      double frequency,
				      NOX::Abstract::Vector& result_real,
				      NOX::Abstract::Vector& result_imag) const
{
  string callingFunction = 
    "LOCA::Bifurcation::HopfBord::AbstractGroup::applyComplex()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  NOX::Abstract::Vector *tmp = input_real.clone(NOX::ShapeCopy);

  // Compute J*y
  finalStatus = applyJacobian(input_real, result_real);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  // Compute B*z
  status = applyMassMatrix(input_imag, *tmp);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute J*y - w*B*z
  result_real.update(-frequency, *tmp, 1.0);

  // Compute J*z
  status = applyJacobian(input_imag, result_imag);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute B*y
  status = applyMassMatrix(input_real, *tmp);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute w*B*y + J*z
  result_imag.update(frequency, *tmp, 1.0);

  delete tmp;

  return finalStatus;
}
