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
#include "LOCA_Bifurcation_HopfBord_AbstractGroup.H"

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::AbstractGroup::applyComplexInverseMulti(
			     NOX::Parameter::List& params,
			     const NOX::Abstract::Vector* const* inputs_real,
			     const NOX::Abstract::Vector* const* inputs_imag,
			     double frequency,
			     NOX::Abstract::Vector** results_real,
			     NOX::Abstract::Vector** results_imag,
			     int nVecs) const
{
  NOX::Abstract::Group::ReturnType res;

  for (int i=0; i<nVecs; i++) {
    res = applyComplexInverse(params, *inputs_real[i], *inputs_imag[i],
			      frequency, *results_real[i], *results_imag[i]);
    if (res != NOX::Abstract::Group::Ok)
      return res;
  }

  return res;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::AbstractGroup::applyComplex(
		                      const NOX::Abstract::Vector& input_real,
				      const NOX::Abstract::Vector& input_imag,
				      double frequency,
				      NOX::Abstract::Vector& result_real,
				      NOX::Abstract::Vector& result_imag) const
{
  NOX::Abstract::Group::ReturnType res;

  NOX::Abstract::Vector *tmp = input_real.clone(NOX::ShapeCopy);

  // Compute J*y
  res = applyJacobian(input_real, result_real);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute B*z
  res = applyMassMatrix(input_imag, *tmp);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute J*y - w*B*z
  result_real.update(-frequency, *tmp, 1.0);

  // Compute J*z
  res = applyJacobian(input_imag, result_imag);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute B*y
  res = applyMassMatrix(input_real, *tmp);
  if (res != NOX::Abstract::Group::Ok)
    return res;

  // Compute w*B*y + J*z
  result_imag.update(frequency, *tmp, 1.0);

  delete tmp;

  return NOX::Abstract::Group::Ok;
}
