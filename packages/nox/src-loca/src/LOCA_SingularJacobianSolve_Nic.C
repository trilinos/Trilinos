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
#include "LOCA_SingularJacobianSolve_Nic.H"
#include "LOCA_ErrorCheck.H"

LOCA::SingularJacobianSolve::Nic::Nic(Teuchos::ParameterList& params)
{
  reset(params);
}

LOCA::SingularJacobianSolve::Nic::Nic(
			  const LOCA::SingularJacobianSolve::Nic& source)
{
}

LOCA::SingularJacobianSolve::Nic::~Nic()
{
}

LOCA::SingularJacobianSolve::Generic*
LOCA::SingularJacobianSolve::Nic::clone() const 
{
  return new Nic(*this);
}

LOCA::SingularJacobianSolve::Generic&
LOCA::SingularJacobianSolve::Nic::operator=(
			  const LOCA::SingularJacobianSolve::Generic& source)
{
  return operator=(dynamic_cast<const LOCA::SingularJacobianSolve::Nic&>(source));
}

LOCA::SingularJacobianSolve::Nic&
LOCA::SingularJacobianSolve::Nic::operator=(
			  const LOCA::SingularJacobianSolve::Nic& source)
{
  return *this;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Nic::reset(Teuchos::ParameterList& params) 
{
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Nic::compute(
				Teuchos::ParameterList& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& input,
			        const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result) 
{
  string callingFunction = 
    "LOCA::SingularJacobianSolve::Nic::compute()";
  NOX::Abstract::Group::ReturnType finalStatus;

  double alpha = approxNullVec.innerProduct(input)
               / approxNullVec.innerProduct(jacApproxNullVec);

  NOX::Abstract::Vector* tmpInput  = input.clone(NOX::DeepCopy);
  tmpInput->update(-alpha, jacApproxNullVec, 1.0);

  finalStatus = grp.applyJacobianInverse(params, *tmpInput, result);
  LOCA::ErrorCheck::checkReturnType(finalStatus, callingFunction);

  delete tmpInput;

  result.update(alpha, approxNullVec, 1.0);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Nic::computeMulti(
				Teuchos::ParameterList& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector*const* inputs,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector** results,
				int nVecs) 
{
  string callingFunction = 
    "LOCA::SingularJacobianSolve::Nic::computeMulti()";
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;

  double denom = approxNullVec.innerProduct(jacApproxNullVec);

  double* alphas = new double[nVecs];
  NOX::Abstract::Vector** tmpInputs  = new NOX::Abstract::Vector*[nVecs];

  for (int i=0; i<nVecs; i++) {
    alphas[i] = approxNullVec.innerProduct(*(inputs[i])) / denom;
    tmpInputs[i] = inputs[i]->clone(NOX::DeepCopy);
    tmpInputs[i]->update(-alphas[i], jacApproxNullVec, 1.0);
  }

  status = grp.applyJacobianInverseMulti(params, tmpInputs, results, nVecs);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  for (int i=0; i<nVecs; i++) {
    results[i]->update(alphas[i], approxNullVec, 1.0);
    delete tmpInputs[i];
  }

  delete [] tmpInputs;
  delete [] alphas;

  return finalStatus;
}
