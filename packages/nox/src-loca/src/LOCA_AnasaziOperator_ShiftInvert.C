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

#include "LOCA_AnasaziOperator_ShiftInvert.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::AnasaziOperator::ShiftInvert::ShiftInvert(
	const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& eigenParams_,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams_,
	const Teuchos::RefCountPtr<LOCA::TimeDependent::AbstractGroup>& grp_)
  : globalData(global_data),
    myLabel("Shift-Invert"),
    eigenParams(eigenParams_),
    solverParams(solverParams_),
    grp(grp_),
    tmp_r(),
    tmp_i(),
    shift(0.0)
{
  string callingFunction = 
    "LOCA::AnasaziOperator::ShiftInvert::ShiftInvert()";

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Get parameters
  shift = eigenParams->get("Shift",0.0);

  // Compute Jacobian matrix
  status = grp->computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Compute mass matrix
  status = grp->computeMassMatrix();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
}

LOCA::AnasaziOperator::ShiftInvert::~ShiftInvert()
{
}

const string&
LOCA::AnasaziOperator::ShiftInvert::label() const
{
  return myLabel;
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::ShiftInvert::apply(
				     const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& output) const
{
  string callingFunction = 
    "LOCA::AnasaziOperator::ShiftInvert::apply()";

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  for (int i=0; i<input.numVectors(); i++) {

    // Allocate temporary vector
    if (tmp_r == Teuchos::null)
      tmp_r = input[i].clone(NOX::ShapeCopy);

    // Compute M*input
    status = grp->applyMassMatrix(input[i], *tmp_r);
    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);

    // Solve (J-omega*M)*output = M*input
    if (shift != 0.0) 
      status = grp->applyShiftedMatrixInverse(*solverParams, *tmp_r, 
					      output[i], 
					      -shift);
    else
      status = grp->applyJacobianInverse(*solverParams, *tmp_r, output[i]);

    finalStatus = 
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status, 
							     finalStatus,
							     callingFunction);
  }

  return finalStatus;
}

void
LOCA::AnasaziOperator::ShiftInvert::transformEigenvalue(double& ev_r, 
							double& ev_i) const
{
  // compute inverse of eigenvalue, then shift
  double mag = ev_r*ev_r + ev_i*ev_i;
  ev_r =  ev_r / mag + shift;
  ev_i = -ev_i / mag;
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::ShiftInvert::rayleighQuotient(
				         const NOX::Abstract::Vector& evec_r,
					 const NOX::Abstract::Vector& evec_i,
					 double& rq_r, double& rq_i) const
{
  string callingFunction = 
    "LOCA::AnasaziOperator::ShiftInvert::rayleighQuotient()";

  // Allocate temporary vectors
  if (tmp_r == Teuchos::null)
    tmp_r = evec_r.clone(NOX::ShapeCopy);
  if (tmp_i == Teuchos::null)
    tmp_i = evec_i.clone(NOX::ShapeCopy);

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure Jacobian is up-to-date
  status = grp->computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Compute z^T J z
  status = grp->applyJacobian(evec_r, *tmp_r);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  
  status = grp->applyJacobian(evec_i, *tmp_i);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  rq_r = evec_r.innerProduct(*tmp_r) + evec_i.innerProduct(*tmp_i);
  rq_i = evec_r.innerProduct(*tmp_i) - evec_i.innerProduct(*tmp_r);

  // Make sure mass matrix is up-to-date
  status = grp->computeMassMatrix();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  // Compute z^T M z
  status = grp->applyMassMatrix(evec_r, *tmp_r);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
  
  status = grp->applyMassMatrix(evec_i, *tmp_i);
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);

  double m_r = evec_r.innerProduct(*tmp_r) + evec_i.innerProduct(*tmp_i);
  double m_i = evec_r.innerProduct(*tmp_i) - evec_i.innerProduct(*tmp_r);
  double m = m_r*m_r + m_i*m_i;

  // Compute z^T J z / z^T M z
  rq_r = (rq_r*m_r + rq_i*m_i) / m;
  rq_i = (rq_i*m_r - rq_r*m_i) / m;

  return finalStatus;
}
