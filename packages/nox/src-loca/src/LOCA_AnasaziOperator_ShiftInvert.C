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

#include "LOCA_AnasaziOperator_ShiftInvert.H"
#include "NOX_Parameter_List.H"
#include "LOCA_ErrorCheck.H"

LOCA::AnasaziOperator::ShiftInvert::ShiftInvert(
					  NOX::Parameter::List& eigenParams_,
					  NOX::Parameter::List& solverParams_,
					  NOX::Abstract::Group& grp_)
  : myLabel("Shift-Invert"),
    tmp_r(NULL),
    tmp_i(NULL)
{
  reset(eigenParams_, solverParams_, grp_);
}

LOCA::AnasaziOperator::ShiftInvert::~ShiftInvert()
{
  if (tmp_r)
    delete tmp_r;
  if (tmp_i)
    delete tmp_i;
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::ShiftInvert::reset(
					  NOX::Parameter::List& eigenParams_,
					  NOX::Parameter::List& solverParams_,
					  NOX::Abstract::Group& grp_)
{
  string callingFunction = 
    "LOCA::AnasaziOperator::ShiftInvert::reset()";

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  eigenParams = &eigenParams_;
  solverParams = &solverParams_;
  grp = dynamic_cast<LOCA::TimeDependent::AbstractGroup*>(&grp_);

  // ensure grp is of the right type
  if (grp == NULL) {
    LOCA::ErrorCheck::throwError(callingFunction,
     "Supplied group is not derived from LOCA::TimeDependent::AbstractGroup!");
  }

  // Get parameters
  shift = eigenParams->getParameter("Shift",0.0);

  if (tmp_r)
    delete tmp_r;
  if (tmp_i)
    delete tmp_i;

  // Compute Jacobian matrix
  status = grp->computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute mass matrix
  status = grp->computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  return finalStatus;
}

const string&
LOCA::AnasaziOperator::ShiftInvert::label() const
{
  return myLabel;
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::ShiftInvert::apply(
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& output) const
{
  string callingFunction = 
    "LOCA::AnasaziOperator::ShiftInvert::apply()";

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Allocate temporary vector
  if (!tmp_r)
    tmp_r = input.clone(NOX::ShapeCopy);

  // Compute M*input
  status = grp->applyMassMatrix(input, *tmp_r);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Solve (J-omega*M)*output = M*input
  if (shift != 0.0) 
    status = grp->applyShiftedMatrixInverse(*solverParams, *tmp_r, output, 
					    -shift);
  else
    status = grp->applyJacobianInverse(*solverParams, *tmp_r, output);

  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
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
  if (!tmp_r)
    tmp_r = evec_r.clone(NOX::ShapeCopy);
  if (!tmp_i)
    tmp_i = evec_i.clone(NOX::ShapeCopy);

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Make sure Jacobian is up-to-date
  status = grp->computeJacobian();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute z^T J z
  status = grp->applyJacobian(evec_r, *tmp_r);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  status = grp->applyJacobian(evec_i, *tmp_i);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  rq_r = evec_r.dot(*tmp_r) + evec_i.dot(*tmp_i);
  rq_i = evec_r.dot(*tmp_i) - evec_i.dot(*tmp_r);

  // Make sure mass matrix is up-to-date
  status = grp->computeMassMatrix();
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  // Compute z^T M z
  status = grp->applyMassMatrix(evec_r, *tmp_r);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);
  
  status = grp->applyMassMatrix(evec_i, *tmp_i);
  finalStatus = 
    LOCA::ErrorCheck::combineAndCheckReturnTypes(status, finalStatus,
						 callingFunction);

  double m_r = evec_r.dot(*tmp_r) + evec_i.dot(*tmp_i);
  double m_i = evec_r.dot(*tmp_i) - evec_i.dot(*tmp_r);
  double m = m_r*m_r + m_i*m_i;

  // Compute z^T J z / z^T M z
  rq_r = (rq_r*m_r + rq_i*m_i) / m;
  rq_i = (rq_i*m_r - rq_r*m_i) / m;

  return finalStatus;
}
