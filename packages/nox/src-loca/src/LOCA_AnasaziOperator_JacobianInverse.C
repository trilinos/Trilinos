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

#include "LOCA_AnasaziOperator_JacobianInverse.H"
#include "LOCA_ErrorCheck.H"

LOCA::AnasaziOperator::JacobianInverse::JacobianInverse(
					  NOX::Parameter::List& eigenParams_,
					  NOX::Parameter::List& solverParams_,
					  NOX::Abstract::Group& grp_)
  : myLabel("Jacobian Inverse"),
    tmp_r(NULL),
    tmp_i(NULL)
{
  reset(eigenParams_, solverParams_, grp_);
}

LOCA::AnasaziOperator::JacobianInverse::~JacobianInverse()
{
  if (tmp_r)
    delete tmp_r;
  if (tmp_i)
    delete tmp_i;
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::JacobianInverse::reset(
					  NOX::Parameter::List& eigenParams_,
					  NOX::Parameter::List& solverParams_,
					  NOX::Abstract::Group& grp_)
{
  eigenParams = &eigenParams_;
  solverParams = &solverParams_;
  grp = &grp_;

  if (tmp_r)
    delete tmp_r;
  if (tmp_i)
    delete tmp_i;

  // make sure Jacobian is up-to-date
  return grp->computeJacobian();
}

const string&
LOCA::AnasaziOperator::JacobianInverse::label() const
{
  return myLabel;
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::JacobianInverse::apply(
					 const NOX::Abstract::Vector& input, 
					 NOX::Abstract::Vector& output) const
{
  return grp->applyJacobianInverse(*solverParams, input, output);
}

void
LOCA::AnasaziOperator::JacobianInverse::transformEigenvalue(double& ev_r, 
							    double& ev_i) const
{
  // compute inverse of eigenvalue
  double mag = ev_r*ev_r + ev_i*ev_i;
  ev_r =  ev_r / mag;
  ev_i = -ev_i / mag;
}

NOX::Abstract::Group::ReturnType 
LOCA::AnasaziOperator::JacobianInverse::rayleighQuotient(
				         const NOX::Abstract::Vector& evec_r,
					 const NOX::Abstract::Vector& evec_i,
					 double& rq_r, double& rq_i) const
{
  string callingFunction = 
    "LOCA::AnasaziOperator::JacobianInverse::rayleighQuotient()";

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

  return finalStatus;
}
