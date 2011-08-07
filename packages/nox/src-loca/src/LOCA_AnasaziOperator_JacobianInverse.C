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

#include "LOCA_AnasaziOperator_JacobianInverse.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "Teuchos_ParameterList.hpp"

LOCA::AnasaziOperator::JacobianInverse::JacobianInverse(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams_,
	const Teuchos::RCP<Teuchos::ParameterList>& solverParams_,
	const Teuchos::RCP<NOX::Abstract::Group>& grp_)
  : globalData(global_data),
    myLabel("Jacobian Inverse"),
    eigenParams(eigenParams_),
    solverParams(solverParams_),
    grp(grp_),
    tmp_r(),
    tmp_i()
{
  string callingFunction = 
    "LOCA::AnasaziOperator::JacobianInverse::JacobianInverse()";

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // make sure Jacobian is up-to-date
  status = grp->computeJacobian();
  finalStatus = 
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
							   callingFunction);
}

LOCA::AnasaziOperator::JacobianInverse::~JacobianInverse()
{
}

const string&
LOCA::AnasaziOperator::JacobianInverse::label() const
{
  return myLabel;
}

void
LOCA::AnasaziOperator::JacobianInverse::apply(
				      const NOX::Abstract::MultiVector& input, 
				      NOX::Abstract::MultiVector& output) const
{
  NOX::Abstract::Group::ReturnType status = 
    grp->applyJacobianInverseMultiVector(*solverParams, input, output);
  globalData->locaErrorCheck->checkReturnType(status, 
		       "LOCA::AnasaziOperator::JacobianInverse::apply()");
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
				         NOX::Abstract::Vector& evec_r,
					 NOX::Abstract::Vector& evec_i,
					 double& rq_r, double& rq_i) const
{
  string callingFunction = 
    "LOCA::AnasaziOperator::JacobianInverse::rayleighQuotient()";

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

  return finalStatus;
}
