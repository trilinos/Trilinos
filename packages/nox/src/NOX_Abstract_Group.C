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

#include "NOX_Abstract_Group.H"

#ifdef HAVE_NOX_MULTIVECS
#include "NOX_Abstract_MultiVector.H"
#endif

NOX::Abstract::Group::ReturnType 
NOX::Abstract::Group::computeJacobian()
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType 
NOX::Abstract::Group::computeGradient()
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType 
NOX::Abstract::Group::computeNewton(NOX::Parameter::List& params)
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType 
NOX::Abstract::Group::applyJacobian(const NOX::Abstract::Vector& input, 
				    NOX::Abstract::Vector& result) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType 
NOX::Abstract::Group::applyJacobianTranspose(const NOX::Abstract::Vector& input, 
					     NOX::Abstract::Vector& result) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType 
NOX::Abstract::Group::applyJacobianInverse(NOX::Parameter::List& params, 
					   const NOX::Abstract::Vector& input, 
					   NOX::Abstract::Vector& result) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType 
NOX::Abstract::Group::applyRightPreconditioning(bool useTranspose,
						NOX::Parameter::List& params, 
						const NOX::Abstract::Vector& input, 
						NOX::Abstract::Vector& result
						) const
{
  return NOX::Abstract::Group::NotDefined;
}

#ifdef HAVE_NOX_MULTIVECS

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianMultiVector(
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;
  
  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobian(input[i], result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus == Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianTransposeMultiVector(
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;
  
  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobianTranspose(input[i], result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus == Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyJacobianInverseMultiVector(
                                    NOX::Parameter::List& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;
  
  for (int i=0; i<input.numVectors(); i++) {
    status = applyJacobianInverse(params, input[i], result[i]);
    
    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus == Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
NOX::Abstract::Group::applyRightPreconditioningMultiVector(
				   bool useTranspose,
				   NOX::Parameter::List& params,
				   const NOX::Abstract::MultiVector& input, 
				   NOX::Abstract::MultiVector& result) const
{
  NOX::Abstract::Group::ReturnType status, finalStatus;
  finalStatus = NOX::Abstract::Group::Ok;
  
  for (int i=0; i<input.numVectors(); i++) {
    status = applyRightPreconditioning(useTranspose, params, input[i],
				       result[i]);

    if (status == NotDefined || status == BadDependency)
      return status;
    else if (status == Failed)
      finalStatus == Failed;
    else if (status == NotConverged && finalStatus != Failed)
      finalStatus = NotConverged;
  }

  return finalStatus;
}

#endif

bool NOX::Abstract::Group::isJacobian() const
{
  return false;
}

bool NOX::Abstract::Group::isGradient() const
{
  return false;
}

bool NOX::Abstract::Group::isNewton() const
{
  return false;
}

NOX::Abstract::Group::ReturnType 
NOX::Abstract::Group::getNormLastLinearSolveResidual(double& residual) const
{
  return NOX::Abstract::Group::NotDefined;
}

