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

#include "NOX_Belos_PreconditionOperator.H"
#include "NOX_Belos_MultiVector.H"
#include "NOX_Abstract_MultiVector.H"
#include "NOX_Parameter_List.H"

NOX::Belos::PreconditionOperator::PreconditionOperator(
			       NOX::Abstract::Group& g,
			       NOX::Parameter::List& preconditionerParameters)
  : grp(g),
    precondParams(preconditionerParameters)
{
}


NOX::Belos::PreconditionOperator::~PreconditionOperator() 
{
}

::Belos::ReturnType
NOX::Belos::PreconditionOperator::Apply(
				 const ::Belos::MultiVec<double>& x, 
				 ::Belos::MultiVec<double>& y, 
				 ::Belos::ETrans trans) const
{
  // Cast x and y to NOX::Belos::MultiVec's
  const NOX::Belos::MultiVector& nox_belos_x = 
    dynamic_cast<const NOX::Belos::MultiVector&>(x);
  NOX::Belos::MultiVector& nox_belos_y = 
    dynamic_cast<NOX::Belos::MultiVector&>(y);

  // Get underlying NOX::Abstract::MultiVector's
  const NOX::Abstract::MultiVector& nox_x = nox_belos_x.getNoxMultiVector();
  NOX::Abstract::MultiVector& nox_y = nox_belos_y.getNoxMultiVector();

  // NOX return type
  NOX::Abstract::Group::ReturnType nox_status;

  bool useTranspose = false;
  if (trans == ::Belos::TRANS)
    useTranspose = true;
  
  nox_status = grp.applyRightPreconditioningMultiVector(useTranspose,
							precondParams, 
							nox_x, 
							nox_y);

  return noxReturnTypeToBelos(nox_status);
}

::Belos::ReturnType
NOX::Belos::PreconditionOperator::ApplyInverse(
				 const ::Belos::MultiVec<double>& x, 
				 ::Belos::MultiVec<double>& y, 
				 ::Belos::ETrans trans) const
{
  return ::Belos::Undefined;
}

::Belos::ReturnType
NOX::Belos::PreconditionOperator::noxReturnTypeToBelos(
			    NOX::Abstract::Group::ReturnType noxStatus) const
{
  if (noxStatus == NOX::Abstract::Group::Ok ||
      noxStatus == NOX::Abstract::Group::NotConverged)
    return ::Belos::Ok;
  else if (noxStatus == NOX::Abstract::Group::NotDefined)
    return ::Belos::Undefined;
  else
    return ::Belos::Error;
}

  
