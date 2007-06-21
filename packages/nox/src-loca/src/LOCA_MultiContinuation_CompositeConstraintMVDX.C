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

#include "LOCA_MultiContinuation_CompositeConstraintMVDX.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::MultiContinuation::CompositeConstraintMVDX::CompositeConstraintMVDX(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const vector< Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterfaceMVDX> >& constraintObjects) :
  LOCA::MultiContinuation::CompositeConstraint(),
  constraintMVDXPtrs(constraintObjects),
  compositeDX()
{
  // Copy constraint object pointers into temporary array of base class
  vector<Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> > tmp(constraintObjects.size());
  for (unsigned int i=0; i<constraintObjects.size(); i++)
    tmp[i] = constraintObjects[i];

  // Initialize parent class
  init(global_data, tmp);

  // Find first constraint that has nonzero DX
  int i=0;
  while (i < numConstraintObjects && constraintPtrs[i]->isDXZero())
    i++;

  // Create new multivector to store all constraint derivatives
  if (i < numConstraintObjects)
    compositeDX = 
      constraintMVDXPtrs[i]->getDX()->clone(totalNumConstraints);
  else
    compositeDX = Teuchos::null;
}

LOCA::MultiContinuation::CompositeConstraintMVDX::CompositeConstraintMVDX(
	      const LOCA::MultiContinuation::CompositeConstraintMVDX& source, 
	      NOX::CopyType type) : 
  LOCA::MultiContinuation::CompositeConstraint(source),
  constraintMVDXPtrs(source.constraintMVDXPtrs),
  compositeDX()
{
  if (source.compositeDX.get() != NULL)
    compositeDX = source.compositeDX->clone(type);
  else
    compositeDX = Teuchos::null;
}

LOCA::MultiContinuation::CompositeConstraintMVDX::~CompositeConstraintMVDX()
{
}

void
LOCA::MultiContinuation::CompositeConstraintMVDX::copy(
		   const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::MultiContinuation::CompositeConstraintMVDX& source = 
    dynamic_cast<const LOCA::MultiContinuation::CompositeConstraintMVDX&>(src);

  if (this != &source) {
    LOCA::MultiContinuation::CompositeConstraint::copy(source);
    constraintMVDXPtrs = source.constraintMVDXPtrs;
    if (compositeDX.get() != NULL && source.compositeDX.get() != NULL)
      *compositeDX = *source.compositeDX;
    else if (source.compositeDX.get() != NULL)
      compositeDX = source.compositeDX->clone(NOX::DeepCopy);
    else
      compositeDX = Teuchos::null;
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::MultiContinuation::CompositeConstraintMVDX::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new CompositeConstraintMVDX(*this, type));
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraintMVDX::computeDX()
{
    string callingFunction = 
    "LOCA::MultiContinuation::CompositeConstraintMVDX::computeConstraints()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;

  if (isValidDX)
    return finalStatus;

  if (isDXZero())
    return finalStatus;

  Teuchos::RCP<NOX::Abstract::MultiVector> dx;
  for (int i=0; i<numConstraintObjects; i++) {

    if (!constraintMVDXPtrs[i]->isDXZero()) {

      // Compute constraint derivative
      status = constraintMVDXPtrs[i]->computeDX();
      finalStatus = 
	globalData->locaErrorCheck->combineAndCheckReturnTypes(
							     status, 
							     finalStatus,
							     callingFunction);
      
      // Copy columns of constraint dervative into composite
      dx = compositeDX->subView(indices[i]);
      *dx = *(constraintMVDXPtrs[i]->getDX());
    }
  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraintMVDX::multiplyDX(
		      double alpha, 
		      const NOX::Abstract::MultiVector& input_x,
	              NOX::Abstract::MultiVector::DenseMatrix& result_p) const
{
  return 
    LOCA::MultiContinuation::ConstraintInterfaceMVDX::multiplyDX(alpha,
								 input_x,
								 result_p);
}

NOX::Abstract::Group::ReturnType
LOCA::MultiContinuation::CompositeConstraintMVDX::addDX(
		              Teuchos::ETransp transb,
			      double alpha, 
	                      const NOX::Abstract::MultiVector::DenseMatrix& b,
			      double beta,
			      NOX::Abstract::MultiVector& result_x) const
{
  return 
    LOCA::MultiContinuation::ConstraintInterfaceMVDX::addDX(transb, alpha, b,
							    beta, result_x);
}

const NOX::Abstract::MultiVector*
LOCA::MultiContinuation::CompositeConstraintMVDX::getDX() const
{
  return compositeDX.get();
}
