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

#include "LOCA_MultiContinuation_ExtendedMultiVector.H" // Class definition
#include "LOCA_MultiContinuation_ExtendedVector.H"  

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
					   const NOX::Abstract::Vector& xVec,
					   int nColumns,
					   int nScalarRows, 
					   NOX::CopyType type) :
  LOCA::Extended::MultiVector(nColumns, 1, nScalarRows)
{
  NOX::Abstract::MultiVector* mv = xVec.createMultiVector(nColumns, 
							  type);
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, mv);
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
				      const NOX::Abstract::MultiVector& xVec, 
				      int nScalarRows) :
  LOCA::Extended::MultiVector(xVec.numVectors(), 1, nScalarRows)
{
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, xVec.clone(NOX::DeepCopy));
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
		       const NOX::Abstract::MultiVector& xVec,
		       const NOX::Abstract::MultiVector::DenseMatrix& params) :
  LOCA::Extended::MultiVector(xVec.numVectors(), 1, params.numRows())
{
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, xVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::getScalars() = params;
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
		  const LOCA::MultiContinuation::ExtendedMultiVector& source, 
		  NOX::CopyType type) :
  LOCA::Extended::MultiVector(source, type)
{
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
		  const LOCA::MultiContinuation::ExtendedMultiVector& source, 
		  int nColumns) :
  LOCA::Extended::MultiVector(source, nColumns)
{
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
		  const LOCA::MultiContinuation::ExtendedMultiVector& source, 
		  vector<int>& index, bool view) :
  LOCA::Extended::MultiVector(source, index, view)
{
}

LOCA::MultiContinuation::ExtendedMultiVector::~ExtendedMultiVector()
{
}

LOCA::Extended::MultiVector& 
LOCA::MultiContinuation::ExtendedMultiVector::operator=(
					 const LOCA::Extended::MultiVector& y)
{
  return 
    operator=(dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(y));
}

NOX::Abstract::MultiVector& 
LOCA::MultiContinuation::ExtendedMultiVector::operator=(
					 const NOX::Abstract::MultiVector& y)
{
 return 
   operator=(dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(y));
}

LOCA::MultiContinuation::ExtendedMultiVector& 
LOCA::MultiContinuation::ExtendedMultiVector::operator=(const 
				  LOCA::MultiContinuation::ExtendedMultiVector& y)
{
  LOCA::Extended::MultiVector::operator=(y);
  return *this;
}

NOX::Abstract::MultiVector* 
LOCA::MultiContinuation::ExtendedMultiVector::clone(NOX::CopyType type) const
{
  return new LOCA::MultiContinuation::ExtendedMultiVector(*this, type);
}

NOX::Abstract::MultiVector* 
LOCA::MultiContinuation::ExtendedMultiVector::clone(int numvecs) const
{
  return new LOCA::MultiContinuation::ExtendedMultiVector(*this, numvecs);
}

NOX::Abstract::MultiVector* 
LOCA::MultiContinuation::ExtendedMultiVector::subCopy(vector<int>& index) const
{
  return new LOCA::MultiContinuation::ExtendedMultiVector(*this, index, false);
}

NOX::Abstract::MultiVector* 
LOCA::MultiContinuation::ExtendedMultiVector::subView(vector<int>& index) const
{
  return new LOCA::MultiContinuation::ExtendedMultiVector(*this, index, true);
}

const NOX::Abstract::MultiVector& 
LOCA::MultiContinuation::ExtendedMultiVector::getXMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

NOX::Abstract::MultiVector& 
LOCA::MultiContinuation::ExtendedMultiVector::getXMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
							    int nColumns,
							    int nScalarRows) :
  LOCA::Extended::MultiVector(nColumns, 1, nScalarRows)
{
}

LOCA::Extended::Vector*
LOCA::MultiContinuation::ExtendedMultiVector::generateVector(
							int nVecs, 
							int nScalarRows) const
{
  return new LOCA::MultiContinuation::ExtendedVector(nScalarRows);
}
