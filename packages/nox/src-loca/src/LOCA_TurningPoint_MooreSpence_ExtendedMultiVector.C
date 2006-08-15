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

#include "LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.H" 
#include "LOCA_TurningPoint_MooreSpence_ExtendedVector.H"  

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& cloneVec,
		    int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 1)
{
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> mv1 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> mv2 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, mv1);
  LOCA::Extended::MultiVector::setMultiVectorPtr(1, mv2);
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
		  const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		  const NOX::Abstract::MultiVector& xVec,
		  const NOX::Abstract::MultiVector& nullVec,
		  const NOX::Abstract::MultiVector::DenseMatrix& bifParams) :
  LOCA::Extended::MultiVector(global_data, xVec.numVectors(), 2, 1)
{
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, xVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::setMultiVectorPtr(1, 
						 nullVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::getScalars()->assign(bifParams);
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
	  const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& source, 
	  NOX::CopyType type) :
  LOCA::Extended::MultiVector(source, type)
{
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
	   const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& source, 
	   int nColumns) :
  LOCA::Extended::MultiVector(source, nColumns)
{
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
	   const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& source, 
	   const vector<int>& index, bool view) :
  LOCA::Extended::MultiVector(source, index, view)
{
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector()
{
}

LOCA::Extended::MultiVector& 
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::operator=(
					 const LOCA::Extended::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(y));
  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::operator=(
					 const NOX::Abstract::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedMultiVector&>(y));
  return *this;
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& 
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::operator=(const 
		      LOCA::TurningPoint::MooreSpence::ExtendedMultiVector& y)
{
  LOCA::Extended::MultiVector::operator=(y);
  return *this;
}

Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::clone(NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(*this, type));
}

Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::clone(int numvecs) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(*this, numvecs));
}

Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::subCopy(
					       const vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(*this, index, false));
}

Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::subView(
					      const vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(*this, index, true));
}

Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getXMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getXMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RefCountPtr<const NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getNullMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

Teuchos::RefCountPtr<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getNullMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
		    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
		    int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 1)
{
}

Teuchos::RefCountPtr<LOCA::Extended::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::generateVector(
							int nVecs, 
							int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedVector(
								 globalData));
}

Teuchos::RefCountPtr<LOCA::TurningPoint::MooreSpence::ExtendedVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getColumn(int i)
{
  return Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedVector>(getVector(i),true);
}

Teuchos::RefCountPtr<const LOCA::TurningPoint::MooreSpence::ExtendedVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getColumn(int i) const
{
  return Teuchos::rcp_dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector>(getVector(i),true);
}
