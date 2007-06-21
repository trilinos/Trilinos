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

#include "LOCA_Hopf_ComplexMultiVector.H" 
#include "LOCA_Hopf_ComplexVector.H"  

LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& cloneVec,
		    int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 0)
{
  Teuchos::RCP<NOX::Abstract::MultiVector> mv1 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv2 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, mv1);
  LOCA::Extended::MultiVector::setMultiVectorPtr(1, mv2);
}

LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data,
		  const NOX::Abstract::MultiVector& realVec,
		  const NOX::Abstract::MultiVector& imagVec) :
  LOCA::Extended::MultiVector(global_data, realVec.numVectors(), 2, 0)
{
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, 
						 realVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::setMultiVectorPtr(1, 
						 imagVec.clone(NOX::DeepCopy));
}

LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(
				 const LOCA::Hopf::ComplexMultiVector& source, 
				 NOX::CopyType type) :
  LOCA::Extended::MultiVector(source, type)
{
}

LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(
				 const LOCA::Hopf::ComplexMultiVector& source, 
				 int nColumns) :
  LOCA::Extended::MultiVector(source, nColumns)
{
}

LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(
				 const LOCA::Hopf::ComplexMultiVector& source, 
				 const vector<int>& index, bool view) :
  LOCA::Extended::MultiVector(source, index, view)
{
}

LOCA::Hopf::ComplexMultiVector::~ComplexMultiVector()
{
}

LOCA::Extended::MultiVector& 
LOCA::Hopf::ComplexMultiVector::operator=(const LOCA::Extended::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::Hopf::ComplexMultiVector&>(y));
  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::Hopf::ComplexMultiVector::operator=(const NOX::Abstract::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::Hopf::ComplexMultiVector&>(y));
  return *this;
}

LOCA::Hopf::ComplexMultiVector& 
LOCA::Hopf::ComplexMultiVector::operator=(
				      const LOCA::Hopf::ComplexMultiVector& y)
{
  LOCA::Extended::MultiVector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::ComplexMultiVector::clone(NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(*this, type));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::ComplexMultiVector::clone(int numvecs) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(*this, numvecs));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::ComplexMultiVector::subCopy(const vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(*this, index, false));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::ComplexMultiVector::subView(const vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(*this, index, true));
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::Hopf::ComplexMultiVector::getRealMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::ComplexMultiVector::getRealMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::Hopf::ComplexMultiVector::getImagMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::ComplexMultiVector::getImagMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

LOCA::Hopf::ComplexMultiVector::ComplexMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 0)
{
}

Teuchos::RCP<LOCA::Extended::Vector>
LOCA::Hopf::ComplexMultiVector::generateVector(int nVecs, 
					       int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::ComplexVector(globalData));
}

Teuchos::RCP<LOCA::Hopf::ComplexVector>
LOCA::Hopf::ComplexMultiVector::getColumn(int i)
{
  return Teuchos::rcp_dynamic_cast<LOCA::Hopf::ComplexVector>(getVector(i),true);
}

Teuchos::RCP<const LOCA::Hopf::ComplexVector>
LOCA::Hopf::ComplexMultiVector::getColumn(int i) const
{
  return Teuchos::rcp_dynamic_cast<const LOCA::Hopf::ComplexVector>(getVector(i),true);
}
