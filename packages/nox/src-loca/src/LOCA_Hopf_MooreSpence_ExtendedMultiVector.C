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

#include "LOCA_Hopf_MooreSpence_ExtendedMultiVector.H" 
#include "LOCA_Hopf_MooreSpence_ExtendedVector.H"  

LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& cloneVec,
		    int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 3, 2)
{
  Teuchos::RCP<NOX::Abstract::MultiVector> mv1 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv2 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv3 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, mv1);
  LOCA::Extended::MultiVector::setMultiVectorPtr(1, mv2);
  LOCA::Extended::MultiVector::setMultiVectorPtr(2, mv3);
}

LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data,
		  const NOX::Abstract::MultiVector& xVec,
		  const NOX::Abstract::MultiVector& realEigenVec,
		  const NOX::Abstract::MultiVector& imagEigenVec,
		  const NOX::Abstract::MultiVector::DenseMatrix& freqs,
		  const NOX::Abstract::MultiVector::DenseMatrix& bifParams) :
  LOCA::Extended::MultiVector(global_data, xVec.numVectors(), 3, 2)
{
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, xVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::setMultiVectorPtr(1, realEigenVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::setMultiVectorPtr(2, imagEigenVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::getScalarRows(1,0)->assign(freqs);
  LOCA::Extended::MultiVector::getScalarRows(1,1)->assign(bifParams);
}

LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
	  const LOCA::Hopf::MooreSpence::ExtendedMultiVector& source, 
	  NOX::CopyType type) :
  LOCA::Extended::MultiVector(source, type)
{
}

LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
	   const LOCA::Hopf::MooreSpence::ExtendedMultiVector& source, 
	   int nColumns) :
  LOCA::Extended::MultiVector(source, nColumns)
{
}

LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
	   const LOCA::Hopf::MooreSpence::ExtendedMultiVector& source, 
	   const vector<int>& index, bool view) :
  LOCA::Extended::MultiVector(source, index, view)
{
}

LOCA::Hopf::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector()
{
}

LOCA::Extended::MultiVector& 
LOCA::Hopf::MooreSpence::ExtendedMultiVector::operator=(
					 const LOCA::Extended::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedMultiVector&>(y));
  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::Hopf::MooreSpence::ExtendedMultiVector::operator=(
					 const NOX::Abstract::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedMultiVector&>(y));
  return *this;
}

LOCA::Hopf::MooreSpence::ExtendedMultiVector& 
LOCA::Hopf::MooreSpence::ExtendedMultiVector::operator=(const 
		      LOCA::Hopf::MooreSpence::ExtendedMultiVector& y)
{
  LOCA::Extended::MultiVector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::clone(NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedMultiVector(*this, type));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::clone(int numvecs) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedMultiVector(*this, numvecs));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::subCopy(
					       const vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedMultiVector(*this, index, false));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::subView(
					      const vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedMultiVector(*this, index, true));
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getXMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getXMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getRealEigenMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getRealEigenMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getImagEigenMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(2);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getImagEigenMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(2);
}

Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getFrequencies() const
{
  return LOCA::Extended::MultiVector::getScalarRows(1, 0);
}

Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getFrequencies()
{
  return LOCA::Extended::MultiVector::getScalarRows(1, 0);
}

Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getBifParams() const
{
  return LOCA::Extended::MultiVector::getScalarRows(1, 1);
}

Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getBifParams()
{
  return LOCA::Extended::MultiVector::getScalarRows(1, 1);
}

LOCA::Hopf::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 3, 2)
{
}

Teuchos::RCP<LOCA::Extended::Vector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::generateVector(
							int nVecs, 
							int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::Hopf::MooreSpence::ExtendedVector(globalData));
}

Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getColumn(int i)
{
  return Teuchos::rcp_dynamic_cast<LOCA::Hopf::MooreSpence::ExtendedVector>(getVector(i),true);
}

Teuchos::RCP<const LOCA::Hopf::MooreSpence::ExtendedVector>
LOCA::Hopf::MooreSpence::ExtendedMultiVector::getColumn(int i) const
{
  return Teuchos::rcp_dynamic_cast<const LOCA::Hopf::MooreSpence::ExtendedVector>(getVector(i),true);
}
