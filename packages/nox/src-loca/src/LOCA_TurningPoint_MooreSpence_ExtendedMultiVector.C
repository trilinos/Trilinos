// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& cloneVec,
		    int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 1)
{
  Teuchos::RCP<NOX::Abstract::MultiVector> mv1 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv2 = 
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, mv1);
  LOCA::Extended::MultiVector::setMultiVectorPtr(1, mv2);
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data,
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
	   const std::vector<int>& index, bool view) :
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

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::clone(NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(*this, type));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::clone(int numvecs) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(*this, numvecs));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::subCopy(
					       const std::vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(*this, index, false));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::subView(
					      const std::vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(*this, index, true));
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getXMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getXMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getNullMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getNullMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 1)
{
}

Teuchos::RCP<LOCA::Extended::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::generateVector(
							int nVecs, 
							int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedVector(
								 globalData));
}

Teuchos::RCP<LOCA::TurningPoint::MooreSpence::ExtendedVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getColumn(int i)
{
  return Teuchos::rcp_dynamic_cast<LOCA::TurningPoint::MooreSpence::ExtendedVector>(getVector(i),true);
}

Teuchos::RCP<const LOCA::TurningPoint::MooreSpence::ExtendedVector>
LOCA::TurningPoint::MooreSpence::ExtendedMultiVector::getColumn(int i) const
{
  return Teuchos::rcp_dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector>(getVector(i),true);
}
