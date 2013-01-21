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

#include "LOCA_MultiContinuation_ExtendedMultiVector.H" // Class definition
#include "LOCA_MultiContinuation_ExtendedVector.H"  

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& xVec,
		    int nColumns,
		    int nScalarRows, 
		    NOX::CopyType type) :
  LOCA::Extended::MultiVector(global_data, nColumns, 1, nScalarRows)
{
  Teuchos::RCP<NOX::Abstract::MultiVector> mv = 
    xVec.createMultiVector(nColumns, type);
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, mv);
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::MultiVector& xVec, 
		    int nScalarRows) :
  LOCA::Extended::MultiVector(global_data, xVec.numVectors(), 1, nScalarRows)
{
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, xVec.clone(NOX::DeepCopy));
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
		     const Teuchos::RCP<LOCA::GlobalData>& global_data,
		     const NOX::Abstract::MultiVector& xVec,
		     const NOX::Abstract::MultiVector::DenseMatrix& params) :
  LOCA::Extended::MultiVector(global_data, xVec.numVectors(), 1, 
			      params.numRows())
{
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, xVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::getScalars()->assign(params);
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
		  const std::vector<int>& index, bool view) :
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
  operator=(dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(y));
  return *this;
}

NOX::Abstract::MultiVector& 
LOCA::MultiContinuation::ExtendedMultiVector::operator=(
					 const NOX::Abstract::MultiVector& y)
{
 operator=(dynamic_cast<const LOCA::MultiContinuation::ExtendedMultiVector&>(y));
 return *this;
}

LOCA::MultiContinuation::ExtendedMultiVector& 
LOCA::MultiContinuation::ExtendedMultiVector::operator=(const 
				  LOCA::MultiContinuation::ExtendedMultiVector& y)
{
  LOCA::Extended::MultiVector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::MultiContinuation::ExtendedMultiVector::clone(NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::MultiContinuation::ExtendedMultiVector(*this, 
								  type));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::MultiContinuation::ExtendedMultiVector::clone(int numvecs) const
{
  return 
    Teuchos::rcp(new LOCA::MultiContinuation::ExtendedMultiVector(*this, 
								  numvecs));
}

Teuchos::RCP<NOX::Abstract::MultiVector> 
LOCA::MultiContinuation::ExtendedMultiVector::subCopy(
					       const std::vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::MultiContinuation::ExtendedMultiVector(*this, 
								  index, 
								  false));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::MultiContinuation::ExtendedMultiVector::subView(
					      const std::vector<int>& index) const
{
  return 
    Teuchos::rcp(new LOCA::MultiContinuation::ExtendedMultiVector(*this, 
								  index, 
								  true));
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::MultiContinuation::ExtendedMultiVector::getXMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::MultiContinuation::ExtendedMultiVector::getXMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

LOCA::MultiContinuation::ExtendedMultiVector::ExtendedMultiVector(
		     const Teuchos::RCP<LOCA::GlobalData>& global_data,
		     int nColumns,
		     int nScalarRows) :
  LOCA::Extended::MultiVector(global_data, nColumns, 1, nScalarRows)
{
}

Teuchos::RCP<LOCA::Extended::Vector>
LOCA::MultiContinuation::ExtendedMultiVector::generateVector(
							int nVecs, 
							int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::MultiContinuation::ExtendedVector(globalData, 
							     nScalarRows));
}
