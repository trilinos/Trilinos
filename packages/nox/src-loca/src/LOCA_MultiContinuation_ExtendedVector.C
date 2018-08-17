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

#include "LOCA_MultiContinuation_ExtendedVector.H"  // Class definition
#include "LOCA_MultiContinuation_ExtendedMultiVector.H"

LOCA::MultiContinuation::ExtendedVector::ExtendedVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& xVec,
            int nScalars) :
  LOCA::Extended::Vector(global_data,1,nScalars)
{
  LOCA::Extended::Vector::setVector(0, xVec);
}

LOCA::MultiContinuation::ExtendedVector::ExtendedVector(
            const LOCA::MultiContinuation::ExtendedVector& source,
            NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}

LOCA::MultiContinuation::ExtendedVector::~ExtendedVector()
{
}

LOCA::Extended::Vector&
LOCA::MultiContinuation::ExtendedVector::operator=(
                         const LOCA::Extended::Vector& y)
{
  operator=(dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y));
  return *this;
}

NOX::Abstract::Vector&
LOCA::MultiContinuation::ExtendedVector::operator=(
                           const NOX::Abstract::Vector& y)
{
 operator=(dynamic_cast<const LOCA::MultiContinuation::ExtendedVector&>(y));
 return *this;
}

LOCA::MultiContinuation::ExtendedVector&
LOCA::MultiContinuation::ExtendedVector::operator=(const
                   LOCA::MultiContinuation::ExtendedVector& y)
{
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::MultiContinuation::ExtendedVector::clone(NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::MultiContinuation::ExtendedVector(*this, type));
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::MultiContinuation::ExtendedVector::getXVec() const
{
  return LOCA::Extended::Vector::getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::MultiContinuation::ExtendedVector::getXVec()
{
  return LOCA::Extended::Vector::getVector(0);
}

LOCA::MultiContinuation::ExtendedVector::ExtendedVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            int nScalars) :
  LOCA::Extended::Vector(global_data,1,nScalars)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::MultiContinuation::ExtendedVector::generateMultiVector(
                            int nColumns,
                            int /* nVectorRows */,
                            int nScalarRows) const
{
  return
    Teuchos::rcp(new LOCA::MultiContinuation::ExtendedMultiVector(globalData,
                                  nColumns,
                                  nScalarRows));
}
