// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
