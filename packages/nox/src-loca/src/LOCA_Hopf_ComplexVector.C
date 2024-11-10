// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Hopf_ComplexVector.H"  // Class definition
#include "LOCA_Hopf_ComplexMultiVector.H"

LOCA::Hopf::ComplexVector::ComplexVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& realVec,
            const NOX::Abstract::Vector& imagVec) :
  LOCA::Extended::Vector(global_data,2,0)
{
  setVector(0, realVec);
  setVector(1, imagVec);
}

LOCA::Hopf::ComplexVector::ComplexVector(
                     const LOCA::Hopf::ComplexVector& source,
                     NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::Hopf::ComplexVector::~ComplexVector()
{
}

NOX::Abstract::Vector&
LOCA::Hopf::ComplexVector::operator=(const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Hopf::ComplexVector&>(y));
}

LOCA::Extended::Vector&
LOCA::Hopf::ComplexVector::operator=(const LOCA::Extended::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::Hopf::ComplexVector&>(y));
}

LOCA::Hopf::ComplexVector&
LOCA::Hopf::ComplexVector::operator=(const LOCA::Hopf::ComplexVector& y)
{
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::clone(NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::Hopf::ComplexVector(*this, type));
}

void
LOCA::Hopf::ComplexVector::setVec(const NOX::Abstract::Vector& realVec,
                  const NOX::Abstract::Vector& imagVec)
{
  setVector(0, realVec);
  setVector(1, imagVec);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::getRealVec() const
{
  return getVector(0);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::getImagVec() const
{
  return getVector(1);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::getRealVec()
{
  return getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Hopf::ComplexVector::getImagVec()
{
  return getVector(1);
}

LOCA::Hopf::ComplexVector::ComplexVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,2,0)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::Hopf::ComplexVector::generateMultiVector(int nColumns,
                           int /* nVectorRows */,
                           int /* nScalarRows */) const
{
  return
    Teuchos::rcp(new LOCA::Hopf::ComplexMultiVector(globalData, nColumns));
}
