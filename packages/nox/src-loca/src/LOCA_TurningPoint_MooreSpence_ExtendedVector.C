// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_TurningPoint_MooreSpence_ExtendedVector.H"  // Class definition
#include "LOCA_TurningPoint_MooreSpence_ExtendedMultiVector.H"

LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& xVec,
            const NOX::Abstract::Vector& nullVec,
            double bifParam) :
  LOCA::Extended::Vector(global_data,2,1)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, bifParam);
}

LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector(
                const LOCA::TurningPoint::MooreSpence::ExtendedVector& source,
        NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::TurningPoint::MooreSpence::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedVector::operator=(
                          const NOX::Abstract::Vector& y)
{
  operator=(dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(y));
  return *this;
}

LOCA::Extended::Vector&
LOCA::TurningPoint::MooreSpence::ExtendedVector::operator=(
                         const LOCA::Extended::Vector& y)
{
  operator=(dynamic_cast<const LOCA::TurningPoint::MooreSpence::ExtendedVector&>(y));
  return *this;
}

LOCA::TurningPoint::MooreSpence::ExtendedVector&
LOCA::TurningPoint::MooreSpence::ExtendedVector::operator=(
                     const LOCA::TurningPoint::MooreSpence::ExtendedVector& y)
{
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::clone(
                            NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedVector(*this,
                                     type));
}

void
LOCA::TurningPoint::MooreSpence::ExtendedVector::setVec(
                    const NOX::Abstract::Vector& xVec,
                    const NOX::Abstract::Vector& nullVec,
                    double bifPar)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, bifPar);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::getXVec() const
{
  return getVector(0);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::getNullVec() const
{
  return getVector(1);
}

double
LOCA::TurningPoint::MooreSpence::ExtendedVector::getBifParam() const
{
  return getScalar(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::getXVec()
{
  return getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::getNullVec()
{
  return getVector(1);
}

double&
LOCA::TurningPoint::MooreSpence::ExtendedVector::getBifParam()
{
  return getScalar(0);
}

LOCA::TurningPoint::MooreSpence::ExtendedVector::ExtendedVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,2,1)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::TurningPoint::MooreSpence::ExtendedVector::generateMultiVector(
                            int nColumns,
                            int /* nVectorRows */,
                            int /* nScalarRows */) const
{
  return
    Teuchos::rcp(new LOCA::TurningPoint::MooreSpence::ExtendedMultiVector(
                                    globalData,
                                    nColumns));
}
