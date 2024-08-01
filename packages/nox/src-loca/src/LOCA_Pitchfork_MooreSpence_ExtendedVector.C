// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Pitchfork_MooreSpence_ExtendedVector.H"  // Class definition
#include "LOCA_Pitchfork_MooreSpence_ExtendedMultiVector.H"

LOCA::Pitchfork::MooreSpence::ExtendedVector::ExtendedVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& xVec,
            const NOX::Abstract::Vector& nullVec,
            double slack,
            double bifParam) :
  LOCA::Extended::Vector(global_data,2,2)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, slack);
  setScalar(1, bifParam);
}

LOCA::Pitchfork::MooreSpence::ExtendedVector::ExtendedVector(
                const LOCA::Pitchfork::MooreSpence::ExtendedVector& source,
        NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::Pitchfork::MooreSpence::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector&
LOCA::Pitchfork::MooreSpence::ExtendedVector::operator=(
                          const NOX::Abstract::Vector& y)
{
  operator=(dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedVector&>(y));
  return *this;
}

LOCA::Extended::Vector&
LOCA::Pitchfork::MooreSpence::ExtendedVector::operator=(
                         const LOCA::Extended::Vector& y)
{
  operator=(dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedVector&>(y));
  return *this;
}

LOCA::Pitchfork::MooreSpence::ExtendedVector&
LOCA::Pitchfork::MooreSpence::ExtendedVector::operator=(
                     const LOCA::Pitchfork::MooreSpence::ExtendedVector& y)
{
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Pitchfork::MooreSpence::ExtendedVector::clone(
                            NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedVector(*this,
                                     type));
}

void
LOCA::Pitchfork::MooreSpence::ExtendedVector::setVec(
                    const NOX::Abstract::Vector& xVec,
                    const NOX::Abstract::Vector& nullVec,
                    double slack,
                    double bifPar)
{
  setVector(0, xVec);
  setVector(1, nullVec);
  setScalar(0, slack);
  setScalar(1, bifPar);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Pitchfork::MooreSpence::ExtendedVector::getXVec() const
{
  return getVector(0);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::Pitchfork::MooreSpence::ExtendedVector::getNullVec() const
{
  return getVector(1);
}

double
LOCA::Pitchfork::MooreSpence::ExtendedVector::getSlack() const
{
  return getScalar(0);
}

double
LOCA::Pitchfork::MooreSpence::ExtendedVector::getBifParam() const
{
  return getScalar(1);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Pitchfork::MooreSpence::ExtendedVector::getXVec()
{
  return getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::Pitchfork::MooreSpence::ExtendedVector::getNullVec()
{
  return getVector(1);
}

double&
LOCA::Pitchfork::MooreSpence::ExtendedVector::getSlack()
{
  return getScalar(0);
}

double&
LOCA::Pitchfork::MooreSpence::ExtendedVector::getBifParam()
{
  return getScalar(1);
}

LOCA::Pitchfork::MooreSpence::ExtendedVector::ExtendedVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,2,2)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedVector::generateMultiVector(
                            int nColumns,
                            int /* nVectorRows */,
                            int /* nScalarRows */) const
{
  return
    Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedMultiVector(
                                    globalData,
                                    nColumns));
}
