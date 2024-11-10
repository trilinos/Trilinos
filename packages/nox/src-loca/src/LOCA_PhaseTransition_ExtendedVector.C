// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_PhaseTransition_ExtendedVector.H"  // Class definition
#include "LOCA_PhaseTransition_ExtendedMultiVector.H"

LOCA::PhaseTransition::ExtendedVector::ExtendedVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& x1Vec,
            const NOX::Abstract::Vector& x2vec,
            double ptp) :
  LOCA::Extended::Vector(global_data,2,1)
{
  setVector(0, x1Vec);
  setVector(1, x2vec);
  setScalar(0, ptp);
}

LOCA::PhaseTransition::ExtendedVector::ExtendedVector(
                const LOCA::PhaseTransition::ExtendedVector& source,
        NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::PhaseTransition::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector&
LOCA::PhaseTransition::ExtendedVector::operator=(
                          const NOX::Abstract::Vector& y)
{
  operator=(dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&>(y));
  return *this;
}

LOCA::Extended::Vector&
LOCA::PhaseTransition::ExtendedVector::operator=(
                         const LOCA::Extended::Vector& y)
{
  operator=(dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&>(y));
  return *this;
}

LOCA::PhaseTransition::ExtendedVector&
LOCA::PhaseTransition::ExtendedVector::operator=(
                     const LOCA::PhaseTransition::ExtendedVector& y)
{
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::clone(
                            NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(*this,
                                     type));
}

void
LOCA::PhaseTransition::ExtendedVector::setVec(
                    const NOX::Abstract::Vector& x1Vec,
                    const NOX::Abstract::Vector& x2vec,
                    double ptp)
{
  setVector(0, x1Vec);
  setVector(1, x2vec);
  setScalar(0, ptp);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::X1() const
{
  return getVector(0);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::X2() const
{
  return getVector(1);
}

double
LOCA::PhaseTransition::ExtendedVector::PTP() const
{
  return getScalar(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::X1()
{
  return getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::X2()
{
  return getVector(1);
}

double&
LOCA::PhaseTransition::ExtendedVector::PTP()
{
  return getScalar(0);
}

LOCA::PhaseTransition::ExtendedVector::ExtendedVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,2,1)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::PhaseTransition::ExtendedVector::generateMultiVector(
                            int nColumns,
                            int /* nVectorRows */,
                            int /* nScalarRows */) const
{
  return
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedMultiVector(
                                    globalData,
                                    nColumns));
}
