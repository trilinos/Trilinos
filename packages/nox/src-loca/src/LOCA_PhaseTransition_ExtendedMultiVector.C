// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_PhaseTransition_ExtendedMultiVector.H"
#include "LOCA_PhaseTransition_ExtendedVector.H"

LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(
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

LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(
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

LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(
      const LOCA::PhaseTransition::ExtendedMultiVector& source,
      NOX::CopyType type) :
  LOCA::Extended::MultiVector(source, type)
{
}

LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(
       const LOCA::PhaseTransition::ExtendedMultiVector& source,
       int nColumns) :
  LOCA::Extended::MultiVector(source, nColumns)
{
}

LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(
       const LOCA::PhaseTransition::ExtendedMultiVector& source,
       const std::vector<int>& index, bool view) :
  LOCA::Extended::MultiVector(source, index, view)
{
}

LOCA::PhaseTransition::ExtendedMultiVector::~ExtendedMultiVector()
{
}

LOCA::Extended::MultiVector&
LOCA::PhaseTransition::ExtendedMultiVector::operator=(
                     const LOCA::Extended::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::PhaseTransition::ExtendedMultiVector&>(y));
  return *this;
}

NOX::Abstract::MultiVector&
LOCA::PhaseTransition::ExtendedMultiVector::operator=(
                     const NOX::Abstract::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::PhaseTransition::ExtendedMultiVector&>(y));
  return *this;
}

LOCA::PhaseTransition::ExtendedMultiVector&
LOCA::PhaseTransition::ExtendedMultiVector::operator=(const
              LOCA::PhaseTransition::ExtendedMultiVector& y)
{
  LOCA::Extended::MultiVector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::PhaseTransition::ExtendedMultiVector::clone(NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedMultiVector(*this, type));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::PhaseTransition::ExtendedMultiVector::clone(int numvecs) const
{
  return
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedMultiVector(*this, numvecs));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::PhaseTransition::ExtendedMultiVector::subCopy(
                           const std::vector<int>& index) const
{
  return
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedMultiVector(*this, index, false));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::PhaseTransition::ExtendedMultiVector::subView(
                          const std::vector<int>& index) const
{
  return
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedMultiVector(*this, index, true));
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::PhaseTransition::ExtendedMultiVector::getXMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::PhaseTransition::ExtendedMultiVector::getXMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::PhaseTransition::ExtendedMultiVector::getNullMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::PhaseTransition::ExtendedMultiVector::getNullMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

LOCA::PhaseTransition::ExtendedMultiVector::ExtendedMultiVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 1)
{
}

Teuchos::RCP<LOCA::Extended::Vector>
LOCA::PhaseTransition::ExtendedMultiVector::generateVector(
                            int /* nVecs */,
                            int /* nScalarRows */) const
{
  return
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(
                                 globalData));
}

Teuchos::RCP<LOCA::PhaseTransition::ExtendedVector>
LOCA::PhaseTransition::ExtendedMultiVector::getColumn(int i)
{
  return Teuchos::rcp_dynamic_cast<LOCA::PhaseTransition::ExtendedVector>(getVector(i),true);
}

Teuchos::RCP<const LOCA::PhaseTransition::ExtendedVector>
LOCA::PhaseTransition::ExtendedMultiVector::getColumn(int i) const
{
  return Teuchos::rcp_dynamic_cast<const LOCA::PhaseTransition::ExtendedVector>(getVector(i),true);
}
