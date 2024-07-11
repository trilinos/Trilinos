// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
                            int /* nVecs */,
                            int nScalarRows) const
{
  return
    Teuchos::rcp(new LOCA::MultiContinuation::ExtendedVector(globalData,
                                 nScalarRows));
}
