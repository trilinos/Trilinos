// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Pitchfork_MooreSpence_ExtendedMultiVector.H"
#include "LOCA_Pitchfork_MooreSpence_ExtendedVector.H"

LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Abstract::Vector& cloneVec,
            int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 2)
{
  Teuchos::RCP<NOX::Abstract::MultiVector> mv1 =
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector> mv2 =
    cloneVec.createMultiVector(nColumns, NOX::ShapeCopy);
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, mv1);
  LOCA::Extended::MultiVector::setMultiVectorPtr(1, mv2);
}

LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
          const Teuchos::RCP<LOCA::GlobalData>& global_data,
          const NOX::Abstract::MultiVector& xVec,
          const NOX::Abstract::MultiVector& nullVec,
          const NOX::Abstract::MultiVector::DenseMatrix& slacks,
          const NOX::Abstract::MultiVector::DenseMatrix& bifParams) :
  LOCA::Extended::MultiVector(global_data, xVec.numVectors(), 2, 2)
{
  LOCA::Extended::MultiVector::setMultiVectorPtr(0, xVec.clone(NOX::DeepCopy));
  LOCA::Extended::MultiVector::setMultiVectorPtr(1,
                         nullVec.clone(NOX::DeepCopy));

  LOCA::Extended::MultiVector::getScalarRows(1,0)->assign(slacks);
  LOCA::Extended::MultiVector::getScalarRows(1,1)->assign(bifParams);
}

LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
      const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& source,
      NOX::CopyType type) :
  LOCA::Extended::MultiVector(source, type)
{
}

LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
       const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& source,
       int nColumns) :
  LOCA::Extended::MultiVector(source, nColumns)
{
}

LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
       const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& source,
       const std::vector<int>& index, bool view) :
  LOCA::Extended::MultiVector(source, index, view)
{
}

LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::~ExtendedMultiVector()
{
}

LOCA::Extended::MultiVector&
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::operator=(
                     const LOCA::Extended::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector&>(y));
  return *this;
}

NOX::Abstract::MultiVector&
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::operator=(
                     const NOX::Abstract::MultiVector& y)
{
  operator=(dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector&>(y));
  return *this;
}

LOCA::Pitchfork::MooreSpence::ExtendedMultiVector&
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::operator=(const
              LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& y)
{
  LOCA::Extended::MultiVector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::clone(NOX::CopyType type) const
{
  return
    Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedMultiVector(*this, type));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::clone(int numvecs) const
{
  return
    Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedMultiVector(*this, numvecs));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::subCopy(
                           const std::vector<int>& index) const
{
  return
    Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedMultiVector(*this, index, false));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::subView(
                          const std::vector<int>& index) const
{
  return
    Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedMultiVector(*this, index, true));
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getXMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getXMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(0);
}

Teuchos::RCP<const NOX::Abstract::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getNullMultiVec() const
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

Teuchos::RCP<NOX::Abstract::MultiVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getNullMultiVec()
{
  return LOCA::Extended::MultiVector::getMultiVector(1);
}

Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getSlacks() const
{
  return LOCA::Extended::MultiVector::getScalarRows(1,0);
}

Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getSlacks()
{
  return LOCA::Extended::MultiVector::getScalarRows(1,0);
}

Teuchos::RCP<const NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getBifParams() const
{
  return LOCA::Extended::MultiVector::getScalarRows(1,1);
}

Teuchos::RCP<NOX::Abstract::MultiVector::DenseMatrix>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getBifParams()
{
  return LOCA::Extended::MultiVector::getScalarRows(1,1);
}

LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::ExtendedMultiVector(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            int nColumns) :
  LOCA::Extended::MultiVector(global_data, nColumns, 2, 2)
{
}

Teuchos::RCP<LOCA::Extended::Vector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::generateVector(
                            int /* nVecs */,
                            int /* nScalarRows */) const
{
  return
    Teuchos::rcp(new LOCA::Pitchfork::MooreSpence::ExtendedVector(
                                 globalData));
}

Teuchos::RCP<LOCA::Pitchfork::MooreSpence::ExtendedVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getColumn(int i)
{
  return Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MooreSpence::ExtendedVector>(getVector(i),true);
}

Teuchos::RCP<const LOCA::Pitchfork::MooreSpence::ExtendedVector>
LOCA::Pitchfork::MooreSpence::ExtendedMultiVector::getColumn(int i) const
{
  return Teuchos::rcp_dynamic_cast<const LOCA::Pitchfork::MooreSpence::ExtendedVector>(getVector(i),true);
}
