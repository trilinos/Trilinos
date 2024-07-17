// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Abstract_Group.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Abstract::Group::Group(
           const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  globalData(global_data)
{
  setDerivUtils(Teuchos::rcp(new LOCA::DerivUtils(globalData)));
}

LOCA::Abstract::Group::Group(
           const Teuchos::RCP<LOCA::GlobalData>& global_data,
           const Teuchos::RCP<LOCA::DerivUtils>& deriv ) :
  globalData(global_data)
{
  setDerivUtils(deriv);
}

LOCA::Abstract::Group::Group(const LOCA::Abstract::Group& source,
                 NOX::CopyType /* type */) :
  globalData(source.globalData)
{
  // We use copy here instead of copy constructors to avoid issues with
  // multiple virtual inheritence
  copy(source);
}


LOCA::Abstract::Group::~Group()
{
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::augmentJacobianForHomotopy(double /* a */, double /* b */)
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::computeShiftedMatrix(double /* alpha */, double /* beta */)
{
  globalData->locaErrorCheck->throwError(
               "LOCA::Abstract::Group::computeShiftedMatrix",
               "Not implemented for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyShiftedMatrix(const NOX::Abstract::Vector& /* input */,
                      NOX::Abstract::Vector& /* result */) const
{
  globalData->locaErrorCheck->throwError(
               "LOCA::Abstract::Group::applyShiftedMatrix",
               "Not implemented for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyShiftedMatrixMultiVector(
                const NOX::Abstract::MultiVector& /* input */,
                NOX::Abstract::MultiVector& /* result */) const
{
  globalData->locaErrorCheck->throwError(
            "LOCA::Abstract::Group::applyShiftedMatrixMultiVector",
            "Not implemented for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyShiftedMatrixInverseMultiVector(
                    Teuchos::ParameterList& /* params */,
                const NOX::Abstract::MultiVector& /* input */,
                NOX::Abstract::MultiVector& /* result */) const
{
  globalData->locaErrorCheck->throwError(
        "LOCA::Abstract::Group::applyShiftedMatrixInverseMultiVector",
        "Not implemented for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::computeSecondShiftedMatrix(double /* alpha */, double /* beta */)
{
  globalData->locaErrorCheck->throwError(
               "LOCA::Abstract::Group::computeSecondShiftedMatrix",
               "Not implemented for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applySecondShiftedMatrix(const NOX::Abstract::Vector& /* input */,
                      NOX::Abstract::Vector& /* result */) const
{
  globalData->locaErrorCheck->throwError(
               "LOCA::Abstract::Group::applySecondShiftedMatrix",
               "Not implemented for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applySecondShiftedMatrixMultiVector(
                const NOX::Abstract::MultiVector& /* input */,
                NOX::Abstract::MultiVector& /* result */) const
{
  globalData->locaErrorCheck->throwError(
            "LOCA::Abstract::Group::applySecondShiftedMatrixMultiVector",
            "Not implemented for group");
  return NOX::Abstract::Group::NotDefined;
}

bool
LOCA::Abstract::Group::isComplex() const
{
  return false;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::computeComplex(double /* frequency */)
{
  globalData->locaErrorCheck->throwError(
                   "LOCA::Abstract::Group::computeComplex",
                   "Method not defined for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplex(const NOX::Abstract::Vector& /* input_real */,
                    const NOX::Abstract::Vector& /* input_imag */,
                    NOX::Abstract::Vector& /* result_real */,
                    NOX::Abstract::Vector& /* result_imag */) const
{
  globalData->locaErrorCheck->throwError("LOCA::Abstract::Group::applyComplex",
                     "Method not defined for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexMultiVector(
                  const NOX::Abstract::MultiVector& /* input_real */,
                  const NOX::Abstract::MultiVector& /* input_imag */,
                  NOX::Abstract::MultiVector& /* result_real */,
                  NOX::Abstract::MultiVector& /* result_imag */) const
{
  globalData->locaErrorCheck->throwError(
                 "LOCA::Abstract::Group::applyComplexMultiVector",
                 "Method not defined for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexInverseMultiVector(
                   Teuchos::ParameterList& /* params */,
                   const NOX::Abstract::MultiVector& /* input_real */,
                   const NOX::Abstract::MultiVector& /* input_imag */,
                   NOX::Abstract::MultiVector& /* result_real */,
                   NOX::Abstract::MultiVector& /* result_imag */) const
{
  globalData->locaErrorCheck->throwError(
                   "LOCA::Abstract::Group::applyComplexInverse",
                   "Method not defined for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexTranspose(
                      const NOX::Abstract::Vector& /* input_real */,
                      const NOX::Abstract::Vector& /* input_imag */,
                      NOX::Abstract::Vector& /* result_real */,
                      NOX::Abstract::Vector& /* result_imag */) const
{
  globalData->locaErrorCheck->throwError(
                 "LOCA::Abstract::Group::applyComplexTranspose",
                 "Method not defined for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexTransposeMultiVector(
                const NOX::Abstract::MultiVector& /* input_real */,
                const NOX::Abstract::MultiVector& /* input_imag */,
                NOX::Abstract::MultiVector& /* result_real */,
                NOX::Abstract::MultiVector& /* result_imag */) const
{
  globalData->locaErrorCheck->throwError(
            "LOCA::Abstract::Group::applyComplexTransposeMultiVector",
            "Method not defined for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexTransposeInverseMultiVector(
                Teuchos::ParameterList& /* params */,
                const NOX::Abstract::MultiVector& /* input_real */,
                const NOX::Abstract::MultiVector& /* input_imag */,
                NOX::Abstract::MultiVector& /* result_real */,
                NOX::Abstract::MultiVector& /* result_imag */) const
{
  globalData->locaErrorCheck->throwError(
          "LOCA::Abstract::Group::applyComplexTransposeInverseMultiVector",
          "Method not defined for group");
  return NOX::Abstract::Group::NotDefined;
}

void
LOCA::Abstract::Group::copy(const NOX::Abstract::Group& src)
{
  const LOCA::Abstract::Group& source =
    dynamic_cast<const LOCA::Abstract::Group&>(src);

  // Copy parent classes
  LOCA::MultiContinuation::FiniteDifferenceGroup::copy(src);
  globalData = source.globalData;
}

void
LOCA::Abstract::Group::setParamsMulti(const std::vector<int>& paramIDs,
        const NOX::Abstract::MultiVector::DenseMatrix& vals)
{
  for (unsigned int i=0; i<paramIDs.size(); i++)
    setParam(paramIDs[i], vals(i,0));
}

void
LOCA::Abstract::Group::notifyCompletedStep()
{
}

NOX::Abstract::Group&
LOCA::Abstract::Group::operator=(const NOX::Abstract::Group& source)
{

  copy(source);
  return *this;
}

double
LOCA::Abstract::Group::computeFreeEnergy()
{
  globalData->locaErrorCheck->throwError(
            "LOCA::Abstract::Group::computeFreeEnergy",
            "Method not defined for group");
  return 0.0;
}

