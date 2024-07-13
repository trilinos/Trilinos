// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_BorderedSolver_ComplexOperator.H"
#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
#include "LOCA_Hopf_ComplexMultiVector.H"
#include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"

LOCA::BorderedSolver::ComplexOperator::
ComplexOperator(const Teuchos::RCP<const LOCA::Hopf::MooreSpence::AbstractGroup>& grp,
        double Omega) :
  grpPtr(grp),
  omega(Omega)
{
}

LOCA::BorderedSolver::ComplexOperator::
~ComplexOperator()
{
}

Teuchos::RCP<const NOX::Abstract::Group>
LOCA::BorderedSolver::ComplexOperator::
getGroup() const
{
  return grpPtr;
}

double
LOCA::BorderedSolver::ComplexOperator::
getFrequency() const
{
  return omega;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::ComplexOperator::
apply(const NOX::Abstract::MultiVector& X,
      NOX::Abstract::MultiVector& Y) const
{
  const LOCA::Hopf::ComplexMultiVector& cX =
    dynamic_cast<const LOCA::Hopf::ComplexMultiVector&>(X);
  LOCA::Hopf::ComplexMultiVector& cY =
    dynamic_cast<LOCA::Hopf::ComplexMultiVector&>(Y);
  return grpPtr->applyComplexMultiVector(*(cX.getRealMultiVec()),
                     *(cX.getImagMultiVec()),
                     *(cY.getRealMultiVec()),
                     *(cY.getImagMultiVec()));
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::ComplexOperator::
applyTranspose(const NOX::Abstract::MultiVector& X,
           NOX::Abstract::MultiVector& Y) const
{
  Teuchos::RCP<const LOCA::Hopf::MinimallyAugmented::AbstractGroup> magrp = Teuchos::rcp_dynamic_cast<const LOCA::Hopf::MinimallyAugmented::AbstractGroup>(grpPtr);
  const LOCA::Hopf::ComplexMultiVector& cX =
    dynamic_cast<const LOCA::Hopf::ComplexMultiVector&>(X);
  LOCA::Hopf::ComplexMultiVector& cY =
    dynamic_cast<LOCA::Hopf::ComplexMultiVector&>(Y);

  if (magrp != Teuchos::null)
    return magrp->applyComplexTransposeMultiVector(*(cX.getRealMultiVec()),
                           *(cX.getImagMultiVec()),
                           *(cY.getRealMultiVec()),
                           *(cY.getImagMultiVec()));
  else
    return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::ComplexOperator::
applyInverse(Teuchos::ParameterList& params,
         const NOX::Abstract::MultiVector& B,
         NOX::Abstract::MultiVector& X) const
{
  const LOCA::Hopf::ComplexMultiVector& cB =
    dynamic_cast<const LOCA::Hopf::ComplexMultiVector&>(B);
  LOCA::Hopf::ComplexMultiVector& cX =
    dynamic_cast<LOCA::Hopf::ComplexMultiVector&>(X);
  return grpPtr->applyComplexInverseMultiVector(params,
                        *(cB.getRealMultiVec()),
                        *(cB.getImagMultiVec()),
                        *(cX.getRealMultiVec()),
                        *(cX.getImagMultiVec()));
}

NOX::Abstract::Group::ReturnType
LOCA::BorderedSolver::ComplexOperator::
applyInverseTranspose(Teuchos::ParameterList& params,
              const NOX::Abstract::MultiVector& B,
              NOX::Abstract::MultiVector& X) const
{
  Teuchos::RCP<const LOCA::Hopf::MinimallyAugmented::AbstractGroup> magrp = Teuchos::rcp_dynamic_cast<const LOCA::Hopf::MinimallyAugmented::AbstractGroup>(grpPtr);
  const LOCA::Hopf::ComplexMultiVector& cB =
    dynamic_cast<const LOCA::Hopf::ComplexMultiVector&>(B);
  LOCA::Hopf::ComplexMultiVector& cX =
    dynamic_cast<LOCA::Hopf::ComplexMultiVector&>(X);

  if (magrp != Teuchos::null)
    return magrp->applyComplexTransposeInverseMultiVector(
                             params,
                             *(cB.getRealMultiVec()),
                             *(cB.getImagMultiVec()),
                             *(cX.getRealMultiVec()),
                             *(cX.getImagMultiVec()));
  else
    return NOX::Abstract::Group::NotDefined;
}
