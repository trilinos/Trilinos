// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

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
