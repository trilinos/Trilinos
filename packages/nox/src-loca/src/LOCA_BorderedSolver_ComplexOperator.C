// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
