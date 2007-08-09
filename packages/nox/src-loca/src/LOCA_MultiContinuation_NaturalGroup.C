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

#include "Teuchos_ParameterList.hpp"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_NaturalConstraint.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& continuationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs)
  : LOCA::MultiContinuation::ExtendedGroup(global_data, topParams,
					   continuationParams,
					   grp, pred, paramIDs)
{
  bool skip_dfdp = continuationParams->get("Skip Parameter Derivative", true);
  Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> cons 
    = Teuchos::rcp(new LOCA::MultiContinuation::NaturalConstraint(
	globalData, Teuchos::rcp(this, false)));
  LOCA::MultiContinuation::ExtendedGroup::setConstraints(cons, skip_dfdp);
}

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
			 const LOCA::MultiContinuation::NaturalGroup& source,
			 NOX::CopyType type)
  : LOCA::MultiContinuation::ExtendedGroup(source, type)
{
  Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::NaturalConstraint>(conGroup->getConstraints())->setNaturalGroup(Teuchos::rcp(this, false));
}


LOCA::MultiContinuation::NaturalGroup::~NaturalGroup() 
{
}

NOX::Abstract::Group&
LOCA::MultiContinuation::NaturalGroup::operator=(
					  const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::MultiContinuation::NaturalGroup::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new LOCA::MultiContinuation::NaturalGroup(*this, type));
}

void
LOCA::MultiContinuation::NaturalGroup::copy(const NOX::Abstract::Group& src) 
{

  // Protect against A = A
  if (this != &src) {
    LOCA::MultiContinuation::ExtendedGroup::copy(src);
  }
}


