// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_NaturalConstraint.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
      const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& continuationParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs)
  : LOCA::MultiContinuation::ExtendedGroup(global_data, topParams,
					   continuationParams,
					   grp, pred, paramIDs)
{
  Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface> cons 
    = Teuchos::rcp(new LOCA::MultiContinuation::NaturalConstraint(
	globalData, Teuchos::rcp(this, false)));
  LOCA::MultiContinuation::ExtendedGroup::setConstraints(cons);
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

Teuchos::RefCountPtr<NOX::Abstract::Group>
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


