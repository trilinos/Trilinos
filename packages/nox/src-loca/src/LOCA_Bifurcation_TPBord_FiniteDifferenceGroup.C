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

#include "LOCA_Bifurcation_TPBord_FiniteDifferenceGroup.H"

LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::FiniteDifferenceGroup(
						     const LOCA::DerivUtils& d)
  : LOCA::Continuation::FiniteDifferenceGroup(d)
{
}

LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::FiniteDifferenceGroup(
               const LOCA::Bifurcation::TPBord::FiniteDifferenceGroup& source, 
	       NOX::CopyType type)
  : LOCA::Continuation::FiniteDifferenceGroup(source, type)
{
}


LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::~FiniteDifferenceGroup() 
{
}

LOCA::Bifurcation::TPBord::FiniteDifferenceGroup&
LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::operator=(
	     const LOCA::Bifurcation::TPBord::FiniteDifferenceGroup& source)
{
  LOCA::Continuation::FiniteDifferenceGroup::operator=(source);

  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::computeDJnDp(
				      const NOX::Abstract::Vector& nullVector,
				      const int param_id, 
				      NOX::Abstract::Vector& result)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDp(*this, nullVector, param_id, result);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::computeDJnDp(
				      const NOX::Abstract::Vector& nullVector,
				      const int param_id, 
				      const NOX::Abstract::Vector& JnVector,
				      NOX::Abstract::Vector& result)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDp(*this, nullVector, param_id, JnVector, result);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::computeDJnDxa(
				      const NOX::Abstract::Vector& nullVector,
				      const NOX::Abstract::Vector& aVector, 
				      NOX::Abstract::Vector& result)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDxa(*this, nullVector, aVector, result);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::computeDJnDxa(
				      const NOX::Abstract::Vector& nullVector,
				      const NOX::Abstract::Vector& aVector, 
				      const NOX::Abstract::Vector& JnVector,
				      NOX::Abstract::Vector& result)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDxa(*this, nullVector, aVector, JnVector, result);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::computeDJnDxa(
				    const NOX::Abstract::Vector& nullVector,
				    const NOX::Abstract::MultiVector& aVector, 
				    NOX::Abstract::MultiVector& result)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDxa(*this, nullVector, aVector, result);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::computeDJnDxa(
				    const NOX::Abstract::Vector& nullVector,
				    const NOX::Abstract::MultiVector& aVector, 
				    const NOX::Abstract::Vector& JnVector,
				    NOX::Abstract::MultiVector& result)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDJnDxa(*this, nullVector, aVector, JnVector, result);
}
