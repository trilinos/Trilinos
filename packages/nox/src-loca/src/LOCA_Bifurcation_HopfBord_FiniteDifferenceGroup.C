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

#include "LOCA_Bifurcation_HopfBord_FiniteDifferenceGroup.H"

LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::FiniteDifferenceGroup(
						 const LOCA::DerivUtils& d)
  : LOCA::Bifurcation::TPBord::FiniteDifferenceGroup(d)
{
}

LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::FiniteDifferenceGroup(
            const LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup& source, 
	    NOX::CopyType type)
  : LOCA::Bifurcation::TPBord::FiniteDifferenceGroup(source, type)
{
}


LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::~FiniteDifferenceGroup() 
{
}

LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup&
LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::operator=(
	    const LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup& source)
{
  LOCA::Bifurcation::TPBord::FiniteDifferenceGroup::operator=(source);

  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::computeDCeDp(
				const NOX::Abstract::Vector& yVector,
				const NOX::Abstract::Vector& zVector,
				double w,
				const int param_id, 
				NOX::Abstract::Vector& result_real,
				NOX::Abstract::Vector& result_imag)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDp(*this, yVector, zVector, w, param_id, result_real, 
		 result_imag);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::computeDCeDp(
				const NOX::Abstract::Vector& yVector,
				const NOX::Abstract::Vector& zVector,
				double w,
				const int param_id, 
				const NOX::Abstract::Vector& Ce_real,
				const NOX::Abstract::Vector& Ce_imag,
				NOX::Abstract::Vector& result_real,
				NOX::Abstract::Vector& result_imag)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDp(*this, yVector, zVector, w, param_id, Ce_real, Ce_imag,
		 result_real, result_imag);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::computeDCeDxa(
				   const NOX::Abstract::Vector& yVector,
				   const NOX::Abstract::Vector& zVector,
				   double w,
				   const NOX::Abstract::Vector& aVector,
				   NOX::Abstract::Vector& result_real,
				   NOX::Abstract::Vector& result_imag)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDxa(*this, yVector, zVector, w, aVector,
		  result_real, result_imag);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::computeDCeDxa(
				   const NOX::Abstract::Vector& yVector,
				   const NOX::Abstract::Vector& zVector,
				   double w,
				   const NOX::Abstract::Vector& aVector,
				   const NOX::Abstract::Vector& Ce_real,
				   const NOX::Abstract::Vector& Ce_imag,
				   NOX::Abstract::Vector& result_real,
				   NOX::Abstract::Vector& result_imag)
{
  return LOCA::Continuation::FiniteDifferenceGroup::derivPtr->
    computeDCeDxa(*this, yVector, zVector, w, aVector, Ce_real, Ce_imag,
		  result_real, result_imag);
}
