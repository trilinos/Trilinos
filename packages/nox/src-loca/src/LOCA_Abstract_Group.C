// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "LOCA_Abstract_Group.H"
#include "NOX_Parameter_List.H"

LOCA::Abstract::Group::Group(const LOCA::DerivUtils& d)
  : LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup(d),
    LOCA::Bifurcation::TPBord::SingularSolveGroup()
{
}

LOCA::Abstract::Group::Group(NOX::Parameter::List& params, 
			     const LOCA::DerivUtils& d)
  : LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup(d),
    LOCA::Bifurcation::TPBord::SingularSolveGroup(params)
{
}

LOCA::Abstract::Group::Group(const LOCA::Abstract::Group& source, 
			     NOX::CopyType type)
  : LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup(source, type),
    LOCA::Bifurcation::TPBord::SingularSolveGroup(source, type)
{
}


LOCA::Abstract::Group::~Group() 
{
}

LOCA::Continuation::AbstractGroup&
LOCA::Abstract::Group::operator=(
			    const LOCA::Continuation::AbstractGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::Continuation::FiniteDifferenceGroup&
LOCA::Abstract::Group::operator=(
		      const LOCA::Continuation::FiniteDifferenceGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::Bifurcation::TPBord::AbstractGroup&
LOCA::Abstract::Group::operator=(
		     const LOCA::Bifurcation::TPBord::AbstractGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::Bifurcation::TPBord::FiniteDifferenceGroup&
LOCA::Abstract::Group::operator=(
	      const LOCA::Bifurcation::TPBord::FiniteDifferenceGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::Bifurcation::TPBord::SingularSolveGroup&
LOCA::Abstract::Group::operator=(
		   const LOCA::Bifurcation::TPBord::SingularSolveGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::TimeDependent::AbstractGroup&
LOCA::Abstract::Group::operator=(
			    const LOCA::TimeDependent::AbstractGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::Bifurcation::HopfBord::AbstractGroup&
LOCA::Abstract::Group::operator=(
	       const LOCA::Bifurcation::HopfBord::AbstractGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup&
LOCA::Abstract::Group::operator=(
	       const LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::Homotopy::AbstractGroup&
LOCA::Abstract::Group::operator=(
			    const LOCA::Homotopy::AbstractGroup& source)
{
  return operator=(dynamic_cast<const LOCA::Abstract::Group&>(source));
}

LOCA::Abstract::Group&
LOCA::Abstract::Group::operator=(const LOCA::Abstract::Group& source)
{

  // Copy parent classes
  LOCA::Bifurcation::HopfBord::FiniteDifferenceGroup::operator=(source);
  LOCA::Bifurcation::TPBord::SingularSolveGroup::operator=(source);
  
  return *this;
}

NOX::Abstract::Group::ReturnType 
LOCA::Abstract::Group::augmentJacobianForHomotopy(double conParamValue)
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::computeMassMatrix()
{
  errorCheck.throwError("LOCA::Abstract::Group::computeMassMatrix",
			"No mass matrix defined for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyMassMatrix(const NOX::Abstract::Vector& input,
				       NOX::Abstract::Vector& result) const
{
  errorCheck.throwError("LOCA::Abstract::Group::applyMassMatrix",
			       "No mass matrix defined for group");
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::Abstract::Group::applyComplexInverse(
			       NOX::Parameter::List& params,
			       const NOX::Abstract::Vector& input_real,
			       const NOX::Abstract::Vector& input_imag,
			       double frequency,
			       NOX::Abstract::Vector& result_real,
			       NOX::Abstract::Vector& result_imag) const
{
  errorCheck.throwError("LOCA::Abstract::Group::applyComplexInverse",
			"No mass matrix defined for group");
  return NOX::Abstract::Group::NotDefined;
}

