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

#include "NOX_Parameter_List.H"
#include "LOCA_Bifurcation_TPBord_SingularSolveGroup.H"

LOCA::Bifurcation::TPBord::SingularSolveGroup::SingularSolveGroup(
					   const NOX::Parameter::List& params)
  : singularSolveManager(params)
{
}

LOCA::Bifurcation::TPBord::SingularSolveGroup::SingularSolveGroup(
						 NOX::Parameter::List& params)
  : singularSolveManager(params.sublist("Singular Solve"))
{
}

LOCA::Bifurcation::TPBord::SingularSolveGroup::SingularSolveGroup(
                 const LOCA::Bifurcation::TPBord::SingularSolveGroup& source, 
		 NOX::CopyType type)
  : singularSolveManager(source.singularSolveManager)
{
}


LOCA::Bifurcation::TPBord::SingularSolveGroup::~SingularSolveGroup() 
{
}

LOCA::Bifurcation::TPBord::SingularSolveGroup&
LOCA::Bifurcation::TPBord::SingularSolveGroup::operator=(
		  const LOCA::Bifurcation::TPBord::SingularSolveGroup& source)
{
  singularSolveManager = source.singularSolveManager;

  return *this;
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::SingularSolveGroup::applySingularJacobianInverse(
                              NOX::Parameter::List& params,
			      const NOX::Abstract::Vector& input,
			      const NOX::Abstract::Vector& approxNullVec,
                              const NOX::Abstract::Vector& jacApproxNullVec,
			      NOX::Abstract::Vector& result)
{
  return singularSolveManager.compute(params, *this, input, approxNullVec, 
				      jacApproxNullVec, result);
}

NOX::Abstract::Group::ReturnType
LOCA::Bifurcation::TPBord::SingularSolveGroup::applySingularJacobianInverseMulti(
                          NOX::Parameter::List& params,
                          const NOX::Abstract::Vector*const* inputs,
                          const NOX::Abstract::Vector& approxNullVec,
                          const NOX::Abstract::Vector& jacApproxNullVec,
                          NOX::Abstract::Vector** results, int nVecs)
{
  return singularSolveManager.computeMulti(params, *this, inputs, 
					   approxNullVec, jacApproxNullVec, 
					   results, nVecs);
}
