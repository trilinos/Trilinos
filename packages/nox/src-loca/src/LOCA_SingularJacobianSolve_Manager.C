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

#include "LOCA_Continuation_AbstractGroup.H"
#include "NOX_Parameter_List.H"
#include "LOCA_SingularJacobianSolve_Manager.H"
#include "LOCA_SingularJacobianSolve_Default.H"
#include "LOCA_SingularJacobianSolve_Nic.H"
#include "LOCA_SingularJacobianSolve_NicDay.H"
#include "LOCA_SingularJacobianSolve_ItRef.H"

LOCA::SingularJacobianSolve::Manager::Manager(NOX::Parameter::List& params) :
  method(),
  singularSolverPtr(NULL)
{
  reset(params);
}

LOCA::SingularJacobianSolve::Manager::Manager(const NOX::Parameter::List& params) :
  method(),
  singularSolverPtr(NULL)
{
  NOX::Parameter::List p(params);
  reset(p);
}

LOCA::SingularJacobianSolve::Manager::Manager(
			  const LOCA::SingularJacobianSolve::Manager& source) 
  : method(source.method),
    singularSolverPtr(source.singularSolverPtr->clone())
{
}

LOCA::SingularJacobianSolve::Manager::~Manager()
{
  delete singularSolverPtr;
}

LOCA::SingularJacobianSolve::Generic*
LOCA::SingularJacobianSolve::Manager::clone() const 
{
  return new Manager(*this);
}

LOCA::SingularJacobianSolve::Generic&
LOCA::SingularJacobianSolve::Manager::operator=(
			  const LOCA::SingularJacobianSolve::Generic& source)
{
  return operator=(dynamic_cast<const LOCA::SingularJacobianSolve::Manager&>(source));
}

LOCA::SingularJacobianSolve::Manager&
LOCA::SingularJacobianSolve::Manager::operator=(
			  const LOCA::SingularJacobianSolve::Manager& source)
{
  // protect against A = A
  if (this != &source) {
    delete singularSolverPtr;

    method = source.method;
    singularSolverPtr = source.singularSolverPtr->clone();
  }

  return *this;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Manager::reset(NOX::Parameter::List& params) 
{
  string newmethod = params.getParameter("Method", "Default");

  if (method != newmethod) {
    delete singularSolverPtr;

    method = newmethod;

    if (method == "Default")
      singularSolverPtr = new LOCA::SingularJacobianSolve::Default(params);
    else if (method == "Nic")
      singularSolverPtr = new LOCA::SingularJacobianSolve::Nic(params);
    else if (method == "Nic-Day")
      singularSolverPtr = new LOCA::SingularJacobianSolve::NicDay(params);
    else if (method == "Iterative Refinement")
      singularSolverPtr = new LOCA::SingularJacobianSolve::ItRef(params);
    else {
      cerr << "LOCA::SingularJacobianSolve::Manager::reset() - invalid choice (" 
	   << method << ") for singular solve method " << endl;
      return NOX::Abstract::Group::Failed;
    }
  }

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Manager::compute(
				NOX::Parameter::List& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& input,
			        const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result) 
{
  if (singularSolverPtr == NULL) {
    cerr << "LOCA::SingularJacobianSolve::Manager::compute - Null pointer error" << endl;
    return NOX::Abstract::Group::Failed;
  }

  cout << "\n\tCalling singular solver with method: " << method << endl;

  return singularSolverPtr->compute(params, grp, input, approxNullVec,
				    jacApproxNullVec, result);
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Manager::computeMulti(
				NOX::Parameter::List& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector*const* inputs,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector** results,
				int nVecs) 
{
  if (singularSolverPtr == NULL) {
    cerr << "LOCA::SingularJacobianSolve::Manager::compute - Null pointer error" << endl;
    return NOX::Abstract::Group::Failed;
  }

  cout << "\n\tCalling singular solver with method: " << method << endl;

  return singularSolverPtr->computeMulti(params, grp, inputs, approxNullVec,
					 jacApproxNullVec, results, nVecs);
}

const string&
LOCA::SingularJacobianSolve::Manager::getMethod() const 
{
  return method;
}
