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

#include "NOX_Solver_Manager.H"	// class definition
#include "NOX_Utils.H"		// for static function doPrint

// Header files for different solvers
#include "NOX_Solver_LineSearchBased.H"	 // LineSearch method
#include "NOX_Solver_TrustRegionBased.H" // Trust region method
#ifdef WITH_PRERELEASE
#include "NOX_Solver_InexactTrustRegionBased.H" // Inexact Trust region method
#include "NOX_Solver_TensorBased.H"  // Tensor method
#include "NOX_Solver_TensorBasedTest.H"  // Tensor-Krylov method
#endif

NOX::Solver::Manager::Manager(Abstract::Group& grp, StatusTest::Generic &t, Parameter::List& p) :
  method(""),
  ptr(NULL)
{
  reset(grp, t, p);
}

NOX::Solver::Manager::Manager() :
  method(""),
  ptr(NULL)
{
}

NOX::Solver::Manager::~Manager()
{
  delete ptr;
}

bool NOX::Solver::Manager::reset(Abstract::Group& grp, 
				 StatusTest::Generic& tests, 
				 Parameter::List& params)
{
  string newmethod = params.getParameter("Nonlinear Solver", "Line Search Based");

  if ((method == newmethod) && (ptr != NULL))
  {
    return ptr->reset(grp, tests, params);
  }
  else 
  {
    method = newmethod;

    delete ptr;
    ptr = NULL;
    
    if ((method == "Newton") || (method == "Line Search")) // deprecated
    {	
      deprecated(method, "Line Search Based");
      ptr = new LineSearchBased(grp, tests, params);
    } 
    else if (method == "Line Search Based") 
    {
      ptr = new LineSearchBased(grp, tests, params);
    } 
    else if (method == "Trust Region")  // deprecated
    {
      deprecated(method, "Trust Region Based");
      ptr = new TrustRegionBased(grp, tests, params);
    } 
    else if (method == "Trust Region Based") 
    {
      ptr = new TrustRegionBased(grp, tests, params);
    } 
#ifdef WITH_PRERELEASE
    else if (method == "Inexact Trust Region Based") 
    {
      ptr = new InexactTrustRegionBased(grp, tests, params);
    } 
    else if (method == "Tensor Based") 
    {
      ptr = new TensorBased(grp, tests, params);
    } 
    else if (method == "Tensor-Krylov Based") 
    {
      ptr = new TensorBasedTest(grp, tests, params);
    } 
#endif
    else 
    {
      cout << "ERROR: NOX::Solver::Manager::reset - Invalid solver choice " << method << endl;
      throw "NOX Error";
    }

    if (ptr == NULL) 
    {
      cerr << "NOX::Solver::Manager::reset - Null pointer error" << endl;
      return false;
    }

    return true;
  }
}

bool NOX::Solver::Manager::reset(Abstract::Group& grp, 
				 StatusTest::Generic& tests)
{
  return ptr->reset(grp, tests);
}

// PRIVATE
void NOX::Solver::Manager::deprecated(const string& oldName, const string& newName) const
{
  cout << "Warning: NOX::Solver::Manager::reset - " 
       << "Nonlinear Solver choice \"" << oldName << "\" is deprecated.\n"
       << "                                       " 
       << "Use \"" << newName << "\" instead." 
       << endl;
}

NOX::StatusTest::StatusType NOX::Solver::Manager::getStatus()
{
  checkNullPtr("getStatus");
  return ptr->getStatus();
}

NOX::StatusTest::StatusType NOX::Solver::Manager::iterate()
{
  checkNullPtr("iterate");
  return ptr->iterate();
}

NOX::StatusTest::StatusType NOX::Solver::Manager::solve()
{
  checkNullPtr("solve");
  return ptr->solve();
}

const NOX::Abstract::Group& NOX::Solver::Manager::getSolutionGroup() const
{
  checkNullPtr("getSolutionGroup");
  return ptr->getSolutionGroup();
}

const NOX::Abstract::Group& NOX::Solver::Manager::getPreviousSolutionGroup() const
{
  checkNullPtr("getPreviousSolutionGroup");
  return ptr->getPreviousSolutionGroup();
}

int NOX::Solver::Manager::getNumIterations() const
{
  if (ptr == NULL)
    return 0;

  return ptr->getNumIterations();
}

const NOX::Parameter::List& NOX::Solver::Manager::getParameterList() const
{
  checkNullPtr("getParameterList");
  return ptr->getParameterList();
}

// PRIVATE
void NOX::Solver::Manager::checkNullPtr(const string& fname) const
{
  if (ptr == NULL) 
  {
    cout << "NOX::Solver::Manager::" << fname << " - Null pointer error" << endl;
    throw "NOX Error";
  }
}

