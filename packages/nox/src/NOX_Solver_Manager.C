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

#include "NOX_Solver_Manager.H"	// class definition
#include "NOX_Utils.H"		// for static function doPrint

// Header files for different solvers
#include "NOX_Solver_LineSearchBased.H"	 // LineSearch method
#include "NOX_Solver_TrustRegionBased.H" // Trust region method
#include "NOX_Solver_NonlinearCG.H"      // Nonlinear Conjugate Gradient method

using namespace NOX;
using namespace NOX::Solver;

Manager::Manager(Abstract::Group& grp, StatusTest::Generic &t, const Parameter::List& p) :
  method(""),
  ptr(NULL)
{
  reset(grp, t, p);
}

Manager::~Manager()
{
  delete ptr;
}

bool Manager::reset(Abstract::Group& grp, StatusTest::Generic& tests, const Parameter::List& params)
{
  string newmethod = params.getParameter("Nonlinear Solver", "Newton");

  if (method != newmethod) {
    
    method = newmethod;

    delete ptr;
    ptr = NULL;
    
    if (method == "Newton") {
      ptr = new LineSearchBased(grp, tests, params);
    } 
    else if (method == "Line Search") {
      ptr = new LineSearchBased(grp, tests, params);
    } 
    else if (method == "NonlinearCG") {
      ptr = new NonlinearCG(grp, tests, params);
    } 
    else if (method == "Trust Region") {
      ptr = new TrustRegionBased(grp, tests, params);
    } 
    else {
      cout << "ERROR: NOX::Solver::Manager - Invalid solver choice" << endl;
      throw "NOX Error";
    }

    if (ptr == NULL) {
      cerr << "NOX::Solver::Manager::reset - Null pointer error" << endl;
      return false;
    }

    return true;
  }
  else {

    if (ptr == NULL) {
      cerr << "NOX::Solver::Manager::reset - Null pointer error" << endl;
      return false;
    }

    return ptr->reset(grp, tests, params);
  }
}

NOX::StatusTest::StatusType Manager::getStatus()
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::getStatus - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->getStatus();
}

NOX::StatusTest::StatusType Manager::iterate()
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::iterate - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->iterate();
}

NOX::StatusTest::StatusType Manager::solve()
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::solve - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->solve();
}

const Abstract::Group& Manager::getSolutionGroup() const
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::getSolutionGroup - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->getSolutionGroup();
}

const Abstract::Group& Manager::getPreviousSolutionGroup() const
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::getPreviousSolutionGroup - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->getPreviousSolutionGroup();
}

int Manager::getNumIterations() const
{
  if (ptr == NULL)
    return 0;

  return ptr->getNumIterations();
}

const Parameter::List& Manager::getOutputParameters() const
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::getOutputParameters - Null pointer error" << endl;
    throw "NOX Error";
  }
    
  return ptr->getOutputParameters();
}




