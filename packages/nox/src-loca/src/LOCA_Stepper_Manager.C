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

#include "LOCA_Stepper_Manager.H"	// class definition

//NOX includes
#include "NOX_Parameter_List.H"

// LOCA Includes
#include "LOCA_Utils.H"              // for static function doPrint
#include "LOCA_Solver_Generic.H"

// Header files for different stepper algorithms
#include "LOCA_Stepper_ZeroOrder.H"
//#include "LOCA_Stepper_FirstOrder.H"
//#include "LOCA_Stepper_ArcLength.H"

using namespace LOCA;
using namespace LOCA::Stepper;

Manager::Manager(Solver::Generic& s) :
  Generic(s), // ctor for base class LOCA::Stepper::Generic
  method(""),
  ptr(0)
{
  reset(s);
}

Manager::~Manager()
{
  delete ptr;
}

bool Manager::reset(Solver::Generic& s)
{
  // reset the generic base class members
  resetGenericMembers(s);

  NOX::Parameter::List& p = s.getParameterList().sublist("Stepper");

  // check to see if the method has changed.  defaults to zero order cont
  string newmethod = p.getParameter("Stepper Method", "Zero Order");

  if (method != newmethod) {
    
    method = newmethod;

    delete ptr;
    ptr = 0;
    
    if (method == "Zero Order") {
      ptr = new ZeroOrder(s);
      cout << "Constructing a \"Zero Order\" method!" << endl;
    } 
    else if (method == "First Order") {
      //ptr = new FirstOrder(s);
      cout << "Constructing a \"First Order\" method!" << endl;
    } 
    else if (method == "Arc-length") {
      //ptr = new ArcLength(s);
      cout << "Constructing an \"Arc-Length\" method!" << endl;
    } 
    else {
      cout << "ERROR: LOCA::Stepper::Manager::reset() - \"Stepper "
	   << "Method\" parameter: \"" << method << "\" is invalid!" << endl;
      throw "LOCA Error";
    }

    // Make sure we allocated the object
    if (ptr == NULL) {
      cout << "ERROR: LOCA::Stepper::Manager::reset() - Null pointer error" 
	   << endl;
      throw "LOCA Error";
    }

    return true;
  }
  else {

    if (ptr == NULL) {
      cout << "ERROR: LOCA::Solver::Manager::reset() - Null pointer error" 
	   << endl;
      return false;
    }

    return ptr->reset(s);
  }
}

bool Manager::resetFromFailedStep()
{
  if (ptr == NULL) {
    cout << "ERROR: LOCA::Solver::Manager::resetFromFailedStep() - Null pointer error" 
	 << endl;
    throw "LOCA Error";
  }
  return ptr->resetFromFailedStep();
}

StatusType Manager::getStatus()
{
  if (ptr == NULL) {
    cout << "ERROR: LOCA::Solver::Manager::getStatus() - Null pointer error" 
	 << endl;
    throw "LOCA Error";
  }

  return ptr->getStatus();
}

StatusType Manager::solve()
{
  if (ptr == NULL) {
    cout << "ERROR: LOCA::Solver::Manager::solve() - Null pointer error" 
	 << endl;
    throw "LOCA Error";
  }

  return ptr->solve();
}

StatusType Manager::step()
{
  if (ptr == NULL) {
    cout << "ERROR: LOCA::Solver::Manager::step() - Null pointer error" 
	 << endl;
    throw "LOCA Error";
  }

  return ptr->step();
}

const Abstract::Group& Manager::getSolutionGroup() const
{
  if (ptr == NULL) {
    cout << "ERROR: LOCA::Solver::Manager::getSolutionGroup() - Null "
	 << "pointer error" << endl;
    throw "LOCA Error";
  }

  return ptr->getSolutionGroup();
}

const Abstract::Group& Manager::getPreviousSolutionGroup() const
{
  if (ptr == NULL) {
    cout << "ERROR: LOCA::Solver::Manager::getPreviousSolutionGroup() - Null "
	 << "pointer error" << endl;
    throw "LOCA Error";
  }

  return ptr->getPreviousSolutionGroup();
}

int Manager::getNumContinuationSteps() const
{
  if (ptr == NULL) {
    cout << "ERROR: LOCA::Solver::Manager::getPreviousSolutionGroup() - Null "
	 << "pointer error" << endl;
    throw "LOCA Error";
  }

  return ptr->getNumContinuationSteps();
}

const NOX::Parameter::List& Manager::getParameterList() const
{
  if (ptr == NULL) {
    cout << "ERROR: LOCA::Solver::Manager::getPreviousSolutionGroup() - Null "
	 << "pointer error" << endl;
    throw "LOCA Error";
  }
    
  return ptr->getParameterList();
}
