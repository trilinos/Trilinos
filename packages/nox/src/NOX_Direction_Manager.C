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

#include "NOX_Direction_Manager.H" // class definition
#include "NOX_Utils.H"
#include "NOX_Parameter_List.H"
#include "NOX_Direction_Newton.H"
#include "NOX_Direction_SteepestDescent.H"
#include "NOX_Parameter_DirectionConstructor.H"

#ifdef WITH_PRERELEASE
#include "NOX_Direction_NonlinearCG.H"
#include "NOX_Direction_Tensor.H"
#include "NOX_Direction_ModifiedNewton.H"
#include "NOX_Direction_QuasiNewton.H"
#include "NOX_Direction_Broyden.H"
#endif

NOX::Direction::Manager::Manager(const NOX::Utils& u) :
  utils(u),
  method(""),
  ptr(NULL)
{
}

NOX::Direction::Manager::Manager(const NOX::Utils& u, NOX::Parameter::List& params) :
  utils(u),
  method(""),
  ptr(NULL)
{
  reset(params);
}

NOX::Direction::Manager::~Manager()
{
  delete ptr;
}

bool NOX::Direction::Manager::reset(NOX::Parameter::List& params)
{
  string newmethod = params.getParameter("Method", "Newton");

  // If the method has not changeed, just call reset on the method.
  if (method == newmethod) 
  {
    return ptr->reset(params);
  }

  // Otherwise, construct a whole new direction
  method = newmethod;
  
  delete ptr;
  ptr = NULL;
  
  if (method == "Newton")
    ptr = new Newton(utils, params);
  else if (method == "Steepest Descent")
    ptr = new SteepestDescent(utils, params);
#ifdef WITH_PRERELEASE
  else if (method == "NonlinearCG")
    ptr = new NonlinearCG(utils, params);
  else if (method == "Tensor")
    ptr = new Tensor(utils, params);
  else if (method == "Modified-Newton")
    ptr = new ModifiedNewton(utils, params);
  else if (method == "Quasi-Newton")
    ptr = new QuasiNewton(utils, params);
  else if (method == "Broyden")
    ptr = new Broyden(utils, params);
#endif
  else if (method == "User Defined")
  {
    // Check that the corresponding Direction parameter exists
    if (!params.isParameterArbitrary("User Defined Constructor"))
    {
      printWarning("reset", "No \"User Defined Constructor\" specified");
      return false;
    }
    
    // Extract the Arbitrary Parameter
    const NOX::Parameter::Arbitrary& ap = params.getArbitraryParameter("User Defined Constructor");
    
    // Dynamically cast the Arbitrary Parameter to a DirectionConstructor Parameter
    const NOX::Parameter::DirectionConstructor* dcPtr = 
      dynamic_cast<const NOX::Parameter::DirectionConstructor*>(&ap);
    
    // Check that the cast was successful
    if (dcPtr == NULL)
    {
      printWarning("reset", "Cannot do dynamic cast from Arbitrary to DirectionConstructor");
      return false;
    }
    
    // Create a new direction from the DirectionConstructor object
    ptr = dcPtr->newDirection(utils, params);
    
    // Check that the creation was successful
    if (ptr == NULL) 
    {
      printWarning("reset", "DirectionConstructor object failed to create new direction");
      return false;
    }
  }
  else 
  {
    printWarning("reset", "invalid choice (" + method + ") for direction method");
    return false;
  }

  return (ptr != NULL);
}


bool NOX::Direction::Manager::compute(Abstract::Vector& dir, Abstract::Group& grp, 
		      const Solver::Generic& solver) 
{
  if (ptr == NULL) 
  {
    if (utils.doPrint(NOX::Utils::Warning)) 
      cout << "Calling NOX::Direction::Manager::compute on uninitialized direction" << endl;
    return false;
  }

  return ptr->compute(dir, grp, solver);
}

bool NOX::Direction::Manager::compute(Abstract::Vector& dir, Abstract::Group& grp, 
		      const Solver::LineSearchBased& solver) 
{
  if (ptr == NULL) 
  {
    if (utils.doPrint(NOX::Utils::Warning)) 
      cout << "Calling NOX::Direction::Manager::compute on uninitialized direction" << endl;
    return false;
  }

  return ptr->compute(dir, grp, solver);
}

void NOX::Direction::Manager::printWarning(const string& name, const string& warning)
{
  if (utils.doPrint(NOX::Utils::Warning)) 
    cout << "Calling NOX::Direction::Manager::" << name << " - " << warning << endl;
}
