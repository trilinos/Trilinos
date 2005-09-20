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

#include "NOX_LineSearch_Manager.H" // class definition

#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"
#include "NOX_Parameter_LineSearchConstructor.H"

// All the different line searches
#include "NOX_LineSearch_FullStep.H"
#include "NOX_LineSearch_Backtrack.H"
#include "NOX_LineSearch_Polynomial.H"
#include "NOX_LineSearch_MoreThuente.H"
#ifdef WITH_PRERELEASE
#include "NOX_LineSearch_NonlinearCG.H"
#endif

NOX::LineSearch::Manager::Manager(const NOX::Utils& u, NOX::Parameter::List& params) :
  utils(u),
  method(""),
  ptr(NULL)
{
  reset(params);
}

NOX::LineSearch::Manager::~Manager()
{
  delete ptr;
}

bool NOX::LineSearch::Manager::reset(Parameter::List& params)
{
   string newmethod = params.getParameter("Method", "Full Step");

  // If the method has not changeed, just call reset on the method.
  if (method == newmethod) 
  {
    return ptr->reset(params);
  }

  method = newmethod;
  delete ptr;
  ptr = NULL;
    
  if (method == "Full Step")
    ptr = new FullStep(params);
  else if (method == "Backtrack")
    ptr = new Backtrack(utils, params);
  else if (method == "Polynomial")
    ptr = new Polynomial(utils, params);
  else if (method == "More'-Thuente")
    ptr = new MoreThuente(utils, params);
#ifdef WITH_PRERELEASE
  else if (method == "NonlinearCG")
    ptr = new NonlinearCG(utils, params);
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
    const NOX::Parameter::Arbitrary& ap = 
      params.getArbitraryParameter("User Defined Constructor");
    
    // Dynamically cast the Arbitrary Parameter to a LineSearchConstructor Parameter
    const NOX::Parameter::LineSearchConstructor* lscPtr = 
      dynamic_cast<const NOX::Parameter::LineSearchConstructor*>(&ap);
    
    // Check that the cast was successful
    if (lscPtr == NULL)
    {
      printWarning("reset", "Cannot do dynamic cast from Arbitrary to LineSearchConstructor");
      return false;
    }
    
    // Create a new direction from the LineSearchConstructor object
    ptr = lscPtr->newLineSearch(utils, params);
    
    // Check that the creation was successful
    if (ptr == NULL) 
    {
      printWarning("reset", "LineSearchConstructor object failed to create new direction");
      return false;
    }
  }
  else 
  {
    printWarning("reset", "invalid choice (" + method + ") for linesearch method");
    return false;
  }

  return (ptr != NULL);
}

bool NOX::LineSearch::Manager::compute(Abstract::Group& newgrp, double& step, 
		      const Abstract::Vector& dir,
		      const Solver::Generic& s) 
{
  return ptr->compute(newgrp, step, dir, s);
}


void NOX::LineSearch::Manager::printWarning(const string& name, const string& warning)
{
  if (utils.isPrintType(NOX::Utils::Warning)) 
    utils.out() << "Calling NOX::LineSearch::Manager::" << name << " - " << warning << endl;
}
