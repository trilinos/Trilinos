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

#include "NOX_Direction_Manager.H" // class definition

// All the different direction methods 
#include "NOX_Direction_Newton.H"
#include "NOX_Direction_NonlinearCG.H"
#include "NOX_Direction_SteepestDescent.H"
#include "NOX_Direction_Tensor.H"

#include "NOX_Utils.H"
#include "NOX_Parameter_List.H"

using namespace NOX;
using namespace NOX::Direction;

Manager::Manager() :
  method(""),
  ptr(NULL)
{
}

Manager::Manager(Parameter::List& params) :
  method(""),
  ptr(NULL)
{
  reset(params);
}

Manager::~Manager()
{
  delete ptr;
}

bool Manager::reset(Parameter::List& params)
{
  string newmethod = params.getParameter("Method", "Newton");

  if (method != newmethod) {
    
    method = newmethod;
    
    delete ptr;
    
    if (method == "Newton")
      ptr = new Newton(params);
    else if (method == "NonlinearCG")
      ptr = new NonlinearCG(params);
    else if (method == "Steepest Descent")
      ptr = new SteepestDescent(params);
    else if (method == "Tensor")
      ptr = new Tensor(params);
    else {
      ptr = NULL;
      if (Utils::doPrint(NOX::Utils::Warning)) {
	cerr << "NOX::Direction::Manager::reset() - invalid choice (" 
	     << method << ") for direction method " << endl;
      }
      return false;
    }
  }

  return ptr->reset(params);
}

bool Manager::compute(Abstract::Vector& dir, Abstract::Group& grp, 
			 const Solver::Generic& solver) 
{
  if (ptr == NULL) {
    if (Utils::doPrint(NOX::Utils::Warning)) 
      cout << "Calling NOX::Direction::Manager::compute on uninitialized direction" << endl;
    return false;
  }

  return ptr->compute(dir, grp, solver);
}

