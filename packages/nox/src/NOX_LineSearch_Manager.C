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

#include "NOX_LineSearch_Manager.H" // class definition

#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

// All the different line searches
#include "NOX_LineSearch_FullStep.H"
#include "NOX_LineSearch_Backtrack.H"
#include "NOX_LineSearch_Polynomial.H"
#include "NOX_LineSearch_MoreThuente.H"
#include "NOX_LineSearch_Secant.H"
#ifdef WITH_PRERELEASE
#include "NOX_LineSearch_Quadratic.H"
#include "NOX_LineSearch_MoreThuente2.H"
#endif

using namespace NOX;
using namespace NOX::LineSearch;

Manager::Manager(const NOX::Utils& u, NOX::Parameter::List& params) :
  utils(u),
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
   string newmethod = params.getParameter("Method", "Full Step");

  if (method != newmethod) {
    
    method = newmethod;
    
    delete ptr;
    
    if (method == "Full Step")
      ptr = new FullStep(params);
    else if (method == "Backtrack")
      ptr = new Backtrack(utils, params);
    else if (method == "Polynomial")
      ptr = new Polynomial(utils, params);
    else if (method == "More'-Thuente")
      ptr = new MoreThuente(utils, params);
    else if (method == "Secant")
      ptr = new Secant(utils, params);
#ifdef WITH_PRERELEASE
    else if (method == "Quadratic")
      ptr = new Quadratic(utils, params);
    else if (method == "More'-Thuente2")
      ptr = new MoreThuente2(utils, params);
#endif
    else {
      ptr = NULL;
      cout << "ERROR: NOX::LineSearch::Manager - invalid choice \"" 
	   << method << "\" for line search method " << endl;
      throw "NOX Error";
    }
  }

  return ptr->reset(params);
}

bool Manager::compute(Abstract::Group& newgrp, double& step, 
		      const Abstract::Vector& dir,
		      const Solver::Generic& s) 
{
  return ptr->compute(newgrp, step, dir, s);
}


