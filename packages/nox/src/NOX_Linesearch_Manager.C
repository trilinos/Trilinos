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

#include "NOX_Linesearch_Manager.H" // class definition

#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Parameter_List.H"
#include "NOX_Utils.H"

// All the different line searches
#include "NOX_Linesearch_FullStep.H"
#include "NOX_Linesearch_Backtrack.H"
#include "NOX_Linesearch_Polynomial.H"
#include "NOX_Linesearch_MoreThuente.H"

using namespace NOX;
using namespace NOX::Linesearch;

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
   string newmethod = params.getParameter("Method", "Full Step");

  if (method != newmethod) {
    
    method = newmethod;
    
    delete ptr;
    
    if (method == "Full Step")
      ptr = new FullStep(params);
    else if ((method == "Interval Halving") // deprecated
	     || (method == "Backtrack"))
      ptr = new Backtrack(params);
    else if (method == "Polynomial")
      ptr = new Polynomial(params);
    else if (method == "More'-Thuente")
      ptr = new MoreThuente(params);
    else {
      ptr = NULL;
      cout << "ERROR: NOX::Linesearch::Manager - invalid choice \"" 
	   << method << "\" for line search method " << endl;
      throw "NOX Error";
    }
  }

  return ptr->reset(params);
}

bool Manager::operator()(Abstract::Group& newgrp, double& step, 
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) 
{
  return ptr->operator()(newgrp, step, oldgrp, dir);
}


