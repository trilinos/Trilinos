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

#include "NOX_Linesearch_Manager.H" // base class

#include <string>

// Different line searches
#include "NOX_Linesearch_Fullstep.H"
#include "NOX_Linesearch_Halving.H"
#include "NOX_Linesearch_Polynomial.H"
#include "NOX_Linesearch_MoreThunte.H"

using namespace NOX;
using namespace NOX::Linesearch;

Manager::Manager(const Parameter::List& params) :
  ptr(NULL),
  method("")
{
  if (!isNewMethod(params))
    method = "Full Step";		// default line search

  newPtr(params);
}

Manager::~Manager()
{
  delete ptr;
}

void Manager::reset(const Parameter::List& params)
{
  if (isNewMethod(params))
    newPtr(params);
  
  ptr->reset(params);
}

bool Manager::operator()(Abstract::Group& newgrp, double& step, 
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) 
{
  return ptr->operator()(newgrp, step, oldgrp, dir);
}

// private
bool Manager::isNewMethod(const Parameter::List& params) 
{
  //! If a new method isn't specified, than nothing changed.
  if (!params.isParameter("Method"))
    return false;
    
  //! Check to see if the method name has changed or not
  string newmethod = params.getParameter("Method", "bogus");
  if (method == newmethod) 
    return false;

  // Method name has changed! Copy it and return true!
  method = newmethod;
  return true;
}

//private
void Manager::newPtr(const Parameter::List& params)
{
  delete ptr;

  if (method == "Full Step")
    ptr = new FullStep(params);
  else if (method == "Interval Halving")
    ptr = new Halving(params);
  else if (method == "Polynomial")
    ptr = new Polynomial(params);
  else if (method == "More Thunte")
    ptr = new MoreThunte(params);
  else {
    ptr = NULL;
    cout << "ERROR: invalid choice \"" << method << "\" for line search method "
	 << "in Manager constructor" << endl;
    throw 1;
  }
}

