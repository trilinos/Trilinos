// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Linesearch_Manager.H" // base class

#include <string>

// Different line searches
#include "NOX_Linesearch_Fullstep.H"
#include "NOX_Linesearch_Halving.H"

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
			 const Abstract::Group& oldgrp, const Abstract::Vector& dir) const
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
  else {
    ptr = NULL;
    cout << "ERROR: invalid choice \"" << method << "\" for line search method "
	 << "in Manager constructor" << endl;
    throw 1;
  }
}

