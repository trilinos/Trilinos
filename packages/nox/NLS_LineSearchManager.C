// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include <string>
#include "NLS_LineSearchManager.H"
#include "NLS_FullStep.H"
#include "NLS_IntervalHalving.H"

NLS_LineSearchManager::NLS_LineSearchManager(const NLS_ParameterList& params) :
  ptr(NULL),
  method("")
{
  if (!isNewMethod(params))
    method = "FullStep";		// default line search

  newPtr(params);
}

NLS_LineSearchManager::~NLS_LineSearchManager()
{
  delete ptr;
}

void NLS_LineSearchManager::reset(const NLS_ParameterList& params)
{
  if (isNewMethod(params))
    newPtr(params);
  
  ptr->reset(params);
}

bool NLS_LineSearchManager::search(const NLS_Group& oldgrp, const NLS_Vector& dir, 
				   NLS_Group& newgrp) const
{
  return ptr->search(oldgrp, dir, newgrp);
}

// private
bool NLS_LineSearchManager::isNewMethod(const NLS_ParameterList& params) 
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
void NLS_LineSearchManager::newPtr(const NLS_ParameterList& params)
{
  delete ptr;

  if (method == "Full Step")
    ptr = new NLS_FullStep(params);
  else if (method == "Interval Halving")
    ptr = new NLS_IntervalHalving(params);
  else {
    ptr = NULL;
    cout << "ERROR: invalid choice for line search method "
	 << "in NLS_LineSearchManager constructor" << endl;
    throw 1;
  }
}

