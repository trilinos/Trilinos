// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NOX_Solver_Manager.H"	// class definition
#include "NOX_Utils.H"		// for static function doPrint

// Header files for different solvers
#include "NOX_Solver_Newton.H"		// Newton's method
#include "NOX_Solver_NonlinearCG.H"	// Nonlinear Conjugate Gradient method
#include "NOX_Solver_TrustRegion.H"     // Trust region method

using namespace NOX;
using namespace NOX::Solver;

Manager::Manager(Abstract::Group& grp, Status::Test &t, const Parameter::List& p) :
  method(""),
  ptr(NULL)
{
  reset(grp, t, p);
}

Manager::~Manager()
{
  delete ptr;
}

bool Manager::reset(Abstract::Group& grp, Status::Test& tests, const Parameter::List& params)
{
  string newmethod = params.getParameter("Nonlinear Solver", "Newton");

  if (method != newmethod) {
    
    method = newmethod;

    delete ptr;
    
    if (method == "Newton") {
      ptr = new Newton(grp, tests, params);
    } 
    else if (method == "NonlinearCG") {
      ptr = new NonlinearCG(grp, tests, params);
    } 
    else if (method == "Trust Region") {
      ptr = new TrustRegion(grp, tests, params);
    } 
    else {
      cout << "ERROR: NOX::Solver::Manager - invalid choice for nonlinear " 
	   << "solver!" << endl;
      throw "NOX Error";
    }

    return true;
  }
  else 
    return ptr->reset(grp, tests, params);
}

Status::StatusType Manager::getStatus()
{
  return ptr->getStatus();
}

Status::StatusType Manager::iterate()
{
  return ptr->iterate();
}

Status::StatusType Manager::solve()
{
  return ptr->solve();
}

const Abstract::Group& Manager::getSolutionGroup() const
{
  return ptr->getSolutionGroup();
}

const Abstract::Group& Manager::getPreviousSolutionGroup() const
{
  return ptr->getPreviousSolutionGroup();
}

int Manager::getNumIterations() const
{
  return ptr->getNumIterations();
}

const Parameter::List & Manager::getOutputParameters() const
{
  return ptr->getOutputParameters();
}




