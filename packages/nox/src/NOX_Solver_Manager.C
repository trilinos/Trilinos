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


using namespace NOX;
using namespace NOX::Solver;

Manager::Manager(Abstract::Group& initialguess, Status::Test &t,
				     Parameter::List& p) :
  ptr(NULL)
{
  
  Utils::setUtils(p);

  string method = p.getParameter("Nonlinear Solver", "Newton");

  if (Utils::doPrint(Utils::Parameters)) 
    cout << "Nonlinear Solver: " << method << endl; 
  
  if (method == "Newton") {
    ptr = new Newton(initialguess, t, p);
  } 
  else {
    cout << "ERROR: NOX::Solver::Manager - invalid choice for nonlinear " 
	 << "solver!" << endl;
    throw 1;
  }
}

Manager::~Manager()
{
  delete ptr;
}

void Manager::resetInputParameters(Parameter::List& p)
{
  /* NOTE FROM TAMMY: May want to add a hook at some point to be able
     to switch nonlinear solver methods mid-stream. */
  ptr->resetInputParameters(p);
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




