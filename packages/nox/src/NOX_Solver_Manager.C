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
#include "NOX_Solver_LineSearchBased.H"	 // LineSearch method
#include "NOX_Solver_TrustRegionBased.H" // Trust region method

using namespace NOX;
using namespace NOX::Solver;

Manager::Manager(Abstract::Group& grp, StatusTest::Generic &t, const Parameter::List& p) :
  method(""),
  ptr(NULL)
{
  reset(grp, t, p);
}

Manager::~Manager()
{
  delete ptr;
}

bool Manager::reset(Abstract::Group& grp, StatusTest::Generic& tests, const Parameter::List& params)
{
  string newmethod = params.getParameter("Nonlinear Solver", "Newton");

  if (method != newmethod) {
    
    method = newmethod;

    delete ptr;
    ptr = NULL;
    
    if (method == "Newton") {
      ptr = new LineSearchBased(grp, tests, params);
    } 
    else if (method == "Line Search") {
      ptr = new LineSearchBased(grp, tests, params);
    } 
    else if (method == "Trust Region") {
      ptr = new TrustRegionBased(grp, tests, params);
    } 
    else {
      cout << "ERROR: NOX::Solver::Manager - Invalid solver choice" << endl;
      throw "NOX Error";
    }

    if (ptr == NULL) {
      cerr << "NOX::Solver::Manager::reset - Null pointer error" << endl;
      return false;
    }

    return true;
  }
  else {

    if (ptr == NULL) {
      cerr << "NOX::Solver::Manager::reset - Null pointer error" << endl;
      return false;
    }

    return ptr->reset(grp, tests, params);
  }
}

NOX::StatusTest::StatusType Manager::getStatus()
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::getStatus - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->getStatus();
}

NOX::StatusTest::StatusType Manager::iterate()
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::iterate - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->iterate();
}

NOX::StatusTest::StatusType Manager::solve()
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::solve - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->solve();
}

const Abstract::Group& Manager::getSolutionGroup() const
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::getSolutionGroup - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->getSolutionGroup();
}

const Abstract::Group& Manager::getPreviousSolutionGroup() const
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::getPreviousSolutionGroup - Null pointer error" << endl;
    throw "NOX Error";
  }

  return ptr->getPreviousSolutionGroup();
}

int Manager::getNumIterations() const
{
  if (ptr == NULL)
    return 0;

  return ptr->getNumIterations();
}

const Parameter::List& Manager::getOutputParameters() const
{
  if (ptr == NULL) {
    cout << "NOX::Solver::Manager::getOutputParameters - Null pointer error" << endl;
    throw "NOX Error";
  }
    
  return ptr->getOutputParameters();
}




