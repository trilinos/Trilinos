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

NOX::Solver::Manager::Manager(Abstract::Group& grp, StatusTest::Generic &t, Parameter::List& p) :
  method(""),
  ptr(NULL)
{
  reset(grp, t, p);
}

NOX::Solver::Manager::Manager() :
  method(""),
  ptr(NULL)
{
}

NOX::Solver::Manager::~Manager()
{
  delete ptr;
}

bool NOX::Solver::Manager::reset(Abstract::Group& grp, 
				 StatusTest::Generic& tests, 
				 Parameter::List& params)
{
  string newmethod = params.getParameter("Nonlinear Solver", "Line Search Based");

  if ((method == newmethod) && (ptr != NULL))
  {
    return ptr->reset(grp, tests, params);
  }
  else 
  {
    method = newmethod;

    delete ptr;
    ptr = NULL;
    
    if ((method == "Newton") || (method == "Line Search")) // deprecated
    {	
      deprecated(method, "Line Search Based");
      ptr = new LineSearchBased(grp, tests, params);
    } 
    else if (method == "Line Search Based") 
    {
      ptr = new LineSearchBased(grp, tests, params);
    } 
    else if (method == "Trust Region")  // deprecated
    {
      deprecated(method, "Trust Region Based");
      ptr = new TrustRegionBased(grp, tests, params);
    } 
    else if (method == "Trust Region Based") 
    {
      ptr = new TrustRegionBased(grp, tests, params);
    } 
    else 
    {
      cout << "ERROR: NOX::Solver::Manager::reset - Invalid solver choice" << endl;
      throw "NOX Error";
    }

    if (ptr == NULL) 
    {
      cerr << "NOX::Solver::Manager::reset - Null pointer error" << endl;
      return false;
    }

    return true;
  }
}

// PRIVATE
void NOX::Solver::Manager::deprecated(const string& oldName, const string& newName) const
{
  cout << "Warning: NOX::Solver::Manager::reset - " 
       << "Nonlinear Solver choice \"" << oldName << "\" is deprecated.\n"
       << "                                       " 
       << "Use \"" << newName << "\" instead." 
       << endl;
}

NOX::StatusTest::StatusType NOX::Solver::Manager::getStatus()
{
  checkNullPtr("getStatus");
  return ptr->getStatus();
}

NOX::StatusTest::StatusType NOX::Solver::Manager::iterate()
{
  checkNullPtr("iterate");
  return ptr->iterate();
}

NOX::StatusTest::StatusType NOX::Solver::Manager::solve()
{
  checkNullPtr("solve");
  return ptr->solve();
}

const NOX::Abstract::Group& NOX::Solver::Manager::getSolutionGroup() const
{
  checkNullPtr("getSolutionGroup");
  return ptr->getSolutionGroup();
}

const NOX::Abstract::Group& NOX::Solver::Manager::getPreviousSolutionGroup() const
{
  checkNullPtr("getPreviousSolutionGroup");
  return ptr->getPreviousSolutionGroup();
}

int NOX::Solver::Manager::getNumIterations() const
{
  if (ptr == NULL)
    return 0;

  return ptr->getNumIterations();
}

const NOX::Parameter::List& NOX::Solver::Manager::getParameterList() const
{
  checkNullPtr("getParameterList");
  return ptr->getParameterList();
}

// PRIVATE
void NOX::Solver::Manager::checkNullPtr(const string& fname) const
{
  if (ptr == NULL) 
  {
    cout << "NOX::Solver::Manager::" << fname << " - Null pointer error" << endl;
    throw "NOX Error";
  }
}

