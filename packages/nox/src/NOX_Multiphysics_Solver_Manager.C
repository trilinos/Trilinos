//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_Multiphysics_Solver_Manager.H"	// class definition
#include "NOX_Multiphysics_DataExchange_Interface.H"	// class definition
#include "NOX_Utils.H"		// for static function doPrint
#include "Teuchos_RefCountPtr.hpp" // for RefCountPtr 

// Header files for different solvers
#include "NOX_Multiphysics_Solver_FixedPointBased.H"	 // LineSearch method
#ifdef WITH_PRERELEASE
//#include "NOX_Solver_TensorBasedTest.H"  // Tensor-Krylov method
#endif

NOX::Multiphysics::Solver::Manager::Manager(
        const Teuchos::RefCountPtr<vector<NOX::Solver::Manager*> >& solvers, 
        const Teuchos::RefCountPtr<NOX::Multiphysics::DataExchange::Interface>& i, 
	const Teuchos::RefCountPtr<StatusTest::Generic>& t, 
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& p) :
  utils(p->sublist("Printing")),
  method(""),
  cplPtr(NULL)
{
  reset(solvers, i, t, p);
}

NOX::Multiphysics::Solver::Manager::Manager(
        const Teuchos::RefCountPtr<Abstract::Group>& grp, 
	const Teuchos::RefCountPtr<StatusTest::Generic>& t, 
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& p) :
  utils(p->sublist("Printing")),
  method(""),
  cplPtr(NULL)
{
  reset(grp, t, p);
}

NOX::Multiphysics::Solver::Manager::Manager() :
  method(""),
  cplPtr(NULL)
{
}

NOX::Multiphysics::Solver::Manager::~Manager()
{
  delete cplPtr;
}

bool NOX::Multiphysics::Solver::Manager::
reset(const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
      const Teuchos::RefCountPtr<StatusTest::Generic>& tests, 
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& params)
{
  return false;
}

bool NOX::Multiphysics::Solver::Manager::reset(
      const Teuchos::RefCountPtr<vector<NOX::Solver::Manager*> >& solvers, 
      const Teuchos::RefCountPtr<NOX::Multiphysics::DataExchange::Interface>& interface, 
      const Teuchos::RefCountPtr<StatusTest::Generic>& tests, 
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& params)
{
  string newmethod = 
    params->get("Coupling Strategy", "Fixed Point Based");

  if ((method == newmethod) && (cplPtr != NULL))
  {
    return cplPtr->reset(solvers, interface, tests, params);
  }
  else 
  {
    method = newmethod;

    delete cplPtr;
    cplPtr = NULL;
    
    if( method == "Fixed Point Based" ) 
    {	
      cplPtr = new NOX::Multiphysics::Solver::FixedPointBased(solvers, interface, tests, params);
    } 
#ifdef WITH_PRERELEASE
//    else if (method == "Tensor-Krylov Based") 
//    {
//      cplPtr = new NOX::Solver::TensorBasedTest(grp, tests, params);
//    } 
#endif
    else 
    {
      utils.out() << "ERROR: NOX::Multiphysics::Solver::Manager::reset - Invalid solver choice " << method << endl;
      throw "NOX Error";
    }

    if (cplPtr == NULL) 
    {
      utils.err() << "NOX::Multiphysics::Solver::Manager::reset - Null pointer error" << endl;
      return false;
    }

    return true;
  }
}

bool 
NOX::Multiphysics::Solver::Manager::reset(
      const Teuchos::RefCountPtr<Abstract::Group>& grp, 
      const Teuchos::RefCountPtr<StatusTest::Generic>& tests)
{
  return cplPtr->reset(grp, tests);
}

// PRIVATE
void NOX::Multiphysics::Solver::Manager::deprecated(const string& oldName, const string& newName) const
{
  utils.out() << "Warning: NOX::Multiphysics::Solver::Manager::reset - " 
       << "Nonlinear Solver choice \"" << oldName << "\" is deprecated.\n"
       << "                                       " 
       << "Use \"" << newName << "\" instead." 
       << endl;
}

NOX::StatusTest::StatusType NOX::Multiphysics::Solver::Manager::getStatus()
{
  checkNullPtr("getStatus");
  return cplPtr->getStatus();
}

NOX::StatusTest::StatusType NOX::Multiphysics::Solver::Manager::step()
{
  checkNullPtr("step");
  return cplPtr->step();
}

NOX::StatusTest::StatusType NOX::Multiphysics::Solver::Manager::solve()
{
  checkNullPtr("solve");
  return cplPtr->solve();
}

const NOX::Abstract::Group& NOX::Multiphysics::Solver::Manager::getSolutionGroup() const
{
  checkNullPtr("getSolutionGroup");
  return cplPtr->getSolutionGroup();
}

const NOX::Abstract::Group& NOX::Multiphysics::Solver::Manager::getPreviousSolutionGroup() const
{
  checkNullPtr("getPreviousSolutionGroup");
  return cplPtr->getPreviousSolutionGroup();
}

int NOX::Multiphysics::Solver::Manager::getNumIterations() const
{
  if (cplPtr == NULL)
    return 0;

  return cplPtr->getNumIterations();
}

const Teuchos::ParameterList& NOX::Multiphysics::Solver::Manager::getList() const
{
  checkNullPtr("getList");
  return cplPtr->getList();
}

// PRIVATE
void NOX::Multiphysics::Solver::Manager::checkNullPtr(const string& fname) const
{
  if (cplPtr == NULL) 
  {
    utils.out() << "NOX::Multiphysics::Solver::Manager::" << fname << " - Null pointer error" << endl;
    throw "NOX Error";
  }
}

