//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Multiphysics_Solver_Manager.H"	// class definition
#include "NOX_Multiphysics_DataExchange_Interface.H"
#include "NOX_Multiphysics_Solver_FixedPointBased.H"
#include "NOX_Utils.H"
#include "Teuchos_RCP.hpp"


NOX::Multiphysics::Solver::Manager::Manager(
        const Teuchos::RCP<std::vector<Teuchos::RCP<NOX::Solver::Generic> > >& solvers, 
        const Teuchos::RCP<NOX::Multiphysics::DataExchange::Interface>& i, 
	const Teuchos::RCP<NOX::StatusTest::Generic>& t, 
	const Teuchos::RCP<Teuchos::ParameterList>& p) :
  utils(p->sublist("Printing")),
  method(""),
  cplPtr(NULL)
{
  reset(solvers, i, t, p);
}

NOX::Multiphysics::Solver::Manager::Manager(
        const Teuchos::RCP<NOX::Abstract::Group>& grp, 
	const Teuchos::RCP<NOX::StatusTest::Generic>& t, 
	const Teuchos::RCP<Teuchos::ParameterList>& p) :
  utils(p->sublist("Printing")),
  method(""),
  cplPtr(NULL)
{
  
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

bool NOX::Multiphysics::Solver::Manager::reset(
      const Teuchos::RCP<std::vector<Teuchos::RCP<NOX::Solver::Generic> > >& solvers, 
      const Teuchos::RCP<NOX::Multiphysics::DataExchange::Interface>& interface, 
      const Teuchos::RCP<NOX::StatusTest::Generic>& tests, 
      const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  std::string newmethod = 
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
    else 
    {
      utils.out() << "ERROR: NOX::Multiphysics::Solver::Manager::reset - Invalid solver choice " << method << std::endl;
      throw "NOX Error";
    }

    if (cplPtr == NULL) 
    {
      utils.err() << "NOX::Multiphysics::Solver::Manager::reset - Null pointer error" << std::endl;
      return false;
    }

    return true;
  }
}

void
NOX::Multiphysics::Solver::Manager::reset(
      const NOX::Abstract::Vector& initialGuess, 
      const Teuchos::RCP<NOX::StatusTest::Generic>& tests)
{
  cplPtr->reset(initialGuess, tests);
}

void 
NOX::Multiphysics::Solver::Manager::
reset(const Abstract::Vector& initialGuess)
{
  cplPtr->reset(initialGuess);
}

// PRIVATE
void NOX::Multiphysics::Solver::Manager::deprecated(const std::string& oldName, const std::string& newName) const
{
  utils.out() << "Warning: NOX::Multiphysics::Solver::Manager::reset - " 
       << "Nonlinear Solver choice \"" << oldName << "\" is deprecated.\n"
       << "                                       " 
       << "Use \"" << newName << "\" instead." 
       << std::endl;
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

Teuchos::RCP< const NOX::Abstract::Group > NOX::Multiphysics::Solver::Manager::getSolutionGroupPtr() const
{
  checkNullPtr("getSolutionGroupPtr");
  return cplPtr->getSolutionGroupPtr();
}

const NOX::Abstract::Group& NOX::Multiphysics::Solver::Manager::getPreviousSolutionGroup() const
{
  checkNullPtr("getPreviousSolutionGroup");
  return cplPtr->getPreviousSolutionGroup();
}

Teuchos::RCP< const NOX::Abstract::Group > NOX::Multiphysics::Solver::Manager::getPreviousSolutionGroupPtr() const
{
  checkNullPtr("getPreviousSolutionGroupPtr");
  return cplPtr->getPreviousSolutionGroupPtr();
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

Teuchos::RCP< const Teuchos::ParameterList > NOX::Multiphysics::Solver::Manager::getListPtr() const
{
  checkNullPtr("getListPtr");
  return cplPtr->getListPtr();
}

// PRIVATE
void NOX::Multiphysics::Solver::Manager::checkNullPtr(const std::string& fname) const
{
  if (cplPtr == NULL) 
  {
    utils.out() << "NOX::Multiphysics::Solver::Manager::" << fname << " - Null pointer error" << std::endl;
    throw "NOX Error";
  }
}

