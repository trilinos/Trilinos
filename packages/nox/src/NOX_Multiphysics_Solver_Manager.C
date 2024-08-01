// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Multiphysics_Solver_Manager.H"    // class definition
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
        const Teuchos::RCP<NOX::Abstract::Group>& /* grp */,
    const Teuchos::RCP<NOX::StatusTest::Generic>& /* t */,
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
      throw std::runtime_error("NOX Error");
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

void NOX::Multiphysics::Solver::Manager::reset() {}

// PRIVATE
void NOX::Multiphysics::Solver::Manager::deprecated(const std::string& oldName, const std::string& newName) const
{
  utils.out() << "Warning: NOX::Multiphysics::Solver::Manager::reset - "
       << "Nonlinear Solver choice \"" << oldName << "\" is deprecated.\n"
       << "                                       "
       << "Use \"" << newName << "\" instead."
       << std::endl;
}

NOX::StatusTest::StatusType NOX::Multiphysics::Solver::Manager::getStatus() const
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
    throw std::runtime_error("NOX Error");
  }
}

