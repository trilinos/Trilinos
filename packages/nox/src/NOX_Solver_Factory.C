// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "NOX_Solver_Factory.H"    // class definition
#include "NOX_Solver_Generic.H"

// Header files for different solvers
#include "NOX_Solver_LineSearchBased.H"     // LineSearch method
#include "NOX_Solver_TrustRegionBased.H" // Trust region method
#include "NOX_Solver_InexactTrustRegionBased.H" // Inexact Trust region method
#include "NOX_Solver_TensorBased.H"      // Tensor method
#include "NOX_Solver_AndersonAcceleration.H"
#include "NOX_Solver_SingleStep.H"
#ifdef HAVE_NOX_THYRA
#include "NOX_Solver_PseudoTransient.hpp"      // Pseudo-Transient
#endif
#ifdef WITH_PRERELEASE
#include "NOX_Solver_TensorBasedTest.H"  // Tensor-Krylov method
#endif

// ************************************************************************
// ************************************************************************
NOX::Solver::Factory::Factory()
{ }

// ************************************************************************
// ************************************************************************
NOX::Solver::Factory::~Factory()
{ }

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::Solver::Generic> NOX::Solver::Factory::
buildSolver(const Teuchos::RCP<NOX::Abstract::Group>& grp,
        const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
        const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  RCP<NOX::Solver::Generic> solver;

  std::string method = params->get("Nonlinear Solver", "Line Search Based");

  if ((method == "Newton") || (method == "Line Search Based"))
    solver = rcp(new LineSearchBased(grp, tests, params));
  else if (method == "Trust Region Based")
    solver = rcp(new TrustRegionBased(grp, tests, params));
  else if (method == "Inexact Trust Region Based")
    solver = rcp(new InexactTrustRegionBased(grp, tests, params));
  else if (method == "Tensor Based")
    solver = rcp(new TensorBased(grp, tests, params));
  else if (method == "Anderson Accelerated Fixed-Point")
    solver = rcp(new AndersonAcceleration(grp, tests, params));
  else if (method == "Single Step")
    solver = rcp(new SingleStep(grp, params));
#ifdef HAVE_NOX_THYRA
  else if (method == "Pseudo-Transient")
    solver = rcp(new PseudoTransient(grp, tests, params));
#endif
#ifdef WITH_PRERELEASE
  else if (method == "Tensor-Krylov Based")
    solver = rcp(new TensorBasedTest(grp, tests, params));
#endif
  else {
    std::ostringstream msg;
    msg << "Error - NOX::Solver::Manager::buildSolver() - The \"Nonlinear Solver\" parameter \"" << method << "\" is not a valid solver option.  Please fix your parameter list!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
  }

  return solver;
}

// ************************************************************************
// ************************************************************************
// Nonmember function
Teuchos::RCP<NOX::Solver::Generic>
NOX::Solver::buildSolver(const Teuchos::RCP<NOX::Abstract::Group>& grp,
             const Teuchos::RCP<NOX::StatusTest::Generic>& tests,
             const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  NOX::Solver::Factory factory;
  return factory.buildSolver(grp, tests, params);
}

// ************************************************************************
// ************************************************************************

