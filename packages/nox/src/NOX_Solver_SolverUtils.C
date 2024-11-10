// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "NOX_Observer.hpp"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_StatusTest_Generic.H"

#include "NOX_Solver_SolverUtils.H"

// ************************************************************************
// ************************************************************************
NOX::StatusTest::CheckType 
NOX::Solver::parseStatusTestCheckType(Teuchos::ParameterList& p)
{
  return Teuchos::getIntegralValue<NOX::StatusTest::CheckType>(p,"Status Test Check Type");
}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::Observer> 
NOX::Solver::parseObserver(Teuchos::ParameterList& solver_options_list)
{
  Teuchos::RCP<NOX::Observer> o = solver_options_list.get<Teuchos::RCP<NOX::Observer>>("Observer");
  Teuchos::RCP<NOX::Observer> ppo = solver_options_list.get<Teuchos::RCP<NOX::Observer>>("User Defined Pre/Post Operator");

  TEUCHOS_TEST_FOR_EXCEPTION(nonnull(o) && nonnull(ppo),std::runtime_error,
                             "ERROR: NOX::Sovler::parseObserver() - User has registered an \"Observer\" "
                             << "and a \"USer Defined Pre/Post Operator\". Pick one or the other!")

  if (nonnull(o))
    return o;
  if (nonnull(ppo))
    return ppo;

  return Teuchos::rcp(new NOX::Observer);
}

// ************************************************************************
// ************************************************************************
void NOX::Solver::validateSolverOptionsSublist(Teuchos::ParameterList& p)
{
  Teuchos::ParameterList validParams("Valid Params");

  Teuchos::setStringToIntegralParameter<NOX::StatusTest::CheckType>
    ("Status Test Check Type",
     "Minimal",
     "Sets the StatusTest check type.",
     Teuchos::tuple<std::string>("Complete","Minimal","None"),
     Teuchos::tuple<NOX::StatusTest::CheckType>(NOX::StatusTest::Complete,NOX::StatusTest::Minimal,NOX::StatusTest::None),
     &validParams);

  Teuchos::RCP<NOX::Observer> observer;
  validParams.set("Observer",observer);
  // Deprecated old flag
  validParams.set("User Defined Pre/Post Operator",observer);

  Teuchos::RCP<NOX::MeritFunction::Generic> mf;
  validParams.set("User Defined Merit Function",mf);

  Teuchos::setStringToIntegralParameter<int>
    ("Fixed Point Iteration Type",
     "Seidel",
     "Sets iteration type for the fixed point solver.",
     Teuchos::tuple<std::string>("Seidel","Jacobi"),
     &validParams);

  p.validateParametersAndSetDefaults(validParams);  
}

// ************************************************************************
// ************************************************************************

