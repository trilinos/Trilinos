// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "NOX_LineSearch_Factory.H" // class definition

#include "Teuchos_ParameterList.hpp"
#include "NOX_GlobalData.H"

// All the different line searches
#include "NOX_LineSearch_FullStep.H"
#include "NOX_LineSearch_Backtrack.H"
#include "NOX_LineSearch_Polynomial.H"
#include "NOX_LineSearch_MoreThuente.H"
#include "NOX_LineSearch_NonlinearCG.H"
#include "NOX_LineSearch_SafeguardedStep.H"
#include "NOX_LineSearch_UserDefinedFactory.H"
#ifdef HAVE_NOX_THYRA
#include "NOX_LineSearch_SafeguardedDirection.hpp"
#endif

// ************************************************************************
// ************************************************************************
NOX::LineSearch::Factory::Factory()
{

}

// ************************************************************************
// ************************************************************************
NOX::LineSearch::Factory::~Factory()
{

}

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::LineSearch::Generic> NOX::LineSearch::Factory::
buildLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
        Teuchos::ParameterList& params)
{

  Teuchos::RCP<NOX::LineSearch::Generic> line_search;

  std::string method = params.get("Method", "Full Step");

  if (method == "Full Step")
    line_search = Teuchos::rcp(new FullStep(gd, params));
  else if (method == "Backtrack")
    line_search = Teuchos::rcp(new Backtrack(gd, params));
  else if (method == "Polynomial")
    line_search = Teuchos::rcp(new Polynomial(gd, params));
  else if (method == "More'-Thuente")
    line_search = Teuchos::rcp(new MoreThuente(gd, params));
  else if (method == "NonlinearCG")
    line_search = Teuchos::rcp(new NonlinearCG(gd, params));
  else if (method == "Safeguarded Step")
    line_search = Teuchos::rcp(new SafeguardedStep(gd, params));
#ifdef HAVE_NOX_THYRA
  else if (method == "Safeguarded Direction")
    line_search = Teuchos::rcp(new SafeguardedDirection(gd, params));
#endif
  else if (method == "User Defined") {
    using namespace Teuchos;
    if (isParameterType< RCP<NOX::LineSearch::UserDefinedFactory> >
    (params, "User Defined Line Search Factory")) {

      RCP<NOX::LineSearch::UserDefinedFactory> user_factory =
    getParameter< Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> >
    (params, "User Defined Line Search Factory");

      line_search = user_factory->buildLineSearch(gd, params);
    }
    else {
      std::string msg = "Error - NOX::LineSearch::Factory::buildLineSearch() -  a \"User Defined\" line search was chosen for the \"Method\" in the \"Line Search\" sublist, but a Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> object was not found in the parameter list!";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
  }
  else {
    std::string msg = "Error - NOX::LineSearch::Facotry::buildLineSearch() - Invalid choice for \"Method\" in \"Line Search\" sublist!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  return line_search;
}

// ************************************************************************
// ************************************************************************
// nonmember function
Teuchos::RCP<NOX::LineSearch::Generic> NOX::LineSearch::
buildLineSearch(const Teuchos::RCP<NOX::GlobalData>& gd,
        Teuchos::ParameterList& params)
{
  NOX::LineSearch::Factory factory;
  return factory.buildLineSearch(gd, params);
}

// ************************************************************************
// ************************************************************************
