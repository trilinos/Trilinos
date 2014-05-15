// $Id$
// $Source$

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
