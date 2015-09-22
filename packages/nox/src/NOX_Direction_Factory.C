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

#include "NOX_Direction_Factory.H" // class definition

#include "Teuchos_ParameterList.hpp"
#include "NOX_GlobalData.H"
#include "NOX_Common.H"

// All the different directions
#include "NOX_Direction_Newton.H"
#include "NOX_Direction_SteepestDescent.H"
#include "NOX_Direction_NonlinearCG.H"
#include "NOX_Direction_Broyden.H"
#include "NOX_Direction_UserDefinedFactory.H"
#ifdef WITH_PRERELEASE
#include "NOX_Direction_Tensor.H"
#include "NOX_Direction_ModifiedNewton.H"
#include "NOX_Direction_QuasiNewton.H"
#endif


// ************************************************************************
// ************************************************************************
NOX::Direction::Factory::Factory()
{ }

// ************************************************************************
// ************************************************************************
NOX::Direction::Factory::~Factory()
{ }

// ************************************************************************
// ************************************************************************
Teuchos::RCP<NOX::Direction::Generic> NOX::Direction::Factory::
buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
           Teuchos::ParameterList& params)
{

  Teuchos::RCP<NOX::Direction::Generic> direction;

  std::string method = params.get("Method", "Newton");

  if (method == "Newton")
    direction = Teuchos::rcp(new Newton(gd, params));
  else if (method == "Steepest Descent")
    direction = Teuchos::rcp(new SteepestDescent(gd, params));
  else if (method == "NonlinearCG")
    direction = Teuchos::rcp(new NonlinearCG(gd, params));
  else if (method == "Broyden")
    direction = Teuchos::rcp(new Broyden(gd, params));
#ifdef WITH_PRERELEASE
  else if (method == "Tensor")
    direction = Teuchos::rcp(new Tensor(gd, params));
  else if (method == "Modified-Newton")
    direction = Teuchos::rcp(new ModifiedNewton(gd, params));
  else if (method == "Quasi-Newton")
    direction = Teuchos::rcp(new QuasiNewton(gd, params));
#endif
  else if (method == "User Defined") {
    using namespace Teuchos;
    if (isParameterType< RCP<NOX::Direction::UserDefinedFactory> >
    (params, "User Defined Direction Factory")) {

      RCP<NOX::Direction::UserDefinedFactory> user_factory =
    getParameter< Teuchos::RCP<NOX::Direction::UserDefinedFactory> >
    (params, "User Defined Direction Factory");

      direction = user_factory->buildDirection(gd, params);
    }
    else {
      std::string msg = "Error - NOX::Direction::Factory::buildDirection() -  a \"User Defined\" direction was chosen for the \"Method\" in the \"Direction\" sublist, but a Teuchos::RCP<NOX::Direction::UserDefinedFactory> object was not found in the parameter list!";
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
  }
  else {
    std::string msg = "Error - NOX::Direction::Facotry::buildDirection() - Invalid choice for \"Method\" in \"Direction\" sublist!";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
  }

  return direction;
}

// ************************************************************************
// ************************************************************************
// nonmember function
Teuchos::RCP<NOX::Direction::Generic> NOX::Direction::
buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
        Teuchos::ParameterList& params)
{
  NOX::Direction::Factory factory;
  return factory.buildDirection(gd, params);
}

// ************************************************************************
// ************************************************************************
