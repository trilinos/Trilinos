// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
      TEST_FOR_EXCEPTION(true, std::logic_error, msg);
    }
  }
  else {
    std::string msg = "Error - NOX::Direction::Facotry::buildDirection() - Invalid choice for \"Method\" in \"Direction\" sublist!";
    TEST_FOR_EXCEPTION(true, std::logic_error, msg);
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
