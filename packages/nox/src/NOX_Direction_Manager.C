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

#include "NOX_Direction_Manager.H" // class definition
#include "NOX_Utils.H"
#include "NOX_GlobalData.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Direction_Newton.H"
#include "NOX_Direction_SteepestDescent.H"
#include "NOX_Direction_NonlinearCG.H"
#include "NOX_Direction_Broyden.H"

#ifdef WITH_PRERELEASE
#include "NOX_Direction_Tensor.H"
#include "NOX_Direction_ModifiedNewton.H"
#include "NOX_Direction_QuasiNewton.H"
#endif

NOX::Direction::Manager::
Manager(const Teuchos::RefCountPtr<NOX::GlobalData>& gd) :
  utils(gd->getUtils()),
  method("")
{
  
}

NOX::Direction::Manager::
Manager(const Teuchos::RefCountPtr<NOX::GlobalData>& gd, 
	Teuchos::ParameterList& params) :
  method("")
{
  reset(gd, params);
}

NOX::Direction::Manager::~Manager()
{

}

bool NOX::Direction::Manager::
reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
      Teuchos::ParameterList& params)
{
  utils = gd->getUtils();

  string newmethod = params.get("Method", "Newton");

  // If the method has not changeed, just call reset on the method.
  if (method == newmethod) 
  {
    return ptr->reset(gd, params);
  }

  // Otherwise, construct a whole new direction
  method = newmethod;
  
  if (method == "Newton")
    ptr = Teuchos::rcp(new Newton(gd, params));
  else if (method == "Steepest Descent")
    ptr = Teuchos::rcp(new SteepestDescent(gd, params));
  else if (method == "NonlinearCG")
    ptr = Teuchos::rcp(new NonlinearCG(gd, params));
  else if (method == "Broyden")
    ptr = Teuchos::rcp(new Broyden(gd, params));
#ifdef WITH_PRERELEASE
  else if (method == "Tensor")
    ptr = Teuchos::rcp(new Tensor(gd, params));
  else if (method == "Modified-Newton")
    ptr = Teuchos::rcp(new ModifiedNewton(gd, params));
  else if (method == "Quasi-Newton")
    ptr = Teuchos::rcp(new QuasiNewton(gd, params));
#endif
  else if (method == "User Defined") {
    if (params.INVALID_TEMPLATE_QUALIFIER
	isType< Teuchos::RefCountPtr<NOX::Direction::Generic> >
	("User Defined Direction")) {
      
      ptr = params.INVALID_TEMPLATE_QUALIFIER
	get< Teuchos::RefCountPtr<NOX::Direction::Generic> >
	("User Defined Direction");
      ptr->reset(gd, params);
    }
    else {
      this->printWarning("reset", " a \"User Defined\" Direction was chosen for the \"Method\" in the \"Direction\" sublist, but a Teuchos::RefCountPtr<NOX::Direction::Generic> object was not found in the parameter list!");
      throw "NOX Error";
    }
  }
  else {
    this->printWarning("reset", "invalid choice (" + method + 
		       ") for direction method");
    return false;
  }

  return (!Teuchos::is_null(ptr));
}


bool NOX::Direction::Manager::
compute(Abstract::Vector& dir, Abstract::Group& grp, 
	const Solver::Generic& solver) 
{
  if (Teuchos::is_null(ptr)) 
  {
    if (utils->isPrintType(NOX::Utils::Warning)) 
      utils->out() << "Calling NOX::Direction::Manager::compute on uninitialized direction" << endl;
    return false;
  }

  return ptr->compute(dir, grp, solver);
}

bool NOX::Direction::Manager::
compute(Abstract::Vector& dir, Abstract::Group& grp, 
	const Solver::LineSearchBased& solver) 
{
  if (Teuchos::is_null(ptr)) 
  {
    if (utils->isPrintType(NOX::Utils::Warning)) 
      utils->out() << "Calling NOX::Direction::Manager::compute on uninitialized direction" << endl;
    return false;
  }

  return ptr->compute(dir, grp, solver);
}

void NOX::Direction::Manager::
printWarning(const string& name, const string& warning)
{
  if (utils->isPrintType(NOX::Utils::Warning)) 
    utils->out() << "Calling NOX::Direction::Manager::" << name << " - " << warning << endl;
}
