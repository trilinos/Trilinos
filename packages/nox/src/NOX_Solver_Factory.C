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

#include "NOX_Common.H"
#include "NOX_Solver_Factory.H"	// class definition
#include "NOX_Solver_Generic.H"

// Header files for different solvers
#include "NOX_Solver_LineSearchBased.H"	 // LineSearch method
#include "NOX_Solver_TrustRegionBased.H" // Trust region method
#include "NOX_Solver_InexactTrustRegionBased.H" // Inexact Trust region method
#include "NOX_Solver_TensorBased.H"      // Tensor method
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

  string method = params->get("Nonlinear Solver", "Line Search Based");
  
  if ((method == "Newton") || (method == "Line Search Based")) 
    solver = rcp(new LineSearchBased(grp, tests, params));
  else if (method == "Trust Region Based")
    solver = rcp(new TrustRegionBased(grp, tests, params));
  else if (method == "Inexact Trust Region Based") 
    solver = rcp(new InexactTrustRegionBased(grp, tests, params));
  else if (method == "Tensor Based") 
    solver = rcp(new TensorBased(grp, tests, params));
#ifdef WITH_PRERELEASE
  else if (method == "Tensor-Krylov Based")
    solver = rcp(new TensorBasedTest(grp, tests, params));
#endif
  else {
    std::ostringstream msg;
    msg << "Error - NOX::Solver::Manager::buildSolver() - The \"Nonlinear Solver\" parameter \"" << method << "\" is not a valid solver option.  Please fix your parameter list!";
    TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
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

