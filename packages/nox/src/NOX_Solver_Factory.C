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
#include "NOX_Solver_Factory.H"	// class definition
#include "NOX_Solver_Generic.H"

// Header files for different solvers
#include "NOX_Solver_LineSearchBased.H"	 // LineSearch method
#include "NOX_Solver_TrustRegionBased.H" // Trust region method
#include "NOX_Solver_InexactTrustRegionBased.H" // Inexact Trust region method
#include "NOX_Solver_TensorBased.H"      // Tensor method
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

