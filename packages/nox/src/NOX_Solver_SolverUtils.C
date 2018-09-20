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

