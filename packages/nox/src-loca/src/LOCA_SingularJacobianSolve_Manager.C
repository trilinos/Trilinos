// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#include "LOCA_Continuation_AbstractGroup.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_SingularJacobianSolve_Manager.H"
#include "LOCA_SingularJacobianSolve_Default.H"
#include "LOCA_SingularJacobianSolve_Nic.H"
#include "LOCA_SingularJacobianSolve_NicDay.H"
#include "LOCA_SingularJacobianSolve_ItRef.H"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::SingularJacobianSolve::Manager::Manager(Teuchos::ParameterList& params) :
  method(),
  singularSolverPtr(NULL)
{
  reset(params);
}

LOCA::SingularJacobianSolve::Manager::Manager(const Teuchos::ParameterList& params) :
  method(),
  singularSolverPtr(NULL)
{
  Teuchos::ParameterList p(params);
  reset(p);
}

LOCA::SingularJacobianSolve::Manager::Manager(
			  const LOCA::SingularJacobianSolve::Manager& source) 
  : method(source.method),
    singularSolverPtr(source.singularSolverPtr->clone())
{
}

LOCA::SingularJacobianSolve::Manager::~Manager()
{
  delete singularSolverPtr;
}

LOCA::SingularJacobianSolve::Generic*
LOCA::SingularJacobianSolve::Manager::clone() const 
{
  return new Manager(*this);
}

LOCA::SingularJacobianSolve::Generic&
LOCA::SingularJacobianSolve::Manager::operator=(
			  const LOCA::SingularJacobianSolve::Generic& source)
{
  return operator=(dynamic_cast<const LOCA::SingularJacobianSolve::Manager&>(source));
}

LOCA::SingularJacobianSolve::Manager&
LOCA::SingularJacobianSolve::Manager::operator=(
			  const LOCA::SingularJacobianSolve::Manager& source)
{
  // protect against A = A
  if (this != &source) {
    delete singularSolverPtr;

    method = source.method;
    singularSolverPtr = source.singularSolverPtr->clone();
  }

  return *this;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Manager::reset(Teuchos::ParameterList& params) 
{
  std::string newmethod = params.get("Method", "Default");

  if (method != newmethod) {
    delete singularSolverPtr;

    method = newmethod;

    if (method == "Default")
      singularSolverPtr = new LOCA::SingularJacobianSolve::Default(params);
    else if (method == "Nic")
      singularSolverPtr = new LOCA::SingularJacobianSolve::Nic(params);
    else if (method == "Nic-Day")
      singularSolverPtr = new LOCA::SingularJacobianSolve::NicDay(params);
    else if (method == "Iterative Refinement")
      singularSolverPtr = new LOCA::SingularJacobianSolve::ItRef(params);
    else {
      LOCA::ErrorCheck::throwError(
			      "LOCA::SingularJacobianSolve::Manager::reset()",
			      "Invalid choice for singular solve method.");
      return NOX::Abstract::Group::Failed;
    }
  }

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Manager::compute(
				Teuchos::ParameterList& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector& input,
			        const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector& result) 
{
  if (singularSolverPtr == NULL) {
    LOCA::ErrorCheck::throwError(
			 "LOCA::SingularJacobianSolve::Manager::compute()", 
			 "Null pointer error");
    return NOX::Abstract::Group::Failed;
  }

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails))
    std::cout << "\n\tCalling singular solver with method: " << method << std::endl;

  return singularSolverPtr->compute(params, grp, input, approxNullVec,
				    jacApproxNullVec, result);
}

NOX::Abstract::Group::ReturnType 
LOCA::SingularJacobianSolve::Manager::computeMulti(
				Teuchos::ParameterList& params,
				LOCA::Continuation::AbstractGroup& grp,
				const NOX::Abstract::Vector*const* inputs,
				const NOX::Abstract::Vector& approxNullVec,
				const NOX::Abstract::Vector& jacApproxNullVec,
				NOX::Abstract::Vector** results,
				int nVecs) 
{
  if (singularSolverPtr == NULL) {
    LOCA::ErrorCheck::throwError(
		       "LOCA::SingularJacobianSolve::Manager::computeMulti()",
		       "Null pointer error");
    return NOX::Abstract::Group::Failed;
  }

  if (LOCA::Utils::doPrint(LOCA::Utils::StepperDetails))
    std::cout << "\n\tCalling singular solver with method: " << method << std::endl;
  
  return singularSolverPtr->computeMulti(params, grp, inputs, approxNullVec,
					 jacApproxNullVec, results, nVecs);
}

const std::string&
LOCA::SingularJacobianSolve::Manager::getMethod() const 
{
  return method;
}
