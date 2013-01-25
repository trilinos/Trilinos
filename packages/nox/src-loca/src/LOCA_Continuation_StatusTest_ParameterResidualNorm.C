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

#include "LOCA_Continuation_StatusTest_ParameterResidualNorm.H" 
#include "LOCA_Continuation_ExtendedGroup.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

LOCA::Continuation::StatusTest::ParameterResidualNorm::ParameterResidualNorm(
							      double rtol_, 
							      double atol_,
							      double tol_) :
  rtol(rtol_),
  atol(atol_),
  tol(tol_),
  paramResidualNorm(0.0),
  status(NOX::StatusTest::Unconverged)
{
}


LOCA::Continuation::StatusTest::ParameterResidualNorm::~ParameterResidualNorm()
{
}

NOX::StatusTest::StatusType 
LOCA::Continuation::StatusTest::ParameterResidualNorm::checkStatus(
					 const NOX::Solver::Generic& problem)
{
  // Get solution groups from solver
  const NOX::Abstract::Group& soln = problem.getSolutionGroup();

  // Cast soln group to continuation group (for parameter step)
  const LOCA::Continuation::ExtendedGroup* conGroupPtr = 
    dynamic_cast<const LOCA::Continuation::ExtendedGroup*>(&soln);

  // Check that group is a continuation group, return converged if not
  if (conGroupPtr == NULL) {
    paramResidualNorm = 0.0;
    return NOX::StatusTest::Converged;
  }

  // Get residual vector
  const LOCA::Continuation::ExtendedVector& f = 
    dynamic_cast<const LOCA::Continuation::ExtendedVector&>(soln.getF());
  
  paramResidualNorm = 
    fabs(f.getParam()) / (rtol*fabs(conGroupPtr->getStepSize()) + atol);

  if (paramResidualNorm < tol) 
    status = NOX::StatusTest::Converged;
  else
    status = NOX::StatusTest::Unconverged;
  
  return status;
}

NOX::StatusTest::StatusType 
LOCA::Continuation::StatusTest::ParameterResidualNorm::getStatus() const
{
  return status;
}


ostream& 
LOCA::Continuation::StatusTest::ParameterResidualNorm::print(ostream& stream, 
							   int indent) const
{
  for (int j = 0; j < indent; j++)
    stream << ' ';
  stream << status;
  stream << "Continuation Scaled Parameter Residual = " << NOX::Utils::sciformat(paramResidualNorm, 3) << " < " << tol;
  stream << std::endl;

  return stream;
} 

double 
LOCA::Continuation::StatusTest::ParameterResidualNorm::getResidualNorm() const
{
  return paramResidualNorm;
}

double 
LOCA::Continuation::StatusTest::ParameterResidualNorm::getRTOL() const
{
  return rtol;
}

double 
LOCA::Continuation::StatusTest::ParameterResidualNorm::getATOL() const
{
  return atol;
}

double 
LOCA::Continuation::StatusTest::ParameterResidualNorm::getTOL() const
{
  return tol;
}

