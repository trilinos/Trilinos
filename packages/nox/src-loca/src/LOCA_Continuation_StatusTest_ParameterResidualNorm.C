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
  stream << endl;

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

