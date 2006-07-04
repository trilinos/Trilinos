// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_Bifurcation_TPBord_StatusTest_ParameterUpdateNorm.H" 
#include "LOCA_Bifurcation_TPBord_ExtendedGroup.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::ParameterUpdateNorm(double rtol_, double atol_, double tol_) :
  rtol(rtol_),
  atol(atol_),
  tol(tol_),
  paramUpdateNorm(0.0),
  status(NOX::StatusTest::Unconverged)
{
}


LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::~ParameterUpdateNorm()
{
}

NOX::StatusTest::StatusType 
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::checkStatus(
					 const NOX::Solver::Generic& problem)
{
  // Get solution groups from solver
  const NOX::Abstract::Group& soln = problem.getSolutionGroup();
  const NOX::Abstract::Group& oldsoln = problem.getPreviousSolutionGroup();

  // Cast soln group to turning point group
  const LOCA::Bifurcation::TPBord::ExtendedGroup* tpGroupPtr = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedGroup*>(&soln);

  // Check that group is a turning point group, return converged if not
  if (tpGroupPtr == NULL) {
    paramUpdateNorm = 0.0;
    return NOX::StatusTest::Converged;
  }

  // Get solution vectors
  const LOCA::Bifurcation::TPBord::ExtendedVector& x = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(soln.getX());
  const LOCA::Bifurcation::TPBord::ExtendedVector& xold = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(oldsoln.getX());
  
  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid 
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0) 
  {
    paramUpdateNorm = 1.0e+12;
    status = NOX::StatusTest::Unconverged;
    return status;
  } 

  paramUpdateNorm = 
    fabs(x.getBifParam() - xold.getBifParam()) / (rtol*fabs(x.getBifParam()) + atol);

  if (paramUpdateNorm < tol) 
    status = NOX::StatusTest::Converged;
  else
    status = NOX::StatusTest::Unconverged;
  
  return status;
}

NOX::StatusTest::StatusType 
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getStatus() const
{
  return status;
}


ostream& 
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::print(
							    ostream& stream, 
							    int indent) const
{
  for (int j = 0; j < indent; j++)
    stream << ' ';
  stream << status;
  stream << "Turning Point Scaled Parameter Update = " 
	 << NOX::Utils::sciformat(paramUpdateNorm, 3) << " < " << tol;
  stream << endl;

  return stream;
}


double 
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getUpdateNorm() const
{
  return paramUpdateNorm;
}   

double 
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getRTOL() const
{
  return rtol;
}

double 
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getATOL() const
{
  return atol;
}

double 
LOCA::Bifurcation::TPBord::StatusTest::ParameterUpdateNorm::getTOL() const
{
  return tol;
}

