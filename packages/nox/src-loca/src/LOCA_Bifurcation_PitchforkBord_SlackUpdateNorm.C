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

#include "LOCA_Bifurcation_PitchforkBord_SlackUpdateNorm.H" 
#include "LOCA_Bifurcation_PitchforkBord_ExtendedGroup.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::SlackUpdateNorm(double rtol_, double atol_, double tol_) :
  rtol(rtol_),
  atol(atol_),
  tol(tol_),
  slackUpdateNorm(0.0),
  status(NOX::StatusTest::Unconverged)
{
}


LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::~SlackUpdateNorm()
{
}

NOX::StatusTest::StatusType 
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::checkStatus(
					 const NOX::Solver::Generic& problem)
{
  // Get solution groups from solver
  const NOX::Abstract::Group& soln = problem.getSolutionGroup();
  const NOX::Abstract::Group& oldsoln = problem.getPreviousSolutionGroup();

  // Cast soln group to pitchfork group
  const LOCA::Bifurcation::PitchforkBord::ExtendedGroup* pfGroupPtr = 
    dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedGroup*>(&soln);

  // Check that group is a pitchfork group, return converged if not
  if (pfGroupPtr == NULL) {
    slackUpdateNorm = 0.0;
    return NOX::StatusTest::Converged;
  }

  // Get solution vectors
  const LOCA::Bifurcation::PitchforkBord::ExtendedVector& x = 
    dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(soln.getX());
  const LOCA::Bifurcation::PitchforkBord::ExtendedVector& xold = 
    dynamic_cast<const LOCA::Bifurcation::PitchforkBord::ExtendedVector&>(oldsoln.getX());
  
  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid 
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0) 
  {
    slackUpdateNorm = 1.0e+12;
    status = NOX::StatusTest::Unconverged;
    return status;
  } 

  slackUpdateNorm = 
    fabs(x.getSlackVar() - xold.getSlackVar()) / (rtol*fabs(x.getSlackVar()) + atol);

  if (slackUpdateNorm < tol) 
    status = NOX::StatusTest::Converged;
  else
    status = NOX::StatusTest::Unconverged;
  
  return status;
}

NOX::StatusTest::StatusType 
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getStatus() const
{
  return status;
}


ostream& 
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::print(
							    ostream& stream, 
							    int indent) const
{
  for (int j = 0; j < indent; j++)
    stream << ' ';
  stream << status;
  stream << "Pitchfork Scaled Slack Variable Update = " 
	 << NOX::Utils::sciformat(slackUpdateNorm, 3) << " < " << tol;
  stream << endl;

  return stream;
}


double 
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getSlackUpdateNorm() const
{
  return slackUpdateNorm;
}   

double 
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getRTOL() const
{
  return rtol;
}

double 
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getATOL() const
{
  return atol;
}

double 
LOCA::Bifurcation::PitchforkBord::StatusTest::SlackUpdateNorm::getTOL() const
{
  return tol;
}

