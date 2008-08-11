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

#include "LOCA_Bifurcation_TPBord_StatusTest_NullVectorNormWRMS.H" 
#include "LOCA_Bifurcation_TPBord_ExtendedGroup.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::NullVectorNormWRMS(
								double rtol_, 
								double atol_, 
								double tol_) :
  rtol(rtol_),
  atol(atol_),
  tol(tol_),
  normWRMS(0.0),
  status(NOX::StatusTest::Unconverged)
{
}


LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::~NullVectorNormWRMS()
{
}

NOX::StatusTest::StatusType 
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::checkStatus(
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
    normWRMS = 0.0;
    return NOX::StatusTest::Converged;
  }

  // Get solution vectors
  const LOCA::Bifurcation::TPBord::ExtendedVector& x = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(soln.getX());
  const LOCA::Bifurcation::TPBord::ExtendedVector& xold = 
    dynamic_cast<const LOCA::Bifurcation::TPBord::ExtendedVector&>(oldsoln.getX());

  // Get null vectors
  const NOX::Abstract::Vector& y = x.getNullVec();
  const NOX::Abstract::Vector& yold = xold.getNullVec();

  // temporary vectors
  NOX::Abstract::Vector *u = y.clone(NOX::ShapeCopy);
  NOX::Abstract::Vector *v = yold.clone(NOX::ShapeCopy);
  
  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid 
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0) 
  {
    normWRMS = 1.0e+12;
    status = NOX::StatusTest::Unconverged;
    return status;
  } 
  
  // Fill vector with 1's
  u->init(1.0);

  // Compute |y|
  v->abs(y);
  
  // Overwrite u with rtol*|y| + atol
  u->update(rtol, *v, atol);

  // Overwrite v with 1/(rtol*|y| + atol)
  v->reciprocal(*u);

  // Overwrite u with y-yold
  u->update(1.0, y, -1.0, yold, 0.0);

  // Overwrite u with (y-yold)/(rtol*|y| + atol)
  u->scale(*v);

  // Compute sqrt( (y-yold)/(rtol*|y| + atol) ) / sqrt(N)
  normWRMS = u->norm() / sqrt(static_cast<double>(u->length()));

  if (normWRMS < tol) 
    status = NOX::StatusTest::Converged;
  else
    status = NOX::StatusTest::Unconverged;

  delete u;
  delete v;
  
  return status;
}

NOX::StatusTest::StatusType 
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getStatus() const
{
  return status;
}


ostream& 
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::print(
							    ostream& stream, 
							    int indent) const
{
  for (int j = 0; j < indent; j++)
    stream << ' ';
  stream << status;
  stream << "Turning Point Scaled Null Vector Update = " 
	 << NOX::Utils::sciformat(normWRMS, 3) << " < " << tol;
  stream << endl;

  return stream;
}


double 
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getNullVectorNormWRMS() const
{
  return normWRMS;
}   

double 
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getRTOL() const
{
  return rtol;
}

double 
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getATOL() const
{
  return atol;
}

double 
LOCA::Bifurcation::TPBord::StatusTest::NullVectorNormWRMS::getTOL() const
{
  return tol;
}

