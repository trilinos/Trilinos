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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::StatusTest;

NormUpdate::NormUpdate(double tol, Abstract::Vector::NormType ntype, ScaleType stype) :
  status(Unconverged),
  updateVector(0),
  normType(ntype),
  scaleType(stype),
  tolerance(tol)
{
}

NormUpdate::NormUpdate(double tol, ScaleType stype) :
  status(Unconverged),
  updateVector(0),
  normType(NOX::Abstract::Vector::TwoNorm),
  scaleType(stype),
  tolerance(tol)
{
}

NormUpdate::~NormUpdate()
{
  delete updateVector;
}

StatusType NormUpdate::checkStatus(const Solver::Generic& problem)
{
  status = Unconverged;
  const Abstract::Group& grp = problem.getSolutionGroup();
  const Abstract::Vector& oldSoln = problem.getPreviousSolutionGroup().getX();
  const Abstract::Vector& curSoln = problem.getSolutionGroup().getX();

  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid 
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0) {
    status = Unconverged;
    normUpdate = -1.0;
    return status;
  } 

  if (updateVector == 0) 
    updateVector = curSoln.clone();

  *updateVector = curSoln;

  updateVector->update(-1.0, oldSoln, 1.0); 

  int n = grp.getX().length();

  switch (normType) {
    
  case NOX::Abstract::Vector::TwoNorm:
    normUpdate = updateVector->norm();
    if (scaleType == Scaled)
      normUpdate /= sqrt(1.0 * n);
    break;

  default:
    normUpdate = updateVector->norm(normType);
    if (scaleType == Scaled)
      normUpdate /= n;
    break;

  }

  if (normUpdate < tolerance)
    status = Converged;

  return status;
}

StatusType NormUpdate::getStatus() const
{
  return status;
}

ostream& NormUpdate::print(ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Absolute Update-Norm = " << Utils::sciformat(normUpdate, 3) 
	 << " < " << Utils::sciformat(tolerance, 3) << endl;
  return stream;
}

double NOX::StatusTest::NormUpdate::getNormUpdate() const
{
  return normUpdate;
}

double NOX::StatusTest::NormUpdate::getTolerance() const
{
  return tolerance;
}
