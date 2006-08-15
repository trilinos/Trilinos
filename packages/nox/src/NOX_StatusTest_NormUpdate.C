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

#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "NOX_Utils.H"

using namespace NOX;
using namespace NOX::StatusTest;

NormUpdate::NormUpdate(double tol, Abstract::Vector::NormType ntype, ScaleType stype) :
  status(Unevaluated),
  normType(ntype),
  scaleType(stype),
  tolerance(tol),
  normUpdate(0.0)
{
}

NormUpdate::NormUpdate(double tol, ScaleType stype) :
  status(Unevaluated),
  normType(NOX::Abstract::Vector::TwoNorm),
  scaleType(stype),
  tolerance(tol),
  normUpdate(0.0)
{
}

NormUpdate::~NormUpdate()
{

}

StatusType NormUpdate::checkStatus(const Solver::Generic& problem, 
				   NOX::StatusTest::CheckType checkType)
{
  if (checkType == None)
  {
    status = Unevaluated;
    normUpdate = -1.0;
    return status;
  }

  // On the first iteration, the old and current solution are the same so
  // we should return the test as unconverged until there is a valid 
  // old solution (i.e. the number of iterations is greater than zero).
  int niters = problem.getNumIterations();
  if (niters == 0) 
  {
    status = Unconverged;
    normUpdate = -1.0;
    return status;
  } 

  // Check that F exists!
  if (!problem.getSolutionGroup().isF())
  {
    status = Unconverged;
    normUpdate = -1.0;
    return status;
  } 

  const Abstract::Vector& oldSoln = problem.getPreviousSolutionGroup().getX();
  const Abstract::Vector& curSoln = problem.getSolutionGroup().getX();

  if (Teuchos::is_null(updateVectorPtr)) 
    updateVectorPtr = curSoln.clone();

  updateVectorPtr->update(1.0, curSoln, -1.0, oldSoln, 0.0); 

  int n = (scaleType == Scaled) ? updateVectorPtr->length() : 0;

  switch (normType) {
    
  case NOX::Abstract::Vector::TwoNorm:
    normUpdate = updateVectorPtr->norm();
    if (scaleType == Scaled)
      normUpdate /= sqrt(1.0 * n);
    break;

  default:
    normUpdate = updateVectorPtr->norm(normType);
    if (scaleType == Scaled)
      normUpdate /= n;
    break;

  }

  status = (normUpdate < tolerance) ? Converged : Unconverged;
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
