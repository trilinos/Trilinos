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

std::ostream& NormUpdate::print(std::ostream& stream, int indent) const
{
  for (int j = 0; j < indent; j ++)
    stream << ' ';
  stream << status;
  stream << "Absolute Update-Norm = " << Utils::sciformat(normUpdate, 3) 
	 << " < " << Utils::sciformat(tolerance, 3) << std::endl;
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
