// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

  if ( Teuchos::is_null(updateVectorPtr) || (updateVectorPtr->length() != curSoln.length()) )
    updateVectorPtr = curSoln.clone();

  updateVectorPtr->update(1.0, curSoln, -1.0, oldSoln, 0.0);

  NOX::size_type n = (scaleType == Scaled) ? updateVectorPtr->length() : 0;

  switch (normType) {

  case NOX::Abstract::Vector::TwoNorm:
    normUpdate = updateVectorPtr->norm();
    if (scaleType == Scaled)
      normUpdate /= sqrt(1.0 * static_cast<double>(n));
    break;

  default:
    normUpdate = updateVectorPtr->norm(normType);
    if (scaleType == Scaled)
      normUpdate /= static_cast<double>(n);
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
