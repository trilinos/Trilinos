// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Common.H"
#include "NOX_MeritFunction_SumOfSquares.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"

NOX::MeritFunction::SumOfSquares::
SumOfSquares(const Teuchos::RCP<NOX::Utils>& u) :
  meritFunctionName("Sum Of Squares (default): 0.5 * ||F|| * ||F||")
{
  utils = u;
}

NOX::MeritFunction::SumOfSquares::~SumOfSquares()
{
}

double NOX::MeritFunction::SumOfSquares::
computef(const NOX::Abstract::Group& grp) const
{
  if ( !(grp.isF()) ) {
    utils->err()
      << "ERROR: NOX::MeritFunction::SumOfSquares::computef() - "
      << "F has not been computed yet!.  Please call "
      << "computeF() on the group passed into this function."
      << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return (0.5 * grp.getNormF() * grp.getNormF());
}

void NOX::MeritFunction::SumOfSquares::
computeGradient(const NOX::Abstract::Group& grp,
        NOX::Abstract::Vector& result) const
{
  if ( !(grp.isF()) ) {
    utils->err()
      << "ERROR: NOX::MeritFunction::SumOfSquares::computeGradient() - "
      << "F has not been computed yet!.  Please call "
      << "computeF() on the group passed into this function."
      << std::endl;
    throw std::runtime_error("NOX Error");
  }

  if ( !(grp.isJacobian()) ) {
    utils->err()
      << "ERROR: NOX::MeritFunction::SumOfSquares::computeGradient() - "
      << "Jacobian has not been computed yet!.  Please call "
      << "computeJacobian() on the group passed into this function."
      << std::endl;
    throw std::runtime_error("NOX Error");
  }

  NOX::Abstract::Group::ReturnType status =
    grp.applyJacobianTranspose(grp.getF(), result);

  if (status != NOX::Abstract::Group::Ok) {
    utils->err() << "ERROR: NOX::MeritFunction::SumOfSquares::compute"
         << "Gradient - applyJacobianTranspose failed!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return;
}

double NOX::MeritFunction::SumOfSquares::
computeSlope(const NOX::Abstract::Vector& dir,
         const NOX::Abstract::Group& grp) const
{
  if (Teuchos::is_null(tmpVecPtr))
    tmpVecPtr = grp.getF().clone();

  // If the Jacobian is not computed, approximate it with
  // directional derivatives.  dir^T J^T F = F^T Jd
  if (!(grp.isJacobian()))
    return this->computeSlopeWithoutJacobian(dir, grp);
  // If the Jacobian is computed but doesn't support a gradient,
  // employ a different form for the inner product, eg
  // return <v, F> = F' * J * dir = <J'F, dir> = <g, dir>
  else if(!(grp.isGradient()))
    return this->computeSlopeWithoutJacobianTranspose(dir, grp);

  this->computeGradient(grp, *(tmpVecPtr.get()));

  return dir.innerProduct(*(tmpVecPtr.get()));
}

double NOX::MeritFunction::SumOfSquares::
computeQuadraticModel(const NOX::Abstract::Vector& dir,
              const NOX::Abstract::Group& grp) const
{
  if (Teuchos::is_null(tmpVecPtr))
    tmpVecPtr = grp.getF().clone();

  double m = 0.0;

  m = this->computef(grp);

  m += this->computeSlope(dir, grp);

  grp.applyJacobian(dir, *(tmpVecPtr.get()));

  m += 0.5 * tmpVecPtr->innerProduct(*(tmpVecPtr.get()));

  return m;
}

double NOX::MeritFunction::SumOfSquares::
computeSlopeWithoutJacobian(const NOX::Abstract::Vector& dir,
                const NOX::Abstract::Group& grp) const
{
  if (Teuchos::is_null(tmpVecPtr))
    tmpVecPtr = grp.getF().clone(NOX::ShapeCopy);

  if (Teuchos::is_null(tmpGrpPtr))
    tmpGrpPtr = grp.clone(NOX::ShapeCopy);


  // Compute the perturbation parameter
  double lambda = 1.0e-6;
  double denominator = dir.norm();

  // Don't divide by zero
  if (denominator == 0.0)
    denominator = 1.0;

  double eta = lambda * (lambda + grp.getX().norm() / denominator);

  // Don't divide by zero
  if (eta == 0.0)
    eta = 1.0e-6;

  // Perturb the solution vector
  tmpVecPtr->update(eta, dir, 1.0, grp.getX(), 0.0);

  // Compute the new F --> F(x + eta * dir)
  tmpGrpPtr->setX(*(tmpVecPtr.get()));
  tmpGrpPtr->computeF();

  // Compute Js = (F(x + eta * dir) - F(x))/eta
  tmpVecPtr->update(-1.0/eta, grp.getF(), 1.0/eta, tmpGrpPtr->getF(), 0.0);

  return(tmpVecPtr->innerProduct(grp.getF()));
}

double NOX::MeritFunction::SumOfSquares::
computeSlopeWithoutJacobianTranspose(const Abstract::Vector& dir,
                                     const Abstract::Group& grp) const
{
  // Allocate space for vecPtr if necessary
  if (Teuchos::is_null(tmpVecPtr))
    tmpVecPtr = grp.getF().clone(NOX::ShapeCopy);

  // v = J * dir
  NOX::Abstract::Group::ReturnType status = grp.applyJacobian(dir,*tmpVecPtr);

  if (status != NOX::Abstract::Group::Ok)
  {
    utils->out() << "NOX::MeritFunction::SumOfSquares::computeSlopeWithoutJacobianTranspose -  Unable to apply Jacobian!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Check that F exists
  if (!grp.isF())
  {
    utils->out() << "NOX::MeritFunction::SumOfSquares::computeSlopeWithoutJacobianTranspose - Invalid F" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Return <v, F> = F' * J * dir = <J'F, dir> = <g, dir>
  return(tmpVecPtr->innerProduct(grp.getF()));
}

void NOX::MeritFunction::SumOfSquares::
computeQuadraticMinimizer(const NOX::Abstract::Group& grp,
              NOX::Abstract::Vector& result) const
{
  // Clone a temporary vector
  if (Teuchos::is_null(tmpVecPtr))
    tmpVecPtr = grp.getF().clone(NOX::ShapeCopy);

  // Make sure the function and Jacobian have been evaluated
  if ( !(grp.isF()) ) {
    utils->err()
      << "ERROR: NOX::MeritFunction::SumOfSquares::"
      << "computeQuadraticMinimizer() - "
      << "F has not been computed yet!.  Please call "
      << "computeF() on the group passed into this function."
      << std::endl;
    throw std::runtime_error("NOX Error");
  }

  if ( !(grp.isJacobian()) ) {
    utils->err()
      << "ERROR: NOX::MeritFunction::SumOfSquares::"
      << "computeQuadraticMinimizer() - "
      << "Jacobian has not been computed yet!.  Please call "
      << "computeJacobian() on the group passed into this function."
      << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Compute the steepest descent direction = J^T F
  this->computeGradient(grp, result);

  // Compute = J (J^T F)
  NOX::Abstract::Group::ReturnType status =
    grp.applyJacobian(result, *tmpVecPtr);
  if (status != NOX::Abstract::Group::Ok) {
    utils->err()
      << "ERROR: NOX::MeritFunction::SumOfSquares::"
      << "computeQuadraticMinimizer() - grp->applyJacobian() has failed!"
      << std::endl;
    throw std::runtime_error("NOX Error");
  }

  result.scale( -1.0 * result.innerProduct(result) /
        tmpVecPtr->innerProduct(*tmpVecPtr) );

}

const std::string& NOX::MeritFunction::SumOfSquares::name() const
{
  return meritFunctionName;
}
